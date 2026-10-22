#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdint>
#include <exception>
#include <future>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>
#include <R_ext/Print.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

namespace {
// Cancellation is cooperative because R's API, including interrupt checks,
// may only be called by the main thread.  The main thread polls through
// R_ToplevelExec(), which contains R's longjmp, then workers observe this flag
// at bounded row/iteration boundaries and unwind using ordinary C++.
class FitCancelled : public std::exception {
 public:
  const char* what() const noexcept override { return "native GMM fit cancelled"; }
};

class CholeskyFailure : public std::runtime_error {
 public:
  explicit CholeskyFailure(const std::string& message)
    : std::runtime_error(message) {}
};

void check_r_interrupt(void*) {
  R_CheckUserInterrupt();
}

bool r_interrupt_pending() {
  return R_ToplevelExec(check_r_interrupt, nullptr) == FALSE;
}

inline void cancellation_point(const std::atomic<bool>* cancelled) {
  if (cancelled != nullptr && cancelled->load(std::memory_order_relaxed)) {
    throw FitCancelled();
  }
}

struct Data {
  int n, p;
  std::vector<double> x;
  std::vector<double> centered;
  // `x` stays unshifted so seeded K-means++ choices retain the established
  // implementation's exact distance/cumulative-sum path. EM uses `centered`
  // so likelihood and component-centred covariance updates are translation
  // stable.
  std::vector<double> center;
};
struct Fit {
  int k, p, iterations; bool converged;
  double reg, nll;
  double tolerance;
  std::vector<double> means, covs, weights, history;
  std::vector<double> center;
  std::vector<int> labels;
  std::vector<double> initialization_counts;
  std::vector<int> initialization_collapsed;
  int early_wait = 0;
  double early_best = std::numeric_limits<double>::infinity();
  // Diagnostic wall times never participate in numerical decisions.
  double kmeans_seeding_seconds = 0.0;
  double lloyd_seconds = 0.0;
  double initialization_covariance_seconds = 0.0;
  double expectation_seconds = 0.0;
  double covariance_seconds = 0.0;
  double final_scoring_seconds = 0.0;
};

struct ThreadJoiner {
  std::vector<std::thread>& threads;
  ~ThreadJoiner() {
    for (std::thread& thread : threads) {
      if (thread.joinable()) thread.join();
    }
  }
};

enum class JobOutcome { Success, Cancelled, Failed };

struct CompletionEvent {
  int index;
  double elapsed_seconds;
  JobOutcome outcome;
};

// Keep long-running native progress compact and readable without calling the
// R API from worker threads. Durations are diagnostic only and never affect a
// fit, its scheduling, or numerical results.
std::string format_progress_duration(double seconds) {
  if (!std::isfinite(seconds) || seconds < 0.0) return "unknown";
  const long long total = static_cast<long long>(std::llround(seconds));
  char buffer[64];
  if (total >= 86400) {
    std::snprintf(buffer, sizeof(buffer), "%lldd %lldh",
                  total / 86400, (total % 86400) / 3600);
  } else if (total >= 3600) {
    std::snprintf(buffer, sizeof(buffer), "%lldh %lldm",
                  total / 3600, (total % 3600) / 60);
  } else if (total >= 60) {
    std::snprintf(buffer, sizeof(buffer), "%lldm %llds",
                  total / 60, total % 60);
  } else {
    std::snprintf(buffer, sizeof(buffer), "%.1fs", seconds);
  }
  return std::string(buffer);
}

inline double &at(std::vector<double>& a, int cols, int i, int j) {
  return a[static_cast<size_t>(i) * cols + j];
}
inline double atc(const std::vector<double>& a, int cols, int i, int j) {
  return a[static_cast<size_t>(i) * cols + j];
}

// Lower Cholesky. Increasing the ridge is handled by the outer complete-fit
// restart, so this routine never silently changes the fitted distribution.
bool cholesky(const double* covariance, int p, double ridge,
              std::vector<double>& lower) {
  lower.assign(static_cast<size_t>(p) * p, 0.0);
  for (int i = 0; i < p; ++i) for (int j = 0; j <= i; ++j) {
    double value = covariance[static_cast<size_t>(i) * p + j];
    if (i == j) value += ridge;
    for (int q = 0; q < j; ++q)
      value -= atc(lower, p, i, q) * atc(lower, p, j, q);
    if (i == j) {
      if (!(value > 0.0) || !std::isfinite(value)) return false;
      at(lower, p, i, j) = std::sqrt(value);
    } else at(lower, p, i, j) = value / atc(lower, p, j, j);
  }
  return true;
}

std::vector<std::vector<double> > factors(const Fit& f) {
  std::vector<std::vector<double> > out(f.k);
  for (int c = 0; c < f.k; ++c)
    if (!cholesky(&f.covs[static_cast<size_t>(c) * f.p * f.p], f.p, f.reg, out[c]))
      throw CholeskyFailure("Cholesky factorization failed");
  return out;
}

struct ScoringCache {
  // Column-major lower-triangular precision-Cholesky factors.  The rest of
  // the native engine stores observation rows contiguously; a row block is
  // therefore also a column-major p-by-n matrix and can be transformed in
  // place with DTRMM without another layout conversion.
  std::vector<std::vector<double> > precision_lower;
  std::vector<double> log_weights;
  std::vector<double> gaussian_bases;
};

std::vector<double> invert_lower(const std::vector<double>& lower, int p) {
  std::vector<double> inverse(static_cast<size_t>(p)*p,0.0);
  for(int column=0;column<p;++column) {
    for(int row=column;row<p;++row) {
      double value=row==column ? 1.0 : 0.0;
      for(int q=column;q<row;++q)
        value-=atc(lower,p,row,q)*atc(inverse,p,q,column);
      at(inverse,p,row,column)=value/atc(lower,p,row,row);
    }
  }
  return inverse;
}

ScoringCache make_scoring_cache(const Fit& f) {
  ScoringCache out;
  const auto chol=factors(f);
  out.precision_lower.resize(f.k);
  out.log_weights.resize(f.k);
  out.gaussian_bases.resize(f.k);
  const double gaussian_constant =
    f.p * std::log(2.0 * 3.14159265358979323846);
  for (int c=0;c<f.k;++c) {
    double logdet=0.0;
    for (int j=0;j<f.p;++j)
      logdet+=std::log(atc(chol[c],f.p,j,j));
    // PyCave likewise precomputes a precision-Cholesky factor. This removes
    // the dependent triangular solve and scalar divisions from every
    // observation/component score while retaining the same quadratic form.
    const std::vector<double> precision_row_major=
      invert_lower(chol[c],f.p);
    out.precision_lower[c].assign(static_cast<size_t>(f.p)*f.p,0.0);
    for(int row=0;row<f.p;++row) for(int column=0;column<=row;++column)
      out.precision_lower[c][static_cast<size_t>(row)+
        static_cast<size_t>(column)*f.p]=
        precision_row_major[static_cast<size_t>(row)*f.p+column];
    out.log_weights[c]=std::log(f.weights[c]);
    out.gaussian_bases[c]=gaussian_constant+2.0*logdet;
  }
  return out;
}

// Score one deterministic row shard in bounded blocks.  With BLAS configured
// for one thread, native fit threads remain the sole source of parallelism.
// The posterior normalization, likelihood reduction, and sufficient-statistic
// accumulation retain their original observation and component order.
void score_rows(const Data& d, int begin, int end, const Fit& f,
                const ScoringCache& cache,
                std::vector<double>* posterior,
                std::vector<double>* local_nk,
                std::vector<double>* local_sums,
                std::vector<int>* labels,
                double& log_likelihood,
                const std::atomic<bool>* cancelled) {
  const int block_rows=4096;
  std::vector<double> transformed(
    static_cast<size_t>(std::min(block_rows,end-begin))*f.p);
  std::vector<double> scores(
    static_cast<size_t>(std::min(block_rows,end-begin))*f.k);
  const char side='L',uplo='L',transa='N',diag='N';
  const double alpha=1.0;
  for(int block_begin=begin;block_begin<end;block_begin+=block_rows) {
    cancellation_point(cancelled);
    const int rows=std::min(block_rows,end-block_begin);
    for(int c=0;c<f.k;++c) {
      const double* component_mean=&f.means[static_cast<size_t>(c)*f.p];
      for(int i=0;i<rows;++i) {
        const double* observation=&d.centered[
          static_cast<size_t>(block_begin+i)*f.p];
        double* destination=&transformed[static_cast<size_t>(i)*f.p];
        for(int j=0;j<f.p;++j)
          destination[j]=observation[j]-component_mean[j];
      }
      // `transformed` is row-major rows-by-p, which is the same memory layout
      // as a column-major p-by-rows matrix. Apply L^{-1} to every residual.
      F77_CALL(dtrmm)(&side,&uplo,&transa,&diag,&f.p,&rows,&alpha,
                      cache.precision_lower[c].data(),&f.p,
                      transformed.data(),&f.p
                      FCONE FCONE FCONE FCONE);
      for(int i=0;i<rows;++i) {
        const double* value=&transformed[static_cast<size_t>(i)*f.p];
        double quadratic=0.0;
        for(int j=0;j<f.p;++j) quadratic+=value[j]*value[j];
        scores[static_cast<size_t>(i)*f.k+c]=cache.log_weights[c]-
          0.5*(cache.gaussian_bases[c]+quadratic);
      }
    }
    for(int i=0;i<rows;++i) {
      double* values=&scores[static_cast<size_t>(i)*f.k];
      double maximum=-std::numeric_limits<double>::infinity();
      for(int c=0;c<f.k;++c) maximum=std::max(maximum,values[c]);
      double total=0.0;
      for(int c=0;c<f.k;++c) {
        values[c]=std::exp(values[c]-maximum);
        total+=values[c];
      }
      log_likelihood+=maximum+std::log(total);
      int best=0;
      for(int c=0;c<f.k;++c) {
        values[c]/=total;
        if(values[c]>values[best]) best=c;
      }
      const int row=block_begin+i;
      if(labels != nullptr) (*labels)[row]=best;
      if(posterior != nullptr) {
        for(int c=0;c<f.k;++c) {
          (*posterior)[static_cast<size_t>(row)*f.k+c]=values[c];
          (*local_nk)[c]+=values[c];
        }
        const double* observation=&d.centered[static_cast<size_t>(row)*f.p];
        for(int feature=0;feature<f.p;++feature) {
          double* feature_sum=&(*local_sums)[static_cast<size_t>(feature)*f.k];
          const double x=observation[feature];
          for(int c=0;c<f.k;++c) feature_sum[c]+=values[c]*x;
        }
      }
    }
  }
}

void nearest_assignments(const Data& d, const std::vector<double>& means,
                         int k, int threads, std::vector<int>& assignments,
                         const std::atomic<bool>* cancelled) {
  const int nt=std::max(1,std::min(threads,d.n));
  std::vector<std::thread> pool;
  std::vector<std::exception_ptr> errors(nt);
  pool.reserve(nt);
  ThreadJoiner joiner{pool};
  for(int t=0;t<nt;++t) pool.emplace_back([&,t](){
    try {
      const int begin=static_cast<int>((static_cast<int64_t>(d.n)*t)/nt);
      const int end=static_cast<int>((static_cast<int64_t>(d.n)*(t+1))/nt);
      for(int i=begin;i<end;++i) {
        if((i&1023)==0) cancellation_point(cancelled);
        int best=0;
        double best_distance=std::numeric_limits<double>::infinity();
        for(int c=0;c<k;++c) {
          double value=0.0;
          for(int j=0;j<d.p;++j) {
            const double delta=d.x[static_cast<size_t>(i)*d.p+j]-
              means[static_cast<size_t>(c)*d.p+j];
            value+=delta*delta;
          }
          // Strict comparison deliberately preserves first-component tie
          // resolution from the serial implementation.
          if(value<best_distance) { best_distance=value; best=c; }
        }
        assignments[i]=best;
      }
    } catch(...) { errors[t]=std::current_exception(); }
  });
  for(auto& thread:pool) if(thread.joinable()) thread.join();
  for(const auto& error:errors) if(error) std::rethrow_exception(error);
}

// std::uniform_*_distribution is permitted to map an engine's bits
// differently across standard-library implementations. Use the standardized
// mt19937_64 sequence directly so a seed has the same meaning on every
// supported compiler/runtime.
uint64_t bounded_random(std::mt19937_64& rng, uint64_t bound) {
  const uint64_t limit = std::numeric_limits<uint64_t>::max() -
    (std::numeric_limits<uint64_t>::max() % bound);
  uint64_t value;
  do { value = rng(); } while (value >= limit);
  return value % bound;
}

double unit_random(std::mt19937_64& rng) {
  return static_cast<double>(rng() >> 11) * (1.0 / 9007199254740992.0);
}

int covariance_block_rows(int features) {
  // Keep each thread's residual workspace near 1 MiB. The former fixed
  // 16,384-row block occupied 5 MiB at P=40 and 15 MiB at P=120, so several
  // concurrent fits displaced one another from the shared last-level cache.
  // Bounds retain efficient BLAS calls for both very small and large P; the
  // 256-row alignment makes the choice stable and easy to audit.
  const size_t target_bytes=static_cast<size_t>(1)*1024*1024;
  int rows=static_cast<int>(target_bytes/
    (sizeof(double)*static_cast<size_t>(features)));
  rows=std::max(1024,std::min(8192,rows));
  return std::max(1024,(rows/256)*256);
}

std::vector<double> kmeans_pp(const Data& d, int k, uint64_t seed, int threads,
                              double* seeding_seconds,
                              double* lloyd_seconds,
                              const std::atomic<bool>* cancelled = nullptr) {
  const auto seeding_started=std::chrono::steady_clock::now();
  std::mt19937_64 rng(seed); // private RNG: does not consume or mutate R's RNG.
  std::vector<double> means(static_cast<size_t>(k) * d.p);
  const int first = static_cast<int>(bounded_random(rng, static_cast<uint64_t>(d.n)));
  std::copy_n(&d.x[static_cast<size_t>(first) * d.p], d.p, means.begin());
  std::vector<double> distance(d.n, std::numeric_limits<double>::infinity());
  for (int c = 1; c < k; ++c) {
    cancellation_point(cancelled);
    double sum = 0.0;
    for (int i = 0; i < d.n; ++i) {
      if ((i & 1023) == 0) cancellation_point(cancelled);
      double current = 0.0;
      for (int j = 0; j < d.p; ++j) {
        double z = d.x[static_cast<size_t>(i) * d.p + j] - means[static_cast<size_t>(c - 1) * d.p + j];
        current += z * z;
      }
      distance[i] = std::min(distance[i], current); sum += distance[i];
    }
    int selected = c % d.n;
    if (sum > 0.0 && std::isfinite(sum)) {
      double target = unit_random(rng) * sum, cumulative = 0.0;
      for (int i = 0; i < d.n; ++i) if ((cumulative += distance[i]) >= target) { selected = i; break; }
    }
    std::copy_n(&d.x[static_cast<size_t>(selected) * d.p], d.p,
                means.begin() + static_cast<size_t>(c) * d.p);
  }
  if(seeding_seconds != nullptr)
    *seeding_seconds=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-seeding_started).count();
  // Complete the initializer with deterministic Lloyd passes. The final hard
  // assignments are summarized again below rather than treated as soft EM.
  const auto lloyd_started=std::chrono::steady_clock::now();
  std::vector<int> assignment(d.n, 0);
  std::vector<int> next_assignment(d.n,0);
  for (int iteration = 0; iteration < 30; ++iteration) {
    cancellation_point(cancelled);
    std::vector<double> sums(static_cast<size_t>(k) * d.p, 0.0), counts(k, 0.0);
    bool changed = false;
    nearest_assignments(d,means,k,threads,next_assignment,cancelled);
    for (int i = 0; i < d.n; ++i) {
      if ((i & 1023) == 0) cancellation_point(cancelled);
      const int best=next_assignment[i];
      changed = changed || iteration == 0 || assignment[i] != best;
      assignment[i] = best; counts[best] += 1.0;
      for (int j = 0; j < d.p; ++j) sums[static_cast<size_t>(best)*d.p+j] += d.x[static_cast<size_t>(i)*d.p+j];
    }
    for (int c = 0; c < k; ++c) if (counts[c] > 0.0) for (int j = 0; j < d.p; ++j)
      means[static_cast<size_t>(c)*d.p+j] = sums[static_cast<size_t>(c)*d.p+j] / counts[c];
    if (!changed) break;
  }
  if(lloyd_seconds != nullptr)
    *lloyd_seconds=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-lloyd_started).count();
  return means;
}

std::vector<double> covariance_from_hard(
    const Data& d, int k, const std::vector<int>& hard,
    const std::vector<double>& counts, const std::vector<double>& means,
    int threads, const std::atomic<bool>* cancelled) {
  std::vector<double> output(static_cast<size_t>(k)*d.p*d.p,0.0);
  std::atomic<int> next_component(0);
  const int nt=std::max(1,std::min(threads,k));
  const int block_rows=covariance_block_rows(d.p);
  std::vector<std::thread> pool;
  std::vector<std::exception_ptr> errors(nt);
  pool.reserve(nt);
  ThreadJoiner joiner{pool};
  for(int t=0;t<nt;++t) pool.emplace_back([&,t](){
    try {
      std::vector<double> residuals(static_cast<size_t>(block_rows)*d.p);
      std::vector<double> cross(static_cast<size_t>(d.p)*d.p);
      for(;;) {
        const int c=next_component.fetch_add(1);
        if(c>=k) break;
        cancellation_point(cancelled);
        double* target=&output[static_cast<size_t>(c)*d.p*d.p];
        if(counts[c]<1e-6) {
          for(int j=0;j<d.p;++j) target[static_cast<size_t>(j)*d.p+j]=1.0;
          continue;
        }
        std::fill(cross.begin(),cross.end(),0.0);
        int buffered=0;
        bool first=true;
        const double* component_mean=&means[static_cast<size_t>(c)*d.p];
        auto flush=[&](){
          if(buffered==0) return;
          const char uplo='L',trans='N';
          const double alpha=1.0,beta=first ? 0.0 : 1.0;
          F77_CALL(dsyrk)(&uplo,&trans,&d.p,&buffered,&alpha,
                          residuals.data(),&d.p,&beta,cross.data(),&d.p
                          FCONE FCONE);
          first=false;
          buffered=0;
        };
        for(int i=0;i<d.n;++i) {
          if((i&1023)==0) cancellation_point(cancelled);
          if(hard[i]!=c) continue;
          const double* observation=&d.x[static_cast<size_t>(i)*d.p];
          double* destination=&residuals[static_cast<size_t>(buffered)*d.p];
          for(int j=0;j<d.p;++j)
            destination[j]=observation[j]-component_mean[j];
          if(++buffered==block_rows) flush();
        }
        flush();
        for(int row=0;row<d.p;++row) for(int column=0;column<=row;++column) {
          const double value=cross[static_cast<size_t>(row)+
            static_cast<size_t>(column)*d.p]/counts[c];
          target[static_cast<size_t>(row)*d.p+column]=value;
          target[static_cast<size_t>(column)*d.p+row]=value;
        }
      }
    } catch(...) { errors[t]=std::current_exception(); }
  });
  for(auto& thread:pool) if(thread.joinable()) thread.join();
  for(const auto& error:errors) if(error) std::rethrow_exception(error);
  return output;
}

std::vector<double> covariance_from_posterior(
    const Data& d, int k, const std::vector<double>& posterior,
    const std::vector<double>& nk, const std::vector<double>& means,
    const std::vector<bool>& collapsed, int threads,
    const std::atomic<bool>* cancelled) {
  std::vector<double> output(static_cast<size_t>(k)*d.p*d.p,0.0);
  std::atomic<int> next_component(0);
  const int nt=std::max(1,std::min(threads,k));
  const int block_rows=covariance_block_rows(d.p);
  std::vector<std::thread> pool;
  std::vector<std::exception_ptr> errors(nt);
  pool.reserve(nt);
  ThreadJoiner joiner{pool};
  for(int t=0;t<nt;++t) pool.emplace_back([&,t](){
    try {
      std::vector<double> weighted(static_cast<size_t>(block_rows)*d.p);
      std::vector<double> cross(static_cast<size_t>(d.p)*d.p);
      for(;;) {
        const int c=next_component.fetch_add(1);
        if(c>=k) break;
        cancellation_point(cancelled);
        double* target=&output[static_cast<size_t>(c)*d.p*d.p];
        if(collapsed[c]) {
          for(int j=0;j<d.p;++j) target[static_cast<size_t>(j)*d.p+j]=1.0;
          continue;
        }
        std::fill(cross.begin(),cross.end(),0.0);
        bool first=true;
        for(int begin=0;begin<d.n;begin+=block_rows) {
          cancellation_point(cancelled);
          const int rows=std::min(block_rows,d.n-begin);
          for(int i=0;i<rows;++i) {
            const double scale=std::sqrt(
              posterior[static_cast<size_t>(begin+i)*k+c]);
            const double* source=&d.centered[static_cast<size_t>(begin+i)*d.p];
            double* destination=&weighted[static_cast<size_t>(i)*d.p];
            // Form component-centred residuals before the rank-k update.
            // Computing E[xx'] - mean mean' would be algebraically equivalent
            // in exact arithmetic, but catastrophically cancels when a
            // component is far from the cohort-wide centre.
            const double* component_mean=
              &means[static_cast<size_t>(c)*d.p];
            for(int j=0;j<d.p;++j)
              destination[j]=scale*(source[j]-component_mean[j]);
          }
          const char uplo='L',trans='N';
          const double alpha=1.0,beta=first ? 0.0 : 1.0;
          F77_CALL(dsyrk)(&uplo,&trans,&d.p,&rows,&alpha,weighted.data(),
                          &d.p,&beta,cross.data(),&d.p FCONE FCONE);
          first=false;
        }
        for(int a=0;a<d.p;++a) for(int b=0;b<=a;++b) {
          // `cross` is column-major and already contains the weighted outer
          // products of component-centred residuals.
          const double value=cross[static_cast<size_t>(a)+
            static_cast<size_t>(b)*d.p]/nk[c];
          target[static_cast<size_t>(a)*d.p+b]=value;
          target[static_cast<size_t>(b)*d.p+a]=value;
        }
      }
    } catch(...) { errors[t]=std::current_exception(); }
  });
  for(auto& thread:pool) if(thread.joinable()) thread.join();
  for(const auto& error:errors) if(error) std::rethrow_exception(error);
  return output;
}

Fit fit_once(const Data& d, int k, uint64_t seed, int max_iter,
             double tolerance, double reg, int threads_per_fit,
             const std::atomic<bool>* cancelled = nullptr) {
  cancellation_point(cancelled);
  Fit f; f.k = k; f.p = d.p; f.reg = reg; f.iterations = 0; f.converged = false;
  f.tolerance = tolerance;
  f.center = d.center;
  const double collapse_threshold = 1e-6;
  std::vector<double> centers = kmeans_pp(
    d,k,seed,threads_per_fit,&f.kmeans_seeding_seconds,&f.lloyd_seconds,
    cancelled);
  f.means.assign(static_cast<size_t>(k) * d.p, 0.0); f.weights.assign(k, 0.0);
  f.covs.assign(static_cast<size_t>(k) * d.p * d.p, 0.0);
  std::vector<int> hard(d.n, 0); f.initialization_counts.assign(k, 0.0);
  nearest_assignments(d,centers,k,threads_per_fit,hard,cancelled);
  for (int i = 0; i < d.n; ++i) {
    if ((i & 1023) == 0) cancellation_point(cancelled);
    const int best=hard[i];
    f.initialization_counts[best] += 1.0;
    for(int j=0;j<d.p;++j) f.means[static_cast<size_t>(best)*d.p+j] += d.x[static_cast<size_t>(i)*d.p+j];
  }
  for(int c=0;c<k;++c) {
    if(f.initialization_counts[c] < collapse_threshold) {
      f.initialization_collapsed.push_back(c);
      int selected=static_cast<int>((seed+c)%d.n);
      std::copy_n(&d.x[static_cast<size_t>(selected)*d.p],d.p,f.means.begin()+static_cast<size_t>(c)*d.p);
    } else for(int j=0;j<d.p;++j) f.means[static_cast<size_t>(c)*d.p+j] /= f.initialization_counts[c];
    f.weights[c]=std::max(f.initialization_counts[c],collapse_threshold);
  }
  double initial_total=std::accumulate(f.weights.begin(),f.weights.end(),0.0); for(double& w:f.weights) w/=initial_total;
  const auto initialization_covariance_started=
    std::chrono::steady_clock::now();
  f.covs=covariance_from_hard(d,k,hard,f.initialization_counts,f.means,
                              threads_per_fit,cancelled);
  f.initialization_covariance_seconds=std::chrono::duration<double>(
    std::chrono::steady_clock::now()-initialization_covariance_started).count();
  // Initialization above deliberately follows the established unshifted
  // arithmetic. EM then proceeds in a shared centered coordinate system.
  for(int c=0;c<k;++c) for(int j=0;j<d.p;++j)
    f.means[static_cast<size_t>(c)*d.p+j]-=d.center[j];
  double best = std::numeric_limits<double>::infinity(); int wait = 0;
  std::vector<double> nk(k), sums(static_cast<size_t>(k) * d.p);
  std::vector<double> posterior(static_cast<size_t>(d.n)*k);
  const int nt = std::max(1, std::min(threads_per_fit, d.n));
  for (int iteration = 0; iteration < max_iter; ++iteration) {
    cancellation_point(cancelled);
    const auto expectation_started=std::chrono::steady_clock::now();
    const ScoringCache scoring_cache=make_scoring_cache(f);
    std::fill(nk.begin(), nk.end(), 0.0); std::fill(sums.begin(), sums.end(), 0.0);
    double ll = 0.0;
    // Fixed contiguous shards plus shard-order reduction are deterministic for
    // a fixed topology. Worker threads never call the R API.
    std::vector<std::vector<double> > local_nk(nt,std::vector<double>(k,0.0));
    // Feature-major storage makes the component dimension contiguous while
    // each thread writes only to its private deterministic shard.
    std::vector<std::vector<double> > local_sums(
      nt, std::vector<double>(static_cast<size_t>(d.p) * k, 0.0));
    std::vector<double> local_ll(nt,0.0); std::vector<std::thread> pass1;
    std::vector<std::exception_ptr> pass1_errors(nt);
    pass1.reserve(nt); ThreadJoiner pass1_joiner{pass1};
    for(int t=0;t<nt;++t) pass1.emplace_back([&,t](){
      try {
        int begin=static_cast<int>((static_cast<int64_t>(d.n)*t)/nt), end=static_cast<int>((static_cast<int64_t>(d.n)*(t+1))/nt);
        score_rows(d,begin,end,f,scoring_cache,&posterior,&local_nk[t],
                   &local_sums[t],nullptr,local_ll[t],cancelled);
      } catch(...) { pass1_errors[t]=std::current_exception(); }
    });
    for(auto& thread:pass1) if(thread.joinable()) thread.join();
    for(const auto& error:pass1_errors) if(error) std::rethrow_exception(error);
    for(int t=0;t<nt;++t){
      ll+=local_ll[t];
      for(int c=0;c<k;++c) nk[c]+=local_nk[t][c];
      for(int j=0;j<d.p;++j) for(int c=0;c<k;++c)
        sums[static_cast<size_t>(c)*d.p+j] +=
          local_sums[t][static_cast<size_t>(j)*k+c];
    }
    std::vector<double> new_means(sums.size()), new_weights(k);
    std::vector<bool> collapsed(k, false);
    for (int c = 0; c < k; ++c) {
      if (nk[c] < 1e-6) {
        collapsed[c] = true;
        int selected = static_cast<int>((seed + iteration + c) % d.n);
        std::copy_n(&d.centered[static_cast<size_t>(selected) * d.p], d.p,
                    new_means.begin() + static_cast<size_t>(c) * d.p);
        nk[c] = 1e-6;
      } else for (int j = 0; j < d.p; ++j)
        new_means[static_cast<size_t>(c) * d.p + j] =
          sums[static_cast<size_t>(c) * d.p + j] / nk[c];
      new_weights[c] = std::max(nk[c] / d.n, 1e-15);
    }
    double weight_sum = std::accumulate(new_weights.begin(), new_weights.end(), 0.0);
    for (double& w : new_weights) w /= weight_sum;
    f.expectation_seconds+=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-expectation_started).count();
    const auto covariance_started=std::chrono::steady_clock::now();
    std::vector<double> new_covs=covariance_from_posterior(
      d,k,posterior,nk,new_means,collapsed,nt,cancelled);
    f.covariance_seconds+=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-covariance_started).count();
    f.means.swap(new_means); f.weights.swap(new_weights); f.covs.swap(new_covs);
    f.iterations = iteration + 1; f.history.push_back(-ll / d.n);
    if (f.history.back() < best - tolerance) { best = f.history.back(); wait = 0; }
    else if (++wait >= 2) { f.converged = true; break; }
  }
  const auto final_scoring_started=std::chrono::steady_clock::now();
  const auto scoring_cache=make_scoring_cache(f);
  double ll = 0.0;
  cancellation_point(cancelled);
  f.labels.resize(d.n);
  std::vector<double> final_ll(nt,0.0); std::vector<std::thread> scoring;
  std::vector<std::exception_ptr> scoring_errors(nt);
  scoring.reserve(nt); ThreadJoiner scoring_joiner{scoring};
  for(int t=0;t<nt;++t) scoring.emplace_back([&,t](){
    try {
      const int begin=static_cast<int>((static_cast<int64_t>(d.n)*t)/nt);
      const int end=static_cast<int>((static_cast<int64_t>(d.n)*(t+1))/nt);
      score_rows(d,begin,end,f,scoring_cache,nullptr,nullptr,nullptr,&f.labels,
                 final_ll[t],cancelled);
    } catch(...) { scoring_errors[t]=std::current_exception(); }
  });
  for(auto& thread:scoring) if(thread.joinable()) thread.join();
  for(const auto& error:scoring_errors) if(error) std::rethrow_exception(error);
  for(int t=0;t<nt;++t) ll+=final_ll[t];
  f.final_scoring_seconds=std::chrono::duration<double>(
    std::chrono::steady_clock::now()-final_scoring_started).count();
  f.nll = -ll / d.n;
  f.early_wait = wait; f.early_best = best;
  return f;
}

Fit fit_restart(const Data& d, int k, uint64_t seed, int max_iter,
                double tolerance, double requested_reg, int threads_per_fit,
                const std::atomic<bool>* cancelled = nullptr) {
  std::string error;
  for (int attempt = 0; attempt < 8; ++attempt) try {
    cancellation_point(cancelled);
    return fit_once(d, k, seed, max_iter, tolerance,
                    requested_reg * std::pow(10.0, attempt), threads_per_fit,
                    cancelled);
  } catch (const CholeskyFailure& e) { error = e.what(); }
  throw std::runtime_error("GMM failed after 8 regularization attempts: " + error);
}

Data import_matrix(const Rcpp::NumericMatrix& x) {
  Data d; d.n = x.nrow(); d.p = x.ncol();
  d.x.resize(static_cast<size_t>(d.n) * d.p);
  d.centered.resize(static_cast<size_t>(d.n) * d.p);
  d.center.assign(d.p, 0.0);
  for (int j=0;j<d.p;++j) {
    long double total=0.0L;
    for (int i=0;i<d.n;++i) {
      double value=x(i,j);
      if (!std::isfinite(value)) Rcpp::stop("`x` contains non-finite values.");
      total+=static_cast<long double>(value);
    }
    d.center[j]=static_cast<double>(total/static_cast<long double>(d.n));
  }
  for (int i = 0; i < d.n; ++i) for (int j = 0; j < d.p; ++j) {
    d.x[static_cast<size_t>(i) * d.p + j] = x(i,j);
    d.centered[static_cast<size_t>(i) * d.p + j] = x(i,j)-d.center[j];
  }
  return d;
}

Rcpp::List export_fit(const Fit& f, double requested_reg, uint64_t seed,
                      double elapsed_seconds) {
  Rcpp::NumericMatrix means(f.k, f.p);
  for (int c = 0; c < f.k; ++c) {
    for (int j = 0; j < f.p; ++j) {
      means(c, j) = f.means[static_cast<size_t>(c) * f.p + j] + f.center[j];
    }
  }
  Rcpp::List covs(f.k);
  for (int c = 0; c < f.k; ++c) {
    Rcpp::NumericMatrix matrix(f.p, f.p);
    for (int row = 0; row < f.p; ++row) {
      for (int column = 0; column < f.p; ++column) {
        matrix(row, column) = f.covs[
          (static_cast<size_t>(c) * f.p + row) * f.p + column
        ];
      }
    }
    covs[c] = matrix;
  }
  return Rcpp::List::create(Rcpp::_["means"]=means, Rcpp::_["covariances"]=covs,
    Rcpp::_["weights"]=f.weights, Rcpp::_["n_features"]=f.p, Rcpp::_["n_clusters"]=f.k,
    Rcpp::_["nll"]=f.nll, Rcpp::_["iterations"]=f.iterations, Rcpp::_["converged"]=f.converged,
    Rcpp::_["covariance_regularization"]=f.reg,
    Rcpp::_["requested_covariance_regularization"]=requested_reg,
    Rcpp::_["nll_history"]=f.history, Rcpp::_["labels"]=f.labels,
    Rcpp::_["initialization_counts"]=f.initialization_counts,
    Rcpp::_["initialization_collapsed_components"]=f.initialization_collapsed,
    Rcpp::_["early_stopping"]=Rcpp::List::create(Rcpp::_["monitor"]="mean_nll_before_m_step",
      Rcpp::_["min_delta"]=f.tolerance, Rcpp::_["patience"]=2,
      Rcpp::_["wait_count"]=f.early_wait,Rcpp::_["best_nll"]=f.early_best),
    Rcpp::_["seed"]=static_cast<double>(seed),
    Rcpp::_["elapsed_seconds"]=elapsed_seconds,
    Rcpp::_["phase_elapsed_seconds"]=Rcpp::List::create(
      Rcpp::_["kmeans_seeding"]=f.kmeans_seeding_seconds,
      Rcpp::_["lloyd"]=f.lloyd_seconds,
      Rcpp::_["initialization_covariance"]=
        f.initialization_covariance_seconds,
      Rcpp::_["expectation"]=f.expectation_seconds,
      Rcpp::_["covariance"]=f.covariance_seconds,
      Rcpp::_["final_scoring"]=f.final_scoring_seconds));
}
}

extern "C" SEXP C_cpp_gmm_fit(SEXP xs, SEXP ks, SEXP seeds, SEXP maxs, SEXP tols, SEXP regs, SEXP threadss) {
  BEGIN_RCPP
  const auto fit_started=std::chrono::steady_clock::now();
  Data d=import_matrix(Rcpp::as<Rcpp::NumericMatrix>(xs)); int k=Rcpp::as<int>(ks);
  uint64_t seed=static_cast<uint64_t>(Rcpp::as<double>(seeds));
  Fit fit=fit_restart(d,k,seed,Rcpp::as<int>(maxs),Rcpp::as<double>(tols),
                      Rcpp::as<double>(regs),Rcpp::as<int>(threadss));
  const double elapsed_seconds=std::chrono::duration<double>(
    std::chrono::steady_clock::now()-fit_started).count();
  return export_fit(fit,Rcpp::as<double>(regs),seed,elapsed_seconds);
  END_RCPP
}

extern "C" SEXP C_cpp_gmm_fit_many(SEXP xs, SEXP kss, SEXP seedss, SEXP maxs, SEXP tols, SEXP regs, SEXP workerss, SEXP threadss, SEXP progresss) {
  BEGIN_RCPP
  const bool progress=Rcpp::as<bool>(progresss);
  const auto native_started=std::chrono::steady_clock::now();
  Data d=import_matrix(Rcpp::as<Rcpp::NumericMatrix>(xs)); Rcpp::IntegerVector ks(kss); Rcpp::NumericVector seeds(seedss);
  int jobs=ks.size(), workers=std::max(1,Rcpp::as<int>(workerss)), max_iter=Rcpp::as<int>(maxs);
  double tol=Rcpp::as<double>(tols), reg=Rcpp::as<double>(regs); int threads_per_fit=std::max(1,Rcpp::as<int>(threadss));
  std::vector<int> native_ks(jobs);
  std::vector<uint64_t> native_seeds(jobs);
  std::vector<int> job_order(jobs);
  for(int i=0;i<jobs;++i) {
    native_ks[i]=ks[i];
    native_seeds[i]=static_cast<uint64_t>(seeds[i]);
    job_order[i]=i;
  }
  // Independent fits are returned in canonical caller order, but starting the
  // largest K first prevents one expensive fit from becoming a serial tail
  // after the other worker threads have gone idle.
  std::stable_sort(job_order.begin(),job_order.end(),[&](int left,int right){
    return native_ks[left]>native_ks[right];
  });
  if(progress) {
    const double copy_seconds=std::chrono::duration<double>(
      std::chrono::steady_clock::now()-native_started).count();
    Rprintf("CellCharterR C++ GMM: native input copy complete in %.1fs; starting %d fits.\n",
            copy_seconds,jobs);
    R_FlushConsole();
  }
  std::vector<Fit> fits(jobs); std::vector<std::string> errors(jobs);
  std::vector<double> fit_elapsed_seconds(jobs,0.0);
  std::atomic<int> next(0);
  std::atomic<bool> cancelled(false);
  std::atomic<int> workers_finished(0);
  std::mutex progress_mutex;
  std::condition_variable progress_condition;
  std::vector<CompletionEvent> completion_events(jobs);
  int published_events=0;
  auto work=[&]() {
    for (;;) {
      if (cancelled.load(std::memory_order_relaxed)) break;
      int position=next.fetch_add(1);
      if(position>=jobs || cancelled.load(std::memory_order_relaxed)) break;
      int i=job_order[position];
      JobOutcome outcome=JobOutcome::Success;
      const auto fit_started=std::chrono::steady_clock::now();
      try {
        fits[i]=fit_restart(d,native_ks[i],native_seeds[i],max_iter,tol,reg,
                            threads_per_fit,&cancelled);
      } catch(const FitCancelled&) {
        outcome=JobOutcome::Cancelled;
        // Interruption is reported once by the main R thread after every
        // worker has unwound and joined. A peer job can also request the same
        // cooperative cancellation after its first fatal native failure.
      } catch(const std::exception& e) {
        outcome=JobOutcome::Failed;
        errors[i]=e.what();
        cancelled.store(true,std::memory_order_relaxed);
      } catch(...) {
        outcome=JobOutcome::Failed;
        errors[i]="unknown native exception";
        cancelled.store(true,std::memory_order_relaxed);
      }
      const double fit_seconds=std::chrono::duration<double>(
        std::chrono::steady_clock::now()-fit_started).count();
      fit_elapsed_seconds[i]=fit_seconds;
      {
        std::lock_guard<std::mutex> lock(progress_mutex);
        completion_events[published_events++]=CompletionEvent{i,fit_seconds,outcome};
      }
      progress_condition.notify_one();
    }
    workers_finished.fetch_add(1);
    progress_condition.notify_one();
  };
  std::vector<std::thread> pool;
  const int pool_size=std::min(workers,jobs);
  pool.reserve(pool_size);
  ThreadJoiner joiner{pool};
  for(int w=0;w<pool_size;++w) pool.emplace_back(work);

  // Poll even when progress is disabled. Calling R_CheckUserInterrupt() from a
  // worker is unsafe, while blocking the R main thread in join() makes a long
  // native round effectively uninterruptible. R_ToplevelExec() contains the
  // interrupt longjmp; cancellation then follows normal C++ unwinding.
  int reported=0;
  int successful_reported=0;
  bool interrupted=false;
  double completed_fit_seconds=0.0;
  for(;;) {
    std::unique_lock<std::mutex> lock(progress_mutex);
    progress_condition.wait_for(lock,std::chrono::milliseconds(100),[&](){
      return published_events>reported || workers_finished.load()>=pool_size;
    });
    while(reported<published_events) {
      CompletionEvent event=completion_events[reported];
      ++reported;
      if(event.outcome==JobOutcome::Success) {
        ++successful_reported;
        completed_fit_seconds+=event.elapsed_seconds;
      }
      if(progress && event.outcome==JobOutcome::Success) {
        lock.unlock();
        const double round_elapsed=std::chrono::duration<double>(
          std::chrono::steady_clock::now()-native_started).count();
        const double mean_fit_seconds=completed_fit_seconds/successful_reported;
        const int remaining_jobs=jobs-successful_reported;
        // Estimate complete worker-sized batches rather than dividing a tail
        // job by every worker. K-dependent runtimes still make this a rough
        // planning estimate, which is stated explicitly in the output.
        const int remaining_batches=(remaining_jobs+pool_size-1)/pool_size;
        const double estimated_remaining=mean_fit_seconds*remaining_batches;
        const std::string fit_duration=format_progress_duration(event.elapsed_seconds);
        const std::string elapsed_duration=format_progress_duration(round_elapsed);
        const std::string eta_duration=format_progress_duration(estimated_remaining);
        Rprintf(
          "CellCharterR C++ GMM: job %d/%d complete; K=%d; fit=%s; "
          "round elapsed=%s; rough round ETA=%s.\n",
          successful_reported,jobs,native_ks[event.index],fit_duration.c_str(),
          elapsed_duration.c_str(),eta_duration.c_str());
        R_FlushConsole();
        lock.lock();
      }
    }
    const bool all_workers_finished=workers_finished.load()>=pool_size;
    lock.unlock();

    if(!interrupted && r_interrupt_pending()) {
      interrupted=true;
      cancelled.store(true,std::memory_order_relaxed);
      progress_condition.notify_all();
    }
    if(all_workers_finished) break;
  }
  for(auto& thread:pool) if(thread.joinable()) thread.join();
  if(interrupted) Rcpp::stop("Native C++ GMM fitting interrupted by user.");
  for(int i=0;i<jobs;++i) if(!errors[i].empty())
    Rcpp::stop(
      "Native GMM job %d (K=%d) failed; remaining jobs were cancelled: %s",
      i+1,native_ks[i],errors[i]);
  Rcpp::List out(jobs);
  for(int i=0;i<jobs;++i)
    out[i]=export_fit(fits[i],reg,native_seeds[i],fit_elapsed_seconds[i]);
  return out;
  END_RCPP
}
