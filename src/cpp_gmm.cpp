#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <future>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {
struct Data { int n, p; std::vector<double> x; };
struct Fit {
  int k, p, iterations; bool converged;
  double reg, nll;
  double tolerance;
  std::vector<double> means, covs, weights, history;
  std::vector<int> labels;
  std::vector<double> initialization_counts;
  std::vector<int> initialization_collapsed;
  int early_wait = 0;
  double early_best = std::numeric_limits<double>::infinity();
};

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

double component_log_density(const Data& d, int row, const Fit& f, int c,
                             const std::vector<double>& lower,
                             std::vector<double>& solved) {
  double quad = 0.0, logdet = 0.0;
  for (int j = 0; j < f.p; ++j) {
    double v = d.x[static_cast<size_t>(row) * f.p + j] -
      f.means[static_cast<size_t>(c) * f.p + j];
    for (int q = 0; q < j; ++q) v -= atc(lower, f.p, j, q) * solved[q];
    solved[j] = v / atc(lower, f.p, j, j);
    quad += solved[j] * solved[j];
    logdet += std::log(atc(lower, f.p, j, j));
  }
  return std::log(f.weights[c]) - 0.5 *
    (f.p * std::log(2.0 * 3.14159265358979323846) + 2.0 * logdet + quad);
}

std::vector<std::vector<double> > factors(const Fit& f) {
  std::vector<std::vector<double> > out(f.k);
  for (int c = 0; c < f.k; ++c)
    if (!cholesky(&f.covs[static_cast<size_t>(c) * f.p * f.p], f.p, f.reg, out[c]))
      throw std::runtime_error("Cholesky factorization failed");
  return out;
}

double responsibilities(const Data& d, int row, const Fit& f,
                        const std::vector<std::vector<double> >& chol,
                        std::vector<double>& values) {
  double maximum = -std::numeric_limits<double>::infinity();
  thread_local std::vector<double> solved;
  solved.resize(f.p);
  for (int c = 0; c < f.k; ++c) {
    values[c] = component_log_density(d, row, f, c, chol[c], solved);
    maximum = std::max(maximum, values[c]);
  }
  double total = 0.0;
  for (int c = 0; c < f.k; ++c) total += std::exp(values[c] - maximum);
  const double log_norm = maximum + std::log(total);
  for (int c = 0; c < f.k; ++c) values[c] = std::exp(values[c] - log_norm);
  return log_norm;
}

std::vector<double> kmeans_pp(const Data& d, int k, uint64_t seed) {
  std::mt19937_64 rng(seed); // private RNG: does not consume or mutate R's RNG.
  std::uniform_int_distribution<int> first(0, d.n - 1);
  std::vector<double> means(static_cast<size_t>(k) * d.p);
  std::copy_n(&d.x[static_cast<size_t>(first(rng)) * d.p], d.p, means.begin());
  std::vector<double> distance(d.n, std::numeric_limits<double>::infinity());
  for (int c = 1; c < k; ++c) {
    double sum = 0.0;
    for (int i = 0; i < d.n; ++i) {
      double current = 0.0;
      for (int j = 0; j < d.p; ++j) {
        double z = d.x[static_cast<size_t>(i) * d.p + j] - means[static_cast<size_t>(c - 1) * d.p + j];
        current += z * z;
      }
      distance[i] = std::min(distance[i], current); sum += distance[i];
    }
    int selected = c % d.n;
    if (sum > 0.0 && std::isfinite(sum)) {
      std::uniform_real_distribution<double> draw(0.0, sum);
      double target = draw(rng), cumulative = 0.0;
      for (int i = 0; i < d.n; ++i) if ((cumulative += distance[i]) >= target) { selected = i; break; }
    }
    std::copy_n(&d.x[static_cast<size_t>(selected) * d.p], d.p,
                means.begin() + static_cast<size_t>(c) * d.p);
  }
  // Complete the initializer with deterministic Lloyd passes. The final hard
  // assignments are summarized again below rather than treated as soft EM.
  std::vector<int> assignment(d.n, 0);
  for (int iteration = 0; iteration < 30; ++iteration) {
    std::vector<double> sums(static_cast<size_t>(k) * d.p, 0.0), counts(k, 0.0);
    bool changed = false;
    for (int i = 0; i < d.n; ++i) {
      int best = 0; double best_distance = std::numeric_limits<double>::infinity();
      for (int c = 0; c < k; ++c) {
        double value = 0.0; for (int j = 0; j < d.p; ++j) {
          double z = d.x[static_cast<size_t>(i)*d.p+j] - means[static_cast<size_t>(c)*d.p+j]; value += z*z;
        }
        if (value < best_distance) { best_distance = value; best = c; }
      }
      changed = changed || iteration == 0 || assignment[i] != best; assignment[i] = best; counts[best] += 1.0;
      for (int j = 0; j < d.p; ++j) sums[static_cast<size_t>(best)*d.p+j] += d.x[static_cast<size_t>(i)*d.p+j];
    }
    for (int c = 0; c < k; ++c) if (counts[c] > 0.0) for (int j = 0; j < d.p; ++j)
      means[static_cast<size_t>(c)*d.p+j] = sums[static_cast<size_t>(c)*d.p+j] / counts[c];
    if (!changed) break;
  }
  return means;
}

Fit fit_once(const Data& d, int k, uint64_t seed, int max_iter,
             double tolerance, double reg, int threads_per_fit) {
  Fit f; f.k = k; f.p = d.p; f.reg = reg; f.iterations = 0; f.converged = false;
  f.tolerance = tolerance;
  const double collapse_threshold = 1e-6;
  std::vector<double> centers = kmeans_pp(d, k, seed);
  f.means.assign(static_cast<size_t>(k) * d.p, 0.0); f.weights.assign(k, 0.0);
  f.covs.assign(static_cast<size_t>(k) * d.p * d.p, 0.0);
  std::vector<int> hard(d.n, 0); f.initialization_counts.assign(k, 0.0);
  for (int i = 0; i < d.n; ++i) {
    double best_distance = std::numeric_limits<double>::infinity(); int best = 0;
    for (int c = 0; c < k; ++c) { double value=0.0; for(int j=0;j<d.p;++j) {
      double z=d.x[static_cast<size_t>(i)*d.p+j]-centers[static_cast<size_t>(c)*d.p+j]; value+=z*z; }
      if(value<best_distance){best_distance=value;best=c;} }
    hard[i]=best; f.initialization_counts[best] += 1.0;
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
  // Center before multiplying: this is stable for large coordinate offsets.
  for(int i=0;i<d.n;++i) { int c=hard[i]; if(f.initialization_counts[c]<collapse_threshold) continue;
    for(int a=0;a<d.p;++a) for(int b=0;b<=a;++b) f.covs[(static_cast<size_t>(c)*d.p+a)*d.p+b] +=
      (d.x[static_cast<size_t>(i)*d.p+a]-f.means[static_cast<size_t>(c)*d.p+a]) *
      (d.x[static_cast<size_t>(i)*d.p+b]-f.means[static_cast<size_t>(c)*d.p+b]); }
  for(int c=0;c<k;++c) for(int a=0;a<d.p;++a) for(int b=0;b<=a;++b) {
    double value = f.initialization_counts[c] < collapse_threshold ? (a==b ? 1.0 : 0.0) :
      f.covs[(static_cast<size_t>(c)*d.p+a)*d.p+b]/f.initialization_counts[c];
    f.covs[(static_cast<size_t>(c)*d.p+a)*d.p+b]=value;
    f.covs[(static_cast<size_t>(c)*d.p+b)*d.p+a]=value;
  }
  double best = std::numeric_limits<double>::infinity(); int wait = 0;
  std::vector<double> r(k), nk(k), sums(static_cast<size_t>(k) * d.p);
  const int nt = std::max(1, std::min(threads_per_fit, d.n));
  for (int iteration = 0; iteration < max_iter; ++iteration) {
    const std::vector<std::vector<double> > chol = factors(f);
    std::fill(nk.begin(), nk.end(), 0.0); std::fill(sums.begin(), sums.end(), 0.0);
    double ll = 0.0;
    // Fixed contiguous shards plus shard-order reduction are deterministic for
    // a fixed topology. Worker threads never call the R API.
    std::vector<std::vector<double> > local_nk(nt,std::vector<double>(k,0.0));
    std::vector<std::vector<double> > local_sums(nt,std::vector<double>(sums.size(),0.0));
    std::vector<double> local_ll(nt,0.0); std::vector<std::thread> pass1;
    for(int t=0;t<nt;++t) pass1.emplace_back([&,t](){
      int begin=static_cast<int>((static_cast<int64_t>(d.n)*t)/nt), end=static_cast<int>((static_cast<int64_t>(d.n)*(t+1))/nt);
      std::vector<double> lr(k);
      for(int i=begin;i<end;++i){ local_ll[t]+=responsibilities(d,i,f,chol,lr);
        for(int c=0;c<k;++c){ local_nk[t][c]+=lr[c]; for(int j=0;j<d.p;++j) local_sums[t][static_cast<size_t>(c)*d.p+j]+=lr[c]*d.x[static_cast<size_t>(i)*d.p+j]; }}
    });
    for(auto& thread:pass1) thread.join();
    for(int t=0;t<nt;++t){ ll+=local_ll[t]; for(int c=0;c<k;++c) nk[c]+=local_nk[t][c]; for(size_t q=0;q<sums.size();++q) sums[q]+=local_sums[t][q]; }
    std::vector<double> new_means(sums.size()), new_weights(k);
    std::vector<bool> collapsed(k, false);
    for (int c = 0; c < k; ++c) {
      if (nk[c] < 1e-6) {
        collapsed[c] = true;
        int selected = static_cast<int>((seed + iteration + c) % d.n);
        std::copy_n(&d.x[static_cast<size_t>(selected) * d.p], d.p,
                    new_means.begin() + static_cast<size_t>(c) * d.p);
        nk[c] = 1e-6;
      } else for (int j = 0; j < d.p; ++j)
        new_means[static_cast<size_t>(c) * d.p + j] = sums[static_cast<size_t>(c) * d.p + j] / nk[c];
      new_weights[c] = std::max(nk[c] / d.n, 1e-15);
    }
    double weight_sum = std::accumulate(new_weights.begin(), new_weights.end(), 0.0);
    for (double& w : new_weights) w /= weight_sum;
    std::vector<double> new_covs(f.covs.size(), 0.0);
    // Pass two recomputes the same old-parameter responsibilities and forms
    // covariance from (x-new_mean)(x-new_mean)', avoiding E[xx']-mu mu'.
    std::vector<std::vector<double> > local_covs(nt,std::vector<double>(new_covs.size(),0.0)); std::vector<std::thread> pass2;
    for(int t=0;t<nt;++t) pass2.emplace_back([&,t](){
      int begin=static_cast<int>((static_cast<int64_t>(d.n)*t)/nt), end=static_cast<int>((static_cast<int64_t>(d.n)*(t+1))/nt); std::vector<double> lr(k);
      for(int i=begin;i<end;++i){ responsibilities(d,i,f,chol,lr); for(int c=0;c<k;++c) for(int a=0;a<d.p;++a){
        double za=d.x[static_cast<size_t>(i)*d.p+a]-new_means[static_cast<size_t>(c)*d.p+a];
        for(int b=0;b<=a;++b) local_covs[t][(static_cast<size_t>(c)*d.p+a)*d.p+b]+=lr[c]*za*(d.x[static_cast<size_t>(i)*d.p+b]-new_means[static_cast<size_t>(c)*d.p+b]); }}
    });
    for(auto& thread:pass2) thread.join(); for(int t=0;t<nt;++t) for(size_t q=0;q<new_covs.size();++q) new_covs[q]+=local_covs[t][q];
    for (int c = 0; c < k; ++c) for (int a = 0; a < d.p; ++a) for (int b = 0; b <= a; ++b) {
      double value=collapsed[c] ? (a==b ? 1.0 : 0.0) :
        new_covs[(static_cast<size_t>(c)*d.p+a)*d.p+b]/nk[c];
      new_covs[(static_cast<size_t>(c)*d.p+a)*d.p+b]=value;
      new_covs[(static_cast<size_t>(c)*d.p+b)*d.p+a]=value;
    }
    f.means.swap(new_means); f.weights.swap(new_weights); f.covs.swap(new_covs);
    f.iterations = iteration + 1; f.history.push_back(-ll / d.n);
    if (f.history.back() < best - tolerance) { best = f.history.back(); wait = 0; }
    else if (++wait >= 2) { f.converged = true; break; }
  }
  const auto chol = factors(f); double ll = 0.0;
  f.labels.resize(d.n);
  std::vector<double> final_ll(nt,0.0); std::vector<std::thread> scoring;
  for(int t=0;t<nt;++t) scoring.emplace_back([&,t](){ int begin=static_cast<int>((static_cast<int64_t>(d.n)*t)/nt), end=static_cast<int>((static_cast<int64_t>(d.n)*(t+1))/nt); std::vector<double> lr(k);
    for(int i=begin;i<end;++i){ final_ll[t]+=responsibilities(d,i,f,chol,lr); f.labels[i]=static_cast<int>(std::max_element(lr.begin(),lr.end())-lr.begin()); }});
  for(auto& thread:scoring) thread.join(); for(int t=0;t<nt;++t) ll+=final_ll[t];
  f.nll = -ll / d.n;
  f.early_wait = wait; f.early_best = best;
  return f;
}

Fit fit_restart(const Data& d, int k, uint64_t seed, int max_iter,
                double tolerance, double requested_reg, int threads_per_fit) {
  std::string error;
  for (int attempt = 0; attempt < 8; ++attempt) try {
    return fit_once(d, k, seed, max_iter, tolerance, requested_reg * std::pow(10.0, attempt), threads_per_fit);
  } catch (const std::runtime_error& e) { error = e.what(); }
  throw std::runtime_error("GMM failed after 8 regularization attempts: " + error);
}

Data import_matrix(const Rcpp::NumericMatrix& x) {
  Data d; d.n = x.nrow(); d.p = x.ncol(); d.x.resize(static_cast<size_t>(d.n) * d.p);
  for (int i = 0; i < d.n; ++i) for (int j = 0; j < d.p; ++j) {
    double value = x(i, j); if (!std::isfinite(value)) Rcpp::stop("`x` contains non-finite values.");
    d.x[static_cast<size_t>(i) * d.p + j] = value;
  }
  return d;
}

Rcpp::List export_fit(const Fit& f, double requested_reg, uint64_t seed) {
  Rcpp::NumericMatrix means(f.k, f.p);
  for (int c = 0; c < f.k; ++c) for (int j = 0; j < f.p; ++j) means(c,j) = f.means[static_cast<size_t>(c)*f.p+j];
  Rcpp::List covs(f.k);
  for (int c = 0; c < f.k; ++c) { Rcpp::NumericMatrix m(f.p,f.p); for(int a=0;a<f.p;++a) for(int b=0;b<f.p;++b)
    m(a,b)=f.covs[(static_cast<size_t>(c)*f.p+a)*f.p+b]; covs[c]=m; }
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
    Rcpp::_["seed"]=static_cast<double>(seed));
}
}

extern "C" SEXP C_cpp_gmm_fit(SEXP xs, SEXP ks, SEXP seeds, SEXP maxs, SEXP tols, SEXP regs, SEXP threadss) {
  BEGIN_RCPP
  Data d=import_matrix(Rcpp::as<Rcpp::NumericMatrix>(xs)); int k=Rcpp::as<int>(ks);
  uint64_t seed=static_cast<uint64_t>(Rcpp::as<double>(seeds));
  return export_fit(fit_restart(d,k,seed,Rcpp::as<int>(maxs),Rcpp::as<double>(tols),Rcpp::as<double>(regs),Rcpp::as<int>(threadss)),Rcpp::as<double>(regs),seed);
  END_RCPP
}

extern "C" SEXP C_cpp_gmm_fit_many(SEXP xs, SEXP kss, SEXP seedss, SEXP maxs, SEXP tols, SEXP regs, SEXP workerss, SEXP threadss) {
  BEGIN_RCPP
  Data d=import_matrix(Rcpp::as<Rcpp::NumericMatrix>(xs)); Rcpp::IntegerVector ks(kss); Rcpp::NumericVector seeds(seedss);
  int jobs=ks.size(), workers=std::max(1,Rcpp::as<int>(workerss)), max_iter=Rcpp::as<int>(maxs);
  double tol=Rcpp::as<double>(tols), reg=Rcpp::as<double>(regs); int threads_per_fit=std::max(1,Rcpp::as<int>(threadss));
  std::vector<Fit> fits(jobs); std::vector<std::string> errors(jobs); std::atomic<int> next(0);
  auto work=[&]() { for (;;) { int i=next.fetch_add(1); if(i>=jobs) break; try {
    fits[i]=fit_restart(d,ks[i],static_cast<uint64_t>(seeds[i]),max_iter,tol,reg,threads_per_fit);
  } catch(const std::exception& e) { errors[i]=e.what(); } } };
  std::vector<std::thread> pool; for(int w=0;w<std::min(workers,jobs);++w) pool.emplace_back(work);
  for(auto& thread:pool) thread.join();
  for(int i=0;i<jobs;++i) if(!errors[i].empty()) Rcpp::stop("Native GMM job %d failed: %s",i+1,errors[i]);
  Rcpp::List out(jobs); for(int i=0;i<jobs;++i) out[i]=export_fit(fits[i],reg,static_cast<uint64_t>(seeds[i]));
  return out;
  END_RCPP
}
