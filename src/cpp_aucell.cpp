#include <Rcpp.h>
#include <R_ext/Utils.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

class AUCellCancelled : public std::exception {
 public:
  const char* what() const noexcept override { return "native AUCell scoring cancelled"; }
};

void aucell_check_r_interrupt(void*) {
  R_CheckUserInterrupt();
}

bool aucell_interrupt_pending() {
  return R_ToplevelExec(aucell_check_r_interrupt, nullptr) == FALSE;
}

uint64_t splitmix64(uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

uint64_t tie_key(uint64_t seed, uint64_t spot, uint64_t feature) {
  uint64_t value = seed;
  value ^= splitmix64(spot + UINT64_C(0x243f6a8885a308d3));
  value ^= splitmix64(feature + UINT64_C(0x13198a2e03707344));
  return splitmix64(value);
}

struct RankItem {
  double value;
  uint64_t tie;
  int feature;
};

struct RankBefore {
  bool operator()(const RankItem& left, const RankItem& right) const {
    if (left.value != right.value) return left.value > right.value;
    if (left.tie != right.tie) return left.tie < right.tie;
    return left.feature < right.feature;
  }
};

struct AUCellThreadJoiner {
  std::vector<std::thread>& threads;
  ~AUCellThreadJoiner() {
    for (std::thread& thread : threads) {
      if (thread.joinable()) thread.join();
    }
  }
};

}  // namespace

extern "C" SEXP C_cpp_aucell_score(
    SEXP xs,
    SEXP gene_setss,
    SEXP cutoffs,
    SEXP seeds,
    SEXP spot_positionss,
    SEXP threadss,
    SEXP tie_methods) {
  BEGIN_RCPP
  Rcpp::NumericMatrix x(xs);
  Rcpp::List gene_sets(gene_setss);
  Rcpp::IntegerVector spot_positions(spot_positionss);
  const int n_features = x.nrow();
  const int n_spots = x.ncol();
  const int n_programs = gene_sets.size();
  const int cutoff = Rcpp::as<int>(cutoffs);
  const uint64_t seed = static_cast<uint64_t>(Rcpp::as<double>(seeds));
  const int requested_threads = std::max(1, Rcpp::as<int>(threadss));
  const std::string tie_method = Rcpp::as<std::string>(tie_methods);
  const bool random_ties = tie_method == "random";

  if (n_features < 1 || n_spots < 1 || n_programs < 1) {
    Rcpp::stop("Native AUCell input must have features, spots, and programs.");
  }
  if (cutoff < 2 || cutoff > n_features) {
    Rcpp::stop("Native AUCell cutoff is outside the valid feature range.");
  }
  if (spot_positions.size() != n_spots) {
    Rcpp::stop("Native AUCell spot positions do not match the input columns.");
  }
  if (!random_ties && tie_method != "ordinal") {
    Rcpp::stop("Unknown native AUCell tie method.");
  }

  std::vector<uint64_t> native_spot_positions(n_spots);
  for (int spot = 0; spot < n_spots; ++spot) {
    if (spot_positions[spot] == NA_INTEGER || spot_positions[spot] < 1) {
      Rcpp::stop("Native AUCell spot positions must be positive integers.");
    }
    native_spot_positions[spot] =
      static_cast<uint64_t>(spot_positions[spot]);
  }

  std::vector<std::vector<int> > programs_by_feature(n_features);
  std::vector<double> denominators(n_programs);
  for (int program = 0; program < n_programs; ++program) {
    Rcpp::IntegerVector members(gene_sets[program]);
    const int available = members.size();
    if (available < 1) Rcpp::stop("Native AUCell received an empty gene set.");
    const int used = std::min(available, cutoff - 1);
    denominators[program] =
      static_cast<double>(used) * cutoff -
      static_cast<double>(used) * (used + 1) / 2.0;
    for (int member : members) {
      const int feature = member - 1;
      if (feature < 0 || feature >= n_features) {
        Rcpp::stop("Native AUCell gene-set index is outside the assay.");
      }
      programs_by_feature[feature].push_back(program);
    }
  }

  Rcpp::NumericMatrix output(n_spots, n_programs);
  const double* input = REAL(x);
  double* result = REAL(output);
  const int retained = cutoff - 1;
  const int pool_size = std::min(requested_threads, n_spots);
  std::atomic<int> next_spot(0);
  std::atomic<int> workers_finished(0);
  std::atomic<bool> cancelled(false);
  std::mutex completion_mutex;
  std::condition_variable completion_condition;
  std::vector<std::string> errors(pool_size);

  auto work = [&](int worker) {
    try {
      std::vector<RankItem> ranking(n_features);
      std::vector<double> numerators(n_programs);
      RankBefore before;
      for (;;) {
        if (cancelled.load(std::memory_order_relaxed)) throw AUCellCancelled();
        const int spot = next_spot.fetch_add(1, std::memory_order_relaxed);
        if (spot >= n_spots) break;
        const uint64_t global_spot = native_spot_positions[spot];
        const double* values = input + static_cast<size_t>(spot) * n_features;
        for (int feature = 0; feature < n_features; ++feature) {
          if ((feature & 2047) == 0 &&
              cancelled.load(std::memory_order_relaxed)) {
            throw AUCellCancelled();
          }
          ranking[feature] = RankItem{
            values[feature],
            random_ties ? tie_key(seed, global_spot,
                                  static_cast<uint64_t>(feature + 1))
                        : static_cast<uint64_t>(feature),
            feature
          };
        }
        std::nth_element(ranking.begin(), ranking.begin() + retained,
                         ranking.end(), before);
        std::sort(ranking.begin(), ranking.begin() + retained, before);
        std::fill(numerators.begin(), numerators.end(), 0.0);
        for (int position = 0; position < retained; ++position) {
          const int feature = ranking[position].feature;
          const double contribution = cutoff - (position + 1);
          for (int program : programs_by_feature[feature]) {
            numerators[program] += contribution;
          }
        }
        for (int program = 0; program < n_programs; ++program) {
          result[spot + static_cast<size_t>(n_spots) * program] =
            numerators[program] / denominators[program];
        }
      }
    } catch (const AUCellCancelled&) {
      // The main R thread reports interruption after all workers join.
    } catch (const std::exception& error) {
      errors[worker] = error.what();
      cancelled.store(true, std::memory_order_relaxed);
    } catch (...) {
      errors[worker] = "unknown native exception";
      cancelled.store(true, std::memory_order_relaxed);
    }
    workers_finished.fetch_add(1, std::memory_order_relaxed);
    completion_condition.notify_one();
  };

  std::vector<std::thread> pool;
  pool.reserve(pool_size);
  AUCellThreadJoiner joiner{pool};
  try {
    for (int worker = 0; worker < pool_size; ++worker) {
      pool.emplace_back(work, worker);
    }
  } catch (...) {
    cancelled.store(true, std::memory_order_relaxed);
    completion_condition.notify_all();
    throw;
  }

  bool interrupted = false;
  while (workers_finished.load(std::memory_order_relaxed) < pool_size) {
    std::unique_lock<std::mutex> lock(completion_mutex);
    completion_condition.wait_for(lock, std::chrono::milliseconds(100));
    lock.unlock();
    if (!interrupted && aucell_interrupt_pending()) {
      interrupted = true;
      cancelled.store(true, std::memory_order_relaxed);
    }
  }
  for (std::thread& thread : pool) {
    if (thread.joinable()) thread.join();
  }
  if (interrupted) Rcpp::stop("Native C++ AUCell scoring interrupted by user.");
  for (const std::string& error : errors) {
    if (!error.empty()) Rcpp::stop("Native AUCell worker failed: %s", error);
  }
  return output;
  END_RCPP
}
