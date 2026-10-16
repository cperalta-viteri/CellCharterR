.expect_cpp_sigint <- function(fit) {
  marker <- tempfile("cpp-gmm-started-")
  collected <- NULL
  job <- parallel::mcparallel(
    {
      # This handshake tells the parent that the forked child has completed
      # startup and is immediately about to enter the fitted call.
      writeLines("ready", marker)
      force(fit)
    },
    silent = TRUE
  )
  on.exit({
    unlink(marker)
    if (is.null(collected)) {
      try(tools::pskill(job$pid, tools::SIGKILL), silent = TRUE)
      parallel::mccollect(job, wait = FALSE)
    }
  }, add = TRUE)

  handshake_started <- proc.time()[["elapsed"]]
  while (!file.exists(marker) &&
         proc.time()[["elapsed"]] - handshake_started < 10) {
    Sys.sleep(0.01)
  }
  expect_true(file.exists(marker), info = paste(
    "Forked interrupt-test child did not complete its startup handshake."
  ))
  if (!file.exists(marker)) return(NULL)

  # Give the wrapper enough time to cross the native boundary. The expensive
  # fixture remains active substantially longer than this on supported CI.
  Sys.sleep(0.1)
  expect_true(tools::pskill(job$pid, tools::SIGINT))
  interrupted <- proc.time()[["elapsed"]]
  while (is.null(collected) && proc.time()[["elapsed"]] - interrupted < 10) {
    collected <- parallel::mccollect(job, wait = FALSE)
    if (is.null(collected)) Sys.sleep(0.05)
  }

  expect_false(is.null(collected), info = paste(
    "Native workers did not cooperatively stop and join within ten seconds."
  ))
  if (is.null(collected)) return(NULL)
  collected[[1L]]
}

.expect_cpp_interrupt_result <- function(result) {
  expect_s3_class(result, "try-error")
  expect_match(
    as.character(result),
    "Native C\\+\\+ GMM fitting interrupted by user",
    fixed = FALSE
  )
}

test_that("monitored single-fit routing preserves fitted results and metadata", {
  set.seed(1203)
  x <- matrix(stats::rnorm(720L), ncol = 6L)
  expected <- .Call(
    "C_cpp_gmm_fit",
    x,
    4L,
    103,
    12L,
    1e-7,
    1e-6,
    2L,
    PACKAGE = "CellCharterR"
  )
  observed <- cpp_gmm_fit(
    x,
    n_clusters = 4L,
    max_iter = 12L,
    tolerance = 1e-7,
    covariance_regularization = 1e-6,
    seed = 103L,
    num_threads = 2L
  )

  for (field in names(expected)) {
    expect_identical(observed[[field]], expected[[field]], info = field)
  }
  expect_identical(class(observed), c("cpp_gmm", "torch_gmm"))
  expect_identical(
    observed$runtime,
    list(
      device = "cpu",
      backend = "compiled-cpp-double",
      seed = 103L,
      cpu_threads = 2L
    )
  )
})

test_that("compiled multi-fit rounds cooperate with SIGINT when progress is off", {
  skip_on_os("windows")

  # Run the native round in a forked R child so the test runner can deliver a
  # real SIGINT without interrupting its own process. The child starts the same
  # outer and inner native worker topology used by AutoK.
  set.seed(1204)
  x <- matrix(stats::rnorm(150000L * 24L), ncol = 24L)
  result <- .expect_cpp_sigint(
    cpp_gmm_fit_many(
      x,
      n_clusters = 12:15,
      seeds = 101:104,
      workers = 2L,
      threads_per_fit = 2L,
      max_iter = 100L,
      progress = FALSE
    )
  )
  .expect_cpp_interrupt_result(result)
})

test_that("exported compiled single fits cooperate with SIGINT", {
  skip_on_os("windows")

  set.seed(1205)
  x <- matrix(stats::rnorm(150000L * 24L), ncol = 24L)
  result <- .expect_cpp_sigint(
    cpp_gmm_fit(
      x,
      n_clusters = 15L,
      max_iter = 100L,
      seed = 105L,
      num_threads = 2L
    )
  )
  .expect_cpp_interrupt_result(result)
})
