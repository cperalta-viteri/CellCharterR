# Compiled parallel-fit tests exercise correctness and interrupt handling, not
# BLAS thread discovery. Give every test process the production-safe contract
# explicitly so R CMD check is independent of the surrounding CI shell. Tests
# of the rejection path override individual values with withr::local_envvar().
.cellcharterr_test_blas_environment <- c(
  OPENBLAS_NUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  GOTO_NUM_THREADS = "1"
)
do.call(Sys.setenv, as.list(.cellcharterr_test_blas_environment))
