test_that("stage contracts invalidate only material identity changes", {
  implementation <- cellcharter_implementation_digest(
    functions = list(transform = function(x) x + 1)
  )
  make_contract <- function(
      upstream = "a", order = paste(rep("0", 64L), collapse = ""), block_rows = 8L) {
    cellcharter_stage_contract(
      stage = "fixture",
      implementation_digest = implementation,
      upstream = list(data_sha256 = upstream),
      cell_order_sha256 = order,
      parameters = list(block_rows = block_rows, layers = 0:3),
      seeds = list(seed = 7L),
      versions = list(CellCharterR = "0.2.0.9000", reference = "fixture-v1"),
      runtime = list(device = "cpu", num_threads = 1L)
    )
  }

  original <- make_contract()
  expect_identical(make_contract()$contract_id, original$contract_id)
  expect_false(identical(make_contract(upstream = "b")$contract_id, original$contract_id))
  expect_false(identical(
    make_contract(order = paste0("1", substring(original$cell_order_sha256, 2L)))$contract_id,
    original$contract_id
  ))
  expect_false(identical(make_contract(block_rows = 9L)$contract_id, original$contract_id))
})

test_that("implementation regions ignore unrelated source changes", {
  source_path <- tempfile(fileext = ".R")
  writeLines(c(
    "unrelated <- 1L",
    "# CONTRACT BEGIN",
    "material <- function(x) x + 1L",
    "# CONTRACT END",
    "also_unrelated <- 2L"
  ), source_path)
  region <- list(stage = list(
    path = source_path, begin = "# CONTRACT BEGIN", end = "# CONTRACT END"
  ))
  first <- cellcharter_implementation_digest(source_regions = region)

  lines <- readLines(source_path)
  lines[[1L]] <- "unrelated <- 999L"
  writeLines(lines, source_path)
  expect_identical(cellcharter_implementation_digest(source_regions = region), first)

  lines[[3L]] <- "material <- function(x) x + 2L"
  writeLines(lines, source_path)
  expect_false(identical(cellcharter_implementation_digest(source_regions = region), first))
})

test_that("checkpoint reuse checks contracts and artifact content", {
  artifact <- tempfile(fileext = ".txt")
  writeLines("fixed artifact", artifact)
  record_path <- tempfile(fileext = ".rds")
  implementation <- cellcharter_implementation_digest(
    functions = list(identity = base::identity)
  )
  expected <- cellcharter_stage_contract(
    "fixture", implementation,
    upstream = list(input_sha256 = paste(rep("a", 64L), collapse = "")),
    cell_order_sha256 = paste(rep("b", 64L), collapse = ""),
    parameters = list(value = 1L), seeds = list(seed = NA_integer_),
    versions = list(reference = "v1"), runtime = list()
  )
  cellcharter_record_checkpoint(
    expected, c(output = artifact), status = "verified_existing_checkpoint",
    evidence = list(validation = "read-only fixture verification"), path = record_path
  )
  record <- readRDS(record_path)
  expect_identical(record$original_creation_contract, "unknown")
  expect_true(cellcharter_checkpoint_reusable(record_path, expected, c(output = artifact)))

  changed <- expected
  changed$parameters$value <- 2L
  changed$contract_id <- NULL
  changed$contract_id <- digest::digest(changed, algo = "sha256", serialize = TRUE)
  expect_false(cellcharter_checkpoint_reusable(record_path, changed, c(output = artifact)))

  writeLines("changed artifact", artifact)
  expect_false(cellcharter_checkpoint_reusable(record_path, expected, c(output = artifact)))
})

test_that("an allowlisted symbol-only digest migration preserves legacy artifacts", {
  artifact <- tempfile(fileext = ".txt")
  record_path <- tempfile(fileext = ".rds")
  on.exit(unlink(c(artifact, record_path)), add = TRUE)
  writeLines("immutable expensive result", artifact)

  legacy_digest <- cellcharter_implementation_digest(
    functions = list(old_name = function(sce) sce)
  )
  current_digest <- cellcharter_implementation_digest(
    functions = list(new_name = function(spe) spe)
  )
  make_contract <- function(implementation_digest) {
    cellcharter_stage_contract(
      stage = "fixture", implementation_digest = implementation_digest,
      upstream = list(input = "same"),
      cell_order_sha256 = paste(rep("a", 64L), collapse = ""),
      parameters = list(method = "same"), seeds = list(seed = NA_integer_),
      versions = list(CellCharterR = if (identical(
        implementation_digest, legacy_digest
      )) "0.2.0.9000" else "0.99.0", reference = "same"), runtime = list()
    )
  }
  legacy <- make_contract(legacy_digest)
  current <- make_contract(current_digest)
  cellcharter_record_checkpoint(
    legacy, c(output = artifact), path = record_path
  )
  artifact_mtime <- file.info(artifact)$mtime

  adopted <- CellCharterR:::.cc_adopt_legacy_implementation_contract(
    record_path, current, c(output = artifact), legacy_digest,
    allowed_version_values = list(CellCharterR = "0.2.0.9000")
  )
  expect_identical(adopted, legacy)
  expect_identical(file.info(artifact)$mtime, artifact_mtime)
  expect_null(CellCharterR:::.cc_adopt_legacy_implementation_contract(
    record_path, current, c(output = artifact), current_digest,
    allowed_version_values = list(CellCharterR = "0.2.0.9000")
  ))

  changed <- current
  changed$parameters$method <- "different"
  changed$contract_id <- NULL
  changed$contract_id <- digest::digest(
    CellCharterR:::.cc_contract_canonical(changed),
    algo = "sha256", serialize = TRUE
  )
  expect_null(CellCharterR:::.cc_adopt_legacy_implementation_contract(
    record_path, changed, c(output = artifact), legacy_digest,
    allowed_version_values = list(CellCharterR = "0.2.0.9000")
  ))
})
