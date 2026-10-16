workflow_file <- function(...) {
  testthat::test_path("..", "..", ...)
}

collapse_source <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("CODEX workflow keeps primary, standalone K11, and published roles distinct", {
  workflow_path <- workflow_file("scripts", "run_codex_mouse_spleen.R")
  skip_if_not(
    file.exists(workflow_path),
    "repository-only workflow sources are absent from the installed test tree"
  )
  expect_silent(parse(workflow_path))
  workflow <- collapse_source(workflow_path)

  expect_match(workflow, "CellCharterR_mouse_spleen_v2", fixed = TRUE)
  expect_match(workflow, "frozen_result_roots", fixed = TRUE)
  expect_match(workflow, '"tutorial_k11"', fixed = TRUE)
  expect_match(workflow, 'role_id = "r_autok_primary"', fixed = TRUE)
  expect_match(
    workflow, 'role_id = "r_tutorial_k11_sensitivity"', fixed = TRUE
  )
  expect_match(
    workflow, 'role_id = "published_h5ad_reference"', fixed = TRUE
  )
  expect_match(workflow, "tutorial_k11_seed <- 12345L", fixed = TRUE)
  expect_match(workflow, "tutorial_k11_model <- if (identical(gmm_backend, \"cpp\"))", fixed = TRUE)
  expect_match(workflow, "cpp_gmm_fit(", fixed = TRUE)
  expect_match(workflow, "torch_gmm_fit(", fixed = TRUE)
  expect_match(workflow, "n_clusters = 11L", fixed = TRUE)
  expect_match(workflow, "seed = tutorial_k11_seed", fixed = TRUE)
  expect_match(workflow, "tutorial_k11_model_provenance_path", fixed = TRUE)
  expect_match(workflow, "cellcharter_autok_stability_data(", fixed = TRUE)
  expect_match(workflow, "cellcharter_plot_autok_stability(", fixed = TRUE)
  expect_match(workflow, 'k_range = 1:20', fixed = TRUE)
  expect_match(workflow, 'highlight_k = c(4L, 11L)', fixed = TRUE)
  expect_match(workflow, 'stability_k1_k20.rds', fixed = TRUE)
  expect_match(workflow, 'stability_k1_k20.png', fixed = TRUE)
  expect_match(
    workflow, 'reference_call = "Cluster(n_clusters=11, random_state=12345)"',
    fixed = TRUE
  )

  # AutoK still selects the primary K. It must not export its K=11 candidate as
  # the tutorial sensitivity checkpoint.
  expect_match(workflow, "k = autok$best_k", fixed = TRUE)
  expect_false(grepl("k11_row <- autok$best_index", workflow, fixed = TRUE))
  expect_false(grepl("k11_artifact <- readRDS", workflow, fixed = TRUE))
  expect_match(
    workflow, 'sys.source(\n      "scripts/build_second_review_deliverables.R"',
    fixed = TRUE
  )
  expect_match(
    workflow, '"report", "CellCharterR_mouse_spleen_AutoK_K11.html"',
    fixed = TRUE
  )
  expect_match(workflow, "output_file = basename(report_html_path)", fixed = TRUE)

  for (field in c(
    "role_id", "role", "seed", "model_path", "cell_ids",
    "cell_order_sha256", "aggregation_data_id"
  )) {
    expect_match(workflow, field, fixed = TRUE)
  }

  # The ordered-cell identity must exist before the AutoK contract consumes it.
  assignment <- regexpr(
    "cell_order_sha256 <- .cc_observation_order_sha256(colnames(sce))",
    workflow, fixed = TRUE
  )[[1L]]
  autok_contract <- regexpr(
    "autok_contract <- cellcharter_stage_contract(", workflow, fixed = TRUE
  )[[1L]]
  expect_gt(assignment, 0L)
  expect_gt(autok_contract, assignment)
})

test_that("report checkpoint lifecycle creates, reuses, and rejects stale artifacts", {
  temporary <- tempfile("cellcharterr-report-contract-")
  dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  artifacts <- c(
    html = file.path(temporary, "CellCharterR_mouse_spleen_AutoK_K11.html"),
    report_inputs = file.path(temporary, "CellCharterR_primary_K11_manifest.rds")
  )
  writeLines("self-contained report", artifacts[["html"]])
  saveRDS(list(primary_k = 4L, tutorial_k = 11L), artifacts[["report_inputs"]])
  contract <- cellcharter_stage_contract(
    stage = "report",
    implementation_digest = paste(rep("a", 64L), collapse = ""),
    upstream = list(downstream_contract_id = "fixture-downstream"),
    cell_order_sha256 = paste(rep("b", 64L), collapse = ""),
    parameters = list(output_file = basename(artifacts[["html"]])),
    seeds = list(downstream_seed = 0L),
    versions = list(CellCharterR = "test"),
    runtime = list(renderer = "fixture")
  )
  record_path <- file.path(temporary, "report.rds")
  cellcharter_record_checkpoint(
    contract, artifacts, status = "created_under_contract", path = record_path
  )

  expect_true(cellcharter_checkpoint_reusable(record_path, contract, artifacts))
  writeLines("mutated report", artifacts[["html"]])
  expect_false(cellcharter_checkpoint_reusable(record_path, contract, artifacts))
})

test_that("workflow contracts include transitive implementation and reject mutation", {
  workflow_path <- workflow_file("scripts", "run_codex_mouse_spleen.R")
  skip_if_not(
    file.exists(workflow_path),
    "repository-only workflow sources are absent from the installed test tree"
  )
  workflow <- collapse_source(workflow_path)

  for (material_helper in c(
    "verify_legacy_representation = .representation_checkpoint_matches_inputs",
    "read_h5ad_csr_block = .read_h5ad_csr_block",
    "aggregate_rows = .cc_aggregate_rows_block",
    "centered_covariances = .rc_centered_covariances",
    "autok_mirror_stability = .autok_mirror_stability",
    "neighborhood_from_library_stats = .cc_nhood_from_library_stats"
  )) {
    expect_match(workflow, material_helper, fixed = TRUE)
  }
  expect_match(
    workflow,
    "A content-addressed artifact is immutable",
    fixed = TRUE
  )
  expect_false(grepl("data_id_migration", workflow, fixed = TRUE))
  expect_match(workflow, 'reference_pytorch = "1.12.1"', fixed = TRUE)
  expect_match(
    workflow,
    'reference_pycave_commit = "6df56fa15797940a2063f3e70b3a7641f592a9f6"',
    fixed = TRUE
  )
})

test_that("revised deliverables enforce tutorial display and inference semantics", {
  builder_path <- workflow_file("scripts", "build_second_review_deliverables.R")
  skip_if_not(
    file.exists(builder_path),
    "repository-only workflow sources are absent from the installed test tree"
  )
  expect_silent(parse(builder_path))
  builder <- collapse_source(builder_path)
  report <- collapse_source(workflow_file(
    "reports", "CellCharterR_mouse_spleen_AutoK_K11.Rmd"
  ))

  expect_match(builder, "results/CellCharterR_mouse_spleen_v2", fixed = TRUE)
  expect_match(builder, "phase0_dir", fixed = TRUE)
  expect_match(builder, "stability_k1_k20.png", fixed = TRUE)
  expect_match(builder, 'primary, "r_autok_primary"', fixed = TRUE)
  expect_match(
    builder, '"r_tutorial_k11_sensitivity"', fixed = TRUE
  )
  expect_match(builder, '"published_h5ad_reference"', fixed = TRUE)
  expect_match(
    builder,
    "analysis_roles <- list(Primary = primary, TutorialK11 = tutorial_k11)",
    fixed = TRUE
  )
  expect_false(grepl("identical(primary$k, 4L)", builder, fixed = TRUE))
  expect_match(builder, "report_input$cell_type", fixed = TRUE)
  expect_match(builder, "computed_order_sha256", fixed = TRUE)
  expect_match(builder, 'group_key = "cell_type"', fixed = TRUE)
  expect_match(builder, "levels = levels(report_input$cell_type)", fixed = TRUE)
  expect_false(grepl(
    "levels = sort(unique(as.character(report_input$cell_type)))",
    builder, fixed = TRUE
  ))
  expect_match(builder, "size_threshold = 0.58", fixed = TRUE)
  expect_match(builder, "color_threshold = 0.58", fixed = TRUE)
  expect_match(builder, "cluster_labels = TRUE", fixed = TRUE)
  expect_match(builder, "neighborhood_limits = c(-1, 1)", fixed = TRUE)
  expect_match(
    builder,
    'differential_limits = "data range (CellCharter v0.2.0 default)"',
    fixed = TRUE
  )
  expect_match(builder, 'displayed_samples = c("BALBc-1", "MRL-4")', fixed = TRUE)
  expect_match(builder, "significance = 0.05", fixed = TRUE)
  expect_match(
    builder,
    'reverse_display = "effect-only sign reversal; counts and p-values reused"',
    fixed = TRUE
  )
  inference_calls <- gregexpr(
    "\\.cc_diff_nhood_enrichment\\(", builder, perl = TRUE
  )[[1L]]
  expect_equal(sum(inference_calls > 0L), 1L)

  expect_match(report, "primary_k <- manifest$primary_k", fixed = TRUE)
  expect_match(report, "knitr::include_graphics(manifest$autok_stability_plot)", fixed = TRUE)
  expect_match(report, "full `cell_type`", fixed = TRUE)
  expect_match(report, "stats::hclust(", fixed = TRUE)
  expect_match(report, 'method = "complete"', fixed = TRUE)
  expect_match(report, "size_threshold <- 0.58", fixed = TRUE)
  expect_match(report, "color_threshold <- 0.58", fixed = TRUE)
  expect_match(report, "cex = 2.5 * cells$size^0.75", fixed = TRUE)
  expect_false(grepl("sqrt(cells$size)", report, fixed = TRUE))
  expect_match(report, "significance = 0.05", fixed = TRUE)
  expect_match(report, "plot_heatmap <- function(x, title, limits = NULL", fixed = TRUE)
  expect_match(report, "limits <- range(plotted, finite = TRUE)", fixed = TRUE)
  expect_match(report, "MRL-4 (early)", fixed = TRUE)
})

test_that("local pinned notebook evidence matches the workflow contract", {
  notebook <- workflow_file(
    "results", "reference", "cellcharter-v0.2.0", "docs", "notebooks",
    "codex_mouse_spleen.ipynb"
  )
  skip_if_not(file.exists(notebook))
  text <- collapse_source(notebook)

  expect_match(text, "n_clusters=11", fixed = TRUE)
  expect_match(text, "random_state=12345", fixed = TRUE)
  expect_match(text, "group_key='cell_type'", fixed = TRUE)
  expect_match(text, "color_threshold=0.58", fixed = TRUE)
  expect_match(text, "size_threshold=0.58", fixed = TRUE)
  expect_match(text, "n_perms=100", fixed = TRUE)
  expect_match(text, "significance=0.05", fixed = TRUE)
})
