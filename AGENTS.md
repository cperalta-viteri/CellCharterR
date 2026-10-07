# rCellCharter — project instructions

## Mission

Maintain and extend `rCellCharter`, an R-native, disk-backed implementation
benchmarked against the official CellCharter tutorials and version-pinned Python
reference implementation.

Core stages:

1. Read spatial data without materialising the full assay where avoidable.
1.2. Make sure the full workflow is orchestrated with SingleCellExperiment
2. Build sample-local Delaunay graphs and remove long links.
3. Aggregate neighbourhood features across graph hops.
4. Fit batched full-covariance Gaussian-mixture models with R `torch`.
5. Select the number of domains with an R implementation of CellCharter AutoK.
6. Reproduce the mouse-spleen tutorial: produce enrichment, neighbourhood-enrichment, and differential
   neighbourhood-enrichment analyses.

The package must remain usable for very large spatial datasets.

## Canonical references

Read these before changing behaviour or interpreting parity:

- CODEX mouse-spleen tutorial:
  https://cellcharter.readthedocs.io/en/latest/notebooks/codex_mouse_spleen.html

- CellCharter AutoK API:
  https://cellcharter.readthedocs.io/en/latest/generated/cellcharter.tl.ClusterAutoK.html

- CellCharter graph/proximity APIs:
  https://cellcharter.readthedocs.io/en/latest/api/graph.html

- Reference implementation:
  https://github.com/CSOgroup/cellcharter

For exact behavioural comparisons, pin the Python reference version and commit
used for the benchmark. Tutorial text is a guide; the corresponding versioned
Python source is authoritative for implementation details.

## Scientific correctness

- Do not claim exact end-to-end equivalence when the embedding differs.
- For parity with the CODEX tutorial, treat the supplied trVAE embedding as
  fixed input and reproduce the downstream graph, AutoK, GMM, and proximity
  stages in R.
- Do not fix K=11 for the primary R result. First run AutoK over the
  predeclared candidate range using repeated fits, Fowlkes–Mallows stability,
  `max_runs = 10`, and `convergence_tol = 0.01`.
- Use the R-selected K for primary maps, enrichment, and proximity analysis.
- K=11 is permitted only as a clearly labelled tutorial-matched sensitivity
  fit, never as the independent R result.
- For tutorial proximity parity, use the same CellCharter-style analysis:
  BALBc versus MRL, sample as library, `n_perms = 100`, and no methodological
  substitutions unless explicitly requested.

## Storage and scale

- Keep large assays and intermediate matrices disk-backed/HDF5-backed.
- Do not materialise all cells or the complete graph unless memory use has been
  measured and is safe.
- Process independent samples or batches separately when appropriate.
- Use batched tensor operations for the GMM.
- Save durable checkpoints for every expensive stage under `results/`:
  representation, graph, AutoK stability, fitted models, selected labels,
  enrichment matrices, permutation results, and rendered report inputs.
- Never assume an old `.rds` exists. If an input is absent, provide a script
  that regenerates it and document the expected path.

## GPU / R torch

- GPU `torch` is available inside this container, but verify it before any
  large run with `torch::cuda_is_available()`.
- Use `device = "auto"`: CUDA when available, CPU otherwise.
- Log the selected device, torch version, CUDA availability, GPU name when
  available, batch size, seed, and number of CPU threads.
- Keep CPU fallback working and tested.
- Do not add Python or scVI as hidden runtime dependencies for the R package.

## R conventions

- Use explicit namespaces, e.g. `torch::`, `Matrix::`, `DelayedArray::`.
- Write functions with roxygen2 documentation and deterministic seeds.
- Validate arguments early; give actionable errors.
- Keep package code in `R/`, tests in `tests/testthat/`, report sources in
  `reports/`, and executable workflow scripts in `scripts/`.
- Do not rewrite working modules unnecessarily. Make focused changes and
  preserve public function names and existing report content.

## Required tests before handoff

Run the smallest relevant set after each change, then the package checks before
a release:

1. GPU/CPU smoke test for `torch_gmm_fit()` and `torch_gmm_predict()`.
2. Synthetic full-covariance GMM test.
3. Graph and neighbour-aggregation unit tests.
4. AutoK stability/selection test on a small deterministic fixture.
5. Python–R aggregation parity test when the fixture is available.
6. `R CMD check --no-manual .`

For tutorial-scale runs, generate a visible HTML report with all code chunks,
parameters, results, figures, and limitations shown.

## Multi-agent workflow

- Use subagents for independent research, code review, testing, and
  implementation analysis when this can improve speed or quality.
- Keep each subagent's scope narrow and explicit.
- Avoid having multiple agents concurrently modify the same files.
- The primary agent must review and synthesize subagent findings before
  making repository-wide changes.

## Must requirements for a defensible R reproduction of CellCharter

1. Pin the reference implementation

Record the exact CellCharter version and Git commit used for every benchmark, together with the relevant Squidpy, TorchGMM or PyCave, PyTorch, and Python versions. Record all random seeds and known seed behavior. This choice must be explicit because the publication analysis used CellCharter 0.2.0 with PyCave, whereas newer CellCharter releases use TorchGMM and may not behave identically.

2. Use the official fixed trVAE representation for CODEX parity

Export the official mouse-spleen trVAE embedding once, preserve cell identifiers and ordering, store it in a disk-backed format, and record its checksum and provenance. Treat this embedding as immutable input to the R graph, aggregation, AutoK, GMM, enrichment, and proximity stages. Python may be used to generate reference artifacts, but it must not become a hidden runtime dependency of the R package. PCA may remain only as a separately labelled R-native baseline and must not be presented as end-to-end CellCharter parity.

3. Establish graph-construction parity

Create deterministic Python-R fixtures that compare sample-local Delaunay edges, spatial distances, the 99th-percentile cutoff, and the final trimmed adjacency. Resolve and document whether the cutoff is computed independently per sample or globally over the combined block-diagonal graph, following the selected pinned Python reference. Test tied coordinates, duplicated points, degenerate geometries, isolated cells, and sample boundaries. Never create edges between physical samples.

4. Establish neighborhood-aggregation parity

Compare the R and pinned Python outputs for the native representation and exclusive one-, two-, and three-hop neighborhood means. Verify cell order, feature order, sparse-hop semantics, zero-neighbor behavior, numeric precision, and sample isolation. Store a small permanent parity fixture and use tolerances appropriate for the Python float32 reference. Keep tutorial-scale output HDF5-backed and avoid cohort-wide realization.

5. Harden and validate the R Torch GMM

Preserve a batched full-covariance EM implementation with device = "auto" and a working CPU fallback. Add early validation for dimensions, finite values, batch size, iteration count, tolerance, covariance regularization, seeds, and thread counts. Handle empty or collapsed components deterministically, report the likelihood of the final fitted parameters, and make label prediction genuinely bounded-memory without first collecting the complete N-by-K probability matrix. Validate initialization, convergence, covariance regularization, serialization, reproducibility, and label-permutation-invariant agreement against the pinned Python reference. Log the selected device, R Torch version, CUDA availability, GPU name when available, batch size, seed, and CPU thread count.

6. Implement CellCharter AutoK in R

Declare the candidate K range before running the analysis. For every required K, run repeated deterministic GMM fits with max_runs = 10 and convergence_tol = 0.01. Reproduce the pinned CellCharter Fowlkes-Mallows stability calculation between neighboring K solutions, its convergence rule, peak candidates, best-K selection, and selection of the lowest-negative-log-likelihood model for each K. Save every fitted model and label set as it completes, keep large label collections disk-backed, and write durable checkpoints for stability values, convergence history, selected K, and the selected model. Use the R-selected K for all primary results. K = 11 is allowed only as a clearly labelled tutorial-matched sensitivity analysis.

7. Implement the required CellCharter downstream statistics

Implement cell-type enrichment as the CellCharter observed-versus-expected enrichment, including log2 fold change and any requested permutation results; simple within-domain composition is not a substitute. Implement asymmetric, inter-domain neighborhood enrichment using the pinned CellCharter analytical expectation. For the CODEX tutorial, compute neighborhood enrichment separately for BALBc and MRL. Implement differential neighborhood enrichment for BALBc versus MRL using sample as the library unit, n_perms = 100, and permutation of condition assignments at the sample/library level. Save observed, expected, enrichment, differential-effect, permutation-count, and empirical-p-value matrices with complete parameters and seeds.

8. Build a checkpointed, provenance-rich CODEX workflow and report

Provide an executable workflow that regenerates every absent artifact from documented input paths without assuming an old RDS exists. Save durable checkpoints under results/ for the fixed representation, per-sample graphs, HDF5 neighborhood representation, every AutoK fit and label set, stability results, selected model and labels, enrichment matrices, neighborhood-enrichment matrices, differential permutation results, and rendered-report inputs. Record input checksums, reference version and commit, parameters, seeds, package versions, device information, batch size, CPU threads, timings, and measured memory use. Produce a visible HTML report with all code chunks, parameters, results, figures, selected-K primary analyses, the separately labelled K = 11 sensitivity analysis, and explicit limitations.

## Safety

- Work only inside this repository unless the user explicitly approves another path.
- You may access public documentation, tutorials, package repositories, and other
  relevant web resources for read-only research.
- Ask before downloading large files or datasets, especially multi-gigabyte data.
- Do not commit secrets, credentials, raw private data, or generated data caches to Git.


## Second-round independent review workflow

Use this workflow whenever the user asks for a second revision, independent audit,
or final methodological review of the rCellCharter implementation.

### Global rules

- Do not modify code during the review phase.
- The second review must be independent from the implementation reasoning that produced
  the current code. Reviewers must actively try to falsify correctness rather than confirm it.
- Spawn independent subagents only for non-overlapping review scopes.
- Every reviewer must inspect the actual repository code and, when relevant, compare it
  against the official CellCharter tutorial, source implementation, and documented method.
- Public documentation/tutorial/source access is allowed for read-only research.
- Reviewers must cite exact repository files/functions and reference evidence where possible.
- Do not accept "looks reasonable" as evidence of equivalence.
- Do not let reviewers modify files, commit, reformat, or "fix while reviewing."
- The primary agent must wait for all reviewer reports before synthesizing the second-round plan.
- Disagreements between reviewers must be explicitly reconciled before implementation begins.
- Implementation starts only after the user approves the synthesized revision plan.
- During implementation, never let multiple agents concurrently edit the same files or tightly
  coupled functions.
- After implementation, use a fresh reviewer that did not author the changes.

### Phase 0 - Freeze and checkpoint the current implementation

Before spawning reviewers:

1. Confirm the current first-round task is at a safe stopping boundary.
2. Save/checkpoint expensive intermediate outputs required for reproducibility.
3. Record the current implementation status, known limitations, failed tests, and unresolved
   questions in the existing progress/status notes (or create a concise review-status note if
   one does not exist).
4. Inspect `git status` and `git diff` so the exact baseline under review is known.
5. Do not start new implementation work after this point until the second-round synthesis is
   complete and approved.
6. Rename the tool as: CellCharterR. In result rename every relevant file or folder
7. Re-do the analysis of mouse spleen with both, K=4 and K=11 (matching CellCharter tutorial). For     both do:
  - Display just BALBc-1 and MRL4 (early) slides to see the clusters
  - Cluster ~ cell type dot plots with log2 FC. Exaclty as in the original tutorial
  - Proximity analysis:
    - Heatmap displaying neighborhood enrichment for healthy and lupus
    - Heatmap of the differential neighborhood enrichment (MRL vs BALBc) and (BALBc vs MRL)

### Phase 1 - Spawn five independent reviewers in parallel

Spawn all five reviewers at the same time. Reviewers are READ-ONLY.

#### Reviewer 1 - Reference-method equivalence

Audit the current implementation against the official CellCharter tutorial, documented method,
and reference source code. Focus on algorithmic and methodological equivalence.

Check especially:
- preprocessing assumptions and required inputs,
- spatial graph construction expectations,
- neighborhood aggregation semantics,
- AutoK/model-selection logic,
- GMM fitting and prediction semantics,
- defaults/parameters that materially affect results,
- any behavior that has been approximated rather than reproduced.

#### Reviewer 2 - Graph and neighborhood construction

Audit:
- per-sample/per-slide graph isolation,
- Delaunay or other graph construction,
- long-edge filtering,
- boundary handling,
- neighborhood definition and aggregation,
- normalization/scaling applied before or after aggregation, accidental cross-sample links,
- reproducibility and deterministic behavior.

Try to construct edge cases that would expose an incorrect graph or neighborhood implementation.

#### Reviewer 3 - AutoK and GMM statistical correctness

Audit:
- candidate-K search logic,
- initialization,
- covariance parameterization,
- EM E-step and M-step equations,
- covariance regularization,
- log-likelihood computation,
- convergence criteria,
- stopping rules,
- model-selection criterion,
- prediction/posterior assignment,
- empty/degenerate components,
- numerical behavior for very small probabilities,
- whether batching changes the mathematical result.

Do not infer correctness from successful convergence alone.

#### Reviewer 4 - Numerical, GPU, and performance behavior

Audit:
- dtype and device placement,
- hidden CPU fallbacks,
- CPU/GPU behavioral differences,
- deterministic/reproducible seeds,
- batching/chunking semantics,
- memory scaling,
- numerical stability,
- synchronization or transfer bottlenecks,
- whether CUDA execution changes results beyond reasonable floating-point tolerance.

Where feasible, propose CPU-vs-GPU parity tests on a small deterministic fixture.

#### Reviewer 5 - Adversarial test reviewer

Assume the implementation is wrong and try to prove it.

Design tests that would expose:
- silent algorithmic deviations,
- tests that can pass despite incorrect outputs,
- label/permutation mistakes,
- shape/order mismatches,
- sample leakage,
- unstable cluster assignments,
- incorrect edge cases,
- non-deterministic failures,
- numerical pathologies,
- reference comparisons that are too weak to be meaningful.

Prefer small deterministic fixtures with known expected behavior whenever possible.

### Standard output required from every reviewer

Each reviewer must return a structured report containing:

1. Finding title.
2. Severity: BLOCKER / HIGH / MEDIUM / LOW / INFORMATIONAL.
3. Confidence: CONFIRMED / LIKELY / SPECULATIVE.
4. Exact repository file(s) and function(s) involved.
5. What the current implementation does.
6. What the reference method or mathematically correct behavior requires.
7. Evidence supporting the finding.
8. Consequence if left unchanged.
9. A concrete validation test that could confirm or falsify the finding.
10. Recommended action: FIX / TEST FIRST / DOCUMENT / NO CHANGE.

Reviewers must also explicitly report "No issue found" for major areas they checked and found
consistent, so silence is not interpreted as validation.

### Phase 2 - Primary-agent synthesis

Only after ALL five reviewers have finished:

1. Collect all findings into one issue matrix.
2. Merge duplicate findings, preserving independent evidence.
3. For conflicting findings, inspect the relevant code/reference evidence directly and state
   which interpretation is better supported and why.
4. Separate confirmed defects from hypotheses that still require a validation experiment.
5. Define an acceptance test for every BLOCKER and HIGH-severity item.
6. Produce ONE revision plan ordered strictly as:
   a. correctness-critical,
   b. reproducibility/numerical stability,
   c. performance/scalability,
   d. code quality/documentation.
7. Estimate which items touch overlapping files/functions and therefore must be implemented
   sequentially.
8. Do not modify files during synthesis.
9. Present the synthesized plan to the user and WAIT for approval.

### Phase 3 - Correctness implementation

After user approval:

- Implement BLOCKER and HIGH correctness items first.
- Prefer one implementation agent for tightly coupled statistical code.
- A second implementation agent is allowed only when file/function scopes are clearly
  independent.
- After EACH accepted change:
  1. run the targeted acceptance test,
  2. run directly affected unit tests,
  3. checkpoint expensive outputs if regenerated,
  4. inspect `git diff` and `git status`,
  5. stop if the result disagrees with the reference or acceptance criterion.
- Do not start performance optimization until correctness-critical items pass.

### Phase 4 - Reproducibility, GPU, performance, and cleanup

Proceed in this order:

1. reproducibility and deterministic behavior,
2. CPU/GPU numerical parity and stability,
3. memory/performance improvements,
4. code-quality/documentation cleanup.

Do not trade reference equivalence for speed unless the deviation is explicitly documented and
approved by the user.

### Phase 5 - Fresh post-implementation adversarial review

Spawn a NEW reviewer that did not participate in implementation.

Its task is to treat the current repository as an external submission and try to falsify the
claim that rCellCharter reproduces CellCharter faithfully.

The fresh reviewer must report only:
- blocking errors,
- methodological deviations,
- numerical/reproducibility risks,
- missing or weak tests,
- evidence supporting equivalence where it survives attempted falsification.

It must not modify files.

### Phase 6 - Final acceptance

Only after the fresh review is resolved:

- run the complete relevant test suite,
- run package checks appropriate for the repository,
- rerun the smallest meaningful reference-validation workflow,
- verify GPU execution when GPU behavior is part of the implementation,
- produce a final table of:
  - exact matches to reference behavior,
  - deliberate approximations/deviations,
  - unresolved limitations,
  - validation evidence.

Never describe the implementation as a faithful reproduction while a BLOCKER/HIGH methodological
issue remains unresolved.

### Usage-limit / interruption behavior

- Treat every phase boundary as a resumable checkpoint.
- If a usage-limit warning appears, do not begin a new major phase.
- Finish only the current atomic operation needed to leave files/results consistent, then:
  - save/checkpoint expensive outputs,
  - update the review/progress status,
  - record the exact next action,
  - stop cleanly.
- After resumption, first read the status note, inspect repository state, and continue from the
  recorded next action rather than repeating completed review work.

