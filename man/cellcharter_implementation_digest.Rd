\name{cellcharter_implementation_digest}
\alias{cellcharter_implementation_digest}
\alias{cellcharter_stage_contract}
\alias{cellcharter_checkpoint_reusable}
\alias{cellcharter_record_checkpoint}
\title{Content-addressed checkpoint contracts}
\usage{
cellcharter_implementation_digest(functions = list(), source_regions = list())

cellcharter_stage_contract(stage, implementation_digest, upstream,
  cell_order_sha256, parameters, seeds, versions, runtime = list(),
  schema_version = 1L)

cellcharter_checkpoint_reusable(record, expected, artifact_paths)

cellcharter_record_checkpoint(expected, artifact_paths,
  status = "created_under_contract", evidence = list(), path)
}
\arguments{
\item{functions}{Named relevant functions whose formals and bodies are hashed.}
\item{source_regions}{Named source regions delimited by exact begin/end markers.}
\item{stage}{Workflow stage name.}
\item{implementation_digest}{Stage-specific SHA-256 implementation digest.}
\item{upstream}{Named exact upstream artifact hashes or immutable data IDs.}
\item{cell_order_sha256}{Canonical ordered-cell SHA-256.}
\item{parameters}{Named complete material stage parameters.}
\item{seeds}{Named complete stage seeds.}
\item{versions}{Named package and pinned-reference versions.}
\item{runtime}{Named device/runtime fields relevant to reuse.}
\item{schema_version}{Positive contract schema version.}
\item{record}{Checkpoint record or its RDS path.}
\item{expected}{Expected stage contract.}
\item{artifact_paths}{Uniquely named existing artifact paths.}
\item{status}{Creation or verified-existing status.}
\item{evidence}{Named validation evidence for the record.}
\item{path}{Atomic record destination.}
}
\value{The digest, contract, logical reuse decision, or saved record, respectively.}
\description{
Creates deterministic contracts and truthful checkpoint records so resumable
work cannot be selected by file existence or dimensions alone. Verified legacy
records retain an unknown original creation contract.
}
\examples{
implementation <- cellcharter_implementation_digest(
  functions = list(enrichment = cellcharter_enrichment)
)
cells <- digest::digest(c("cell1", "cell2"), algo = "sha256")
contract <- cellcharter_stage_contract(
  stage = "example", implementation_digest = implementation,
  upstream = list(input = "example-v1"), cell_order_sha256 = cells,
  parameters = list(method = "documented"), seeds = list(seed = NA_integer_),
  versions = list(CellCharterR = "0.99.0")
)
artifact <- tempfile(fileext = ".txt")
writeLines("checkpoint", artifact)
record_path <- tempfile(fileext = ".rds")
cellcharter_record_checkpoint(
  contract, c(output = artifact), path = record_path
)
cellcharter_checkpoint_reusable(
  record_path, contract, c(output = artifact)
)
unlink(c(artifact, record_path))
}
