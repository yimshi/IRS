# Rebuild package data. Set IRS_VANDEPUTTE_DATA_DIR to a directory containing
# the three source files to avoid downloading the full Zenodo data archive.

source_dir <- Sys.getenv("IRS_VANDEPUTTE_DATA_DIR", unset = "")
source_files <- c(
  "Vandeputte_2017_16S.tsv",
  "Vandeputte_2017_load.tsv",
  "41586_2017_BFnature24460_MOESM10_ESM.csv"
)

if (!nzchar(source_dir)) {
  archive <- tempfile(fileext = ".zip")
  extraction_dir <- tempfile("irs-vandeputte-")
  dir.create(extraction_dir)
  download.file(
    "https://zenodo.org/records/14243685/files/data.zip?download=1",
    archive,
    mode = "wb",
    method = "libcurl"
  )
  archive_index <- unzip(archive, list = TRUE)$Name
  archive_members <- vapply(source_files, function(file) {
    hits <- archive_index[basename(archive_index) == file]
    if (length(hits) != 1L) {
      stop("Could not locate exactly one archive member named ", file)
    }
    hits
  }, character(1))
  unzip(archive, files = archive_members, exdir = extraction_dir)
  paths <- file.path(extraction_dir, archive_members)
  names(paths) <- source_files
} else {
  paths <- file.path(source_dir, source_files)
  names(paths) <- source_files
}

if (!all(file.exists(paths))) {
  stop("Missing Vandeputte source file(s): ",
       paste(names(paths)[!file.exists(paths)], collapse = ", "))
}

sample_taxa <- read.delim(
  paths[["Vandeputte_2017_16S.tsv"]],
  check.names = FALSE,
  row.names = 1
)
loads <- read.delim(
  paths[["Vandeputte_2017_load.tsv"]],
  check.names = FALSE
)
study_metadata <- read.csv(
  paths[["41586_2017_BFnature24460_MOESM10_ESM.csv"]],
  check.names = FALSE
)

id_map <- setNames(loads$Sample_ID, loads$ENA_ID)
ena_id <- rownames(sample_taxa)
sample_id <- unname(id_map[ena_id])
study_metadata <- study_metadata[
  match(sample_id, study_metadata$Sample),
  ,
  drop = FALSE
]

stopifnot(
  !anyNA(sample_id),
  identical(study_metadata$Sample, sample_id),
  all(study_metadata$`Health status` %in% c("Control", "CD"))
)

counts <- t(as.matrix(sample_taxa))
storage.mode(counts) <- "integer"
colnames(counts) <- sample_id

sample_data <- data.frame(
  sample_id = sample_id,
  ena_id = ena_id,
  group = factor(
    study_metadata$`Health status`,
    levels = c("Control", "CD")
  ),
  microbial_load = unname(setNames(loads$count, loads$Sample_ID)[sample_id]),
  row.names = sample_id,
  check.names = FALSE
)

stopifnot(
  identical(colnames(counts), rownames(sample_data)),
  !anyNA(sample_data$microbial_load),
  all(sample_data$microbial_load > 0),
  identical(dim(counts), c(166L, 95L))
)

vandeputte2017 <- list(
  counts = counts,
  sample_data = sample_data
)

irs_daa_benchmark <- read.csv(
  file.path("data-raw", "irs_daa_benchmark.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
irs_daa_benchmark$analysis_settings <- paste(
  "IRS 0.2.0; DESeq2 1.44.0; edgeR 4.2.2; Maaslin2 1.18.0;",
  "GUniFrac 1.9; LDM 6.0.1; dacomp 1.26; ALDEx2 1.36.0;",
  "ZicoSeq 999 permutations; DACOMP/LDM 1000 permutations;",
  "ALDEx2 128 Monte Carlo samples"
)

dir.create("data", showWarnings = FALSE)
save(vandeputte2017, file = "data/vandeputte2017.rda", compress = "xz", version = 2)
save(irs_daa_benchmark, file = "data/irs_daa_benchmark.rda", compress = "xz", version = 2)
