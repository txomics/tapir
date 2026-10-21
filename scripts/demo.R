library(SummarizedExperiment)
library(plyxp)
library(plyranges)
devtools::load_all("~/proj/bea-campillo/splicelogic")

# simulate a small dataset: 2 genes x 2 isoforms, 3 control + 3 treated samples
set.seed(5)
n_per_group <- 3
coldata <- data.frame(
  condition = factor(rep(c("control", "treated"), each = n_per_group))
)
rownames(coldata) <- paste0("s", 1:6)

rowdata <- data.frame(
  isoform = c("A-201", "A-202", "B-201", "B-202"),
  gene_id = c("A", "A", "B", "B")
)
rownames(rowdata) <- rowdata$isoform

# gene A switches from isoform 201 to 202 upon treatment, gene B is unchanged
mu <- cbind(
  matrix(c(80, 20, 50, 50), nrow = 4, ncol = n_per_group),
  matrix(c(20, 80, 50, 50), nrow = 4, ncol = n_per_group)
)
counts <- matrix(rpois(length(mu), mu), nrow = 4,
                 dimnames = list(rownames(rowdata), rownames(coldata)))

se <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = rowdata,
  colData = coldata
)

xp <- se |>
  new_plyxp()

xp <- xp |>
  group_by(rows(gene_id)) |>
  mutate(
    iso_prop = t(t(counts) / colSums(counts))
  ) |>
  ungroup()

xp
assay(xp, "iso_prop") |> round(2)

# mean and SD of isoform proportions per condition
iso_prop_stats <- xp |>
  group_by(cols(condition)) |>
  summarize(
    prop_mean = rowMeans(iso_prop),
    prop_sd = rowSds(iso_prop),
    cols(.samples = unique(condition))
  )

assay(iso_prop_stats, "prop_mean") |> round(2)
assay(iso_prop_stats, "prop_sd") |> round(3)

# exon ranges for the transcripts, with DTU results
# gene A starts at 101, gene B at 1001; the 202 isoforms (positive
# coefficient) skip the middle exon, the 201 isoforms include it
exon_coords <- data.frame(
  start = c(101, 401, 700),
  end = c(300, 600, 900)
)
coefficient <- c("A-201" = -2, "A-202" = 2, "B-201" = -0.1, "B-202" = 0.1)
pvalue <- c("A-201" = 1e-6, "A-202" = 1e-6, "B-201" = 0.6, "B-202" = 0.6)
gene_exons <- function(offset, gene) {
  tx_id <- rep(paste0(gene, c("-201", "-202")), times = c(3, 2))
  data.frame(
    seqnames = "chr1",
    start = c(exon_coords$start, exon_coords$start[-2]) + offset,
    end = c(exon_coords$end, exon_coords$end[-2]) + offset,
    strand = "+",
    tx_id = tx_id,
    gene_id = gene,
    exon_rank = c(1:3, 1:2),
    coefficient = unname(coefficient[tx_id]),
    pvalue = unname(pvalue[tx_id])
  )
}
exons <- rbind(gene_exons(0, "A"), gene_exons(1000, "B")) |>
  as_granges()

exons

# add isoform proportion mean and SD per condition, looked up by tx_id
exons <- exons |>
  mutate(
    prop_mean_control = assay(iso_prop_stats, "prop_mean")[tx_id, "control"],
    prop_mean_treated = assay(iso_prop_stats, "prop_mean")[tx_id, "treated"],
    prop_sd_control = assay(iso_prop_stats, "prop_sd")[tx_id, "control"],
    prop_sd_treated = assay(iso_prop_stats, "prop_sd")[tx_id, "treated"]
  )

exons

# find skipped exons: exons of negative-coefficient transcripts that
# positive-coefficient transcripts of the same gene skip
skipped <- exons |>
  preprocess(
    coef_col = "coefficient",
    additional_columns = c(
      "pvalue",
      "prop_mean_control", "prop_mean_treated",
      "prop_sd_control", "prop_sd_treated"
    )
  ) |>
  find_skipped_exons()

skipped |>
  select(
    pval = pvalue,
    other_pval = event_pvalue, 
    ctl = prop_mean_control, 
    trt = prop_mean_treated,
    other_ctl = event_prop_mean_control,
    other_trt = event_prop_mean_treated
  )

