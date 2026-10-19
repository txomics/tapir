library(SummarizedExperiment)
library(plyxp)
library(plyranges)
library(splicelogic)

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
assay(xp, "iso_prop")

# mean and SD of isoform proportions per condition
iso_prop_stats <- xp |>
  group_by(cols(condition)) |>
  summarize(
    prop_mean = rowMeans(iso_prop),
    prop_sd = rowSds(iso_prop),
    cols(.samples = unique(condition))
  )

assay(iso_prop_stats)

