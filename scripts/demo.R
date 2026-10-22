library(SummarizedExperiment)
library(plyxp)
library(plyranges)
devtools::load_all("../splicelogic")

# some functions to simulate data and exons
library(here)
source(here("scripts/demo_sims.R"))

# simulate a small dataset: 2 genes x 2 isoforms, 3 control + 3 treated samples
# gene A switches from isoform 201 to 202 upon treatment, gene B is unchanged
set.seed(5)
se <- simulate_se(n_per_group = 3)

xp <- se |>
  new_plyxp()

# compute isoform proportion per sample
xp <- xp |>
  group_by(rows(gene_id)) |>
  mutate(
    iso_prop = t(t(counts) / colSums(counts))
  ) |>
  ungroup()

xp
assay(xp, "iso_prop") |> round(2)

# the true means for the first gene
assay(xp, "mu")[1:2,]

# mean and SD of isoform proportions per condition
iso_stats <- xp |>
  group_by(cols(condition)) |>
  summarize(
    mean = rowMeans(iso_prop),
    sd = rowSds(iso_prop),
    cols(.samples = unique(condition))
  )

# examine the statistics
assay(iso_stats, "mean") |> round(2)
assay(iso_stats, "sd") |> round(3)

# now we could use those statistics while running splicelogic
# as an example, first we create some demo transcripts with exons

# exon ranges for the transcripts, with DTU results
# the 202 isoforms (positive coefficient) skip the middle exon,
# the 201 isoforms include it
exons <- simulate_exons()

exons

# add isoform proportion mean and SD per condition, looked up by tx_id
exons <- exons |>
  mutate(
    mean_ctl = assay(iso_stats, "mean")[tx_id, "control"],
    mean_trt = assay(iso_stats, "mean")[tx_id, "treated"],
    sd_ctl   = assay(iso_stats,   "sd")[tx_id, "control"],
    sd_trt   = assay(iso_stats,   "sd")[tx_id, "treated"]
  )


# find skipped exons: exons of negative-coefficient transcripts that
# positive-coefficient transcripts of the same gene skip.
# first preprocess() marks the key columns for inclusion downstream
exons <- exons |>
  preprocess(
    coef_col = "coefficient",
    additional_columns = c(
      "pvalue",
      "mean_ctl", "mean_trt",
      "sd_ctl", "sd_trt"
    )
  )

# usually we would first subset to significant txps.
# here we include gene B
skipped <- exons |>
  find_skipped_exons()

# make a small table
skipped |>
  select(
    gene_id,
    pvalue,
    up_pvalue = event_pvalue,
    ctl = mean_ctl, 
    trt = mean_trt,
    up_ctl = event_mean_ctl,
    up_trt = event_mean_trt
  )

# interpretation: gene A has two significant txps,
# the one with the skipped exon goes down with treatment (.8 -> .2)
# the one that excludes this exon goes up with treatment (.2 -> .8)
