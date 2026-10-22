# simulation helpers for demo.R

# simulate a small dataset: 2 genes x 2 isoforms, 3 control + 3 treated samples
# gene A switches from isoform 201 to 202 upon treatment, gene B is unchanged
simulate_se <- function(n_per_group = 3) {
  coldata <- data.frame(
    condition = factor(rep(c("control", "treated"), each = n_per_group))
  )
  rownames(coldata) <- paste0("s", seq_len(2 * n_per_group))

  rowdata <- data.frame(
    isoform = c("A-201", "A-202", "B-201", "B-202"),
    gene_id = c("A", "A", "B", "B")
  )
  rownames(rowdata) <- rowdata$isoform

  mu <- cbind(
    matrix(c(80, 20, 50, 50), nrow = 4, ncol = n_per_group),
    matrix(c(20, 80, 50, 50), nrow = 4, ncol = n_per_group)
  )
  dimnames(mu) <- list(rownames(rowdata), rownames(coldata))
  counts <- matrix(rpois(length(mu), mu), nrow = 4, dimnames = dimnames(mu))

  # keep the Poisson means alongside the simulated counts
  SummarizedExperiment(
    assays = list(counts = counts, mu = mu),
    rowData = rowdata,
    colData = coldata
  )
}

# exon ranges for the transcripts, with DTU results
# gene A starts at 101, gene B at 1001; the 202 isoforms (positive
# coefficient) skip the middle exon, the 201 isoforms include it
simulate_exons <- function() {
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
  rbind(gene_exons(0, "A"), gene_exons(1000, "B")) |>
    as_granges()
}
