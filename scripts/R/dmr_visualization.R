library(ggplot2)
library(dplyr)

# --- function to plot one gene region ---
plot_dmr_region <- function(chrom, gene_start, gene_end, dmrs,
                             bedgraph_C, bedgraph_T,
                             gene_name, flank = 1000) {

  region_start <- gene_start - flank
  region_end   <- gene_end   + flank

  # filter bedgraph to region
  C <- bedgraph_C %>%
    filter(chr == chrom, start >= region_start, end <= region_end) %>%
    mutate(condition = "Control", pos = (start + end) / 2)

  T <- bedgraph_T %>%
    filter(chr == chrom, start >= region_start, end <= region_end) %>%
    mutate(condition = "Treatment", pos = (start + end) / 2)

  dat <- bind_rows(C, T)

  # filter DMRs to region
  dmr_region <- dmrs %>%
    filter(chr == chrom, start >= region_start, end <= region_end)

  # build legend data for shading
  legend_df <- data.frame(
    label = c("Gene body", "DMR"),
    fill  = c("#FFB347",   "#E8470A"),
    alpha = c(0.25,         0.35)
  )

  p <- ggplot(dat, aes(x = pos, y = methylation, color = condition)) +
    # gene body shading
    annotate("rect",
             xmin = gene_start, xmax = gene_end,
             ymin = -Inf, ymax = Inf,
             fill = "#FFB347", alpha = 0.25) +
    # DMR shading
    geom_rect(data = dmr_region,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE,
              fill = "#E8470A", alpha = 0.35) +
    # invisible rects for shading legend
    geom_rect(aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf,
                  fill = "Gene body"), color = NA, alpha = 0.25) +
    geom_rect(aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf,
                  fill = "DMR"), color = NA, alpha = 0.35) +
    # smoothed methylation lines
    geom_smooth(method = "loess", span = 0.05,
                se = FALSE, linewidth = 0.8, alpha = 0.85) +
    scale_color_manual(
      name   = NULL,
      values = c("Control" = "#2166AC", "Treatment" = "#D6604D"),
      breaks = c("Control", "Treatment")
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("Gene body" = "#FFB347", "DMR" = "#E8470A"),
      guide  = guide_legend(
        override.aes = list(alpha = c(0.25, 0.35), color = NA)
      )
    ) +
    scale_x_continuous(labels = function(x) paste0(round(x / 1000, 1), " kb")) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title = gene_name,
      x     = paste0(chrom, " position"),
      y     = "Methylation (%)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position    = "bottom",
      legend.box         = "horizontal",
      plot.title         = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )

  return(p)
}

# --- load bedgraphs ---
cat("Loading bedgraphs...\n")

load_bedgraph <- function(path) {
  read.table(path, header = FALSE,
             col.names = c("chr", "start", "end", "methylation"))
}

bg_C <- load_bedgraph("analisis/bam_files/Hlat.C.5mC.bedgraph")
bg_T <- load_bedgraph("analisis/bam_files/Hlat.T.5mC.bedgraph")

# --- load DMRs ---
dmrs <- read.table("analisis/annotation/dmr.5mC.bed", header = FALSE,
                   col.names = c("chr", "start", "end", "length", "nCG",
                                 "meanMethy1", "meanMethy2", "diff", "stat"))

# --- define genes ---
# format: list(chrom, gene_start, gene_end, gene_name)
genes <- list(
  list("Hel_chr12_3", 460734, 468905, "DNMT1"),
  list("Hel_chr4_1",  NA,     NA,     "PAK"),    # add coordinates from dnmt_gene.bed equivalents
  list("Hel_chr1_1",  NA,     NA,     "SIN3A"),
  list("Hel_chr2_1",  NA,     NA,     "DRP1"),
  list("Hel_chr3_1",  NA,     NA,     "MSPS"),
  list("Hel_chr5_1",  NA,     NA,     "CHC")
)

# --- plot DNMT1 (confirmed coordinates) ---
cat("Plotting DNMT1...\n")

p_dnmt <- plot_dmr_region(
  chrom      = "Hel_chr12_3",
  gene_start = 460734,
  gene_end   = 468905,
  dmrs       = dmrs,
  bedgraph_C = bg_C,
  bedgraph_T = bg_T,
  gene_name  = "DNMT1 (Hlat.Hel_chr12_3G19438)",
  flank      = 2000
)

ggsave("analisis/figures/dmr_dnmt1.pdf", p_dnmt, width = 8, height = 4)
ggsave("analisis/figures/dmr_dnmt1.png", p_dnmt, width = 8, height = 4, dpi = 300)

cat("Done. Figures saved to analisis/figures/\n")
cat("For remaining genes, extract coordinates with:\n")
cat("  grep 'GENE_ID' analisis/annotation/dmr.5mC.named_genes.txt\n")
cat("  Then call plot_dmr_region() with those coordinates.\n")
