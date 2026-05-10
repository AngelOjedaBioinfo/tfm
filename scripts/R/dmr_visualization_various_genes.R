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
# chrom, gene_start, gene_end, gene_label, gene_id, flank
genes <- list(
  list("Hel_chr12_3", 460734,   468905,   "DNMT1",  "Hlat.Hel_chr12_3G19438", 2000),
  list("Hel_chr11_3", 4144619,  4153903,  "PAK",    "Hlat.Hel_chr11_3G11048", 2000),
  list("Hel_chr13_5", 14145348, 14182530, "MSPS",   "Hlat.Hel_chr13_5G12433", 2000),
  list("Hel_chr17_2", 6511381,  6533751,  "DRP1",   "Hlat.Hel_chr17_2G18531", 2000),
  list("Hel_chr7_1",  966593, 1008497, "SIN3A",  "Hlat.Hel_chr7_1G4714",   2000),
  list("Hel_chr2_12", 41078,  75121,  "CHC",    "Hlat.Hel_chr2_12G13616", 2000)
)

# --- plot all genes ---
for (g in genes) {
  chrom      <- g[[1]]
  gene_start <- g[[2]]
  gene_end   <- g[[3]]
  label      <- g[[4]]
  gene_id    <- g[[5]]
  flank      <- g[[6]]

  if (is.na(gene_start)) {
    cat("Skipping", label, "— coordinates not yet available\n")
    next
  }

  cat("Plotting", label, "...\n")

  p <- plot_dmr_region(
    chrom      = chrom,
    gene_start = gene_start,
    gene_end   = gene_end,
    dmrs       = dmrs,
    bedgraph_C = bg_C,
    bedgraph_T = bg_T,
    gene_name  = paste0(label, " (", gene_id, ")"),
    flank      = flank
  )

  out_base <- paste0("analisis/figures/dmr_", tolower(label))
  ggsave(paste0(out_base, ".pdf"), p, width = 8, height = 4)
  ggsave(paste0(out_base, ".png"), p, width = 8, height = 4, dpi = 300)
}

cat("Done. Figures saved to analisis/figures/\n")
