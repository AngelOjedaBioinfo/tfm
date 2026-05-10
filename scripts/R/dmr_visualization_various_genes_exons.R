library(ggplot2)
library(dplyr)
library(patchwork)

# --- load GFF3 exons ---
load_exons <- function(gff_path) {
  con <- gzcon(file(gff_path, "rb"))
  lines <- readLines(con)
  close(con)
  lines <- lines[!grepl("^#", lines)]
  parts <- strsplit(lines, "\t")
  parts <- parts[sapply(parts, length) == 9]
  df <- as.data.frame(do.call(rbind, parts), stringsAsFactors = FALSE)
  colnames(df) <- c("chr","source","feature","start","end","score","strand","frame","attr")
  df <- df[df$feature == "exon", ]
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  return(df)
}

cat("Loading GFF3 exons...\n")
exons_all <- load_exons("reference_genome/Hlat.v1.1.annotation.CAT.gff3.gz")

# --- exon track plot ---
plot_exon_track <- function(chrom, gene_start, gene_end, exons_all,
                             dmrs, region_start, region_end) {

  exons <- exons_all %>%
    filter(chr == chrom, start >= region_start, end <= region_end)

  dmr_region <- dmrs %>%
    filter(chr == chrom, start >= region_start, end <= region_end)

  p <- ggplot() +
    annotate("rect",
             xmin = gene_start, xmax = gene_end,
             ymin = 0.2, ymax = 0.8,
             fill = "#FFB347", alpha = 0.25) +
    annotate("segment",
             x = gene_start, xend = gene_end,
             y = 0.5, yend = 0.5,
             color = "#888888", linewidth = 0.5) +
    geom_rect(data = exons,
              aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.8),
              fill = "#4A4A8A", color = NA) +
    geom_rect(data = dmr_region,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE,
              fill = "#E8470A", alpha = 0.35) +
    scale_x_continuous(labels = function(x) paste0(round(x / 1000, 1), " kb"),
                       limits = c(region_start, region_end)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(plot.margin = margin(2, 5, 0, 5))

  return(p)
}

# --- methylation track plot ---
plot_methyl_track <- function(chrom, gene_start, gene_end, dmrs,
                               bedgraph_C, bedgraph_T,
                               region_start, region_end) {

  C <- bedgraph_C %>%
    filter(chr == chrom, start >= region_start, end <= region_end) %>%
    mutate(condition = "Control", pos = (start + end) / 2)

  T <- bedgraph_T %>%
    filter(chr == chrom, start >= region_start, end <= region_end) %>%
    mutate(condition = "Treatment", pos = (start + end) / 2)

  dat <- bind_rows(C, T)

  dmr_region <- dmrs %>%
    filter(chr == chrom, start >= region_start, end <= region_end)

  p <- ggplot(dat, aes(x = pos, y = methylation, color = condition)) +
    annotate("rect",
             xmin = gene_start, xmax = gene_end,
             ymin = -Inf, ymax = Inf,
             fill = "#FFB347", alpha = 0.15) +
    geom_rect(data = dmr_region,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE,
              fill = "#E8470A", alpha = 0.35) +
    geom_rect(aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf,
                  fill = "Gene body"), color = NA, alpha = 0.25) +
    geom_rect(aes(xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = -Inf,
                  fill = "DMR"), color = NA, alpha = 0.35) +
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
      guide  = guide_legend(override.aes = list(alpha = c(0.25, 0.35), color = NA))
    ) +
    scale_x_continuous(labels = function(x) paste0(round(x / 1000, 1), " kb"),
                       limits = c(region_start, region_end)) +
    scale_y_continuous(limits = c(0, 100),
                       labels = function(x) paste0(x, "%")) +
    labs(x = paste0(chrom, " position"), y = "5mC (%)") +
    theme_classic(base_size = 11) +
    theme(
      legend.position    = "bottom",
      legend.box         = "horizontal",
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin        = margin(0, 5, 5, 5)
    )

  return(p)
}

# --- combined plot ---
plot_gene <- function(chrom, gene_start, gene_end, label, gene_id, flank,
                      dmrs, exons_all, bg_C, bg_T) {

  region_start <- gene_start - flank
  region_end   <- gene_end   + flank

  p_exon <- plot_exon_track(chrom, gene_start, gene_end,
                             exons_all, dmrs, region_start, region_end)

  p_meth <- plot_methyl_track(chrom, gene_start, gene_end,
                               dmrs, bg_C, bg_T,
                               region_start, region_end)

  combined <- p_exon / p_meth +
    plot_layout(heights = c(1, 6)) +
    plot_annotation(
      title = paste0(label, " (", gene_id, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )

  return(combined)
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
genes <- list(
  list("Hel_chr12_3", 460734,   468905,   "DNMT1", "Hlat.Hel_chr12_3G19438", 2000),
  list("Hel_chr11_3", 4144619,  4153903,  "PAK",   "Hlat.Hel_chr11_3G11048", 2000),
  list("Hel_chr13_5", 14145348, 14182530, "MSPS",  "Hlat.Hel_chr13_5G12433", 2000),
  list("Hel_chr17_2", 6511381,  6533751,  "DRP1",  "Hlat.Hel_chr17_2G18531", 2000),
  list("Hel_chr7_1",  966593,   1008497,  "SIN3A", "Hlat.Hel_chr7_1G4714",   2000),
  list("Hel_chr2_12", 41078,    75121,    "CHC",   "Hlat.Hel_chr2_12G13616", 2000)
)

# --- plot all genes ---
for (g in genes) {
  chrom      <- g[[1]]
  gene_start <- g[[2]]
  gene_end   <- g[[3]]
  label      <- g[[4]]
  gene_id    <- g[[5]]
  flank      <- g[[6]]

  cat("Plotting", label, "...\n")

  p <- plot_gene(chrom, gene_start, gene_end, label, gene_id, flank,
                 dmrs, exons_all, bg_C, bg_T)

  out_base <- paste0("analisis/figures/dmr_", tolower(label))
  ggsave(paste0(out_base, ".pdf"), p, width = 8, height = 4.5)
  ggsave(paste0(out_base, ".png"), p, width = 8, height = 4.5, dpi = 300)
}

cat("Done. Figures saved to analisis/figures/\n")
