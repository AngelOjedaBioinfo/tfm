library(ggplot2)
library(dplyr)

# --- load windowed methylation ---
load_windows <- function(path, sample) {
  df <- read.table(path, header = FALSE, sep = "\t",
                   col.names = c("chr", "start", "end", "methylation"))
  df$methylation[df$methylation == "."] <- NA
  df$methylation <- as.numeric(df$methylation)
  df$sample <- sample
  df$mid <- (df$start + df$end) / 2
  return(df)
}

cat("Loading windowed data...\n")
C_win <- load_windows("analisis/bam_files/windows_10kb.C.5mC.bed", "Control")
T_win <- load_windows("analisis/bam_files/windows_10kb.T.5mC.bed", "Treatment")

dat <- bind_rows(C_win, T_win)

# --- chromosome number and cumulative position ---
dat <- dat %>%
  mutate(
    chrom_num = as.integer(gsub("Hel_chr([0-9]+)_.*", "\\1", chr)),
    scaffold_idx = as.integer(gsub("Hel_chr[0-9]+_(.*)", "\\1", chr))
  ) %>%
  filter(!is.na(chrom_num), !is.na(methylation))

# compute scaffold offsets within each chromosome
scaffold_sizes <- dat %>%
  group_by(chrom_num, chr, scaffold_idx) %>%
  summarise(scaffold_len = max(end), .groups = "drop") %>%
  arrange(chrom_num, scaffold_idx) %>%
  group_by(chrom_num) %>%
  mutate(offset = cumsum(lag(scaffold_len, default = 0))) %>%
  ungroup() %>%
  select(chr, offset)

dat <- dat %>%
  left_join(scaffold_sizes, by = "chr") %>%
  mutate(cum_pos = (mid + offset) / 1e6)

# filter chromosomes with enough data
chrom_counts <- dat %>%
  group_by(chrom_num, sample) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n >= 50)

valid_chroms <- unique(chrom_counts$chrom_num)
dat <- dat %>% filter(chrom_num %in% valid_chroms)

dat$chrom_label <- factor(paste0("Chr ", dat$chrom_num),
                           levels = paste0("Chr ", sort(unique(dat$chrom_num))))

cat("Chromosomes with data:", length(unique(dat$chrom_num)), "\n")

# --- plot ---
p <- ggplot(dat, aes(x = cum_pos, y = methylation, color = sample)) +
  geom_line(linewidth = 0.3, alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.2,
              se = FALSE, linewidth = 1.0) +
  facet_wrap(~ chrom_label, scales = "free", ncol = 3) +
  scale_color_manual(
    name   = NULL,
    values = c("Control" = "#2166AC", "Treatment" = "#D6604D")
  ) +
  scale_y_continuous(limits = c(0, NA),
                     labels = function(x) paste0(x, "%")) +
  labs(
    x        = "Position (Mb)",
    y        = "Mean 5mC methylation (%)",
    title    = "Methylation landscape across chromosomes",
    subtitle = "10 kb non-overlapping windows - LOESS trend line (span = 0.2)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold", size = 8),
    strip.background   = element_rect(fill = "grey95", color = NA),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8, color = "grey40"),
    axis.text.x        = element_text(size = 7)
  )


ggsave("analisis/figures/methylation_landscape_10kb.pdf",
       p, width = 18, height = 16)
ggsave("analisis/figures/methylation_landscape_10kb.png",
       p, width = 18, height = 16, dpi = 300)

cat("Done. Figure saved to analisis/figures/methylation_landscape_10kb.png\n")
