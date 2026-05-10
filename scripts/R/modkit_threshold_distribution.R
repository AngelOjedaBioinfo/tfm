library(ggplot2)
library(dplyr)

# --- load probability distributions ---
load_probs <- function(path, sample) {
  df <- read.table(path, header = TRUE, sep = "\t",
                   col.names = c("code", "base", "range_start", "range_end",
                                 "count", "frac", "percentile_rank"))
  df$sample <- sample
  return(df)
}

C_probs <- load_probs("analisis/bam_files/prob_distributions/C/probabilities.tsv", "Control")
T_probs <- load_probs("analisis/bam_files/prob_distributions/T/probabilities.tsv", "Treatment")

dat <- bind_rows(C_probs, T_probs)

# keep only 5mC (m) and 5hmC (h)
dat <- dat %>%
  filter(code %in% c("m", "h")) %>%
  mutate(
    modification = case_when(
      code == "m" ~ "5mC",
      code == "h" ~ "5hmC"
    ),
    midpoint = (range_start + range_end) / 2
  )

# threshold line at 0.79
threshold <- 0.79

# --- plot ---
p <- ggplot(dat, aes(x = midpoint, y = frac * 100, color = sample)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_vline(xintercept = threshold,
             linetype = "dashed", color = "#E8470A", linewidth = 0.8) +
  annotate("text", x = threshold + 0.01, y = Inf,
           label = paste0("p10 = ", threshold),
           hjust = 0, vjust = 1.5,
           color = "#E8470A", size = 3.5, fontface = "italic") +
  annotate("rect",
           xmin = -Inf, xmax = threshold,
           ymin = -Inf, ymax = Inf,
           fill = "#E8470A", alpha = 0.05) +
  facet_wrap(~ modification, scales = "free_y", ncol = 2) +
  scale_color_manual(
    name   = NULL,
    values = c("Control" = "#2166AC", "Treatment" = "#D6604D")
  ) +
  scale_x_continuous(
    limits = c(0.33, 1.0),
    breaks = seq(0.4, 1.0, by = 0.1),
    labels = function(x) round(x, 1)
  ) +
  labs(
    x     = "Modification probability",
    y     = "Frequency (%)",
    title = "Empirical distribution of modification probabilities",
    subtitle = "Dashed line: p10 threshold (0.79) — calls to the left are discarded"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "bottom",
    strip.text         = element_text(face = "bold", size = 11),
    strip.background   = element_rect(fill = "grey95", color = NA),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, color = "grey40")
  )

ggsave("analisis/figures/modkit_threshold_distribution.pdf", p, width = 8, height = 4)
ggsave("analisis/figures/modkit_threshold_distribution.png", p, width = 8, height = 4, dpi = 300)

cat("Done. Figure saved to analisis/figures/modkit_threshold_distribution.png\n")
