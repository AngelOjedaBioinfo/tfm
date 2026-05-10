library(ggplot2)

# 1. Volcano plot DMLs 5mC
cat("Cargando DMLs...\n")
dml <- read.table("dml.5mC.txt", header=TRUE, sep="\t")

dml$category <- "Not significant"
dml$category[dml$fdr < 0.05 & dml$diff > 0.1]  <- "Hypomethylated in treatment"
dml$category[dml$fdr < 0.05 & dml$diff < -0.1] <- "Hypermethylated in treatment"

colors <- c(
  "Not significant"             = "grey70",
  "Hypomethylated in treatment" = "#1D9E75",
  "Hypermethylated in treatment" = "#D85A30"
)

p1 <- ggplot(dml, aes(x=diff, y=-log10(fdr), color=category)) +
  geom_point(size=0.4, alpha=0.4) +
  scale_color_manual(values=colors) +
  geom_vline(xintercept=c(-0.1, 0.1), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  labs(x="Methylation difference (Control - Treatment)", y="-log10(FDR)", color="", title="Differentially methylated positions (5mC)") +
  theme_classic(base_size=12) +
  theme(legend.position="bottom")

ggsave("volcano_5mC.pdf", p1, width=8, height=6)
ggsave("volcano_5mC.png", p1, width=8, height=6, dpi=300)
cat("Volcano plot guardado.\n")

# 2. Distribucion de DMRs por cromosoma
cat("Cargando DMRs...\n")
dmr <- read.table("dmr.5mC.txt", header=TRUE, sep="\t")

dmr$chrom <- gsub("Hel_chr", "Chr", gsub("_[0-9]+$", "", dmr$chr))
dmr$direction <- ifelse(dmr$diff.Methy > 0, "Hypomethylated", "Hypermethylated")

dmr_count <- as.data.frame(table(dmr$chrom, dmr$direction))
colnames(dmr_count) <- c("chrom", "direction", "count")
dmr_count$count[dmr_count$direction == "Hypermethylated"] <- -dmr_count$count[dmr_count$direction == "Hypermethylated"]

p2 <- ggplot(dmr_count, aes(x=reorder(chrom, -abs(count)), y=count, fill=direction)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Hypomethylated"="#1D9E75", "Hypermethylated"="#D85A30")) +
  labs(x="Chromosome", y="Number of DMRs", fill="", title="DMR distribution by chromosome (5mC)") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="bottom")

ggsave("dmr_by_chromosome.pdf", p2, width=10, height=6)
ggsave("dmr_by_chromosome.png", p2, width=10, height=6, dpi=300)
cat("Distribucion por cromosoma guardada.\n")

# 3. Distribucion de diferencias en DMRs
dmr$direction <- ifelse(dmr$diff.Methy > 0, "Hypomethylated", "Hypermethylated")

p3 <- ggplot(dmr, aes(x=diff.Methy, fill=direction)) +
  geom_histogram(bins=60, color="white", linewidth=0.2) +
  scale_fill_manual(values=c("Hypomethylated"="#1D9E75", "Hypermethylated"="#D85A30")) +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  labs(x="Methylation difference (Control - Treatment)", y="Number of DMRs", fill="", title="Distribution of methylation differences across DMRs (5mC)") +
  theme_classic(base_size=12) +
  theme(legend.position="bottom")

ggsave("dmr_diff_distribution.pdf", p3, width=8, height=5)
ggsave("dmr_diff_distribution.png", p3, width=8, height=5, dpi=300)
cat("Distribucion guardada.\n")

cat("Figuras generadas: volcano_5mC, dmr_by_chromosome, dmr_diff_distribution\n")