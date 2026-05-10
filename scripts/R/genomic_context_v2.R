library(ggplot2)

count_lines <- function(f) as.numeric(system(paste("wc -l <", f), intern=TRUE))

bam <- "analisis/bam_files/"
ann <- "analisis/annotation/"

# Baseline 5mC
total_5mC  <- count_lines(paste0(bam, "sites.5mC.C.bed"))
cds_5mC    <- count_lines(paste0(bam, "sites.5mC.C.in_cds.dedup.bed"))
utr_5mC    <- count_lines(paste0(bam, "sites.5mC.C.in_utr.dedup.bed"))
intron_5mC <- count_lines(paste0(bam, "sites.5mC.C.in_introns.bed"))
inter_5mC  <- count_lines(paste0(bam, "sites.5mC.C.intergenic.bed"))

# DMRs 5mC
total_dmr5mC  <- count_lines(paste0(ann, "dmr.5mC.bed"))
cds_dmr5mC    <- count_lines(paste0(ann, "dmr.5mC.in_cds.bed"))
utr_dmr5mC    <- count_lines(paste0(ann, "dmr.5mC.in_utr.bed"))
intron_dmr5mC <- count_lines(paste0(ann, "dmr.5mC.in_introns.bed"))
inter_dmr5mC  <- count_lines(paste0(ann, "dmr.5mC.intergenic.bed"))

# DMRs 5hmC
total_dmr5hmC  <- count_lines(paste0(ann, "dmr.5hmC.bed"))
cds_dmr5hmC    <- count_lines(paste0(ann, "dmr.5hmC.in_cds.bed"))
utr_dmr5hmC    <- count_lines(paste0(ann, "dmr.5hmC.in_utr.bed"))
intron_dmr5hmC <- count_lines(paste0(ann, "dmr.5hmC.in_introns.bed"))
inter_dmr5hmC  <- count_lines(paste0(ann, "dmr.5hmC.intergenic.bed"))

context_data <- data.frame(
  context = rep(c("CDS", "UTR", "Intron", "Intergenic"), 3),
  type = c(rep("Baseline CpG", 4),
           rep("5mC DMRs", 4),
           rep("5hmC DMRs", 4)),
  percentage = c(
    cds_5mC/total_5mC*100,    utr_5mC/total_5mC*100,
    intron_5mC/total_5mC*100, inter_5mC/total_5mC*100,
    cds_dmr5mC/total_dmr5mC*100,    utr_dmr5mC/total_dmr5mC*100,
    intron_dmr5mC/total_dmr5mC*100, inter_dmr5mC/total_dmr5mC*100,
    cds_dmr5hmC/total_dmr5hmC*100,    utr_dmr5hmC/total_dmr5hmC*100,
    intron_dmr5hmC/total_dmr5hmC*100, inter_dmr5hmC/total_dmr5hmC*100
  ),
  n = c(
    cds_5mC,    utr_5mC,    intron_5mC, inter_5mC,
    cds_dmr5mC, utr_dmr5mC, intron_dmr5mC, inter_dmr5mC,
    cds_dmr5hmC, utr_dmr5hmC, intron_dmr5hmC, inter_dmr5hmC
  )
)

context_data$context <- factor(context_data$context,
                                levels = c("CDS", "UTR", "Intron", "Intergenic"))
context_data$type <- factor(context_data$type,
                             levels = c("Baseline CpG", "5mC DMRs", "5hmC DMRs"))

print(context_data)

p <- ggplot(context_data, aes(x = type, y = percentage, fill = context)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(percentage, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 3.5, color = "#3D3D3A", fontface = "bold") +
  scale_fill_manual(values = c(
    "CDS"        = "#E8470A",
    "UTR"        = "#C0392B",
    "Intron"     = "#FFB347",
    "Intergenic" = "#D3D1C7"
  ), limits = rev) +
  labs(
    x    = "",
    y    = "Percentage (%)",
    fill = "Genomic context",
    title = "CpG site and DMR distribution by genomic context"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("analisis/figures/genomic_context.pdf", p, width = 7, height = 5)
ggsave("analisis/figures/genomic_context.png", p, width = 7, height = 5, dpi = 300)
