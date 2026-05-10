library(ggplot2)
library(ggpubr)

mean_meth <- function(bed_file, feature_file, invert=FALSE) {
  flag <- if(invert) "-v" else "-u"
  cmd <- paste("bedtools intersect -a", bed_file, "-b", feature_file, flag,
               "| awk '{sum+=$11; n++} END {print sum/n}'")
  as.numeric(system(cmd, intern=TRUE))
}

calc_intron <- function(bed_file) {
  tmp <- tempfile()
  system(paste("bedtools intersect -a", bed_file,
               "-b", paste0(ann, "genes.gff3"), "-u >", tmp))
  cmd <- paste("bedtools intersect -a", tmp,
               "-b", paste0(ann, "exons.gff3"),
               "-v | awk '{sum+=$11; n++} END {print sum/n}'")
  as.numeric(system(cmd, intern=TRUE))
}

bam <- "analisis/bam_files/"
ann <- "analisis/annotation/"

cat("Calculating methylation by context...\n")

values <- list(
  C_5mC_exon    = mean_meth(paste0(bam,"Hlat.C.5mC.cov10.bed"),  paste0(ann,"exons.gff3")),
  C_5mC_intron  = calc_intron(paste0(bam,"Hlat.C.5mC.cov10.bed")),
  C_5mC_inter   = mean_meth(paste0(bam,"Hlat.C.5mC.cov10.bed"),  paste0(ann,"genes.gff3"), invert=TRUE),
  C_5hmC_exon   = mean_meth(paste0(bam,"Hlat.C.5hmC.cov10.bed"), paste0(ann,"exons.gff3")),
  C_5hmC_intron = calc_intron(paste0(bam,"Hlat.C.5hmC.cov10.bed")),
  C_5hmC_inter  = mean_meth(paste0(bam,"Hlat.C.5hmC.cov10.bed"), paste0(ann,"genes.gff3"), invert=TRUE),
  T_5mC_exon    = mean_meth(paste0(bam,"Hlat.T.5mC.cov10.bed"),  paste0(ann,"exons.gff3")),
  T_5mC_intron  = calc_intron(paste0(bam,"Hlat.T.5mC.cov10.bed")),
  T_5mC_inter   = mean_meth(paste0(bam,"Hlat.T.5mC.cov10.bed"),  paste0(ann,"genes.gff3"), invert=TRUE),
  T_5hmC_exon   = mean_meth(paste0(bam,"Hlat.T.5hmC.cov10.bed"), paste0(ann,"exons.gff3")),
  T_5hmC_intron = calc_intron(paste0(bam,"Hlat.T.5hmC.cov10.bed")),
  T_5hmC_inter  = mean_meth(paste0(bam,"Hlat.T.5hmC.cov10.bed"), paste0(ann,"genes.gff3"), invert=TRUE)
)

cat("Values calculated:\n")
print(unlist(values))

data <- data.frame(
  context      = rep(c("Exon","Intron","Intergenic"), 4),
  condition    = c(rep("Control",3), rep("Control",3), rep("Treatment",3), rep("Treatment",3)),
  modification = c(rep("5mC",3), rep("5hmC",3), rep("5mC",3), rep("5hmC",3)),
  methylation  = c(
    values$C_5mC_exon,  values$C_5mC_intron,  values$C_5mC_inter,
    values$C_5hmC_exon, values$C_5hmC_intron, values$C_5hmC_inter,
    values$T_5mC_exon,  values$T_5mC_intron,  values$T_5mC_inter,
    values$T_5hmC_exon, values$T_5hmC_intron, values$T_5hmC_inter
  )
)

data$context      <- factor(data$context,      levels=c("Exon","Intron","Intergenic"))
data$condition    <- factor(data$condition,     levels=c("Control","Treatment"))
data$modification <- factor(data$modification, levels=c("5mC","5hmC"))

p_5mC <- ggplot(data[data$modification=="5mC",],
                aes(x=context, y=methylation, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  geom_text(aes(label=paste0(round(methylation, 2), "%")),
            position=position_dodge(width=0.7),
            vjust=-0.4, size=3, fontface="bold", color="#3D3D3A") +
  scale_fill_manual(values=c("Control"="#D3D1C7", "Treatment"="#E8470A")) +
  scale_y_continuous(limits=c(0, 2.8)) +
  labs(x="Genomic context", y="Mean methylation (%)",
       fill="", title="5mC") +
  theme_classic(base_size=12) +
  theme(legend.position="bottom",
        plot.title=element_text(face="bold", hjust=0.5),
        plot.background=element_rect(fill="white", color=NA))

p_5hmC <- ggplot(data[data$modification=="5hmC",],
                 aes(x=context, y=methylation, fill=condition)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  geom_text(aes(label=paste0(round(methylation, 2), "%")),
            position=position_dodge(width=0.7),
            vjust=-0.4, size=3, fontface="bold", color="#3D3D3A") +
  scale_fill_manual(values=c("Control"="#D3D1C7", "Treatment"="#E8470A")) +
  scale_y_continuous(limits=c(0, 0.28)) +
  labs(x="Genomic context", y="Mean methylation (%)",
       fill="", title="5hmC") +
  theme_classic(base_size=12) +
  theme(legend.position="bottom",
        plot.title=element_text(face="bold", hjust=0.5),
        plot.background=element_rect(fill="white", color=NA))

p_combined <- ggarrange(p_5mC, p_5hmC,
                        ncol=2, nrow=1,
                        common.legend=TRUE,
                        legend="bottom",
                        labels=c("A", "B"))

p_combined <- annotate_figure(p_combined,
  top=text_grob("Mean methylation by genomic context",
                face="bold", size=13, color="black"))

ggsave("analisis/figures/methylation_by_context.pdf",
       p_combined, width=10, height=5, bg="white")
ggsave("analisis/figures/methylation_by_context.png",
       p_combined, width=10, height=5, dpi=300, bg="white")
cat("Done.\n")