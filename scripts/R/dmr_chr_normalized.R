library(ggplot2)

fai <- read.table("~/tfm/reference_genome/Heliconius_erato_lativitta_v1_-_scaffolds.fa.fai",
                  col.names=c("scaffold","length","offset","linebases","linewidth"))
fai$chr <- gsub("_[0-9]+$", "", fai$scaffold)
chr_lengths <- aggregate(length ~ chr, data=fai, FUN=sum)
chr_lengths$length_mb <- chr_lengths$length / 1e6

dmr_5mC  <- read.table("analisis/annotation/dmr.5mC.txt",  header=TRUE, sep="\t")
dmr_5hmC <- read.table("analisis/annotation/dmr.5hmC.txt", header=TRUE, sep="\t")

dmr_5mC$chr  <- gsub("_[0-9]+$", "", dmr_5mC$chr)
dmr_5hmC$chr <- gsub("_[0-9]+$", "", dmr_5hmC$chr)

dmr_5mC$direction  <- ifelse(dmr_5mC$diff.Methy  > 0, "Hypomethylated", "Hypermethylated")
dmr_5hmC$direction <- ifelse(dmr_5hmC$diff.Methy > 0, "Hypomethylated", "Hypermethylated")

count_dir <- function(dmr_df, lengths_df) {
  hypo  <- as.data.frame(table(dmr_df$chr[dmr_df$direction=="Hypomethylated"]))
  hyper <- as.data.frame(table(dmr_df$chr[dmr_df$direction=="Hypermethylated"]))
  colnames(hypo)  <- c("chr", "hypo")
  colnames(hyper) <- c("chr", "hyper")
  data <- merge(lengths_df, hypo,  by="chr", all.x=TRUE)
  data <- merge(data,       hyper, by="chr", all.x=TRUE)
  data[is.na(data)] <- 0
  data$hypo_per_mb  <- -data$hypo  / data$length_mb
  data$hyper_per_mb <-  data$hyper / data$length_mb
  data$chr_label <- gsub("Hel_chr", "Chr", data$chr)
  data
}

data_5mC  <- count_dir(dmr_5mC,  chr_lengths)
data_5hmC <- count_dir(dmr_5hmC, chr_lengths)

order_chr <- data_5mC$chr_label[order(-(data_5mC$hyper_per_mb + abs(data_5mC$hypo_per_mb)))]
data_5mC$chr_label  <- factor(data_5mC$chr_label,  levels=order_chr)
data_5hmC$chr_label <- factor(data_5hmC$chr_label, levels=order_chr)

make_panel <- function(data, title, ymax) {
  ggplot(data) +
    geom_bar(aes(x=chr_label, y=hyper_per_mb),
             stat="identity", fill="#E8470A") +
    geom_bar(aes(x=chr_label, y=hypo_per_mb),
             stat="identity", fill="#D3D1C7") +
    geom_hline(yintercept=0, color="grey40", linewidth=0.5) +
    scale_y_continuous(limits=c(-ymax, ymax),
                       labels=function(x) abs(round(x, 1))) +
    annotate("text", x=Inf, y=ymax*0.85, hjust=1.1, vjust=1.5,
             size=3.5, color="#E8470A", fontface="bold",
             label="Hypermethylated") +
    annotate("text", x=Inf, y=-ymax*0.85, hjust=1.1, vjust=-0.5,
             size=3.5, color="#888780", fontface="bold",
             label="Hypomethylated") +
    labs(x="Chromosome", y="DMRs per Mb", title=title) +
    theme_classic(base_size=12) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          plot.title=element_text(face="bold", hjust=0.5))
}

p_5mC  <- make_panel(data_5mC,  "5mC DMR density by chromosome",  15)
p_5hmC <- make_panel(data_5hmC, "5hmC DMR density by chromosome",  3)

ggsave("analisis/figures/dmr_chr_normalized_5mC.pdf",
       p_5mC, width=10, height=6, bg="white")
ggsave("analisis/figures/dmr_chr_normalized_5mC.png",
       p_5mC, width=10, height=6, dpi=300, bg="white")

ggsave("analisis/figures/dmr_chr_normalized_5hmC.pdf",
       p_5hmC, width=10, height=6, bg="white")
ggsave("analisis/figures/dmr_chr_normalized_5hmC.png",
       p_5hmC, width=10, height=6, dpi=300, bg="white")

cat("Done.\n")