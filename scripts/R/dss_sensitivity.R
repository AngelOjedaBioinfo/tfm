library(DSS)

C.5mC  <- read.table("bam_files/dss.5mC.C.txt",  header=FALSE, col.names=c("chr","pos","N","X"))
T.5mC  <- read.table("bam_files/dss.5mC.T.txt",  header=FALSE, col.names=c("chr","pos","N","X"))
C.5hmC <- read.table("bam_files/dss.5hmC.C.txt", header=FALSE, col.names=c("chr","pos","N","X"))
T.5hmC <- read.table("bam_files/dss.5hmC.T.txt", header=FALSE, col.names=c("chr","pos","N","X"))

BSobj.5mC  <- makeBSseqData(list(C.5mC,  T.5mC),  sampleNames=c("Control","Treatment"))
BSobj.5hmC <- makeBSseqData(list(C.5hmC, T.5hmC), sampleNames=c("Control","Treatment"))

cat("DML test equal.disp=TRUE...\n")
dml.5mC.ed  <- DMLtest(BSobj.5mC,  group1="Control", group2="Treatment", equal.disp=TRUE)
dml.5hmC.ed <- DMLtest(BSobj.5hmC, group1="Control", group2="Treatment", equal.disp=TRUE)

dmr.5mC.ed  <- callDMR(dml.5mC.ed,  p.threshold=0.05)
dmr.5hmC.ed <- callDMR(dml.5hmC.ed, p.threshold=0.05)

dmr.5mC.orig  <- read.table("annotation/dmr.5mC.txt",  header=TRUE)
dmr.5hmC.orig <- read.table("annotation/dmr.5hmC.txt", header=TRUE)

cat("\n=== COMPARACION smoothing vs equal.disp ===\n")
cat("DMRs 5mC  smoothing:", nrow(dmr.5mC.orig),  "| equal.disp:", nrow(dmr.5mC.ed),  "\n")
cat("DMRs 5hmC smoothing:", nrow(dmr.5hmC.orig), "| equal.disp:", nrow(dmr.5hmC.ed), "\n")

write.table(dmr.5mC.ed,  "sensitivity_dmr.5mC.equaldisp.txt",  sep="\t", quote=FALSE, row.names=FALSE)
write.table(dmr.5hmC.ed, "sensitivity_dmr.5hmC.equaldisp.txt", sep="\t", quote=FALSE, row.names=FALSE)
