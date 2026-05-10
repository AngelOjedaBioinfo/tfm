library(DSS)

# Cargar datos
cat("Cargando datos 5mC...\n")
C.5mC <- read.table("bam_files/dss.5mC.C.txt", header=FALSE,
                     col.names=c("chr","pos","N","X"))
T.5mC <- read.table("bam_files/dss.5mC.T.txt", header=FALSE,
                     col.names=c("chr","pos","N","X"))

cat("Cargando datos 5hmC...\n")
C.5hmC <- read.table("bam_files/dss.5hmC.C.txt", header=FALSE,
                      col.names=c("chr","pos","N","X"))
T.5hmC <- read.table("bam_files/dss.5hmC.T.txt", header=FALSE,
                      col.names=c("chr","pos","N","X"))

# Crear objetos BSseq
cat("Creando objetos BSseq...\n")
BSobj.5mC  <- makeBSseqData(list(C.5mC, T.5mC),
                             sampleNames=c("Control","Treatment"))
BSobj.5hmC <- makeBSseqData(list(C.5hmC, T.5hmC),
                             sampleNames=c("Control","Treatment"))

# Test diferencial por sitio (DML)
cat("Corriendo DML test 5mC...\n")
dml.5mC <- DMLtest(BSobj.5mC,
                   group1="Control",
                   group2="Treatment",
                   smoothing=TRUE)

cat("Corriendo DML test 5hmC...\n")
dml.5hmC <- DMLtest(BSobj.5hmC,
                    group1="Control",
                    group2="Treatment",
                    smoothing=TRUE)

# Llamar DMRs
cat("Llamando DMRs 5mC...\n")
dmr.5mC <- callDMR(dml.5mC, p.threshold=0.05)

cat("Llamando DMRs 5hmC...\n")
dmr.5hmC <- callDMR(dml.5hmC, p.threshold=0.05)

# Guardar resultados
cat("Guardando resultados...\n")
write.table(dml.5mC,  "dml.5mC.txt",  sep="\t", quote=FALSE, row.names=FALSE)
write.table(dml.5hmC, "dml.5hmC.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(dmr.5mC,  "dmr.5mC.txt",  sep="\t", quote=FALSE, row.names=FALSE)
write.table(dmr.5hmC, "dmr.5hmC.txt", sep="\t", quote=FALSE, row.names=FALSE)

cat("Analisis completado.\n")
cat("DMLs 5mC: ", nrow(dml.5mC), "\n")
cat("DMLs 5hmC:", nrow(dml.5hmC), "\n")
cat("DMRs 5mC: ", nrow(dmr.5mC), "\n")
cat("DMRs 5hmC:", nrow(dmr.5hmC), "\n")
