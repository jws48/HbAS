library("DESeq2")
library("GenomicFeatures")
library("pcaExplorer")
library("pheatmap")
library("genefilter")
library("vsn")
library("RColorBrewer")
library("gplots")
library("limma")
library("edgeR")
library("Glimma")
library("org.Pf.plasmo.db")
library("ggplot2")

data_dir <- "/Users/joe/R Projects/inVitro_Salmon/data"
cvsfile <- file.path(data_dir, "sampleTable.csv")
sampleTable <- read.csv(cvsfile, row.names = 1)
files <- file.path(data_dir, "Pfal", sampleTable$sampleName, "quant.sf")
file.exists(files)

gtffile <- file.path(data_dir, "Plasmodium_falciparum.ASM276v2.45.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
colnames(txi$counts) <- rownames(sampleTable)
colnames(txi$abundance) <- rownames(sampleTable)
colnames(txi$length) <- rownames(sampleTable)
coldata <- sampleTable
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~Genotype)
dds$group <- factor(paste0(dds$Genotype, dds$HPI, dds$Strain))
design(dds) <- ~ group
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("Stage", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Stage, shape=Genotype)) +
  geom_point(aes(fill=Stage, shape=Genotype), size=3, color="black") +
  scale_shape_manual(values = c(21, 24)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

pcaData <- plotPCA(vsd, intgroup=c("Stage", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Stage, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

countfile <- file.path(data_dir, "counts_dropVSA_salmon.csv")
countdata <- read.csv(countfile, row.names = 1)
countdata[1:128] <- lapply(countdata[1:128], as.integer)

(AS3vAA3_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS33D7", "AA33D7")))
(AS6vAA6_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS63D7", "AA63D7")))
(AS9vAA9_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS93D7", "AA93D7")))
(AS12vAA12_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS123D7", "AA123D7")))
(AS15vAA15_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS153D7", "AA153D7")))
(AS18vAA18_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS183D7", "AA183D7")))
(AS21vAA21_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS213D7", "AA213D7")))
(AS24vAA24_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS243D7", "AA243D7")))
(AS27vAA27_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS273D7", "AA273D7")))
(AS30vAA30_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS303D7", "AA303D7")))
(AS33vAA33_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS333D7", "AA333D7")))
(AS36vAA36_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS363D7", "AA363D7")))
(AS39vAA39_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS393D7", "AA393D7")))
(AS42vAA42_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS423D7", "AA423D7")))
(AS45vAA45_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS453D7", "AA453D7")))
(AS48vAA48_3D7 <- results(dds,alpha=0.05, contrast = c("group", "AS483D7", "AA483D7")))
keysdds <- rownames(AS3vAA3_3D7)
annosdds <- select(org.Pf.plasmo.db, keys = keysdds, columns =c("GENENAME"), keytype = "ORF")
AS3vAA3_3D7$genes <- annosdds
AS6vAA6_3D7$genes <- annosdds
AS9vAA9_3D7$genes <- annosdds
AS12vAA12_3D7$genes <- annosdds
AS15vAA15_3D7$genes <- annosdds
AS18vAA18_3D7$genes <- annosdds
AS21vAA21_3D7$genes <- annosdds
AS24vAA24_3D7$genes <- annosdds
AS27vAA27_3D7$genes <- annosdds
AS30vAA30_3D7$genes <- annosdds
AS33vAA33_3D7$genes <- annosdds
AS36vAA36_3D7$genes <- annosdds
AS39vAA39_3D7$genes <- annosdds
AS42vAA42_3D7$genes <- annosdds
AS45vAA45_3D7$genes <- annosdds
AS48vAA48_3D7$genes <- annosdds

write.csv(AS3vAA3_3D7, file = "ASvAA_3hpi__3D7_DESeq2_results.csv")
write.csv(AS6vAA6_3D7, file = "ASvAA_6hpi__3D7_DESeq2_results.csv")
write.csv(AS9vAA9_3D7, file = "ASvAA_9hpi__3D7_DESeq2_results.csv")
write.csv(AS12vAA12_3D7, file = "ASvAA_12hpi__3D7_DESeq2_results.csv")
write.csv(AS15vAA15_3D7, file = "ASvAA_15hpi__3D7_DESeq2_results.csv")
write.csv(AS18vAA18_3D7, file = "ASvAA_18hpi__3D7_DESeq2_results.csv")
write.csv(AS21vAA21_3D7, file = "ASvAA_21hpi__3D7_DESeq2_results.csv")
write.csv(AS24vAA24_3D7, file = "ASvAA_24hpi__3D7_DESeq2_results.csv")
write.csv(AS27vAA27_3D7, file = "ASvAA_27hpi__3D7_DESeq2_results.csv")
write.csv(AS30vAA30_3D7, file = "ASvAA_30hpi__3D7_DESeq2_results.csv")
write.csv(AS33vAA33_3D7, file = "ASvAA_33hpi__3D7_DESeq2_results.csv")
write.csv(AS36vAA36_3D7, file = "ASvAA_36hpi__3D7_DESeq2_results.csv")
write.csv(AS39vAA39_3D7, file = "ASvAA_39hpi__3D7_DESeq2_results.csv")
write.csv(AS42vAA42_3D7, file = "ASvAA_42hpi__3D7_DESeq2_results.csv")
write.csv(AS45vAA45_3D7, file = "ASvAA_45hpi__3D7_DESeq2_results.csv")
write.csv(AS48vAA48_3D7, file = "ASvAA_48hpi__3D7_DESeq2_results.csv")

(AS3vAA3_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS3FUP", "AA3FUP")))
(AS6vAA6_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS6FUP", "AA6FUP")))
(AS9vAA9_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS9FUP", "AA9FUP")))
(AS12vAA12_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS12FUP", "AA12FUP")))
(AS15vAA15_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS15FUP", "AA15FUP")))
(AS18vAA18_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS18FUP", "AA18FUP")))
(AS21vAA21_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS21FUP", "AA21FUP")))
(AS24vAA24_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS24FUP", "AA24FUP")))
(AS27vAA27_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS27FUP", "AA27FUP")))
(AS30vAA30_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS30FUP", "AA30FUP")))
(AS33vAA33_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS33FUP", "AA33FUP")))
(AS36vAA36_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS36FUP", "AA36FUP")))
(AS39vAA39_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS39FUP", "AA39FUP")))
(AS42vAA42_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS42FUP", "AA42FUP")))
(AS45vAA45_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS45FUP", "AA45FUP")))
(AS48vAA48_FUP <- results(dds, alpha = 0.05, contrast = c("group", "AS48FUP", "AA48FUP")))



keysdds <- rownames(AS3vAA3_FUP)
annosdds <- select(org.Pf.plasmo.db, keys = keysdds, columns =c("GENENAME"), keytype = "ORF")
AS3vAA3_FUP$genes <- annosdds
AS6vAA6_FUP$genes <- annosdds
AS9vAA9_FUP$genes <- annosdds
AS12vAA12_FUP$genes <- annosdds
AS15vAA15_FUP$genes <- annosdds
AS18vAA18_FUP$genes <- annosdds
AS21vAA21_FUP$genes <- annosdds
AS24vAA24_FUP$genes <- annosdds
AS27vAA27_FUP$genes <- annosdds
AS30vAA30_FUP$genes <- annosdds
AS33vAA33_FUP$genes <- annosdds
AS36vAA36_FUP$genes <- annosdds
AS39vAA39_FUP$genes <- annosdds
AS42vAA42_FUP$genes <- annosdds
AS45vAA45_FUP$genes <- annosdds
AS48vAA48_FUP$genes <- annosdds

write.csv(AS3vAA3_FUP, file = "ASvAA_3hpi__FUP_DESeq2_results.csv")
write.csv(AS6vAA6_FUP, file = "ASvAA_6hpi__FUP_DESeq2_results.csv")
write.csv(AS9vAA9_FUP, file = "ASvAA_9hpi__FUP_DESeq2_results.csv")
write.csv(AS12vAA12_FUP, file = "ASvAA_12hpi__FUP_DESeq2_results.csv")
write.csv(AS15vAA15_FUP, file = "ASvAA_15hpi__FUP_DESeq2_results.csv")
write.csv(AS18vAA18_FUP, file = "ASvAA_18hpi__FUP_DESeq2_results.csv")
write.csv(AS21vAA21_FUP, file = "ASvAA_21hpi__FUP_DESeq2_results.csv")
write.csv(AS24vAA24_FUP, file = "ASvAA_24hpi__FUP_DESeq2_results.csv")
write.csv(AS27vAA27_FUP, file = "ASvAA_27hpi__FUP_DESeq2_results.csv")
write.csv(AS30vAA30_FUP, file = "ASvAA_30hpi__FUP_DESeq2_results.csv")
write.csv(AS33vAA33_FUP, file = "ASvAA_33hpi__FUP_DESeq2_results.csv")
write.csv(AS36vAA36_FUP, file = "ASvAA_36hpi__FUP_DESeq2_results.csv")
write.csv(AS39vAA39_FUP, file = "ASvAA_39hpi__FUP_DESeq2_results.csv")
write.csv(AS42vAA42_FUP, file = "ASvAA_42hpi__FUP_DESeq2_results.csv")
write.csv(AS45vAA45_FUP, file = "ASvAA_45hpi__FUP_DESeq2_results.csv")
write.csv(AS48vAA48_FUP, file = "ASvAA_48hpi__FUP_DESeq2_results.csv")


txi_limma <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
library(limma)
colnames(txi_limma$counts) <- rownames(sampleTable)
colnames(txi_limma$abundance) <- rownames(sampleTable)
colnames(txi_limma$length) <- rownames(sampleTable)
y <- DGEList(txi_limma$counts)
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)
group <- paste(sampleTable$Genotype, sampleTable$HPI, sampleTable$Strain, sep = ".")
group <- factor(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
v <- voom(y,design,plot = TRUE)
fit <- lmFit(v)
fit <- lmFit(v)

cont.matrix <- makeContrasts(
ASvAA_3hpi_3D7 = AS.3.3D7 - AA.3.3D7, #1
ASvAA_3hpi_FUP = AS.3.FUP - AA.3.FUP, #2
FUPv3D7_3hpi_AA = AA.3.FUP - AA.3.3D7, #3
FUPv3D7_3hpi_AS = AS.3.FUP - AS.3.3D7, #4
AA_6v3hpi_3D7 = AA.6.3D7 - AA.3.3D7, #5
AS_6v3hpi_3D7 = AS.6.3D7 - AS.3.3D7, #6
AA_6v3hpi_FUP = AA.6.FUP - AA.3.FUP, #7
AS_6v3hpi_FUP = AS.6.FUP - AS.3.FUP, #8
ASvAA_6hpi_3D7 = AS.6.3D7 - AA.6.3D7, #9
ASvAA_6hpi_FUP = AS.6.FUP - AA.6.FUP, #10
FUPv3D7_6hpi_AA = AA.6.FUP - AA.6.3D7, #11
FUPv3D7_6hpi_AS = AS.6.FUP - AS.6.3D7, #12
AA_9v6hpi_3D7 = AA.9.3D7 - AA.6.3D7, #13
AS_9v6hpi_3D7 = AS.9.3D7 - AS.6.3D7, #14
AA_9v6hpi_FUP = AA.9.FUP - AA.6.FUP, #15
AS_9v6hpi_FUP = AS.9.FUP - AS.6.FUP, #16
ASvAA_9hpi_3D7 = AS.9.3D7 - AA.9.3D7, #17
ASvAA_9hpi_FUP = AS.9.FUP - AA.9.FUP, #18
FUPv3D7_9hpi_AA = AA.9.FUP - AA.9.3D7, #19
FUPv3D7_9hpi_AS = AS.9.FUP - AS.9.3D7, #20
AA_12v9hpi_3D7 = AA.12.3D7 - AA.9.3D7, #21
AS_12v9hpi_3D7 = AS.12.3D7 - AS.9.3D7, #22
AA_12v9hpi_FUP = AA.12.FUP - AA.9.FUP, #23
AS_12v9hpi_FUP = AS.12.FUP - AS.9.FUP, #24
ASvAA_12hpi_3D7 = AS.12.3D7 - AA.12.3D7, #25
ASvAA_12hpi_FUP = AS.12.FUP - AA.12.FUP, #26
FUPv3D7_12hpi_AA = AA.12.FUP - AA.12.3D7, #27
FUPv3D7_12hpi_AS = AS.12.FUP - AS.12.3D7, #28
AA_15v12hpi_3D7 = AA.15.3D7 - AA.12.3D7, #29
AS_15v12hpi_3D7 = AS.15.3D7 - AS.12.3D7, #30
AA_15v12hpi_FUP = AA.15.FUP - AA.12.FUP, #31
AS_15v12hpi_FUP = AS.15.FUP - AS.12.FUP, #32
ASvAA_15hpi_3D7 = AS.15.3D7 - AA.15.3D7, #33
ASvAA_15hpi_FUP = AS.15.FUP - AA.15.FUP, #34
FUPv3D7_15hpi_AA = AA.15.FUP - AA.15.3D7, #35
FUPv3D7_15hpi_AS = AS.15.FUP - AS.15.3D7, #36
AA_18v15hpi_3D7 = AA.18.3D7 - AA.15.3D7, #37
AS_18v15hpi_3D7 = AS.18.3D7 - AS.15.3D7, #38
AA_18v15hpi_FUP = AA.18.FUP - AA.15.FUP, #39
AS_18v15hpi_FUP = AS.18.FUP - AS.15.FUP, #40
ASvAA_18hpi_3D7 = AS.18.3D7 - AA.18.3D7, #41
ASvAA_18hpi_FUP = AS.18.FUP - AA.18.FUP, #42
FUPv3D7_18hpi_AA = AA.18.FUP - AA.18.3D7, #43
FUPv3D7_18hpi_AS = AS.18.FUP - AS.18.3D7, #44
AA_21v18hpi_3D7 = AA.21.3D7 - AA.18.3D7, #45
AS_21v18hpi_3D7 = AS.21.3D7 - AS.18.3D7, #46
AA_21v18hpi_FUP = AA.21.FUP - AA.18.FUP, #47
AS_21v18hpi_FUP = AS.21.FUP - AS.18.FUP, #48
ASvAA_21hpi_3D7 = AS.21.3D7 - AA.21.3D7, #49
ASvAA_21hpi_FUP = AS.21.FUP - AA.21.FUP, #50
FUPv3D7_21hpi_AA = AA.21.FUP - AA.21.3D7, #51
FUPv3D7_21hpi_AS = AS.21.FUP - AS.21.3D7, #52
AA_24v21hpi_3D7 = AA.24.3D7 - AA.21.3D7,
AS_24v21hpi_3D7 = AS.24.3D7 - AS.21.3D7,
AA_24v21hpi_FUP = AA.24.FUP - AA.21.FUP,
AS_24v21hpi_FUP = AS.24.FUP - AS.21.FUP,
ASvAA_24hpi_3D7 = AS.24.3D7 - AA.24.3D7,
ASvAA_24hpi_FUP = AS.24.FUP - AA.24.FUP,
FUPv3D7_24hpi_AA = AA.24.FUP - AA.24.3D7,
FUPv3D7_24hpi_AS = AS.24.FUP - AS.24.3D7,
AA_27v24hpi_3D7 = AA.27.3D7 - AA.24.3D7,
AS_27v24hpi_3D7 = AS.27.3D7 - AS.24.3D7,
AA_27v24hpi_FUP = AA.27.FUP - AA.24.FUP,
AS_27v24hpi_FUP = AS.27.FUP - AS.24.FUP,
ASvAA_27hpi_3D7 = AS.27.3D7 - AA.27.3D7,
ASvAA_27hpi_FUP = AS.27.FUP - AA.27.FUP,
FUPv3D7_27hpi_AA = AA.27.FUP - AA.27.3D7,
FUPv3D7_27hpi_AS = AS.27.FUP - AS.27.3D7,
AA_30v27hpi_3D7 = AA.30.3D7 - AA.27.3D7,
AS_30v27hpi_3D7 = AS.30.3D7 - AS.27.3D7,
AA_30v27hpi_FUP = AA.30.FUP - AA.27.FUP,
AS_30v27hpi_FUP = AS.30.FUP - AS.27.FUP,
ASvAA_30hpi_3D7 = AS.30.3D7 - AA.30.3D7,
ASvAA_30hpi_FUP = AS.30.FUP - AA.30.FUP,
FUPv3D7_30hpi_AA = AA.30.FUP - AA.30.3D7,
FUPv3D7_30hpi_AS = AS.30.FUP - AS.30.3D7,
AA_33v30hpi_3D7 = AA.33.3D7 - AA.30.3D7,
AS_33v30hpi_3D7 = AS.33.3D7 - AS.30.3D7,
AA_33v30hpi_FUP = AA.33.FUP - AA.30.FUP,
AS_33v30hpi_FUP = AS.33.FUP - AS.30.FUP,
ASvAA_33hpi_3D7 = AS.33.3D7 - AA.33.3D7,
ASvAA_33hpi_FUP = AS.33.FUP - AA.33.FUP,
FUPv3D7_33hpi_AA = AA.33.FUP - AA.33.3D7,
FUPv3D7_33hpi_AS = AS.33.FUP - AS.33.3D7,
AA_36v33hpi_3D7 = AA.36.3D7 - AA.33.3D7,
AS_36v33hpi_3D7 = AS.36.3D7 - AS.33.3D7,
AA_36v33hpi_FUP = AA.36.FUP - AA.33.FUP,
AS_36v33hpi_FUP = AS.36.FUP - AS.33.FUP,
ASvAA_36hpi_3D7 = AS.36.3D7 - AA.36.3D7,
ASvAA_36hpi_FUP = AS.36.FUP - AA.36.FUP,
FUPv3D7_36hpi_AA = AA.36.FUP - AA.36.3D7,
FUPv3D7_36hpi_AS = AS.36.FUP - AS.36.3D7,
AA_39v36hpi_3D7 = AA.39.3D7 - AA.36.3D7,
AS_39v36hpi_3D7 = AS.39.3D7 - AS.36.3D7,
AA_39v36hpi_FUP = AA.39.FUP - AA.36.FUP,
AS_39v36hpi_FUP = AS.39.FUP - AS.36.FUP,
ASvAA_39hpi_3D7 = AS.39.3D7 - AA.39.3D7,
ASvAA_39hpi_FUP = AS.39.FUP - AA.39.FUP,
FUPv3D7_39hpi_AA = AA.39.FUP - AA.39.3D7,
FUPv3D7_39hpi_AS = AS.39.FUP - AS.39.3D7,
AA_42v39hpi_3D7 = AA.42.3D7 - AA.39.3D7,
AS_42v39hpi_3D7 = AS.42.3D7 - AS.39.3D7,
AA_42v39hpi_FUP = AA.42.FUP - AA.39.FUP,
AS_42v39hpi_FUP = AS.42.FUP - AS.39.FUP,
ASvAA_42hpi_3D7 = AS.42.3D7 - AA.42.3D7,
ASvAA_42hpi_FUP = AS.42.FUP - AA.42.FUP,
FUPv3D7_42hpi_AA = AA.42.FUP - AA.42.3D7,
FUPv3D7_42hpi_AS = AS.42.FUP - AS.42.3D7,
AA_45v42hpi_3D7 = AA.45.3D7 - AA.42.3D7,
AS_45v42hpi_3D7 = AS.45.3D7 - AS.42.3D7,
AA_45v42hpi_FUP = AA.45.FUP - AA.42.FUP,
AS_45v42hpi_FUP = AS.45.FUP - AS.42.FUP,
ASvAA_45hpi_3D7 = AS.45.3D7 - AA.45.3D7,
ASvAA_45hpi_FUP = AS.45.FUP - AA.45.FUP,
FUPv3D7_45hpi_AA = AA.45.FUP - AA.45.3D7,
FUPv3D7_45hpi_AS = AS.45.FUP - AS.45.3D7,
AA_48v45hpi_3D7 = AA.48.3D7 - AA.45.3D7,
AS_48v45hpi_3D7 = AS.48.3D7 - AS.45.3D7,
AA_48v45hpi_FUP = AA.48.FUP - AA.45.FUP,
AS_48v45hpi_FUP = AS.48.FUP - AS.45.FUP,
ASvAA_48hpi_3D7 = AS.48.3D7 - AA.48.3D7,
ASvAA_48hpi_FUP = AS.48.FUP - AA.48.FUP,
FUPv3D7_48hpi_AA = AA.48.FUP - AA.48.3D7,
FUPv3D7_48hpi_AS = AS.48.FUP - AS.48.3D7,
AA_48v3hpi_3D7 = AA.48.3D7 - AA.3.3D7,
AS_48v3hpi_3D7 = AS.48.3D7 - AS.3.3D7,
AA_48v3hpi_FUP = AA.48.FUP - AA.3.FUP,
AS_48v3hpi_FUP = AS.48.FUP - AS.3.FUP, levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
keys <- rownames(fit.cont)
annos <- select(org.Pf.plasmo.db, keys = keys, columns =c("GENENAME"), keytype = "ORF")
fit.cont$genes <- annos
ASvAA_3hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_3hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_3hpi_3D7.res, file = "ASvAA_3hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_6hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_6hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_6hpi_3D7.res, file = "ASvAA_6hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_9hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_9hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_9hpi_3D7.res, file = "ASvAA_9hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_12hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_12hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_12hpi_3D7.res, file = "ASvAA_12hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_15hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_15hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_15hpi_3D7.res, file = "ASvAA_15hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_18hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_18hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_18hpi_3D7.res, file = "ASvAA_18hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_21hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_21hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_21hpi_3D7.res, file = "ASvAA_21hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_24hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_24hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_24hpi_3D7.res, file = "ASvAA_24hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_27hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_27hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_27hpi_3D7.res, file = "ASvAA_27hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_30hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_30hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_30hpi_3D7.res, file = "ASvAA_30hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_33hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_33hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_33hpi_3D7.res, file = "ASvAA_33hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_36hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_36hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_36hpi_3D7.res, file = "ASvAA_36hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_39hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_39hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_39hpi_3D7.res, file = "ASvAA_39hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_42hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_42hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_42hpi_3D7.res, file = "ASvAA_42hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_45hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_45hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_45hpi_3D7.res, file = "ASvAA_45hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_48hpi_3D7.res <- topTable(fit.cont, coef = "ASvAA_48hpi_3D7",sort.by = "p", n="Inf")
write.csv(ASvAA_48hpi_3D7.res, file = "ASvAA_48hpi_3D7_limma_results.csv",row.names = FALSE)
ASvAA_3hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_3hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_3hpi_FUP.res, file = "ASvAA_3hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_6hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_6hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_6hpi_FUP.res, file = "ASvAA_6hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_9hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_9hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_9hpi_FUP.res, file = "ASvAA_9hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_12hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_12hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_12hpi_FUP.res, file = "ASvAA_12hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_15hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_15hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_15hpi_FUP.res, file = "ASvAA_15hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_18hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_18hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_18hpi_FUP.res, file = "ASvAA_18hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_21hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_21hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_21hpi_FUP.res, file = "ASvAA_21hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_24hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_24hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_24hpi_FUP.res, file = "ASvAA_24hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_27hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_27hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_27hpi_FUP.res, file = "ASvAA_27hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_30hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_30hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_30hpi_FUP.res, file = "ASvAA_30hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_33hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_33hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_33hpi_FUP.res, file = "ASvAA_33hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_36hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_36hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_36hpi_FUP.res, file = "ASvAA_36hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_39hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_39hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_39hpi_FUP.res, file = "ASvAA_39hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_42hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_42hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_42hpi_FUP.res, file = "ASvAA_42hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_45hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_45hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_45hpi_FUP.res, file = "ASvAA_45hpi_FUP_limma_results.csv",row.names = FALSE)
ASvAA_48hpi_FUP.res <- topTable(fit.cont, coef = "ASvAA_48hpi_FUP",sort.by = "p", n="Inf")
write.csv(ASvAA_48hpi_FUP.res, file = "ASvAA_48hpi_FUP_limma_results.csv",row.names = FALSE)
glXYPlot(x = fit.cont$coefficients[,1], y = fit.cont$lods[,1], xlab = "logFC", ylab = "B", main = "ASvAA_3hpi_3D7", counts = y$counts, groups = group, status = summa.fit[,1], anno = fit.cont$genes, side.main = "ORF", folder = "output/xy/ASvAA_3hpi_3D7")
glXYPlot(x = fit.cont$coefficients[,2], y = fit.cont$lods[,2], xlab = "logFC", ylab = "B", main = "ASvAA_3hpi_FUP", counts = y$counts, groups = group, status = summa.fit[,2], anno = fit.cont$genes, side.main = "ORF", folder = "output/xy/ASvAA_3hpi_FUP")
glXYPlot(x = fit.cont$coefficients[,9], y = fit.cont$lods[,9], xlab = "logFC", ylab = "B", main = "ASvAA_6hpi_3D7", counts = y$counts, groups = group, status = summa.fit[,9], anno = fit.cont$genes, side.main = "ORF", folder = "output/xy/ASvAA_6hpi_3D7")
glXYPlot(x = fit.cont$coefficients[,10], y = fit.cont$lods[,10], xlab = "logFC", ylab = "B", main = "ASvAA_6hpi_FUP", counts = y$counts, groups = group, status = summa.fit[,10], anno = fit.cont$genes, side.main = "ORF", folder = "output/xy/ASvAA_6hpi_FUP")


su_rpkm <- read.csv("/Users/joe/R Projects/maliSalmon/data/rpkm_7SexualAndAsexualLifeStages_suetal.csv", header = TRUE, sep = ",")
our_log_rpkm = rpkm(txi_limma$counts, txi_limma$length)
our_log_rpkm <- log2(1+our_log_rpkm)
su_rpkm$Ookinete <- NULL
rownames(su_rpkm) <- su_rpkm$ID
su_rpkm$ID <- NULL
su_log_rpkm <- log2(1 + su_rpkm)

findMix <- function(Y, X) {
  X[is.na(X)] <- t(replicate(ncol(X), apply(X, 1, mean, na.rm = T)))[is.na(X)]
  Rinv <- solve(chol(t(X) %*% X))
  C <- cbind(rep(1, ncol(X)), diag(ncol(X)))
  b <- c(1, rep(0, ncol(X)))
  d <- t(Y) %*% X
  qp <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  sol <- qp$solution
  sol[sol < 1e-10] <- 0
  return(sol)
}

inter <- intersect(rownames(our_log_rpkm), rownames(su_log_rpkm))
O <- our_log_rpkm[rownames(our_log_rpkm) %in% inter, ]
O <- O[order(rownames(O)), ]
S <- su_log_rpkm[rownames(su_log_rpkm) %in% inter, ]
S <- S[order(rownames(S)), ]
ourPlotData <- data.frame()
for (i in 1:ncol(O)) { 
  mix <- findMix(O[, i], as.matrix(S))
  ourPlotData <- rbind(ourPlotData, data.frame(sample = rep(colnames(O)[i], ncol(S)), 
                                               stage = colnames(S), proportion = mix)) }

# Organise the results.
ourPlotData$stage <- gsub("Gametocyte.*", "Gametocyte", ourPlotData$stage)
ourPlotData <- aggregate(proportion ~ sample + stage, data = ourPlotData, FUN = sum)
ourPlotData <- within(ourPlotData, stage <- factor(stage, levels = c("Ring", "Early.Trophozoite", "Late.Trophozoite", "Schizont", "Gametocyte")))
ourPlotData$phenotype <- ifelse(substring(ourPlotData$sample, 1, 2) == "AS", 
                                "Sickle-trait", "Normal")

gg <- ggplot(ourPlotData, aes(x = factor(sample), y = proportion, fill = factor(phenotype))) +
  geom_bar(stat = "identity")
gg <- gg + scale_fill_manual(values = c('Normal' = "#054fb4", 'Sickle-trait' = "#d3160a"))
gg <- gg + facet_wrap(~stage, ncol = 1)
gg <- gg + theme_bw()
gg <- gg + theme(axis.text.x = element_text(size = 12, angle = 90), axis.text.y = element_text(size = 12,
      angle = 0), axis.title = element_text(size = 14, face = "bold"), strip.text.x = element_text(size = 16,
                                                                                                                                                                                            face = "bold"))
gg <- gg + labs(x = "Sample", y = "Proportion", fill = "Stage")
gg <- gg + theme(legend.text = element_text(size = 14))
gg <- gg + theme(legend.key.size = unit(0.25, "in"))
gg <- gg + theme(legend.title = element_text(size = 16, face = "bold"))
gg <- gg + guides(fill = guide_legend(title = "Phenotype"))
gg