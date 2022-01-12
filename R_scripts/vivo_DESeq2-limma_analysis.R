library("DESeq2")
library("GenomicFeatures")
data_dir <- "/Users/joe/R Projects/maliSalmon/data"
csvfile <- file.path(data_dir, "sampleTable_update.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
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
dds$group <- factor(paste0(dds$Genotype, dds$Stage))
design(dds) <- ~ group
design(dds)
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
  coord_fixed() +
  scale_fill_manual(values=c("red2", "#0080FF"))

pcaData <- plotPCA(vsd, intgroup=c("Stage", "Genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Stage, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_color_manual(values=c("red2", "#0080FF"))


resultsNames(dds)
(ASvAA_Ring <- results(dds, contrast = c("group", "ASRing", "AARing"), alpha = 0.1))
(ASvAA_Troph <- results(dds, contrast = c("group", "ASTrophozoite", "AATrophozoite"), alpha = 0.1))
library(org.Pf.plasmo.db)
keysdds <- rownames(ASvAA_Ring)
annosdds <- select(org.Pf.plasmo.db, keys = keysdds, columns =c("GENENAME"), keytype = "ORF")
ASvAA_Ring$genes <- annosdds
ASvAA_Troph$genes <- annosdds
write.csv(ASvAA_Ring, file = "ASvAA_Ring_mali_DESeq2_results.csv")
write.csv(ASvAA_Troph, file = "ASvAA_Troph_mali_DESeq2_results.csv")

library("limma")
library("edgeR")
library("Glimma")

txi_limma <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
colnames(txi_limma$counts) <- rownames(sampleTable)
colnames(txi_limma$abundance) <- rownames(sampleTable)
colnames(txi_limma$length) <- rownames(sampleTable)
y <- DGEList(txi_limma$counts)
keep <- filterByExpr(y)
y <- y[keep, ]
y <- calcNormFactors(y)

group <- paste(sampleTable$Genotype, sampleTable$Stage, sep = ".")
labels <- paste(sampleTable$sampleName, sampleTable$Genotype, sampleTable$Stage)
group <- factor(group)
glMDSPlot(y, labels=labels, groups=group, folder="Glimma_plots/mds")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

v <- voom(y,design,plot = TRUE)
fit <- lmFit(v)

cont.matrix <- makeContrasts(
ASvAA_rings = AS.Ring - AA.Ring, #1
ASvAA_trophs = AS.Trophozoite - AA.Trophozoite, #2
ASvAA_mix = AS.Mixed - AA.Mixed, #3
levels = design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)

library(ggplot2)

top <- topTable(fit.cont, coef = 2, p.value = 0.1, sort.by = "p", number = Inf, adjust.method = "BH", confint = TRUE)
hist <- data.frame(topTable(fit.cont, coef = 2, number = Inf))
gg <- ggplot(hist, aes(x = P.Value)) + geom_histogram(binwidth = 0.01)
gg <- gg + theme_bw()
gg <- gg + theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12, angle = 0), axis.title = element_text(size = 12, face = "bold"))
gg <- gg + labs(x = "P value", y = "Count")
gg


keys <- rownames(fit.cont)
annos <- select(org.Pf.plasmo.db, keys = keys, columns =c("GENENAME"), keytype = "ORF")
fit.cont$genes <- annos
ASvAA_rings.res <- topTable(fit.cont, coef = "ASvAA_rings",sort.by = "p", n="Inf")
write.csv(ASvAA_rings.res, file = "ASvAA_rings_limma_results.csv",row.names = FALSE)
ASvAA_trophs.res <- topTable(fit.cont, coef = "ASvAA_trophs",sort.by = "p", n="Inf")
write.csv(ASvAA_trophs.res, file = "ASvAA_trophs_limma_results.csv",row.names = FALSE)

glXYPlot(x = fit.cont$coefficients[,1], y = fit.cont$lods[,1], xlab = "logFC", ylab = "B", main = "ASvAA_Rings", counts = y$counts, groups = group, status = summa.fit[,1], anno = fit.cont$genes, side.main = "ORF", folder = "Glimma_plots/xy/ASvAA_Rings")
glXYPlot(x = fit.cont$coefficients[,2], y = fit.cont$lods[,2], xlab = "logFC", ylab = "B", main = "ASvAA_Trophs", counts = y$counts, groups = group, status = summa.fit[,2], anno = fit.cont$genes, side.main = "ORF", folder = "Glimma_plots/xy/ASvAA_Trophs")


library(reshape2)


bplot <- melt(colSums(y$counts))
colnames(bplot) <- c("Read.Counts")
bplot$Sample <- rownames(bplot)
gg <- ggplot(bplot, aes(x = factor(Sample), y = Read.Counts)) + geom_bar(stat = "identity")
gg <- gg + theme_bw() + scale_y_sqrt(breaks = c(0, 10000, 1e+06, 1e+07, 1e+08))
gg <- gg + theme(axis.text.x = element_text(size = 12, angle = 90), axis.text.y = element_text(size = 12,
angle = 0), axis.title = element_text(size = 14, face = "bold"))
gg <- gg + labs(x = "Sample", y = "Read Count")
gg <- gg + geom_hline(aes(yintercept = 1e+06), col = "red")
gg

library(quadprog)
library(ruv)

colors <- c("#054fb4", "#d3160a")

su_rpkm <- read.csv("/Users/joe/R Projects/maliSalmon/data/rpkm_7SexualAndAsexualLifeStages_suetal.csv", header = TRUE, sep = ",")
our_log_rpkm = rpkm(txi$counts, txi_limma$length)
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
ourPlotData <- within(ourPlotData, stage <- factor(stage, levels = c("Ring", 
                "Early.Trophozoite", "Late.Trophozoite", "Schizont", "Gametocyte")))
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

plotVoomRLE <- function(E, colours) {
mn <- apply(E, 1, median)
rle <- data.frame(sweep(E, MARGIN = 1, STATS = mn, FUN = "-"))
boxplot(rle, col = colours, outline = FALSE, las = 2, ylim = c(-7, 7))
abline(h = 0, col = "black")
}
plotVoomRLE(v$E, colors[categories])

library(EnhancedVolcano)

EnhancedVolcano(ASvAA_Ring, lab = rownames(ASvAA_Ring), x = 'log2FoldChange', y = 'pvalue', FCcutoff = log2(1.5), pointSize = 3.5, labSize = 3.75)

EnhancedVolcano(ASvAA_Troph, lab = rownames(ASvAA_Troph), x = 'log2FoldChange', y = 'pvalue', FCcutoff = log2(1.5), pCutoff = 0.001, xlim= c(-6,6), pointSize = 3.5, labSize = 3.75)


