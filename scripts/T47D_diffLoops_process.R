## Perform differential loops analysis with DESeq2 on T47D -isDroso false loops

library(DESeq2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(pheatmap)
library(apeglm)

# Load extracted loops counts ---------------------------------------------

loopCounts <- readRDS("data/raw/hic/loops/YAPP_hic_diff_loopCounts_noDroso.rds")

## add a loop_count column
loopCounts$loop_name <- glue("loop_{1:length(loopCounts)}")

# Create a matrix for countData -------------------------------------------

m <- mcols(loopCounts)[, grep("*inter.*", colnames(mcols(loopCounts)))] %>% 
  as.matrix()

# Construct colData/metadata ----------------------------------------------

## String split the colnames 

colData <- as.data.frame(do.call(rbind, strsplit(colnames(m), "_")), stringsAsFactors = T)
rownames(colData) <- colnames(m)
colnames(colData) <- c("Project", "Cell_Type", "Treatment", "Replicate")
colData <- colData[, c(1:4)]

## Correct replicate number
colData <- colData |> 
  mutate(Replicate = str_replace(Replicate, "4", "1")) |> 
  mutate(Replicate = str_replace(Replicate, "5", "2")) |> 
  mutate(Replicate = str_replace(Replicate, "6", "3"))

## Make sure sample names match
all(colnames(m) == rownames(colData))

# Run DESeq2 --------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = m,
                              colData = colData,
                              design = ~ Replicate + Treatment)

## disable DESeq's default normalization 
sizeFactors(dds) <- rep(1, ncol(dds))

## Hypothesis testing with Wald with `betaPrior = F`
dds <- DESeq(dds)

## outputs
res <- results(dds)
resultsNames(dds)
summary(res)
res <- lfcShrink(dds, coef="Treatment_sorbitol_vs_control", type= "apeglm")

summary(res)
pdf(file = "plots/YAPP_diffLoops_sorb_noDroso_MA_hic.pdf")
plotMA(res, ylim=c(-4,4), main = "Differential Loop Analysis",
       ylab = "LFC",
       xlab = "mean of norm. counts")
dev.off()

# Concatenate loopCounts and DESeqResults ---------------------------------
mcols(loopCounts) <- cbind(mcols(loopCounts), res)
diff_loopCounts <- loopCounts

## save as .rds
saveRDS(diff_loopCounts, file = "data/raw/hic/loops/YAPP_hic_diff_loopCounts_noDroso.rds")
sessionInfo()
