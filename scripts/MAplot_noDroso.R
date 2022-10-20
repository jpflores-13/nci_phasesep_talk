## Perform differential loops analysis with DESeq2

library(DESeq2)
library(ggplot2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(pheatmap)
library(apeglm)
library(RColorBrewer)
library(plotgardener)

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
plotMA(res)


# manipulate dataset for custom MA plot -----------------------------------

res_gg <- res |>
  data.frame(res) |>
  select(baseMean, log2FoldChange, padj) |> 
  mutate(isDE = case_when(
    log2FoldChange > 0 & padj < 0.1 ~ "TRUE - upreg",
    log2FoldChange < 0 & padj < 0.1 ~ "TRUE - downreg",
    padj > 0.1  ~ "FALSE",
    is.character("NA") ~ "FALSE"
  ))


# custom MA plot ----------------------------------------------------------

(ggplot <- res_gg |> 
   ggplot(aes(x = baseMean, y = log2FoldChange, color = isDE)) +
   geom_point(alpha= 1, size = 0.5) +
   geom_hline(yintercept = 0,
              lwd = 2,
              color = "grey40") +
   scale_color_manual(values = c("TRUE - upreg" = "#DC3220",
                                 "TRUE - downreg" = "#005AB5",
                                 "FALSE" = "grey60")) +
   ylim(c(-4,4)) +
   scale_x_log10(breaks=c(1,2,5,10,20,50,100)) +
   labs(y = "LFC",
        x = "mean of normalized counts") +
   theme_classic() +
   theme(legend.position = "NONE",
         axis.text.x=element_text(size=15),
         axis.text.y = element_text(size = 15),
         axis.title.x = element_text(size=15),
         axis.title.y = element_text(size =15)))

ggsave(plot = ggplot,
       "plots/MAplot_noDroso.pdf",
       width = 4,
       height = 3,
       units = "in",
       device = "pdf")
