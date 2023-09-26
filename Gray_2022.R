library(tidyverse)
library(conflicted)
library(GenomicFeatures)
library(DESeq2)
library(ggrepel)
library(tools)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
df_transcripts <- read_tsv("GRCH38.p14_biomart.tsv")
tx2gene <- df_transcripts %>% 
  rename("TXNAME"="Transcript stable ID version",
         "GENEID"="Gene name") %>%
  select(TXNAME,GENEID) %>% unique %>%
  drop_na()
samples <- list.files("tx")
abundance_files <- setNames(file.path("tx",samples,"abundance.tsv"),samples)
txi <- tximport::tximport(abundance_files,type="kallisto",tx2gene=tx2gene)
df_samples <- tibble("Sample_ID"=samples) %>%
  mutate(Line=if_else(str_detect(Sample_ID,"WT"),"WT","KO")) %>%
  mutate(Rep=str_replace(Sample_ID,"\\..*","")) %>%
  mutate(Treatment=if_else(str_detect(Sample_ID,".N"),"Normoxia","Hypoxia")) %>%
  mutate(Line=factor(Line,levels=c("WT","KO")),
         Treatment=factor(Treatment,levels=c("Normoxia","Hypoxia"))) %>%
  write_tsv("Gray_2022_sample_info.tsv")
sampleSheet <- df_samples  %>%
  column_to_rownames(var="Sample_ID")
df_kallisto <- samples %>% map(function(smp){
  read_tsv(file.path("tx",smp,"abundance.tsv"),show_col_types = F) %>%
    mutate(Sample_ID=smp)
}) %>% bind_rows
sampleSheet <- sampleSheet[names(abundance_files),]
counts <- round(txi$counts)
mode(counts) <- "integer"
colData <- df_samples
class(colData) <- "data.frame"
# dds <- DESeqDataSetFromTximport(txi,df_samples,
#                                 ~Line + Treatment + Line:Treatment ) %>%
#        DESeq(.) %>%
#        estimateSizeFactors(.)
samples <- colData
save(list = c("counts","samples"),
     file = "hypoxia.RData")
resaveRdaFiles("hypoxia.RData")
wt <- which(colData$Line == "WT")
ko <- which(colData$Line == "KO")
colData_wt <- colData[wt,]
colData_ko <- colData[ko,]
counts_wt  <- counts[,wt]
counts_ko  <- counts[,ko]
dat_wt <- DESeqDataSetFromMatrix(counts_wt,colData_wt,~Treatment)
dat_ko <- DESeqDataSetFromMatrix(counts_ko,colData_ko,~Treatment)
deseq_wt <- DESeq(dat_wt)
deseq_ko <- DESeq(dat_ko)
res_wt <- results(deseq_wt,alpha = 0.05,
                  contrast = list("Treatment_Hypoxia_vs_Normoxia"))
res_ko <- results(deseq_ko,alpha = 0.05,
                  contrast = list("Treatment_Hypoxia_vs_Normoxia"))
#res_wt<-lfcShrink(deseq_wt,coef="Treatment_Hypoxia_vs_Normoxia",type="ashr")
#res_ko<-lfcShrink(deseq_ko,coef="Treatment_Hypoxia_vs_Normoxia",type="ashr")
pdat <- data.frame(wt = res_wt$log2FoldChange,
                   ko = res_ko$log2FoldChange)
library(cowplot)
rows <- with(pdat,which(abs(wt) > 1 | abs(ko) > 1))
plot(res_wt$log2FoldChange,res_ko$log2FoldChange,pch = 20)
coef(lm(res_ko$log2FoldChange[rows] ~ res_wt$log2FoldChange[rows]))
abline(a = -0.1326,b = 0.4986,col = "magenta",lty = "dashed")
subset(res_wt,log2FoldChange > 10)
