## -----------------------------------------------------------------------------
library(tidyverse)
library(conflicted)
library(GenomicFeatures)
library(DESeq2)
library(ggrepel)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)


## -----------------------------------------------------------------------------
df_transcripts <- read_tsv("index/GRCH38.p14_biomart.tsv")
tx2gene <- df_transcripts %>% 
  rename("TXNAME"="Transcript stable ID version",
         "GENEID"="Gene name") %>%
  select(TXNAME,GENEID) %>% unique %>%
  drop_na()
samples <- list.files("aligned")
abundance_files <- setNames(file.path("aligned",samples,"abundance.tsv"),samples)
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
  read_tsv(file.path("aligned",smp,"abundance.tsv"),show_col_types = F) %>%
    mutate(Sample_ID=smp)
}) %>% bind_rows
sampleSheet <- sampleSheet[names(abundance_files),]
dds <- DESeqDataSetFromTximport(txi, df_samples, ~ Line + Treatment + Line:Treatment )%>%
    DESeq(.) %>% estimateSizeFactors(.)
df <- results(dds,tidy=T,alpha=0.05,contrast=list(c("Treatment_Hypoxia_vs_Normoxia")))
p <- ggplot(df,
       aes(log2FoldChange,-log10(padj)))+
  geom_point(alpha=0.3)+
  geom_text_repel(data=df%>%filter(row=="EGLN3"),aes(label=row))
p

