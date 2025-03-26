library(tidyverse)
library(DreamAI)
library(PCAtools)
library(sva)
library(proDA)
library(uwot)
library(stringr)
library(msImpute)
library(MSnSet.utils)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(cowplot)
library(data.table)
#devtools::install_version("Matrix", version = "1.6.1.1")
library(Matrix)
library(edgeR)
library(glmGamPoi)
library(sctransform)
library(gprofiler2)
library(pheatmap)
library(EnhancedVolcano)
library(reshape2)
library(RColorBrewer)

getwd()
setwd("C:/Users/dawa726/OneDrive - PNNL/Desktop/nanoSPINS_Pranav")
load("../../../repos/nanoSPINS_data_scripts/nanoSPINS_proteome.RData")
#scRNA-seq data analysis

annota <- read.delim("masterkey.txt")

#annota %>% select(Cell.type, Condition, gene_count) %>% filter(!Cell.type == "NoCell") %>% #filter(Condition == "nanoPOTS") %>% group_by(Cell.type) %>%   ggplot() +  aes(x = reorder(Cell.type, gene_count), y = gene_count) +  geom_violin(aes(colour = Cell.type, fill = Cell.type), outlier.size = 0, width = 0.5) + stat_summary(fun.y="mean",color="black", geom = "point") + #geom_jitter(aes(colour = Cell.type), width = 0.25) + scale_y_continuous(limits = c(0,10000), breaks = seq(0,10000, by = 2000)) + ylab("Number of Genes Identified (n)") + xlab("cell_type") + ggtitle("scRNA-seq") + theme_minimal(base_size = 14)
#ggsave("number_of_genes_identified_scRNAseq.png", width = 5, height = 4, bg = "white")

annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Cell.type != "NoCell") %>%
  group_by(Condition, Cell.type) %>%
  ggplot() +
  aes(x = fct_reorder(Condition, gene_count), 
      y = gene_count) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(aes(colour = Cell.type), width = 0.25) +
  stat_summary(fun="mean",color= "NA", geom = "point") +
  geom_hline(yintercept = 3022.365, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 3040.979, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 6135.792, linetype = "dashed", color = "purple", linewidth = 1) +
  geom_hline(yintercept = 6867.556, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  #scale_fill_manual(values = c(C10 = "green3", SVEC = "cornflowerblue", NoCell = "red")) +
  scale_y_continuous(limits = c(0,16000), 
                     breaks = seq(0,16000, by = 2000)) +
  ylab("Number of Genes Identified (n)") + 
  xlab("Condition") + 
  theme_minimal(base_size = 18)
#ggsave("number_of_genes_identified_in_each_condition_scRNAseq.png", width = 5, height = 4, bg = "white")

mean_genes_identified <- annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Cell.type != "NoCell") %>%
  #filter(Condition == "nanoPOTS") %>%
  group_by(Condition, Cell.type) %>%
  summarise(mean(gene_count))

annota %>% 
  select(Cell.Barcode, Cell.type, Condition,Total, 
         UMI_count, gene_count, Ratio_mapped) %>% 
  filter(Cell.type != "NoCell") %>% 
  ggplot() + aes(x = UMI_count, y = gene_count, 
                 colour = Cell.type) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 200000), breaks = seq(0, 200000, by = 50000)) +
  scale_y_continuous(limits = c(0, 15000), breaks = seq(0, 15000, by = 5000)) + theme_minimal(base_size = 14)
#ggsave("C:/Users/dawa726/UMI_count&gene_count.png", width = 6, height = 4, bg = "white")

annota %>% 
  select(Cell.Barcode, Cell.type, Condition,Total, 
         UMI_count, gene_count, Ratio_mapped) %>% 
  filter(Cell.type != "NoCell") %>% 
  ggplot() + aes(x = UMI_count, y = Ratio_mapped, 
                 colour = Cell.type) + 
  geom_point() + 
  scale_x_continuous(limits = c(0, 200000), breaks = seq(0, 200000, by = 50000)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25)) + theme_minimal(base_size = 14)
#ggsave("C:/Users/dawa726/UMI_count&ratio_mapped.png", width = 6, height = 4, bg = "white")

convert <- read.delim("convert_names.txt")
protein_samples <- colnames(protein_input_impute)
RNA_BCs_test <- annota %>%
  filter(Annotation %in% protein_samples) %>% 
  filter(gene_count >= 1000) %>% 
  filter(Ratio_mapped >= 50) %>%
  arrange(Cell.Barcode)

RNA_selected_samples_new_test <- annota %>%
  filter(Cell.type != "NoCell") %>%
  filter(gene_count >= 1000) %>%
  filter(Ratio_mapped >= 50) %>%
  filter(!Annotation == "C10_2_8_4")
RNA_selected_samples_new_test_sorted <- RNA_selected_samples_new_test %>%
  arrange(Cell.Barcode)

RNA_selected_samples_new_test_sorted %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Cell.type != "NoCell") %>%
  filter(Condition == "nanoPOTS") %>%
  group_by(Condition, Cell.type) %>%
  summarise(mean(gene_count))

RNA_selected_samples_new_test_sorted %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(!Cell.type == "NoCell") %>%
  filter(Condition == "nanoPOTS") %>%
  group_by(Cell.type) %>% 
  ggplot() +
  aes(x = reorder(Cell.type, -gene_count), 
      y = gene_count) +
  geom_violin(aes(colour = Cell.type, fill = Cell.type), width = 0.5) + 
  stat_summary(fun.y="mean",color="black", geom = "point") +
  #geom_jitter(aes(colour = Cell.type), width = 0.25) +
  scale_y_continuous(limits = c(0,15000), 
                     breaks = seq(0,15000, by = 3000)) +
  ylab("Number of Genes Identified (n)") + 
  xlab("cell_type") + 
  ggtitle("scRNA-seq") + facet_wrap(~ Condition) +
  theme_minimal(base_size = 18)
#ggsave("C:/Users/dawa726/number_of_genes_identified_scRNAseq_new.png", width = 5, height = 4, bg = "white")

RNA_selected_samples_new_test_sorted %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Cell.type != "NoCell") %>%
  #filter(Condition == "nanoPOTS") %>%
  group_by(Condition, Cell.type) %>%
  summarise(mean(gene_count)) %>% 
  view()

RNA_input_all_new <- read.delim("./SB22_12_PNNL_21.umicount.inex.all.tsv", 
                        row.names = "Gene") %>% dplyr::select(one_of(RNA_selected_samples_new_test_sorted$Cell.Barcode))

RNA_input_all_new2 <- RNA_input_all_new %>% 
  rownames_to_column(var = "gene_id") %>%
  inner_join(., convert) %>% 
  mutate(Avg = rowSums(.[2:278]/277)) %>% 
  group_by(gene_name) %>% 
  slice_max(Avg, n = 1, with_ties = F ) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select( -gene_id, -Avg, -gene_id.1)

RNA_input_all_new2.2 <- RNA_input_all_new2 %>%
  select(sort(names(RNA_input_all_new2)))

all(colnames(RNA_input_all_new2.2) == RNA_selected_samples_new_test_sorted$Cell.Barcode)

colnames(RNA_input_all_new2.2) <- RNA_selected_samples_new_test_sorted$Annotation

RNA_input_all_new2.2.1 <- RNA_input_all_new2.2 %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS_C10",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS_SVEC", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  filter(!reads == 0) %>%
  group_by(Gene, Group) %>% 
  add_count(name = "n") %>% 
  select(Gene, SampleID, Group, reads, method, n) %>%
  ungroup() %>%
  group_by(Gene) %>% 
  add_count(name = "n2") %>% 
  mutate(n3 = n2 - n) %>% 
  mutate(Impute = case_when(reads > 1 & n >= 5 & Group == "SVEC" & n3 >= 5 ~ "Keep", 
                            reads > 1 & n >= 5 & Group == "C10" & n3 >= 5 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>% 
  ungroup() %>%
  select(Gene, SampleID, Group, method, reads) %>%
  ungroup() %>%
  select(Gene, SampleID, reads) %>%
  pivot_wider(names_from = SampleID, values_from = reads) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  column_to_rownames("Gene")

dge_object_1 <- DGEList(counts = RNA_input_all_new2.2.1, 
                        remove.zeros = TRUE)
dge_object_1 <- calcNormFactors(dge_object_1)
normalized_counts_1 <- cpm(dge_object_1, 
                           normalized.lib.sizes = TRUE)
df_normalized_counts_1 <- as.data.frame(normalized_counts_1)
df_normalized_counts_1 <- replace(df_normalized_counts_1, df_normalized_counts_1 == 0, NA)
df_normalized_counts_1 <- log2(df_normalized_counts_1)

length(which(grepl("SVEC_[A-Z]", colnames(df_normalized_counts_1))))
length(which(grepl("C10_[A-Z]", colnames(df_normalized_counts_1))))

median_CV_RNA_1 <- df_normalized_counts_1 %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS_C10",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS_SVEC", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  group_by(Gene, method, Group) %>%
  filter(!is.na(normalized_reads)) %>%
  mutate(raw_reads = 2^normalized_reads) %>%
  mutate(Avg = mean(raw_reads)) %>%
  mutate(CV = (sd(raw_reads)/(Avg))) %>%
  distinct(Gene, CV, method, Group) %>%
  ungroup() %>%
  select(Group, method, CV) %>%
  group_by(Group, method) %>%
  filter(!is.na(CV)) %>%
  summarise(median(CV))

df_normalized_counts_1_CVs <- df_normalized_counts_1 %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted")) %>%
  group_by(Gene, method, Group) %>%
  filter(!is.na(normalized_reads)) %>%
  mutate(raw_reads = 2^normalized_reads) %>%
  mutate(Avg = mean(raw_reads)) %>%
  mutate(CV = (sd(raw_reads)/(Avg))) %>%
  distinct(Gene, CV, method, Group) %>%
  filter(!is.na(CV))

df_normalized_counts_1_CVs %>%
  filter(method == "directsorted") %>%
  ggplot()+
  aes(x = Group, y = CV, fill = Group)+
  geom_violin(drop = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               color = "black")+
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5, by = 0.5))+
  ylab("Coefficient of Variation")+
  xlab("")+
  #facet_wrap(~ method) +
  theme_minimal(base_size = 18)
#ggsave("C:/Users/dawa726/median_CVs_1_DS.png", width = 6, height = 4, bg = "white")

df_normalized_reads_nanoSPINS_1 <- df_normalized_counts_1 %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10",
                           grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC",
                           grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                           grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS",
                            str_detect(Group, "directsorted") ~ "directsorting")) %>%
  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10",
                              str_detect(SampleID, "SVEC") ~ "SVEC")) %>% 
  filter(method == "nanoPOTS") %>%
  #filter(SampleID %in% protein_samples) %>%
  select(gene, SampleID, normalized_reads) %>%
  pivot_wider(names_from = SampleID, values_from = normalized_reads) %>%
  column_to_rownames(var= "gene")

sorted_colnames_RNA_1_NS <- sort(colnames(df_normalized_reads_nanoSPINS_1))
df_normalized_reads_nanoSPINS_1 <- df_normalized_reads_nanoSPINS_1[, sorted_colnames_RNA_1_NS]
df_normalized_reads_nanoSPINS_1 <- replace(df_normalized_reads_nanoSPINS_1, is.na(df_normalized_reads_nanoSPINS_1), 0)
RNA_BCs_NS <- annota %>%
  #filter(Annotation %in% protein_samples) %>% 
  filter(gene_count >= 1000) %>% 
  filter(Ratio_mapped >= 50) %>%
  filter(Type == "nanoPOTS") %>%
  arrange(Cell.Barcode)
RNA_BCs_NS <- RNA_BCs_NS[order(RNA_BCs_NS$Annotation),]
pca_RNA_1_NS <- PCAtools::pca(df_normalized_reads_nanoSPINS_1, 
                           scale = F, 
                           center = T)

PCAtools::biplot(pca_RNA_1_NS, x = "PC2", y =  "PC4", 
                 lab = RNA_BCs_NS$Cell.type,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, 
                 ellipseLevel = 0.90, 
                 ylim = c(-200, 200),  
                 xlim = c(-200, 200)
                 )
#ggsave("C:/Users/dawa726/PCA_scRNAseq_1_PC2_PC4.png", width = 10, height = 6, bg = "white")

RNA_BCs_DS <- annota %>%
  #filter(Annotation %in% protein_samples) %>% 
  filter(gene_count >= 1000) %>% 
  filter(Ratio_mapped >= 50) %>%
  filter(Type == "Direct") %>%
  filter(!Annotation == "C10_2_8_4") %>%
  arrange(Cell.Barcode)
RNA_BCs_DS <- RNA_BCs_DS[order(RNA_BCs_DS$Annotation),]

df_normalized_reads_directsorted_1 <- df_normalized_counts_1 %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10",
                           grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC",
                           grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                           grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS",
                            str_detect(Group, "directsorted") ~ "directsorting")) %>%
  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10",
                              str_detect(SampleID, "SVEC") ~ "SVEC")) %>% 
  filter(method == "directsorting") %>%
  filter(SampleID %in% RNA_BCs_DS$Annotation) %>%
  select(gene, SampleID, normalized_reads) %>%
  pivot_wider(names_from = SampleID, values_from = normalized_reads) %>%
  column_to_rownames(var= "gene")

sorted_colnames_RNA_DS_1 <- sort(colnames(df_normalized_reads_directsorted_1))
df_normalized_reads_directsorted_1 <- df_normalized_reads_directsorted_1[, sorted_colnames_RNA_DS_1]
df_normalized_reads_directsorted_1 <- replace(df_normalized_reads_directsorted_1, is.na(df_normalized_reads_directsorted_1), 0)

pca_RNA_DS_1 <- PCAtools::pca(df_normalized_reads_directsorted_1, 
                              scale = F, 
                              center = T)

PCAtools::biplot(pca_RNA_DS_1, x = "PC2", y =  "PC3", 
                 lab = RNA_BCs_DS$Cell.type,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, 
                 ellipseLevel = 0.90, 
                 ylim = c(-250, 250),  
                 xlim = c(-250, 250)
)
#ggsave("C:/Users/dawa726/PCA_scRNAseq_1_DS_PC2_PC3.png", width = 10, height = 6, bg = "white")

unique_genes_scRNAseq_NS_1 <- df_normalized_reads_nanoSPINS_1 %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10",
                           grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC")) %>%
  mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS")) %>%
  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10",
                              str_detect(SampleID, "SVEC") ~ "SVEC")) %>%
  filter(normalized_reads != 0) %>%
  select(gene, celltype) %>% 
  distinct()
#write.csv(unique_genes_scRNAseq_NS_1, file = "C:/Users/dawa726/unique_genes_scRNAseq_1.csv")

unique_genes_scRNAseq_DS_1 <- df_normalized_counts_1 %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                           grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10",
                              str_detect(SampleID, "SVEC") ~ "SVEC")) %>%
  filter(normalized_reads != 0) %>%
  select(gene, celltype) %>% 
  distinct()
#write.csv(unique_genes_scRNAseq_DS_1, file = "C:/Users/dawa726/unique_genes_scRNAseq_DS1.csv")

correlation_matrix <- cor(as.data.frame(lapply(df_normalized_counts_1, as.numeric)), method = "pearson", use = "pairwise.complete.obs")
colnames(correlation_matrix) <- ifelse(
  grepl("C10_[A-Z]", colnames(correlation_matrix)), 
  paste0("nanoSPINS_", colnames(correlation_matrix)), 
  ifelse(
    grepl("SVEC_[A-Z]", colnames(correlation_matrix)), 
    paste0("nanoSPINS_", colnames(correlation_matrix)),
    ifelse(
      grepl("SVEC_[1-9]", colnames(correlation_matrix)), 
      paste0("Direct_sorting_", colnames(correlation_matrix)), 
      ifelse(
        grepl("C10_[1-9]", colnames(correlation_matrix)), 
        paste0("Direct_sorting_", colnames(correlation_matrix)), 
        colnames(correlation_matrix)
      )
    )
  )
)

rownames(correlation_matrix) <- ifelse(
  grepl("C10_[A-Z]", rownames(correlation_matrix)), 
  paste0("nanoSPINS_", rownames(correlation_matrix)), 
  ifelse(
    grepl("SVEC_[A-Z]", rownames(correlation_matrix)), 
    paste0("nanoSPINS_", rownames(correlation_matrix)),
    ifelse(
      grepl("SVEC_[1-9]", rownames(correlation_matrix)), 
      paste0("Direct_sorting_", rownames(correlation_matrix)), 
      ifelse(
        grepl("C10_[1-9]", rownames(correlation_matrix)), 
        paste0("Direct_sorting_", rownames(correlation_matrix)), 
        rownames(correlation_matrix)
      )
    )
  )
)

hc <- hclust(dist(correlation_matrix))
reordered_correlation_matrix <- correlation_matrix[hc$order, hc$order]
png("correlation_plot_scRNAseq_all_031725.png", width = 3000, height = 3000, bg = "white", res = 300)
corrplot::corrplot(correlation_matrix, 
                   method = 'shade', 
                   diag = FALSE, 
                   tl.cex = 0.2, 
                   order = 'alphabet', 
                   col.lim = c(0, 1)) %>% 
  corrplot::corrRect(c(1,65,133,203,278), 
                     lwd = 2, col = "black")
dev.off()

nanoSPINS_C10_RNA_matrix <- correlation_matrix[grep("nanoSPINS_C10", rownames(correlation_matrix)), grep("nanoSPINS_C10", colnames(correlation_matrix))]
average_correlation_C10 <- mean(nanoSPINS_C10_RNA_matrix[upper.tri(nanoSPINS_C10_RNA_matrix) | lower.tri(nanoSPINS_C10_RNA_matrix)])
print(average_correlation_C10)

nanoSPINS_SVEC_RNA_matrix <- correlation_matrix[grep("nanoSPINS_SVEC", rownames(correlation_matrix)), grep("nanoSPINS_SVEC", colnames(correlation_matrix))]
average_correlation_SVEC <- mean(nanoSPINS_SVEC_RNA_matrix[upper.tri(nanoSPINS_SVEC_RNA_matrix) | lower.tri(nanoSPINS_SVEC_RNA_matrix)])
print(average_correlation_SVEC)

Direct_sorting_C10_RNA_matrix <- correlation_matrix[grep("Direct_sorting_C10", rownames(correlation_matrix)), grep("Direct_sorting_C10", colnames(correlation_matrix))]
average_correlation_C10_DS <- mean(Direct_sorting_C10_RNA_matrix[upper.tri(Direct_sorting_C10_RNA_matrix) | lower.tri(Direct_sorting_C10_RNA_matrix)])
print(average_correlation_C10_DS)

Direct_sorting_SVEC_RNA_matrix <- correlation_matrix[grep("Direct_sorting_SVEC", rownames(correlation_matrix)), grep("Direct_sorting_SVEC", colnames(correlation_matrix))]
average_correlation_SVEC_DS <- mean(Direct_sorting_SVEC_RNA_matrix[upper.tri(Direct_sorting_SVEC_RNA_matrix) | lower.tri(Direct_sorting_SVEC_RNA_matrix)])
print(average_correlation_SVEC_DS)

SVEC_RNA_matrix <- correlation_matrix[grep("SVEC", rownames(correlation_matrix)), grep("SVEC", colnames(correlation_matrix))]
average_correlation_SVEC <- mean(SVEC_RNA_matrix[upper.tri(SVEC_RNA_matrix) | lower.tri(SVEC_RNA_matrix)])
print(average_correlation_SVEC)

C10_RNA_matrix <- correlation_matrix[grep("C10", rownames(correlation_matrix)), grep("C10", colnames(correlation_matrix))]
average_correlation_C10 <- mean(C10_RNA_matrix[upper.tri(C10_RNA_matrix) | lower.tri(C10_RNA_matrix)])
print(average_correlation_C10)