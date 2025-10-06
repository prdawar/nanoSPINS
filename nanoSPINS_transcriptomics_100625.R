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
library(Matrix)
library(edgeR)
library(glmGamPoi)
library(sctransform)
library(gprofiler2)
library(pheatmap)
library(EnhancedVolcano)
library(reshape2)
library(RColorBrewer)
library(VennDiagram)

getwd()
setwd("C:/Users/dawa726/OneDrive - PNNL/Desktop/nanoSPINS_Pranav")
annota <- read.delim("masterkey.txt")

#For supplementary Figure 5A
annota %>%
  select(Cell.type, Condition, gene_count, Ratio_mapped) %>%
  mutate(color_group = case_when(
    Condition == "nanoPOTS" & 
      gene_count >= 1500 & Ratio_mapped >= 50 ~ "nanoPOTS_QC_pass",
    Condition == "Direct_Sorting" & 
      gene_count >= 1500 & Ratio_mapped >= 50 ~ "Direct_Sorting_QC_pass",
    Cell.type == "NoCell" ~ "Blank",
    TRUE ~ "QC_fail"
  )) %>%
  filter(str_detect(color_group, "_QC_pass$")) %>%
  group_by(color_group) %>%
  summarise(mean(gene_count))

annota %>%
  select(Cell.type, Condition, gene_count, Ratio_mapped) %>%
  mutate(color_group = case_when(
    Condition == "nanoPOTS" & 
      gene_count >= 1500 & Ratio_mapped >= 40 ~ "nanoPOTS_QC_pass",
    Condition == "Direct_Sorting" & 
      gene_count >= 1500 & Ratio_mapped >= 50 ~ "Direct_Sorting_QC_pass",
    Cell.type == "NoCell" ~ "Blank",
    TRUE ~ "QC_fail"
  )) %>% 
  ggplot(aes(x = fct_reorder(Condition, gene_count),
             y = gene_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(colour = color_group), width = 0.25, alpha = 0.75, size = 2) +
  geom_hline(yintercept = 3646, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 7003, linetype = "dashed", color = "purple", linewidth = 1) +
  scale_colour_manual(name = "QC Status",
    values = c(
      "nanoPOTS_QC_pass" = "#F8766D",
      "Direct_Sorting_QC_pass" = "#00BFC4",
      "Blank" = "#808000",
      "QC_fail" = "#708090")) +
  scale_y_continuous(limits = c(0, 13500),
    breaks = seq(0, 13500, by = 1500)) +
  ylab("Number of Genes Identified (n)") +
  xlab(NULL) +
  theme_minimal(base_size = 18)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
#ggsave("gene_count_072725.png", width = 8, height = 6, bg = "white")

convert <- read.delim("convert_names.txt")
#protein_input_impute <- read.csv(file = "protein_input_impute.csv")
protein_samples <- colnames(protein_input_impute)
RNA_NS_QCpass <- annota %>%
  filter(str_detect(Annotation, "SVEC_[A-C]|C10_[A-C]|Control_[A-C]")) %>%
  filter(gene_count >= 1500) %>% 
  filter(Ratio_mapped >= 50) %>% 
  filter(Annotation != "Control_C5") %>%
  arrange(Cell.Barcode)

RNA_NS_QCpass %>% 
  select(Cell.type, UMI_count) %>% 
  group_by(Cell.type) %>% 
  summarise(mean = mean(UMI_count))

RNA_NS_QCpass %>% 
  select(Cell.type, gene_count) %>% 
  group_by(Cell.type) %>% 
  summarise(mean = mean(gene_count))

#Figure 4D
RNA_NS_QCpass %>% 
  select(Cell.type, Condition, gene_count) %>%
  mutate(Cell.type = factor(Cell.type, levels = c("C10", "SVEC"))) %>%
  ggplot() +
  aes(x = Cell.type, 
      y = gene_count) +
  geom_boxplot(aes(colour = Cell.type), 
               outlier.size = -1, 
               width = 0.45, na.rm = TRUE) +
  geom_jitter(aes(colour = Cell.type), 
              width = 0.2, na.rm = TRUE, 
              size = 3, alpha = 0.75) +
  stat_summary(fun = mean,
               geom = "point", 
               color = "black", 
               size = 2) +
  scale_y_continuous(limits = c(0,8000), 
                     breaks = seq(0,8000, by = 1000)) +
  scale_colour_manual(values = c("C10" = "#F8766D", 
                                 "SVEC" = "#00BFC4"))+
  ylab("Number of Genes Identified (n)")+
  xlab(NULL)+ 
  ggtitle(NULL)+
  theme_minimal(base_size = 18)
#ggsave("gene_IDs_QCpass_NS_072725.png", width = 6, height = 4, bg = "white")

RNA_all_QC_pass <- annota %>%
  filter(Cell.type != "NoCell") %>%
  filter(gene_count >= 1500) %>%
  filter(Ratio_mapped >= 50) %>%
  filter(!Annotation == "C10_2_8_4") %>%
  arrange(Cell.Barcode)

RNA_input_all_new <- read.delim("./SB22_12_PNNL_21.umicount.inex.all.tsv", 
                                row.names = "Gene") %>% 
  dplyr::select(one_of(RNA_all_QC_pass$Cell.Barcode))

RNA_input_all_new <- RNA_input_all_new %>% 
  rownames_to_column(var = "gene_id") %>%
  inner_join(., convert) %>% 
  mutate(Avg = rowMeans(select(., 2:272))) %>% 
  group_by(gene_name) %>% 
  slice_max(Avg, n = 1, with_ties = F ) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select( -gene_id, -Avg, -gene_id.1)

RNA_input_all_new <- RNA_input_all_new %>%
  select(sort(names(RNA_input_all_new)))

all(colnames(RNA_input_all_new) == RNA_all_QC_pass$Cell.Barcode)

colnames(RNA_input_all_new) <- RNA_all_QC_pass$Annotation

RNA_input_all_new.1 <- RNA_input_all_new %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS_C10",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS_SVEC", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  filter(reads > 0) %>%
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

dge_object <- DGEList(counts = RNA_input_all_new.1, 
                        remove.zeros = TRUE)
dge_object <- calcNormFactors(dge_object)
normalized_counts <- cpm(dge_object, 
                           normalized.lib.sizes = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)
df_normalized_counts <- replace(df_normalized_counts, df_normalized_counts == 0, NA)
df_normalized_counts <- log2(df_normalized_counts)

#length(which(grepl("SVEC_[A-Z]", colnames(df_normalized_counts))))
#length(which(grepl("C10_[A-Z]", colnames(df_normalized_counts))))

median_CV_RNA <- df_normalized_counts %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC"),
                           Group = factor(Group, levels = c("C10", "SVEC"))) %>%
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

df_normalized_counts_CVs <- df_normalized_counts %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, 
               names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC"),
                           Group = factor(Group, levels = c("C10", "SVEC", "Blank"))) %>%
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

#Supplementary Figure 5F
df_normalized_counts_CVs %>%
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
#ggsave("CVs_distribution_DS.png", width = 6, height = 4, bg = "white")

#Figure 4E
df_normalized_counts_CVs %>%
  filter(method == "nanoSPINS") %>% 
  ggplot()+
  aes(x = Group, y = CV, fill = Group) +
  geom_violin(drop = FALSE) +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 1))+
  ylab("Coefficient of Variation")+
  xlab("")+
  #facet_wrap(~ method) +
  theme_minimal(base_size = 18)
#ggsave("CVs_distribution_NS.png", width = 6, height = 4, bg = "white")

df_normalized_reads_nanoSPINS <- df_normalized_counts %>% 
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
  filter(!is.na(normalized_reads)) %>%
  select(gene, SampleID, normalized_reads) %>%
  pivot_wider(names_from = SampleID, values_from = normalized_reads) %>%
  column_to_rownames(var= "gene")

sorted_colnames_RNA_NS <- sort(colnames(df_normalized_reads_nanoSPINS))
df_normalized_reads_nanoSPINS <- df_normalized_reads_nanoSPINS[, sorted_colnames_RNA_NS]
df_normalized_reads_nanoSPINS <- replace(df_normalized_reads_nanoSPINS, is.na(df_normalized_reads_nanoSPINS), 0)
RNA_BCs_NS <- RNA_all_QC_pass %>%
  filter(Type == "nanoPOTS") %>%
  arrange(Cell.Barcode)
RNA_BCs_NS <- RNA_BCs_NS[order(RNA_BCs_NS$Annotation),]
pca_RNA_NS <- PCAtools::pca(df_normalized_reads_nanoSPINS, 
                              scale = F, 
                              center = T)
pca_RNA_NS$metadata <- RNA_BCs_NS

#Figure 4F
PCAtools::biplot(pca_RNA_NS, x = "PC2", y =  "PC4", 
                 lab = RNA_BCs_NS$Cell.type,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, 
                 ellipseLevel = 0.90, 
                 ylim = c(-100, 200),  
                 xlim = c(-125, 125)
)
#ggsave("scRNAseq_PCA_NS_072725.png", width = 8, height = 4, bg = "white")

nanoSPINS_c10 <- df_normalized_counts %>% 
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
  filter(Group == "nanoPOTS_C10") %>%
  filter(!is.na(normalized_reads)) %>%
  select(gene, Group) %>%
  unique()

nanoSPINS_SVEC <- df_normalized_counts %>% 
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
  filter(Group == "nanoPOTS_SVEC") %>%
  filter(!is.na(normalized_reads)) %>%
  select(gene, Group) %>%
  unique()

#Supplementary Figure 5B
venn.diagram(
  x = list(
    nanoSPINS_c10 = nanoSPINS_c10$gene,
    nanoSPINS_SVEC = nanoSPINS_SVEC$gene
  ),
  category.names = c("nanoSPINS_c10", "nanoSPINS_SVEC"),
  height = 8,
  width = 8,
  units = "in",
  resolution = 300,
  disable.logging = TRUE,
  col = c("#F8766D", "#00BFC4"),
  fill = c(alpha("#f8766d", 0.3), alpha("#00BFC4", 0.4)),
  cex = 1.5,
  ext.text = TRUE,
  cat.cex = 1.5,
  ext.line.lwd = 0.000001,
  print.mode = 'raw',
  filename = "venn_compare2.png",
  output = TRUE,
  cat.dist = c(0.05, 0.05)
)

#Supplementary Figure 5D
RNA_all_QC_pass %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Condition == "Direct_Sorting") %>%
  mutate(Cell.type = factor(Cell.type, levels = c("C10", "SVEC"))) %>%
  ggplot() +
  aes(x = Cell.type, 
      y = gene_count) +
  geom_boxplot(aes(colour = Cell.type), 
               outlier.size = -1, 
               width = 0.45, na.rm = TRUE) +
  geom_jitter(aes(colour = Cell.type), 
              width = 0.2, na.rm = TRUE, 
              size = 3, alpha = 0.75) +
  stat_summary(fun = mean,
               geom = "point", 
               color = "black", 
               size = 2) +
  scale_y_continuous(limits = c(0,13500), 
                     breaks = seq(0,13500, by = 1500)) +
  scale_colour_manual(values = c("C10" = "#F8766D", 
                                 "SVEC" = "#00BFC4"))+
  ylab("Number of Genes Identified (n)")+
  xlab(NULL)+ 
  ggtitle(NULL)+
  theme_minimal(base_size = 18)
#ggsave("gene_ids_DS.png", width = 6, height = 4, bg = "white")

RNA_all_QC_pass %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Condition == "Direct_Sorting") %>%
  mutate(Cell.type = factor(Cell.type, levels = c("C10", "SVEC"))) %>%
  group_by(Cell.type) %>%
  summarise(mean(gene_count))

QC_DS_genes <- df_normalized_counts %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10",
                           grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC",
                           grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                           grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS",
                            str_detect(Group, "directsorted") ~ "directsorting")) %>%
  filter(method == "directsorting") %>%
  filter(!is.na(normalized_reads)) %>%
  select(gene) %>% 
  unique()

QC_NS_genes <- df_normalized_counts %>% 
  rownames_to_column(var = "gene") %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10",
                           grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC",
                           grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                           grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS",
                            str_detect(Group, "directsorted") ~ "directsorting")) %>%
  filter(method == "nanoPOTS") %>%
  filter(!is.na(normalized_reads)) %>%
  select(gene) %>% 
  unique()

#Supplementary Figure 5E
venn.diagram(
  x = list(
    QC_DS_genes = QC_DS_genes$gene,
    QC_NS_genes = QC_NS_genes$gene
  ),
  category.names = c("QC_DS_genes", "QC_NS_genes"),
  height = 15,
  width = 15,
  units = "in",
  resolution = 300,
  disable.logging = TRUE,
  col = c("#4169E1", "#87CEFA"),
  fill = c(alpha("#4169E1", 0.3), alpha("#87CEFA", 0.4)),
  cex = 2,
  ext.text = TRUE,
  cat.cex = 2,
  ext.line.lwd = 0.000001,
  print.mode = 'raw',
  filename = "venn_compare3.png",
  output = TRUE,
  cat.dist = c(0.05, 0.05)
)

correlation_matrix_scRNA <- cor(as.data.frame(lapply(df_normalized_counts, as.numeric)), method = "pearson", use = "pairwise.complete.obs")
colnames(correlation_matrix_scRNA) <- ifelse(
  grepl("C10_[A-Z]", colnames(correlation_matrix_scRNA)), 
  paste0("nanoSPINS_", colnames(correlation_matrix_scRNA)), 
  ifelse(
    grepl("SVEC_[A-Z]", colnames(correlation_matrix_scRNA)), 
    paste0("nanoSPINS_", colnames(correlation_matrix_scRNA)),
    ifelse(
      grepl("SVEC_[1-9]", colnames(correlation_matrix_scRNA)), 
      paste0("Direct_sorting_", colnames(correlation_matrix_scRNA)), 
      ifelse(
        grepl("C10_[1-9]", colnames(correlation_matrix_scRNA)), 
        paste0("Direct_sorting_", colnames(correlation_matrix_scRNA)), 
        colnames(correlation_matrix_scRNA)
      )
    )
  )
)

rownames(correlation_matrix_scRNA) <- ifelse(
  grepl("C10_[A-Z]", rownames(correlation_matrix_scRNA)), 
  paste0("nanoSPINS_", rownames(correlation_matrix_scRNA)), 
  ifelse(
    grepl("SVEC_[A-Z]", rownames(correlation_matrix_scRNA)), 
    paste0("nanoSPINS_", rownames(correlation_matrix_scRNA)),
    ifelse(
      grepl("SVEC_[1-9]", rownames(correlation_matrix_scRNA)), 
      paste0("Direct_sorting_", rownames(correlation_matrix_scRNA)), 
      ifelse(
        grepl("C10_[1-9]", rownames(correlation_matrix_scRNA)), 
        paste0("Direct_sorting_", rownames(correlation_matrix_scRNA)), 
        rownames(correlation_matrix_scRNA)
      )
    )
  )
)

hc_scRNA <- hclust(dist(correlation_matrix_scRNA))
reordered_correlation_matrix_scRNA <- correlation_matrix_scRNA[hc_scRNA$order, hc_scRNA$order]

#Supplementary Figure 6
png("correlation_plot_scRNAseq_all_072825.png", width = 3000, height = 3000, bg = "white", res = 300)
corrplot::corrplot(correlation_matrix_scRNA, 
                   method = 'shade', 
                   diag = FALSE, 
                   tl.cex = 0.2, 
                   order = 'alphabet', 
                   col.lim = c(0, 1)) %>% 
  corrplot::corrRect(c(1,65,133,198,271), 
                     lwd = 2, col = "black")
dev.off()

nanoSPINS_C10_RNA_matrix <- correlation_matrix_scRNA[grep("nanoSPINS_C10", rownames(correlation_matrix_scRNA)), grep("nanoSPINS_C10", colnames(correlation_matrix_scRNA))]
average_correlation_C10_NS <- mean(nanoSPINS_C10_RNA_matrix[upper.tri(nanoSPINS_C10_RNA_matrix) | lower.tri(nanoSPINS_C10_RNA_matrix)])
#print(average_correlation_C10_NS)

nanoSPINS_SVEC_RNA_matrix <- correlation_matrix_scRNA[grep("nanoSPINS_SVEC", rownames(correlation_matrix_scRNA)), grep("nanoSPINS_SVEC", colnames(correlation_matrix_scRNA))]
average_correlation_SVEC_NS <- mean(nanoSPINS_SVEC_RNA_matrix[upper.tri(nanoSPINS_SVEC_RNA_matrix) | lower.tri(nanoSPINS_SVEC_RNA_matrix)])
#print(average_correlation_SVEC_NS)

Direct_sorting_C10_RNA_matrix <- correlation_matrix_scRNA[grep("Direct_sorting_C10", rownames(correlation_matrix_scRNA)), grep("Direct_sorting_C10", colnames(correlation_matrix_scRNA))]
average_correlation_C10_DS <- mean(Direct_sorting_C10_RNA_matrix[upper.tri(Direct_sorting_C10_RNA_matrix) | lower.tri(Direct_sorting_C10_RNA_matrix)])
#print(average_correlation_C10_DS)

Direct_sorting_SVEC_RNA_matrix <- correlation_matrix_scRNA[grep("Direct_sorting_SVEC", rownames(correlation_matrix_scRNA)), grep("Direct_sorting_SVEC", colnames(correlation_matrix_scRNA))]
average_correlation_SVEC_DS <- mean(Direct_sorting_SVEC_RNA_matrix[upper.tri(Direct_sorting_SVEC_RNA_matrix) | lower.tri(Direct_sorting_SVEC_RNA_matrix)])
#print(average_correlation_SVEC_DS)

SVEC_RNA_matrix <- correlation_matrix_scRNA[grep("SVEC", rownames(correlation_matrix_scRNA)), grep("SVEC", colnames(correlation_matrix_scRNA))]
average_correlation_SVEC <- mean(SVEC_RNA_matrix[upper.tri(SVEC_RNA_matrix) | lower.tri(SVEC_RNA_matrix)])
#print(average_correlation_SVEC)

C10_RNA_matrix <- correlation_matrix_scRNA[grep("C10", rownames(correlation_matrix_scRNA)), grep("C10", colnames(correlation_matrix_scRNA))]
average_correlation_C10 <- mean(C10_RNA_matrix[upper.tri(C10_RNA_matrix) | lower.tri(C10_RNA_matrix)])
#print(average_correlation_C10)

#Figure 5A
venn.diagram(
  x = list(
    prot = protein_ids$Gene,
    QC_NS_genes = QC_NS_genes$gene
  ),
  category.names = c("protein_ids", "QC_NS_genes"),
  height = 15,
  width = 15,
  units = "in",
  resolution = 300,
  disable.logging = TRUE,
  col = c("#4B0082", "#006400"),
  fill = c(alpha("#4B0082", 0.3), alpha("#006400", 0.4)),
  cex = 2,
  ext.text = TRUE,
  cat.cex = 2,
  ext.line.lwd = 0.000001,
  print.mode = 'raw',
  filename = "venn_compare4_new.png",
  output = TRUE,
  cat.dist = c(0.05, 0.05)
)

#For mitochondrial reads percentage and differential expression/abundance analysis in both modalities.
scRNAseq_all <- RNA_input_all_new.1 %>%
  as.sparse()
scRNA_all_seuratObj <- CreateSeuratObject(counts = scRNAseq_all, 
                                          min.cells = 0, 
                                          min.features = 0)
DefaultAssay(scRNA_all_seuratObj)
scRNA_all_seuratObj <- PercentageFeatureSet(scRNA_all_seuratObj, 
                                            pattern = "^mt-", 
                                            col.name = "percent.mt")
RNA_all_QC_pass <- RNA_all_QC_pass[match(colnames(scRNAseq_all), 
                                         RNA_all_QC_pass$Annotation), ]
scRNA_all_seuratObj$celltype <- RNA_all_QC_pass$Cell.type
scRNA_all_seuratObj$condition <- RNA_all_QC_pass$Condition
scRNA_all_seuratObj_metadata <- scRNA_all_seuratObj@meta.data
scRNA_all_seuratObj_metadata$condition <- factor(scRNA_all_seuratObj_metadata$condition,
                                                 levels = c("nanoPOTS", "Direct_Sorting"))

#Subset the Seurat object based on barcodes in RNA_BCs_NS$Annotation
scRNA_NS_seuratObj <- subset(scRNA_all_seuratObj, 
                             cells = RNA_BCs_NS$Annotation)
RNA_NS_QCpass <- RNA_NS_QCpass[match(colnames(scRNA_NS_seuratObj), 
                                     RNA_NS_QCpass$Annotation), ]
scRNA_NS_seuratObj$celltype <- RNA_NS_QCpass$Cell.type
scRNA_NS_seuratObj_metadata <- scRNA_NS_seuratObj@meta.data
scRNA_NS_seuratObj_metadata$celltype <- factor(scRNA_NS_seuratObj_metadata$celltype,
                                               levels = c("C10", "SVEC"))
#Figure 4C
ggplot(scRNA_NS_seuratObj_metadata, 
       aes(x = celltype, y = percent.mt)) +
  geom_boxplot(aes(colour = celltype), 
               outlier.size = -1, 
               width = 0.45, na.rm = TRUE) +
  geom_jitter(aes(colour = celltype), 
              width = 0.2, na.rm = TRUE, 
              size = 3, alpha = 0.75) +
  stat_summary(fun = median,
               geom = "point", 
               color = "black", 
               size = 2) +
  scale_y_continuous(limits = c(0,20), 
                     breaks = seq(0,20, by = 4)) +
  scale_colour_manual(values = c("C10" = "#F8766D", 
                                 "SVEC" = "#00BFC4"))+
  ylab("Percent Mitochondrial Reads")+
  xlab(NULL)+ 
  ggtitle(NULL)+
  theme_minimal(base_size = 18)
#ggsave("mitochondrial_reads_NS.png", bg = "white", width = 6, height = 4)

Assays(scRNA_NS_seuratObj)
DefaultAssay(scRNA_NS_seuratObj)
scRNA_NS_seuratObj <- SCTransform(scRNA_NS_seuratObj, 
                                  vst.flavor = "v2", 
                                  verbose = FALSE)
scRNA_NS_seuratObj <- RunPCA(scRNA_NS_seuratObj, 
                             verbose = FALSE)
VizDimLoadings(scRNA_NS_seuratObj, dims = c(1:10), 
               reduction = "pca")
DimPlot(scRNA_NS_seuratObj, reduction = "pca", 
        dims = c(1, 3), 
        pt.size = 2, 
        alpha = 0.9, 
        group.by = "celltype")
DimHeatmap(scRNA_NS_seuratObj, dims = 1, balanced = TRUE)
DimHeatmap(scRNA_NS_seuratObj, dims = 1:10, balanced = TRUE, assays = 'SCT')
ElbowPlot(scRNA_NS_seuratObj)
scRNA_NS_seuratObj <- FindNeighbors(scRNA_NS_seuratObj, 
                                    dims = 1:5)
scRNA_NS_seuratObj <- FindClusters(scRNA_NS_seuratObj, 
                                   resolution = 1, 
                                   verbose = FALSE)
scRNA_NS_seuratObj <- RunUMAP(scRNA_NS_seuratObj, 
                              dims = 1:10, n.neighbors = 20)
DimPlot(scRNA_NS_seuratObj, 
        label = TRUE, 
        reduction = "umap",
        group.by = "celltype", 
        pt.size = 4, 
        alpha = 0.9)
scRNA_NS_seuratObj_C10markers_RNA <- FindMarkers(scRNA_NS_seuratObj, 
                                                 assay = "SCT", 
                                                 ident.1 = "C10", 
                                                 ident.2 = "SVEC",
                                                 group.by = "celltype",
                                                 only.pos = FALSE,
                                                 test.use = "wilcox_limma" ,
                                                 # min.diff.pct = 0.2, 
                                                 logfc.threshold = 0)

#protein_input_batch_corrected.1 <- read.csv(file = "protein_input_batch_corrected.1.csv", row.names = 1)
#meta_prot <- read.csv(file = "meta_prot.csv")
prot_matrix <- protein_input_batch_corrected.1 %>%
  as.matrix()
group_factor <- factor(meta_prot$Group)
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)
fit <- lmFit(prot_matrix, design)
contrast_matrix <- makeContrasts(C10_vs_SVEC = C10 - SVEC, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
diff_abundance_test <- topTable(fit2, 
                                coef = "C10_vs_SVEC", 
                                number = Inf, 
                                adjust.method = "BH") %>%
  rownames_to_column(var = "gene")

base_plot <- EnhancedVolcano(diff_abundance_test,
                             lab = rep("", nrow(diff_abundance_test)),
                             x = 'logFC',
                             y = 'P.Value',
                             pCutoff = 0.0224,
                             FCcutoff = 0.5,
                             pointSize = 4.0,
                             labSize = 0,
                             drawConnectors = FALSE,  # disable default connectors
                             xlim = c(-2, 2),
                             ylim = c(0, 45),
                             cutoffLineWidth = 1,
                             cutoffLineCol = "darkgrey",
                             title = "scProt nanoSPINS",
                             subtitle = NULL) +
  theme_minimal(base_size = 18)

label_df <- subset(diff_abundance_test, gene %in% c("Col1a1", "Col3a1", "Fn1", "H2-K1", "H2-D1"))
label_df$label_color <- ifelse(label_df$gene %in% c("Col1a1", "Col3a1", "Fn1"),
                               "black", "purple")
label_layer <- geom_text_repel(
  data = label_df,
  aes(x = logFC,
      y = -log10(P.Value),
      label = gene),
  color = label_df$label_color,
  size = 5,
  segment.color = "black",
  segment.size = 1,
  min.segment.length = 0,
  force = 2,                 
  max.overlaps = Inf,
  box.padding = 0.6,         
  point.padding = 0.5,       
  nudge_y = 5,
  nudge_x = 0.2,
  show.legend = FALSE
)
EVP1_Prot <- base_plot + label_layer
EVP1_Prot #Figure 5C
#ggsave("scProt_nanoSPINS_EVP_073125.png", bg = "white", width = 10, height = 7)

#Figure 5C
scRNA_sig <- scRNA_NS_seuratObj_C10markers_RNA %>%
  filter(abs(avg_log2FC) >= 1 & p_val_adj < 0.05) %>%
  rownames_to_column(var = "gene") %>%
  select(gene)
scProt_sig <- diff_abundance_test %>%
  filter(abs(logFC) >= 0.5 & adj.P.Val < 0.05) %>%
  select(gene)
common_sig <- intersect(scRNA_sig$gene, scProt_sig$gene)
base_plot_1 <- EnhancedVolcano(scRNA_NS_seuratObj_C10markers_RNA,
                             lab = rep("", nrow(scRNA_NS_seuratObj_C10markers_RNA)),
                             x = 'avg_log2FC', 
                             y = 'p_val',
                             pCutoff = 6.985574e-06, 
                             FCcutoff = 1, 
                             pointSize = 4.0,
                             labSize = 0,
                             drawConnectors = FALSE,
                             xlim = c(-10,10), 
                             ylim = c(0, 30),
                             cutoffLineWidth = 1,
                             cutoffLineCol = "darkgrey",
                             title = "scRNA nanoSPINS",
                             subtitle = NULL) +
  theme_minimal(base_size = 18)
label_df_1 <- subset(scRNA_NS_seuratObj_C10markers_RNA, rownames(scRNA_NS_seuratObj_C10markers_RNA) %in% c("Col1a1", "Col3a1", "Fn1", "H2-K1", "H2-D1"))
label_df_1$gene <- rownames(label_df_1)
label_df_1$label_color <- ifelse(label_df_1$gene %in% c("Col1a1", "Col3a1", "Fn1"),
                                 "black", "purple")
label_layer_1 <- geom_text_repel(
  data = label_df_1,
  aes(x = avg_log2FC,
      y = -log10(p_val),
      label = gene),
  color = label_df_1$label_color,
  size = 5,
  segment.color = "black",
  segment.size = 1,
  min.segment.length = 0,
  force = 2,                 
  max.overlaps = Inf,
  box.padding = 0.6,         
  point.padding = 0.5,       
  nudge_y = 5,
  nudge_x = 0.2,
  show.legend = FALSE
)
EVP1_scRNAseq <- base_plot_1 + label_layer_1
EVP1_scRNAseq #Figure 5D
#ggsave("scRNAseq_nanoSPINS_EVP_073125.png", bg = "white", width = 10, height = 7)

#Figure 5B
venn.diagram(
  x = list(
    scProt_sig = scProt_sig$gene,
    scRNA_sig = scRNA_sig$gene
  ),
  category.names = c("scProt_sig", "scRNA_sig"),
  height = 15,
  width = 15,
  units = "in",
  resolution = 300,
  disable.logging = TRUE,
  col = c("#4B0082", "#006400"),
  fill = c(alpha("#4B0082", 0.3), alpha("#006400", 0.4)),
  cex = 2,
  ext.text = TRUE,
  cat.cex = 2,
  ext.line.lwd = 0.000001,
  print.mode = 'raw',
  filename = "venn_compare6_new.png",
  output = TRUE,
  cat.dist = c(0.05, 0.05)
)
