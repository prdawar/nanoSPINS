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
#scRNA-seq data analysis
annota <- read.delim("masterkey.txt")
annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(!Cell.type == "NoCell") %>%
  filter(Condition == "nanoPOTS") %>%
  group_by(Cell.type) %>% 
  ggplot() +
  aes(x = reorder(Cell.type, gene_count), 
      y = gene_count) +
  geom_violin(aes(colour = Cell.type, fill = Cell.type), 
               outlier.size = 0, width = 0.5) + 
  stat_summary(fun.y="mean",color="black", geom = "point") +
  #geom_jitter(aes(colour = Cell.type), width = 0.25) +
  scale_y_continuous(limits = c(0,10000), 
                     breaks = seq(0,10000, by = 2000)) +
  ylab("Number of Genes Identified (n)") + 
  xlab("cell_type") + 
  ggtitle("scRNA-seq") +
  theme_minimal(base_size = 14)
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
  #scale_fill_manual(values = c(C10 = "green3", SVEC = "cornflowerblue", NoCell = "red")) +
  scale_y_continuous(limits = c(0,16000), 
                     breaks = seq(0,16000, by = 2000)) +
  ylab("Number of Genes Identified (n)") + 
  xlab("Condition") + 
  theme_minimal(base_size = 14)
#ggsave("number_of_genes_identified_in_each_condition_scRNAseq.png", width = 5, height = 4, bg = "white")

mean_genes_identified <- annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Cell.type != "NoCell") %>%
  #filter(Condition == "nanoPOTS") %>%
  group_by(Condition, Cell.type) %>%
  summarise(mean(gene_count))

convert <- read.delim("convert_names.txt")
#TMT_batch <- read.csv("22BJ4_prot_combat.na.csv", row.names = 1) 

protein_samples <- colnames(protein_input_impute)

RNA_BCs <- annota %>%
  filter(Annotation %in% protein_samples) %>% 
  filter(gene_count >= 1000) %>% 
  arrange(Cell.Barcode)
#write.csv(RNA_BCs, file = "RNA_sample_annotations.csv")

RNA_selected_samples_new <- annota %>%
  filter(Cell.type != "NoCell") %>%
  filter(gene_count >= 1000)
#write.csv(RNA_selected_samples_new, file = "RNA_sample_annotations.csv")
RNA_spins <- read.delim("SB22_12_PNNL_21.umicount.inex.all.tsv", 
                    row.names = "Gene") %>% 
  dplyr::select(one_of(RNA_selected_samples_new$Cell.Barcode)) 

RNA_spins2 <- RNA_spins %>% 
  rownames_to_column(var = "gene_id") %>%
  inner_join(., convert) %>% 
  mutate(Avg = rowSums(.[2:303]/302)) %>% 
  group_by(gene_name) %>% 
  slice_max(Avg, n = 1, with_ties = F ) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select( -gene_id, -Avg)
#rownames_to_column(var = "gene_id") %>% pivot_longer(!gene_id, names_to = "SampleID", values_to = "raw_counts") %>% group_by(SampleID) %>% summarise(Count = sum(raw_counts))

dge_object <- DGEList(counts = RNA_spins2, 
                      remove.zeros = TRUE)
dge_object <- calcNormFactors(dge_object)
normalized_counts <- cpm(dge_object, 
                         normalized.lib.sizes = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)
df_normalized_counts <- replace(df_normalized_counts, df_normalized_counts == 0, NA)
df_normalized_counts <- log2(df_normalized_counts)
#write.csv(df_normalized_counts, file = "df_normalized_counts.csv")
df_normalized_counts <- read.csv(file = "./df_normalized_counts.csv", row.names = 1)
 
all_identified_genes <- df_normalized_counts %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS_C10",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS_SVEC", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Group) %>% 
  add_count(name = "n") %>% 
  group_by(Gene) %>% 
  add_count(name = "n2") %>% 
  mutate(n3 = n2 - n) %>% 
  mutate(Impute = case_when(n >= 5 & Group == "SVEC"  & n3 >= 5 ~ "Keep",
                            n >= 5 & Group == "C10" & n3 >= 5 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>% 
  ungroup() %>%
  select(Gene, Group, method) %>%
  unique()
#write.csv(all_identified_genes, file = "all_identified_genes_cutoff5.csv")


df_normalized_counts %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted")) %>%
  filter(!is.na(normalized_reads)) %>%
  filter(normalized_reads >= 1) %>%
  group_by(Gene, Group) %>% 
  add_count(name = "n") %>% 
  group_by(Gene) %>% 
  add_count(name = "n2") %>% 
  mutate(n3 = n2 - n) %>% 
  mutate(Impute = case_when(n >= 5 & Group == "SVEC"  & n3 >= 5 ~ "Keep",
                            n >= 5 & Group == "C10" & n3 >= 5 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>%
  group_by(Gene, method, Group) %>% 
  filter(!is.na(normalized_reads)) %>%
  mutate(logReads = 2^normalized_reads) %>%
  mutate(Avg = mean(logReads)) %>%
  mutate(CV = (sd(logReads)/(Avg))) %>%
  distinct(Gene, CV, method, Group) %>%
  #filter(!is.na(CV)) %>%
  ggplot() +
  aes(x = reorder(Group, CV), y = CV, fill = Group) +
  ylab("Coefficient of Variation") +
  xlab("") +  
  ggtitle("") +
  scale_y_continuous(limits = c(0,5)) +
  geom_violin() + 
  facet_wrap(~method) +
  #geom_jitter() +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  theme_minimal(base_size = 14)
#ggsave("Coefficient_of_variation_scRNAseq.png", width = 6, height = 4, bg =  "white")

median_CVs_RNA <- df_normalized_counts %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted")) %>%
  filter(!is.na(normalized_reads)) %>%
  group_by(Gene, Group) %>% 
  add_count(name = "n") %>% 
  group_by(Gene) %>% 
  add_count(name = "n2") %>% 
  mutate(n3 = n2 - n) %>% 
  mutate(Impute = case_when(n >= 5 & Group == "SVEC"  & n3 >= 5 ~ "Keep",
                            n >= 5 & Group == "C10" & n3 >= 5 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>%
  group_by(Gene, method, Group) %>% 
  filter(!is.na(normalized_reads)) %>%
  filter(normalized_reads >= 1) %>%
  mutate(logReads = 2^normalized_reads) %>%
  mutate(Avg = mean(logReads)) %>%
  mutate(CV = (sd(logReads)/(Avg))) %>%
  distinct(Gene, CV, method, Group) %>%
  ungroup() %>%
  select(Group, method, CV) %>%
  group_by(Group, method) %>%
  filter(!is.na(CV)) %>%
  summarise(median(CV))

#all_identified_genes<- df_normalized_counts %>% rownames_to_column(var = "Gene") %>%   pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10", grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC", grepl("C10_[1-9]", SampleID) ~ "directsorted_C10", grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>% mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS",  str_detect(Group, "directsorted") ~ "directsorting")) %>% mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC")) %>% group_by(Gene, method, celltype) %>% filter(!is.na(normalized_reads)) %>%  filter(method == "nanoPOTS") %>%   filter(normalized_reads >= 1) %>%  ungroup() %>% select(Gene) %>% unique()

#C10_SVEC_genes<- df_normalized_counts %>% rownames_to_column(var = "Gene") %>% pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%  mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10", grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC", grepl("C10_[1-9]", SampleID) ~ "directsorted_C10", grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>% mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS", str_detect(Group, "directsorted") ~ "directsorting")) %>% mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC")) %>% group_by(Gene, method, celltype) %>% filter(!is.na(normalized_reads)) %>% filter(method == "nanoPOTS") %>% filter(normalized_reads >= 0.5) %>%  ungroup() %>% select(Gene, celltype) %>%  unique()
#write.csv(C10_SVEC_genes, file = "C10_SVEC_genes.csv")

#Directsorting_vs_nanoSPINS_genes <- df_normalized_counts %>% rownames_to_column(var = "Gene") %>%  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>% mutate(Group = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoPOTS_C10", grepl("SVEC_[A-Z]", SampleID) ~ "nanoPOTS_SVEC", grepl("C10_[1-9]", SampleID) ~ "directsorted_C10",  grepl("SVEC_[1-9]", SampleID) ~ "directsorted_SVEC")) %>% mutate(method = case_when(str_detect(Group, "nanoPOTS") ~ "nanoPOTS",  str_detect(Group, "directsorted") ~ "directsorting")) %>%  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC")) %>% group_by(Gene, method, celltype) %>% filter(!is.na(normalized_reads)) %>% filter(normalized_reads >= 1) %>%  ungroup() %>% select(Gene, method) %>%  unique()
#write.csv(Directsorting_vs_nanoSPINS_genes, file = "Directsorting_vs_nanoSPINS_genes.csv")

df_normalized_counts_wide <- df_normalized_counts %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "normalized_reads") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(method = case_when(grepl("C10_[A-Z]", SampleID) ~ "nanoSPINS",
                            grepl("SVEC_[A-Z]", SampleID) ~ "nanoSPINS", 
                            grepl("C10_[1-9]", SampleID) ~ "directsorted",
                            grepl("SVEC_[1-9]", SampleID) ~ "directsorted")) %>%
  filter(!is.na(normalized_reads)) %>%
  filter(normalized_reads >= 1) %>%
  group_by(Gene, Group) %>% 
  add_count(name = "n") %>% 
  group_by(Gene) %>% 
  add_count(name = "n2") %>% 
  mutate(n3 = n2 - n) %>% 
  mutate(Impute = case_when(n >= 5 & Group == "SVEC"  & n3 >= 5 ~ "Keep",
                            n >= 5 & Group == "C10" & n3 >= 5 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>%
  ungroup() %>%
  select(Gene, SampleID, normalized_reads) %>%
  pivot_wider(names_from = SampleID, values_from = normalized_reads)

df_normalized_counts_wide <- replace(df_normalized_counts_wide, is.na(df_normalized_counts_wide), 0)

df_normalized_reads_nanoSPINS <- df_normalized_counts_wide %>%
  pivot_longer(!Gene, names_to = "SampleID", 
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
  filter(SampleID %in% protein_samples) %>%
  select(Gene, SampleID, normalized_reads) %>%
  pivot_wider(names_from = SampleID, values_from = normalized_reads) %>%
  column_to_rownames(var= "Gene")

sorted_colnames_RNA <- sort(colnames(df_normalized_reads_nanoSPINS))
df_normalized_reads_nanoSPINS <- df_normalized_reads_nanoSPINS[, sorted_colnames_RNA]

RNA_BCs <- RNA_BCs[order(RNA_BCs$Annotation),]

pca_RNA <- PCAtools::pca(df_normalized_reads_nanoSPINS, 
                         scale = F, 
                         center = T)

PCAtools::biplot(pca_RNA, x = "PC2", y =  "PC4", 
                 lab = RNA_BCs$Cell.type,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, ellipseLevel = 0.90)
#ggsave("scRNASeq_PCA_plot.png", width = 10, height = 6, bg = "white")

sorted_colnames_RNA_new <- sort(colnames(df_normalized_counts_wide))
df_normalized_counts_wide <- df_normalized_counts_wide[, sorted_colnames_RNA_new] %>% column_to_rownames(var = "Gene")

RNA_selected_samples_new <- RNA_selected_samples_new[order(RNA_selected_samples_new$Annotation),]

pca_RNA_all <- PCAtools::pca(df_normalized_counts_wide, 
                             scale = F, 
                             center = T)
PCAtools::biplot(pca_RNA_all, x = "PC1", y =  "PC2", 
                 lab = RNA_selected_samples_new$Condition,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, ellipseLevel = 0.90)
#ggsave("scRNASeq_PCA_plot_nanoSPINSvsDirectsorting.png", width = 10, height =6, bg = "white")

df_normalized_reads_direct_sorting <- df_normalized_counts_wide %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", 
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
  #filter(SampleID %in% protein_samples) %>%
  select(Gene, SampleID, normalized_reads) %>%
  pivot_wider(names_from = SampleID, values_from = normalized_reads) %>%
  column_to_rownames(var= "Gene")

sorted_colnames_DS <- sort(colnames(df_normalized_reads_direct_sorting))
df_normalized_reads_direct_sorting <- df_normalized_reads_direct_sorting[, sorted_colnames_DS]

RNA_DS <- RNA_selected_samples_new %>% 
  filter(Condition == "Direct_Sorting")
RNA_DS <- RNA_DS[order(RNA_DS$Annotation),]
pca_RNA_DS <- PCAtools::pca(df_normalized_reads_direct_sorting, 
                             scale = F, 
                             center = T)
PCAtools::findElbowPoint(df_normalized_reads_direct_sorting)
PCAtools::biplot(pca_RNA_DS, x = "PC2", y =  "PC3", 
                 lab = RNA_DS$Cell.type,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, ellipseLevel = 0.90)
#ggsave("scRNASeq_PCA_plot_Directsorting.png", width = 10, height = 6, bg = "white")

elbow_DS <- findElbowPoint(pca_RNA_DS$variance)
horn_DS <- parallelPCA(df_normalized_reads_direct_sorting)
screeplot(pca_RNA_DS,
          components = getComponents(pca_RNA_DS, 1:10),
          vline = c(horn_DS$n, elbow_DS), ylim = c(0,25)) + 
  geom_label(aes(x = horn_DS$n, y = 20,
                 label = 'Horn\'s', vjust = 1, size = 8)) +
  geom_label(aes(x = elbow_DS, y = 20,
                 label = 'Elbow method', vjust = 1, size = 8))

#correlation analysis at library/sample level
protein_input_750wide_norm_new <- protein_input_750wide_norm %>% 
  setNames(paste0(names(.), "_Prot")) %>% rownames_to_column(var = "gene") 
df_normalized_reads_nanoSPINS_new <- df_normalized_reads_nanoSPINS %>%
  setNames(paste0(names(.), "_RNA")) %>%
  rownames_to_column(var = "gene")
combined_new <- full_join(protein_input_750wide_norm_new, df_normalized_reads_nanoSPINS_new) %>% column_to_rownames(var = "gene")

correlation_matrix_new <- cor(as.data.frame(lapply(combined_new, as.numeric)), method = "pearson", use = "pairwise.complete.obs")

#correlation_matrix_C10_new <- correlation_matrix_new[grep("^C10", rownames(correlation_matrix_new)), grep("^C10", colnames(correlation_matrix_new))]
#hc_C10_new <- hclust(dist(correlation_matrix_C10_new))
#reordered_correlation_matrix_C10_new <- correlation_matrix_C10_new[hc_C10_new$order, hc_C10_new$order]
#pheatmap(reordered_correlation_matrix_C10_new, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", color = colorRampPalette(c("white", "darkblue"))(100), main = "C10 Correlation heatmap",  fontsize_row = 15, fontsize_col = 15, cellwidth = 17, cellheight = 17, fontsize = 15, border_color = NA, filename = "./pheatmap_C10_new.png")

#correlation_matrix_SVEC_new <- correlation_matrix_new[grep("^SVEC", rownames(correlation_matrix_new)), grep("^SVEC", colnames(correlation_matrix_new))]
#hc_SVEC_new <- hclust(dist(correlation_matrix_SVEC_new))
#reordered_correlation_matrix_SVEC_new <- correlation_matrix_SVEC_new[hc_SVEC_new$order, hc_SVEC_new$order]
#pheatmap(reordered_correlation_matrix_SVEC_new, clustering_distance_rows = "correlation",  clustering_distance_cols = "correlation", color = colorRampPalette(c("white", "darkred"))(100), main = "SVEC Correlation heatmap", fontsize_row = 15, fontsize_col = 15, cellwidth = 17, cellheight = 17, fontsize = 15, border_color = NA, filename = "./pheatmap_SVEC_new.png")

correlation_matrix_new2 <- cor(as.data.frame(lapply(df_normalized_counts, as.numeric)), method = "pearson", use = "pairwise.complete.obs")
colnames(correlation_matrix_new2) <- ifelse(
  grepl("C10_[A-Z]", colnames(correlation_matrix_new2)), 
  paste0("nanoSPINS_", colnames(correlation_matrix_new2)), 
  ifelse(
    grepl("SVEC_[A-Z]", colnames(correlation_matrix_new2)), 
    paste0("nanoSPINS_", colnames(correlation_matrix_new2)),
    ifelse(
      grepl("SVEC_[1-9]", colnames(correlation_matrix_new2)), 
      paste0("Direct_sorting_", colnames(correlation_matrix_new2)), 
      ifelse(
        grepl("C10_[1-9]", colnames(correlation_matrix_new2)), 
        paste0("Direct_sorting_", colnames(correlation_matrix_new2)), 
        colnames(correlation_matrix_new2)
      )
    )
  )
)

rownames(correlation_matrix_new2) <- ifelse(
  grepl("C10_[A-Z]", rownames(correlation_matrix_new2)), 
  paste0("nanoSPINS_", rownames(correlation_matrix_new2)), 
  ifelse(
    grepl("SVEC_[A-Z]", rownames(correlation_matrix_new2)), 
    paste0("nanoSPINS_", rownames(correlation_matrix_new2)),
    ifelse(
      grepl("SVEC_[1-9]", rownames(correlation_matrix_new2)), 
      paste0("Direct_sorting_", rownames(correlation_matrix_new2)), 
      ifelse(
        grepl("C10_[1-9]", rownames(correlation_matrix_new2)), 
        paste0("Direct_sorting_", rownames(correlation_matrix_new2)), 
        rownames(correlation_matrix_new2)
      )
    )
  )
)

hc <- hclust(dist(correlation_matrix_new2))
reordered_correlation_matrix_new2 <- correlation_matrix_new2[hc$order, hc$order]
#pheatmap(reordered_correlation_matrix_new2, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", color = colorRampPalette(c("white", "darkblue"))(100), main = "C10 Correlation heatmap", fontsize_row = 15, fontsize_col = 15, cellwidth = 17, cellheight = 17, fontsize = 15, border_color = NA, filename = "./pheatmap_SPINS_vs_DS.png")

png("correlation_plot_scRNAseq_all.png", width = 2500, height = 2500, bg = "white", res = 300)
corrplot::corrplot(correlation_matrix_new2, method = 'shade', diag = FALSE, tl.cex = 0.2, order = 'alphabet') %>% corrplot::corrRect(c(1,67,135,217,303), lwd = 2, col = "red")
dev.off()

nanoSPINS_C10_RNA_matrix <- correlation_matrix_new2[grep("nanoSPINS_C10", rownames(correlation_matrix_new2)), grep("nanoSPINS_C10", colnames(correlation_matrix_new2))]
average_correlation_C10 <- mean(nanoSPINS_C10_RNA_matrix[upper.tri(nanoSPINS_C10_RNA_matrix) | lower.tri(nanoSPINS_C10_RNA_matrix)])
print(average_correlation_C10)

nanoSPINS_SVEC_RNA_matrix <- correlation_matrix_new2[grep("nanoSPINS_SVEC", rownames(correlation_matrix_new2)), grep("nanoSPINS_SVEC", colnames(correlation_matrix_new2))]
average_correlation_SVEC <- mean(nanoSPINS_SVEC_RNA_matrix[upper.tri(nanoSPINS_SVEC_RNA_matrix) | lower.tri(nanoSPINS_SVEC_RNA_matrix)])
print(average_correlation_SVEC)

Direct_sorting_C10_RNA_matrix <- correlation_matrix_new2[grep("Direct_sorting_C10", rownames(correlation_matrix_new2)), grep("Direct_sorting_C10", colnames(correlation_matrix_new2))]
average_correlation_C10_DS <- mean(Direct_sorting_C10_RNA_matrix[upper.tri(Direct_sorting_C10_RNA_matrix) | lower.tri(Direct_sorting_C10_RNA_matrix)])
print(average_correlation_C10_DS)

Direct_sorting_SVEC_RNA_matrix <- correlation_matrix_new2[grep("Direct_sorting_SVEC", rownames(correlation_matrix_new2)), grep("Direct_sorting_SVEC", colnames(correlation_matrix_new2))]
average_correlation_SVEC_DS <- mean(Direct_sorting_SVEC_RNA_matrix[upper.tri(Direct_sorting_SVEC_RNA_matrix) | lower.tri(Direct_sorting_SVEC_RNA_matrix)])
print(average_correlation_SVEC_DS)

#grep("nanoSPINS", rownames(cor_matrix)) and grep("nanoSPINS", colnames(cor_matrix)) identify rows and columns that contain "nanoSPINS" in their names.
#cor_matrix[ , ] subsets the correlation matrix based on the selected rows and columns.
#mean(nanoSPINS_matrix[upper.tri(nanoSPINS_matrix) | lower.tri(nanoSPINS_matrix)]) calculates the average of the off-diagonal elements, which represent pairwise correlations without the self-correlation on the diagonal.

SVEC_RNA_matrix <- correlation_matrix_new2[grep("SVEC", rownames(correlation_matrix_new2)), grep("SVEC", colnames(correlation_matrix_new2))]
average_correlation_SVEC_all <- mean(SVEC_RNA_matrix[upper.tri(SVEC_RNA_matrix) | lower.tri(SVEC_RNA_matrix)])
print(average_correlation_SVEC_all)

C10_RNA_matrix <- correlation_matrix_new2[grep("C10", rownames(correlation_matrix_new2)), grep("C10", colnames(correlation_matrix_new2))]
average_correlation_C10_all <- mean(C10_RNA_matrix[upper.tri(C10_RNA_matrix) | lower.tri(C10_RNA_matrix)])
print(average_correlation_C10_all)

#png("cross_madality_correlation_plot.png", width = 2500, height = 2500, bg = "white", res = 300)
corrplot::corrplot(correlation_matrix_new, method = 'shade', diag = FALSE, order = 'original', tl.cex = 0.2) %>% corrplot::corrRect(c(1,70,161,218,301), lwd = 2, col = "red")
dev.off()

correlation_matrix_new_df <- as.data.frame(correlation_matrix_new) %>% rownames_to_column(var = "SampleID1") %>% pivot_longer(!SampleID1, names_to = "SampleID2", values_to = "corr_value")

correlation_matrix_new_df %>%  
  mutate(
  SampleID1_1 = str_extract(SampleID1, "C10_[A-Z][0-9]_[0-9]|SVEC_[A-Z][0-9]_[0-9]"),
  SampleID2_2 = str_extract(SampleID2, "C10_[A-Z][0-9]_[0-9]|SVEC_[A-Z][0-9]_[0-9]")) %>%
  filter(str_detect(SampleID1, "Prot") & str_detect(SampleID2, "RNA") &
           SampleID1_1 == SampleID2_2) %>% 
  select(SampleID1, SampleID2, corr_value) %>%
  summarise(mean(corr_value))

###combined_new has NA and zero... make all NA and then re-run!!!
df_normalized_reads_nanoSPINS_new_0toNA <-  df_normalized_reads_nanoSPINS_new %>% 
  mutate(across(everything(), ~ replace(.x, .x == 0, NA)))
combined_new_0toNA <- full_join(protein_input_750wide_norm_new, df_normalized_reads_nanoSPINS_new_0toNA) %>% 
  column_to_rownames(var = "gene")
correlation_matrix_new_0toNA <- cor(as.data.frame(lapply(combined_new, as.numeric)), method = "pearson", use = "pairwise.complete.obs")

correlation_matrix_new_0toNAdf <- as.data.frame(correlation_matrix_new_0toNA) %>% 
  rownames_to_column(var = "SampleID1") %>% 
  pivot_longer(!SampleID1, names_to = "SampleID2", values_to = "corr_value")

correlation_matrix_new_0toNAdf %>%  mutate(
  SampleID1_1 = str_extract(SampleID1, "C10_[A-Z][0-9]_[0-9]|SVEC_[A-Z][0-9]_[0-9]"),
  SampleID2_2 = str_extract(SampleID2, "C10_[A-Z][0-9]_[0-9]|SVEC_[A-Z][0-9]_[0-9]")) %>%
  filter(str_detect(SampleID1, "Prot") & str_detect(SampleID2, "RNA") &
           SampleID1_1 == SampleID2_2) %>% 
  select(SampleID1, SampleID2, corr_value) %>% 
  summarise(mean(corr_value))
#No impact of 0's or NA's on the correlation values.