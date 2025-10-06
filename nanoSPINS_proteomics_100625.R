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

getwd()
setwd("C:/Users/dawa726/OneDrive - PNNL/Desktop/nanoSPINS_Pranav")

#QC filtering -> normalization -> batch correction followed by identification depth plots -> coefficients of variation -> PCA -> differential abundance volcano plot (one for each modality) 

filter_prot <- read.csv("./filter.csv")
peptide_input <- read.delim("./abundance_peptide_GN.tsv") %>%
  select(contains("C10")| contains("SVEC")|
           contains("Control") | contains ("Index") | contains("Index_2") | contains("Peptide")) %>%
  filter(!Index_2 %in% filter_prot$Index) %>%
  dplyr::select(-c(Index, Index_2)) %>%
  pivot_longer(!Peptide, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC",
    grepl("Control", SampleID) ~ "Blank"),
  Group = factor(Group, levels = c("C10", "SVEC", "Blank")))

protein_input <- read.delim("./abundance_protein_GN.tsv") %>%
  select(contains("C10")| contains("SVEC")| contains("Control") | 
           contains ("Index") | contains("Gene")) %>%
  filter(!Index %in% filter_prot$Index) %>%
  dplyr::select(-c(Index)) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC",
    grepl("Control", SampleID) ~ "Blank"),
  Group = factor(Group, levels = c("C10", "SVEC", "Blank")))

#Supplementary Figure 4B
protein_input %>%
  filter(!is.na(Intensity)) %>%
  group_by(SampleID) %>%  
  add_count(name = "n") %>% 
  select(SampleID, Group, n) %>% 
  distinct(SampleID, Group, n) %>%
  ggplot()+
  aes(x = Group, y = n)+
  geom_boxplot(aes(colour = Group), 
               outlier.size = -1, 
               width = 0.45, na.rm = TRUE) +
  geom_jitter(aes(colour = Group), 
              width = 0.2, na.rm = TRUE, 
              size = 3, alpha = 0.75) +
  stat_summary(fun = mean,
               geom = "point", 
               color = "black", 
               size = 2) +
  scale_y_continuous(limits = c(0,1500), 
                     breaks = seq(0,1500, by = 250)) +
  scale_colour_manual(values = c("Blank" = "darkgreen", 
                                 "C10" = "#F8766D", 
                                 "SVEC" = "#00BFC4"))+
  ylab("Number of Proteins Identified (n)")+
  xlab(NULL)+ 
  ggtitle(NULL)+
  theme_minimal(base_size = 18)
#ggsave("Protein_IDs_beforeQC_072325.png", width = 8, height = 5, bg = "white")

overall_protein_identifications <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  select(Gene) %>%
  unique() %>%
  summarize(n = n())

protein_input_long750 <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  mutate(Channel = paste(Group, str_sub(SampleID, -1, -1), 
                         sep = "_")) %>%
  group_by(SampleID) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(Channel != "C10_2", 
         n >= 750, 
         Group != "Blank")

protein_input_750wide <- protein_input_long750 %>% 
  dplyr::select(Gene, SampleID, Intensity) %>%
  pivot_wider(names_from = SampleID, 
              values_from = Intensity) %>%
  column_to_rownames(var= "Gene")

protein_input_750wide_norm <- as.data.frame(median_normalization(as.matrix(protein_input_750wide)))
sorted_colnames <- sort(colnames(protein_input_750wide_norm))
protein_input_750wide_norm <- protein_input_750wide_norm[, sorted_colnames]

meta_prot <- data.frame(colnames(protein_input_750wide_norm)) %>%
  mutate(SampleID = `colnames.protein_input_750wide_norm.`) %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  mutate(Channel = paste(Group, str_sub(SampleID,-1,-1), sep = "_"))

protein_input_batch_corrected <- ComBat.NA(protein_input_750wide_norm, 
                                           meta_prot$Batch, 
                                           par.prior = TRUE,
                                           mean.only = FALSE,
                                           prior.plots = FALSE)

protein_input_batch_corrected <- as.data.frame(protein_input_batch_corrected$`corrected data`)

protein_input_long750_bc <- protein_input_batch_corrected %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  filter(!is.na(Intensity))

mean_proteins_overall <- protein_input_long750_bc %>%
  filter(!is.na(Intensity)) %>%
  select(Gene, SampleID) %>%
  distinct() %>%
  group_by(SampleID) %>%
  summarise(n = n()) %>%
  summarise(mean = mean(n))

mean_peptides_overall <- peptide_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  filter(SampleID %in% protein_input_long750_bc$SampleID) %>%
  select(Peptide, SampleID) %>%
  distinct() %>%
  group_by(SampleID) %>%
  summarise(n = n()) %>%
  summarise(mean = mean(n))

mean_proteins_cellgroup <- protein_input_long750_bc %>%
  select(Gene, SampleID, Group) %>%
  distinct() %>%
  group_by(SampleID, Group) %>%
  summarise(n = n()) %>%
  group_by(Group) %>%
  summarise(mean = mean(n))

mean_peptides_cellgroup <- peptide_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  filter(SampleID %in% protein_input_long750_bc$SampleID) %>%
  select(Peptide, SampleID, Group) %>%
  distinct() %>%
  group_by(SampleID, Group) %>%
  summarise(n = n()) %>%
  group_by(Group) %>%
  summarise(mean = mean(n))

#Figure 4A
protein_input_long750_bc %>%
  group_by(SampleID) %>%  
  add_count(name = "n") %>% 
  select(SampleID, Group, n) %>% 
  distinct(SampleID, Group, n) %>%
  ggplot()+
  aes(x = Group, y = n)+
  geom_boxplot(aes(colour = Group), 
               outlier.size = -1, 
               width = 0.45, 
               na.rm = TRUE) +
  geom_jitter(aes(colour = Group), 
              width = 0.2, 
              na.rm = TRUE, 
              size = 3, 
              alpha = 0.75) +
  stat_summary(fun = mean,
               geom = "point", 
               color = "black", 
               size = 2) +
  scale_y_continuous(limits = c(0,1500), 
                     breaks = seq(0,1500, 
                                  by = 250)) +
  scale_colour_manual(values = c("C10" = "#F8766D", "SVEC" = "#00BFC4"))+
  ylab("Number of Proteins Identified (n)")+
  xlab(NULL)+ 
  ggtitle(NULL)+
  theme_minimal(base_size = 18)
#ggsave("Protein_IDs_per_celltype_072325.png", width = 6, height = 5, bg = "white")

#length(unique(protein_input_long750_bc$SampleID)[grepl("C10", unique(protein_input_long750$SampleID))])
#length(unique(protein_input_long750_bc$SampleID)[grepl("SVEC", unique(protein_input_long750$SampleID))])

#Figure 4B
protein_input_long750_bc %>% 
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Group) %>%
  mutate(Intensity = 2^Intensity,
         Avg = mean(Intensity),
         CV = sd(Intensity) / Avg) %>% 
  filter(!is.na(CV)) %>%
  distinct(Gene, CV, Group) %>%
  ggplot(aes(x = Group, y = CV, fill = Group)) +
  geom_violin(trim = FALSE, scale = "width") + 
  stat_summary(fun = median,
               geom = "point",
               color = "black") +
  ylab("Coefficient of Variation") +
  xlab("") +
  scale_y_continuous(limits = c(0, 3), 
                     breaks = seq(0, 3, by = 1)) +
  theme_minimal(base_size = 18) +
  ggtitle(NULL)
#ggsave("CVs_scProt_072325.png", width = 6, height = 5, bg = "white")

median_CVs <- protein_input_long750_bc %>%
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Group) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(Avg = mean(Intensity)) %>%
  mutate(CV = (sd(Intensity)/(Avg))) %>%
  filter(!is.na(CV)) %>%
  distinct(Gene, CV, Group) %>% 
  ungroup() %>%
  group_by(Group) %>%
  summarise(median = median(CV))

#Supplementary Figure 4A
protein_stats <- protein_input_long750_bc %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  dplyr::select(SampleID, n) %>% 
  distinct() %>%
  ungroup() %>%
  summarize(mean_n = mean(n),
            sd_n = sd(n))
peptide_stats <- peptide_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  filter(SampleID %in% protein_input_long750_bc$SampleID) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  dplyr::select(SampleID, n) %>% 
  distinct() %>%
  ungroup() %>%
  summarize(
    mean_n = mean(n),  
    sd_n = sd(n))
combined_stats <- bind_rows(
  protein_stats %>% mutate(type = "Protein"),
  peptide_stats %>% mutate(type = "Peptide")) %>% 
  mutate(type = factor(type, levels = c("Peptide", "Protein")))
ggplot(combined_stats, aes(x = type, y = mean_n, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5) +
  geom_errorbar(aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n), 
                width = 0.2, 
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 7000), breaks = seq(0, 7000, 1000))+
  scale_fill_manual(values = c("Protein" = "#4682B4", "Peptide" = "#87CEFA")) +
  ggtitle(NULL) +
  xlab(NULL) +
  ylab("mean identifications (n)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none")
#ggsave("Mean_Protein_and_Peptide_Identifications.png", width = 8, height = 5, bg = "white")

#Supplementary Figure 4C
correlation_matrix <- cor(as.data.frame(lapply(protein_input_batch_corrected, as.numeric)), method = "pearson", use = "pairwise.complete.obs")
hc <- hclust(dist(correlation_matrix))
reordered_correlation_matrix <- correlation_matrix[hc$order, hc$order]
png("correlation_plot_scProt_all_072325.png", width = 3000, height = 3000, bg = "white", res = 300)
corrplot::corrplot(correlation_matrix, 
                   method = 'shade', 
                   diag = FALSE, 
                   tl.cex = 0.3, 
                   order = 'alphabet', 
                   col.lim = c(0, 1)) 
dev.off()

nanoSPINS_C10_Protein_matrix <- correlation_matrix[grep("C10", rownames(correlation_matrix)), grep("C10", colnames(correlation_matrix))]
average_correlation_C10 <- mean(nanoSPINS_C10_Protein_matrix[upper.tri(nanoSPINS_C10_Protein_matrix) | lower.tri(nanoSPINS_C10_Protein_matrix)])
#print(average_correlation_C10)
nanoSPINS_SVEC_Protein_matrix <- correlation_matrix[grep("SVEC", rownames(correlation_matrix)), grep("SVEC", colnames(correlation_matrix))]
average_correlation_SVEC <- mean(nanoSPINS_SVEC_Protein_matrix[upper.tri(nanoSPINS_SVEC_Protein_matrix) | lower.tri(nanoSPINS_SVEC_Protein_matrix)])
#print(average_correlation_SVEC)
average_correlation_overall <- mean(correlation_matrix[upper.tri(correlation_matrix) | lower.tri(correlation_matrix)])
#print(average_correlation_overall)

protein_input_batch_corrected.1 <- protein_input_batch_corrected %>% 
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Group) %>% 
  add_count(name = "n") %>% 
  group_by(Gene) %>% 
  add_count(name = "n2") %>% 
  mutate(n3 = n2 - n) %>% 
  mutate(Impute = case_when(n >= 45 & Group == "SVEC"  & n3 >= 34 ~ "Keep",
                            n >= 34 & Group == "C10" & n3 >= 45 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>% 
  ungroup() %>% 
  dplyr::distinct(Gene, SampleID, .keep_all = TRUE) %>% 
  dplyr::select(Gene, SampleID, Intensity) %>%
  pivot_wider(names_from = SampleID, values_from = Intensity) %>%
  column_to_rownames(var= "Gene")
protein_input_batch_corrected.1 <- protein_input_batch_corrected.1[, sorted_colnames]
protein_input_impute <-  DreamAI(protein_input_batch_corrected.1, 
                                 k = 50, 
                                 maxiter_MF = 10, 
                                 ntree = 100,
                                 maxnodes = NULL,
                                 maxiter_ADMIN = 30, 
                                 tol = 10^(-2),
                                 gamma_ADMIN = NA, 
                                 gamma = 50, 
                                 CV = FALSE,
                                 fillmethod = "row_mean" , 
                                 maxiter_RegImpute = 10,
                                 conv_nrmse = 1e-06, 
                                 iter_SpectroFM = 40, 
                                 method = "KNN",
                                 out = c("KNN"))

protein_input_impute <- as.data.frame(protein_input_impute$KNN)
protein_input_impute <- protein_input_impute[, sorted_colnames]
#write.csv(protein_input_impute, file = "protein_input_impute.csv")

#Figure 4C
pca_protein_input_impute <- PCAtools::pca(protein_input_impute, scale = T, center = T)
pca_protein_input_impute$metadata <- meta_prot
PCAtools::biplot(pca_protein_input_impute, x = "PC2", y =  "PC3", 
                 lab = 'Group',
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = 'right', 
                 colby = 'Group',
                 #encircle = TRUE, 
                 ellipse = TRUE,
                 ellipseLevel = 0.95)
#ggsave("scProteomics_PCA_plot_072325.png", width = 8, height = 4, bg = "white")

protein_ids <- protein_input_long750 %>%
  filter(!is.na(Intensity)) %>%
  select(Gene) %>% 
  unique()
