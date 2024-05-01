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

getwd()
setwd("C:/Users/dawa726/OneDrive - PNNL/Desktop/nanoSPINS_Pranav/nanoSPINS_Pranav/")

# scProteomics data analysis
filter_prot <- read.csv("./filter.csv")

TMT_long <- read.delim("abundance_protein_GN_new.tsv") %>%
  select(contains("C10")| contains("SVEC")|
           contains("Control") | contains ("Index") | contains ("Gene")) %>%
  filter(!Index %in% filter_prot$Index) %>%
  dplyr::select(-Index) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC",
                           grepl("Control", SampleID) ~ "Blank"))

TMT_long %>%
  filter(!is.na(Intensity)) %>%
  #filter(Group != "Blank") %>%
  #  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  #  mutate(Batch = Batch[,2]) %>%
  group_by(SampleID) %>% ## Can add Batch for studying batch effects 
  add_count(name = "n") %>%
  ggplot()+
  aes(x = fct_reorder(Group, -n), y = n, fill = Group)+
  geom_boxplot()+ #facet_wrap(~Batch)+
  scale_y_continuous(limits = c(0,1500))+
  ylab("Number of Proteins Identified")+
  xlab("cell_type")+
  geom_text(data = . %>% filter(Group == "SVEC" & n < 750), aes(label = SampleID), hjust = -0.2, vjust = 0.2, size = 2, check_overlap = TRUE) +
  geom_text(data = . %>% filter(Group == "C10" & n < 750), aes(label = SampleID), hjust = -0.2, vjust = 0.2, size = 2, check_overlap = TRUE) +
  theme_minimal(base_size = 18)
#ggsave("number_of_proteins_identified.png", width = 8, height = 6, bg = "white")

TMT_long %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  #  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  #  mutate(Batch = Batch[,2]) %>%
  group_by(SampleID) %>% ## Can add Batch for studying batch effects 
  add_count(name = "n") %>% 
  select(SampleID, Group, n) %>% 
  distinct(SampleID, Group, n) %>%
  ggplot()+
  aes(x = fct_reorder(Group, -n), y = n)+
  geom_boxplot(aes(colour = Group), outlier.size = -1) + geom_jitter(aes(colour = Group), width = 0.25) +
  scale_y_continuous(limits = c(0,1500), breaks = seq(0,1500, by = 250))+
  ylab("Number of Proteins Identified (n)")+
  xlab("cell_type")+ ggtitle("scProteomics") +
  #geom_text(data = . %>% filter(Group == "SVEC" & n < 750), aes(label = SampleID), hjust = -0.2, vjust = 0.2, size = 3, check_overlap = TRUE) +
  #geom_text(data = . %>% filter(Group == "C10" & n < 750), aes(label = SampleID), hjust = -0.2, vjust = 0.2, size = 3, check_overlap = TRUE) +
  theme_minimal(base_size = 20)
#ggsave("../nanoSPINS_Pranav/02102024/number_of_proteins_identified_new.png", width = 8, height = 6, bg = "white")

TMT_long %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  group_by(Group, Gene) %>%
  summarize(Med = median(Intensity)) %>%
  ggplot()+
  aes(x = Med, fill = Group)+
  geom_histogram(bins = 100, alpha = 0.5)+
  xlab("Log2(Median Protein Intensity)")+
  ylab("Count (n)")+
  theme_minimal(base_size = 18)
#ggsave("protein_intensities_per_group.png", width = 8, height = 6, bg = "white")

TMT_long %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  ggplot()+
  aes(x = Group, y = Intensity, fill = Group)+
  geom_violin()+ #geom_jitter() +
  theme_minimal(base_size = 18)+
  xlab("")+
  ylab("log2(Protein Intensity)")
#ggsave("sample_intensities.png", width = 8, height = 6, bg = "white")

TMT_long %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  group_by(Group, Gene) %>%
  summarize(mean = mean(Intensity)) %>%
  ggplot()+ 
  aes(x = reorder(Group, -mean), y = mean, fill = Group)+
  geom_violin(aes(colour = Group)) +
  theme_minimal(base_size = 18)+
  xlab("Group")+
  ylab("log2(Protein Intensity)") +
  scale_y_continuous(limits = c(0,30), breaks = seq(0,30, by = 10))
#ggsave("../nanoSPINS_Pranav/02102024/Mean_protein_intensities.png", width = 8, height = 6, bg = "white")

TMT_wide <- TMT_long %>%
  filter(Group != "Blank") %>%
  dplyr::distinct(Gene, SampleID, .keep_all = TRUE) %>% 
  dplyr::select(Gene, SampleID, Intensity) %>%
  pivot_wider(names_from = SampleID, values_from = Intensity) %>%
  column_to_rownames(var= "Gene")

TMT_wide_750 <- TMT_wide %>% rownames_to_column(var = "Gene") %>% 
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  filter(!is.na(Intensity)) %>%
  mutate(Channel = paste(Group, str_sub(SampleID,-1,-1), sep = "_")) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  filter(Channel != "C10_2") %>% 
  filter(n >= 750) %>% 
  ungroup() %>%
  dplyr::select(Gene, SampleID, Intensity) %>%
  pivot_wider(names_from = SampleID, values_from = Intensity) %>%
  column_to_rownames(var= "Gene")

TMT_wide_750 %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Batch, Group) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(Avg = mean(Intensity)) %>%
  mutate(CV = (sd(Intensity)/(Avg))) %>%
  filter(!is.na(CV)) %>%
  distinct(Gene, CV, Group) %>% 
  #filter(Batch %in% c("A6", "A7", "C7", "C8")) %>%
  ggplot()+
  aes(x = reorder(Group, -CV), y = CV, fill = Group)+
  ylab("Coefficient of Variation")+
  xlab("")+  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 0.5)) +
  geom_violin() + 
  #geom_jitter() +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  #stat_summary(fun = "median", geom = "point", color = "blue") +
  #facet_wrap(~Batch) + 
  theme_minimal(base_size = 16) + ggtitle("scProteomics (normalized; non-imputed)")

TMT_norm <- as.data.frame(median_normalization(as.matrix(TMT_wide_750)))

meta <- data.frame(colnames(TMT_norm)) %>%
  mutate(SampleID = `colnames.TMT_norm.`) %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  #grepl("Blank", SampleID) ~ "Blank"))
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  mutate(Channel = paste(Group, str_sub(SampleID,-1,-1), sep = "_"))

TMT_batch_filter750 <- ComBat.NA(TMT_norm, meta$Batch, par.prior = TRUE,
                                 mean.only = FALSE,
                                 prior.plots = FALSE)

TMT_batch_filter750 <- as.data.frame(TMT_batch_filter750$`corrected data`)

TMT_batch_wide <- TMT_batch_filter750 %>% 
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

TMT_impute <-  DreamAI(TMT_batch_wide, k = 10, maxiter_MF = 10, ntree = 100, 
                       maxnodes = NULL, maxiter_ADMIN = 30, tol = 10^(-2), 
                       gamma_ADMIN = NA, gamma = 50, CV = FALSE, 
                       fillmethod = "row_mean", maxiter_RegImpute = 10, 
                       conv_nrmse = 1e-06, iter_SpectroFM = 40, method = "KNN",
                       out = c("KNN"))

TMT_impute <- as.data.frame(TMT_impute$KNN)

TMT_impute %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  ggplot()+
  aes(x = SampleID, y = Intensity)+
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, face = "bold", colour = "black"))

TMT_impute %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Batch, Group) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(Avg = mean(Intensity)) %>%
  mutate(CV = (sd(Intensity)/(Avg))) %>%
  filter(!is.na(CV)) %>%
  distinct(Gene, CV, Group) %>% 
  #filter(Batch %in% c("A6", "A7", "C7", "C8")) %>%
  ggplot()+
  aes(x = reorder(Group, -CV), y = CV, fill = Group)+
  ylab("Coefficient of Variation")+
  xlab("")+  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 0.5)) +
  geom_violin() + 
  #geom_jitter() +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  #stat_summary(fun = "median", geom = "point", color = "blue") +
  #facet_wrap(~Batch) + 
  theme_minimal(base_size = 16) + ggtitle("scProteomics (normalized; imputed)")

pca_test <- PCAtools::pca(TMT_impute, scale = T, center = T)
pca_test$metadata <- meta
PCAtools::biplot(pca_test, x = "PC2", y =  "PC3", 
                 lab = 'Group',
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = 'right', colby = 'Group',
                 #encircle = TRUE, 
                 ellipse = TRUE,
                 ellipseLevel = 0.9)

scores <- as.data.frame(pca_test$rotated)
set.seed(888)
check <- uwot::umap(scores[, c(2,3,4,5)], 
                    n_neighbors = 25, 
                    n_epochs = 1000,
                    bandwidth = 1.1, 
                    ret_nn = T,
                    metric = "euclidean"
                    #scale = "none"
                    #init  = "pca"
                    # repulsion_strength = 1
)
df.plot <- as.data.frame(check$embedding) %>%
  rownames_to_column(var = "SampleID") %>%
  inner_join(.,meta) %>%
  mutate(Group = case_when(Group == "C10" ~ "C10",
                           TRUE ~ "SVEC"))
ggplot(df.plot, aes(x = V1, y = V2, color = Group)) + 
  geom_point(size = 3, alpha = 0.5)+
  stat_ellipse()+
  theme_minimal(base_size = 22)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white"))+
  scale_color_manual(values = c("blue", "red"))

TMT_impute_long <- TMT_impute %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity")

#write.csv(TMT_impute, "22BJ4_prot_combat.na_02282024.csv")

TMT_batch_filter750 %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  #filter(., Group == "C10") %>%
  group_by(Gene, Group) %>% 
  summarize(Avg = mean(Intensity)) %>%
  pivot_wider(names_from = Group, values_from = Avg) %>%
  mutate(Intensity_Difference = abs(C10 - SVEC)) %>%
  distinct(Gene, .keep_all = TRUE) %>% 
  ggplot(aes(x = C10, y = SVEC)) +
  geom_point() + 
  scale_x_continuous(limits = c(10,30),  
                     breaks = seq(10,30, by = 2.5)) + 
  scale_y_continuous(limits = c(10,30),  
                     breaks = seq(10,30, by = 2.5)) +
  labs(title = "Distribution of protein intensities",
       x = "C10", 
       y = "SVEC") + 
  theme_minimal(base_size = 15)

TMT_impute %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC"
  )) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  group_by(Gene, Group) %>% 
  summarize(Avg = mean(Intensity)) %>%
  pivot_wider(names_from = Group, values_from = Avg) %>%
  mutate(Intensity_Difference = abs(C10 - SVEC)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  ggplot(aes(x = C10, y = SVEC, label = ifelse(Intensity_Difference >= 1, 
                                               as.character(Gene), ""))) +
  geom_point() +
  geom_text(hjust = -0.1, vjust = -0.5, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "C10", y = "SVEC", title = "Protein intensity distribution between C10 and SVEC") + 
  scale_x_continuous(limits = c(10,30),  
                     breaks = seq(10,30, by = 2.5)) + 
  scale_y_continuous(limits = c(10,30),  
                     breaks = seq(10,30, by = 2.5)) +
  theme_minimal(base_size = 15)

#ggsave("Distribution_of_protein_intensities_imputed.png", width = 15, height = 10, bg = "white")

TMT_batch_filter750 %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  filter(!is.na(Intensity)) %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC"
  )) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  group_by(Gene, Group) %>% 
  summarize(Avg = mean(Intensity)) %>%
  pivot_wider(names_from = Group, values_from = Avg) %>%
  mutate(Intensity_Difference = abs(C10 - SVEC)) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  ggplot(aes(x = C10, y = SVEC, label = ifelse(Intensity_Difference >= 1, 
                                               as.character(Gene), ""))) +
  geom_point() +
  geom_text(hjust = -0.1, vjust = -0.5, size = 3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "C10", y = "SVEC", title = "Protein intensity distribution between C10 and SVEC") + 
  scale_x_continuous(limits = c(10,30),  
                     breaks = seq(10,30, by = 2.5)) + 
  scale_y_continuous(limits = c(10,30),  
                     breaks = seq(10,30, by = 2.5)) +
  theme_minimal(base_size = 15)

#ggsave("Distribution_of_protein_intensities.png", width = 15, height = 10, bg = "white")

# scRNA-seq data analysis
annota <- read.delim("masterkey.txt")
annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  #filter(Cell.type != "NoCell") %>%
  group_by(Condition) %>% 
  ggplot() +
  aes(x = fct_reorder(Condition, gene_count), 
      y = gene_count, fill = Cell.type) +
  geom_boxplot()+ 
  scale_fill_manual(values = c(C10 = "green3", SVEC = "cornflowerblue", NoCell = "red")) +
  scale_y_continuous(limits = c(0,14000), breaks = seq(0,14000, by = 2000))+
  theme_minimal(base_size = 15) + 
  ylab("Number of Genes Identified")+ 
  xlab("Condition")
#facet_wrap(~Batch)+
#ggsave("Number_of_Genes_Identified_scRNAseq.png", width = 10, height = 8, bg = "white")

annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(!Cell.type == "NoCell") %>%
  filter(Condition == "nanoPOTS") %>%
  group_by(Cell.type) %>% 
  ggplot() +
  aes(x = reorder(Cell.type, -gene_count), 
      y = gene_count) +
  geom_boxplot(aes(colour = Cell.type), 
               outlier.size = -1) + 
  geom_jitter(aes(colour = Cell.type), width = 0.25) + 
  #scale_fill_manual(values = c(C10 = "green3", SVEC = "cornflowerblue", NoCell = "red")) +
  scale_y_continuous(limits = c(0,8000), 
                     breaks = seq(0,8000, by = 1000)) +
  ylab("Number of Genes Identified (n)") + 
  xlab("cell_type") + 
  ggtitle("scRNA-seq") +
  theme_minimal(base_size = 20)
#ggsave("../nanoSPINS_Pranav/02102024/number_of_genes_identified.png", width = 8, height = 6, bg = "white")

annota %>% 
  select(Cell.type, Condition, gene_count) %>%
  filter(Cell.type != "NoCell") %>%
  group_by(Condition, Cell.type) %>% 
  ggplot() +
  aes(x = fct_reorder(Condition, gene_count), 
      y = gene_count) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(aes(colour = Cell.type), width = 0.20) +
  #scale_fill_manual(values = c(C10 = "green3", SVEC = "cornflowerblue", NoCell = "red")) +
  scale_y_continuous(limits = c(0,14000), 
                     breaks = seq(0,14000, by = 1000)) +
  ylab("Number of Genes Identified (n)") + 
  xlab("Condition") + 
  theme_minimal(base_size = 20)
#ggsave("../nanoSPINS_Pranav/02102024/number_of_genes_identified_in_each_condition.png", width = 8, height = 6, bg = "white")

convert <- read.delim("convert_names.txt")
#TMT_batch <- read.csv("22BJ4_prot_combat.na.csv", row.names = 1) 
TMT_impute <- TMT_impute %>% select(order(colnames(.))) 
TMTnames <- colnames(TMT_impute)

RNA_BCs <- annota %>%
  filter(Annotation %in% TMTnames) %>% 
  filter(gene_count > 1500) %>% 
  arrange(Cell.Barcode)

RNA_selected_samples_new <- annota %>%
  filter(Cell.type != "NoCell") %>%
  filter(gene_count >= 1500)

test1 <- read.delim("SB22_12_PNNL_21.umicount.inex.all.tsv", 
                    row.names = "Gene") %>% 
  dplyr::select(one_of(RNA_selected_samples_new$Cell.Barcode)) 

test2 <- test1 %>% 
  rownames_to_column(var = "gene_id") %>%
  inner_join(., convert) %>% mutate(Avg = rowSums(.[2:294])/293) %>% 
  group_by(gene_name) %>% slice_max(Avg, n = 1, with_ties = F ) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select( -gene_id, -Avg)
#rownames_to_column(var = "gene_id") %>% pivot_longer(!gene_id, names_to = "SampleID", values_to = "raw_counts") %>% group_by(SampleID) %>% summarise(Count = sum(raw_counts))

dge_object <- DGEList(counts = test2, remove.zeros = TRUE)
dge_object <- calcNormFactors(dge_object)
normalized_counts <- cpm(dge_object, normalized.lib.sizes = TRUE)
df_normalized_counts <- as.data.frame(normalized_counts)
test3 <- replace(df_normalized_counts, df_normalized_counts == 0, NA)
test4 <- log2(test3)
#write.csv(test4, file = "test4.csv"). I saved this file and manually changed the column names from cell.barcode to Annotation names. I'll return to this step and write a code to do this.

test4_meta <- RNA_selected_samples_new %>% select(Cell.Barcode, Annotation, Condition)

test5 <- read.csv(file = "test4.csv", row.names = 1)

test5 %>%
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
  group_by(Gene, method, celltype) %>% 
  filter(!is.na(normalized_reads)) %>%
  mutate(logReads = 2^normalized_reads) %>%
  mutate(Avg = mean(logReads)) %>%
  mutate(CV = (sd(logReads)/(Avg))) %>%
  distinct(Gene, CV, method, celltype) %>% 
  #filter(!is.na(CV)) %>%
  ggplot()+
  aes(x = reorder(celltype, CV), y = CV, fill = celltype)+
  ylab("Coefficient of Variation")+
  #xlab("")+  
  ggtitle("") +
  scale_y_continuous(limits = c(0,5)) +
  geom_violin() + 
  facet_wrap(~method) +
  #geom_jitter() +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  theme_minimal(base_size = 16)
#ggsave("Coefficient_of_variation_scRNAseq.png", width = 10, height = 8, bg =  "white")

#Seurat with Raw data
#test6 <- test2
#colnames(test6) <- colnames(test5)
#scRNAseq_rawdata <- test6 %>% dplyr::select(one_of(RNA_BCs$Annotation)) %>% as.sparse()
#scProt <- TMT_impute
#scProt <- scProt %>% dplyr::select(one_of(RNA_BCs$Annotation)) %>% as.sparse()
#all.equal(colnames(scRNAseq_rawdata), colnames(scProt))
#C10SVEC_rawdata_new <- CreateSeuratObject(counts = scRNAseq_rawdata)
#prot_assay <- CreateAssayObject(counts = scProt)
#C10SVEC_rawdata_new[["Prot"]] <- prot_assay
#Assays(C10SVEC_rawdata_new)
#DefaultAssay(C10SVEC_rawdata_new)
#C10SVEC_rawdata_new <- PercentageFeatureSet(C10SVEC_rawdata_new, pattern = "^mt-", col.name = "percent.mt")
#C10SVEC_rawdata_new$celltype <- RNA_BCs$Cell.type
#VlnPlot(C10SVEC_rawdata_new, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "celltype", pt.size = 1)
#ggsave("features_scRNAseq.png", width = 8, height = 6)
#FeatureScatter(C10SVEC_rawdata_new, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "celltype")
#FeatureScatter(C10SVEC_rawdata_new, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "celltype")
#C10SVEC_rawdata_new <- NormalizeData(C10SVEC_rawdata_new)
#C10SVEC_rawdata_new <- FindVariableFeatures(C10SVEC_rawdata_new, selection.method = "vst")
#top10new <- head(VariableFeatures(C10SVEC_rawdata_new), 10)
#VariableFeaturePlot(C10SVEC_rawdata_new)
#C10SVEC_rawdata_new <- ScaleData(C10SVEC_rawdata_new)
#C10SVEC_rawdata_new <- RunPCA(C10SVEC_rawdata_new, verbose = FALSE)
#print(C10SVEC_rawdata_new[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(C10SVEC_rawdata_new, dims = 1:2, reduction = "pca")
#DimPlot(C10SVEC_rawdata_new, reduction = "pca") + NoLegend()
#DimHeatmap(C10SVEC_rawdata_new, dims = 1, balanced = TRUE)
#DimHeatmap(C10SVEC_rawdata_new, dims = 1:10, balanced = TRUE)
#ElbowPlot(C10SVEC_rawdata_new)
#ggsave("elbowplot_scRNAseq_normalprotocol.png", width = 8, height = 6, bg = "white")
#DefaultAssay(C10SVEC_rawdata_new)
#C10SVEC_rawdata_new <- FindNeighbors(C10SVEC_rawdata_new, dims = 1:30)
#C10SVEC_rawdata_new <- FindClusters(C10SVEC_rawdata_new,resolution = 1, verbose = FALSE)
#C10SVEC_rawdata_new <- RunUMAP(C10SVEC_rawdata_new, dims = 1:5)
#DimPlot(C10SVEC_rawdata_new, label = TRUE, reduction = "umap", group.by = "celltype")
#ggsave("UMAP_scRNAseq_normalprotocol.png", width = 8, height = 6)

#DefaultAssay(C10SVEC_rawdata_new) <- 'Prot'
#VariableFeatures(C10SVEC_rawdata_new) <- rownames(C10SVEC_rawdata_new[["Prot"]])
#C10SVEC_rawdata_new <-  ScaleData(C10SVEC_rawdata_new) %>% RunPCA(reduction.name = 'apca')
#C10SVEC_rawdata_new <- FindMultiModalNeighbors(C10SVEC_rawdata_new, reduction.list = list("pca", "apca"), dims.list = list(1:4, 1:10), knn.range = 50, smooth = FALSE, return.intermediate = TRUE, modality.weight.name = c("RNA.weight", "Prot.weight"))

#C10SVEC_rawdata_new <- RunUMAP(C10SVEC_rawdata_new, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#C10SVEC_rawdata_new <- FindClusters(C10SVEC_rawdata_new, graph.name = "wsnn", algorithm = 1, resolution = 1, verbose = FALSE)
#p1 <- DimPlot(C10SVEC_rawdata_new, group.by = "celltype", reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
#p1
#C10SVEC_rawdata_new <- RunUMAP(C10SVEC_rawdata_new, reduction = 'pca', dims = 1:5, assay = 'RNA',reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
#C10SVEC_rawdata_new <- RunUMAP(C10SVEC_rawdata_new, reduction = 'apca', dims = 1:10, assay = 'Prot', reduction.name = 'prot.umap', reduction.key = 'protUMAP_')
#p3 <- DimPlot(C10SVEC_rawdata_new, reduction = 'rna.umap', group.by = "celltype", label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
#p4 <- DimPlot(C10SVEC_rawdata_new, group.by = "celltype", reduction = 'prot.umap', label = TRUE, repel = TRUE, label.size = 5) + NoLegend()
#p1+p3+p4
#ggsave("all_umaps_normalprocessing.png", width = 8, height = 6, bg = "white")

#Idents(C10SVEC_rawdata_new) <- C10SVEC_rawdata_new$celltype
#C10SVEC_rawdata_new.markers_RNA <- FindAllMarkers(C10SVEC_rawdata_new, assay = "RNA",  only.pos = TRUE, test.use = "wilcox_limma",# min.diff.pct = 0.2, logfc.threshold = 1)

#top10_RNA_normalprocessing <- C10SVEC_rawdata_new.markers_RNA %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
#top10_RNA_normalprocessing_log2fc <- C10SVEC_rawdata_new.markers_RNA %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#colorpalette <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
#DoHeatmap(C10SVEC_rawdata_new, group.by = "celltype", assay  ="RNA", features = as.character(top10_RNA_normalprocessing$gene), group.colors = c("gold", "seagreen"), label = TRUE, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(colorpalette)) + guides(color="none")
#ggsave("scRNAseq_normalprocessing_markers.png", width = 8, height = 6)

#C10SVEC_rawdata_new.markers_prot <- FindAllMarkers(C10SVEC_rawdata_new, assay = "Prot", only.pos = TRUE, test.use = "wilcox_limma" , # min.diff.pct = 0.2, logfc.threshold = 1)

#top10_Prot_normalprocessing_avgl2fc <- C10SVEC_rawdata_new.markers_prot %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
#top10_Prot_normalprocessing_padj <- C10SVEC_rawdata_new.markers_prot %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
#DoHeatmap(C10SVEC_rawdata_new, group.by = "celltype", assay  ="Prot", features = top10_Prot_normalprocessing_padj$gene, group.colors = c("gold", "seagreen"), label = TRUE, draw.lines = FALSE) + scale_fill_gradientn(colours = rev(colorpalette)) + guides(color="none")
#ggsave("scProt_normalprocessing_markers.png", width = 8, height = 6)

#scTRansform with raw data
test6 <- test2
colnames(test6) <- colnames(test5)
scRNAseq_rawdata <- test6 %>% 
  dplyr::select(one_of(RNA_BCs$Annotation)) %>% 
  as.sparse()
C10SVEC_rawdata <- CreateSeuratObject(counts = scRNAseq_rawdata)
scProt <- TMT_impute
scProt <- scProt %>% dplyr::select(one_of(RNA_BCs$Annotation)) %>% as.sparse()
all.equal(colnames(scRNAseq_rawdata), colnames(scProt))
C10SVEC_rawdata <- CreateSeuratObject(counts = scRNAseq_rawdata)
prot_assay <- CreateAssayObject(counts = scProt)
C10SVEC_rawdata[["Prot"]] <- prot_assay
Assays(C10SVEC_rawdata)
DefaultAssay(C10SVEC_rawdata)
C10SVEC_rawdata <- PercentageFeatureSet(C10SVEC_rawdata, pattern = "^mt-", col.name = "percent.mt")
C10SVEC_rawdata$celltype <- RNA_BCs$Cell.type
VlnPlot(C10SVEC_rawdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "celltype", pt.size = 1)
C10SVEC_rawdata <- SCTransform(C10SVEC_rawdata, vst.flavor = "v2", verbose = FALSE)
C10SVEC_rawdata <- RunPCA(C10SVEC_rawdata, verbose = FALSE)
VizDimLoadings(C10SVEC_rawdata, dims = 1:2, reduction = "pca")
DimPlot(C10SVEC_rawdata, reduction = "pca") + NoLegend()
DimHeatmap(C10SVEC_rawdata, dims = 1, balanced = TRUE)
DimHeatmap(C10SVEC_rawdata, dims = 1:20, balanced = TRUE, assays = 'SCT')
ElbowPlot(C10SVEC_rawdata)
#ggsave("elbowplot_scRNAseq_rawdata_scT.png", width = 8, height = 6, bg = "white")
C10SVEC_rawdata <- FindNeighbors(C10SVEC_rawdata, dims = 1:30) ##why dims = 1:4?
C10SVEC_rawdata <- FindClusters(C10SVEC_rawdata,resolution = 1, verbose = FALSE)
C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, dims = 1:5)
DimPlot(C10SVEC_rawdata, 
        label = TRUE, 
        reduction = "umap",
        group.by = "celltype")
#ggsave("UMAP_scRNAseq_rawdata_scT.png", width = 8, height = 6)
Assays(C10SVEC_rawdata)
DefaultAssay(C10SVEC_rawdata) <- 'Prot'
VariableFeatures(C10SVEC_rawdata) <- rownames(C10SVEC_rawdata[["Prot"]])
C10SVEC_rawdata <-  ScaleData(C10SVEC_rawdata) %>% RunPCA(reduction.name = 'apca')
C10SVEC_rawdata <- FindMultiModalNeighbors(C10SVEC_rawdata, 
                                           reduction.list = list("pca", "apca"),
                                           dims.list = list(1:4, 1:10),
                                           knn.range = 50, 
                                           smooth = FALSE,
                                           return.intermediate = TRUE,
                                           modality.weight.name = c("RNA.weight", "Prot.weight"))

C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           nn.name = "weighted.nn", 
                           reduction.name = "wnn.umap", 
                           reduction.key = "wnnUMAP_"
)
C10SVEC_rawdata <- FindClusters(C10SVEC_rawdata, 
                                graph.name = "wsnn", 
                                algorithm = 1, 
                                resolution = 1, 
                                verbose = FALSE
)
p1 <- DimPlot(C10SVEC_rawdata, 
              group.by = "celltype",
              reduction = 'wnn.umap',
              label = TRUE, 
              repel = TRUE,
              label.size = 5) + NoLegend()
p1

C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           reduction = 'pca', 
                           dims = 1:5, 
                           assay = 'RNA',
                           reduction.name = 'rna.umap', 
                           reduction.key = 'rnaUMAP_'
)
C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata,
                           reduction = 'apca',
                           dims = 1:10,
                           assay = 'Prot',
                           reduction.name = 'prot.umap',
                           reduction.key = 'protUMAP_'
)
p3 <- DimPlot(C10SVEC_rawdata, reduction = 'rna.umap', 
              group.by = "celltype",
              label = TRUE, 
              repel = TRUE, 
              label.size = 5) + NoLegend()
p4 <- DimPlot(C10SVEC_rawdata, group.by = "celltype",
              reduction = 'prot.umap', label = TRUE, 
              repel = TRUE, label.size = 5) + NoLegend()
p1+p3+p4
#ggsave("all_umaps_scTranform.png", width = 8, height = 6, bg = "white")

Idents(C10SVEC_rawdata) <- C10SVEC_rawdata$celltype
C10SVEC_rawdata.markers_RNA <- FindAllMarkers(C10SVEC_rawdata, 
                                              assay = "SCT",
                                              only.pos = TRUE,
                                              test.use = "wilcox_limma" ,
                                              # min.diff.pct = 0.2, 
                                              logfc.threshold = 1)

top10_RNA_scT <- C10SVEC_rawdata.markers_RNA %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
top10_RNA_scT_log2fc <- C10SVEC_rawdata.markers_RNA %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
colorpalette <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
heat1 <- DoHeatmap(C10SVEC_rawdata, 
                   group.by = "celltype", 
                   assay  ="SCT", 
                   features = as.character(top10_RNA_scT$gene), 
                   group.colors = c("gold", "seagreen"), 
                   label = TRUE, 
                   draw.lines = FALSE) + 
  scale_fill_gradientn(colours = rev(colorpalette)) + 
  guides(color="none")
#ggsave("scRNAseq_scT_markers.png", width = 8, height = 6)

C10SVEC_rawdata.markers_prot <- FindAllMarkers(C10SVEC_rawdata, 
                                               assay = "Prot",
                                               only.pos = TRUE,
                                               test.use = "wilcox_limma" ,
                                               # min.diff.pct = 0.2, 
                                               logfc.threshold = 1)

top10_Prot_scT_avgl2fc <- C10SVEC_rawdata.markers_prot %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10_Prot_scT_padj <- C10SVEC_rawdata.markers_prot %>% group_by(cluster) %>% top_n(n = -5, wt = p_val_adj)
heat2 <- DoHeatmap(C10SVEC_rawdata, 
                   group.by = "celltype",
                   assay  ="Prot", 
                   features = top10_Prot_scT_padj$gene, 
                   group.colors = c("gold", "seagreen"), 
                   label = TRUE, draw.lines = FALSE) + 
  scale_fill_gradientn(colours = rev(colorpalette)) + 
  guides(color="none")
heat1 + heat2 + theme_minimal(base_size = 16)
#ggsave("scRNAseq_scProt_scT_markers.png", width = 16, height = 8)

#plotting PCAs for Seurat normal processing and ScTranform approach
#normalized_scRNAseq_normalprocedure <- as.data.frame(as.matrix(GetAssayData(C10SVEC_rawdata_new, assay= 'RNA', slot = 'data')))
normalized_scRNAseq_scT <- as.data.frame(as.matrix(GetAssayData(C10SVEC_rawdata, assay= 'SCT', slot = 'data')))
#For CV
normalized_scRNAseq_scT %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10",
                              str_detect(SampleID, "SVEC") ~ "SVEC")) %>%
  group_by(Gene, celltype) %>% 
  filter(!is.na(normalized_reads)) %>%
  mutate(logReads = 2^normalized_reads) %>%
  mutate(Avg = mean(logReads)) %>%
  mutate(CV = (sd(logReads)/(Avg))) %>%
  distinct(Gene, CV, celltype) %>% 
  #filter(!is.na(CV)) %>%
  ggplot()+
  aes(x = reorder(celltype, -CV), y = CV, fill = celltype)+
  ylab("Coefficient of Variation")+
  xlab("")+  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 0.5)) +
  geom_violin() + 
  #geom_jitter() +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  theme_minimal(base_size = 16) + ggtitle("scRNA-seq (normalized; non-imputed)")
meta_RNA_pca <- RNA_selected_samples_new %>% filter(Condition != "Direct_Sorting") %>% select(Cell.Barcode,Annotation,Cell.type)
#pca_RNA_test <- df_normalized_counts %>% dplyr::select(one_of(meta_RNA_pca$Cell.Barcode)) %>% log10(.)
#pca_RNA_test[pca_RNA_test == -Inf] <- 0
#pca_RNA <- PCAtools::pca(pca_RNA_test, scale = F, center = T)
#PCAtools::biplot(pca_RNA, x = "PC2", y =  "PC3", lab = meta_RNA_pca$Cell.type, showLoadings = FALSE, labSize = 0,   drawConnectors = FALSE, legendPosition = "right", #encircle = TRUE, ellipse = TRUE, ellipseLevel = 0.85)
#ggsave("PCA_scRNAseq_normalized_edgeR.png", width = 8, height = 6, bg = "white")
#pca_RNA <- PCAtools::pca(normalized_scRNAseq_normalprocedure, scale = F, center = T)
#PCAtools::biplot(pca_RNA, x = "PC2", y =  "PC3", lab = RNA_BCs$Cell.type, showLoadings = FALSE, labSize = 0, drawConnectors = FALSE, legendPosition = "right", #encircle = TRUE, ellipse = TRUE, ellipseLevel = 0.85)
#ggsave("PCA_scRNAseq_normalized_Seurat_NP.png", width = 8, height = 6, bg = "white")

pca_RNA <- PCAtools::pca(normalized_scRNAseq_scT, scale = F, center = T)
PCAtools::biplot(pca_RNA, x = "PC2", y =  "PC3", 
                 lab = RNA_BCs$Cell.type,
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = "right",
                 #encircle = TRUE, 
                 ellipse = TRUE, ellipseLevel = 0.85)

#ggsave("PCA_scRNAseq_normalized_Seurat_SCT.png", width = 8, height = 6, bg = "white")

# Ontology enrichment analysis
DefaultAssay(C10SVEC_rawdata) <- "SCT"
all.markers_SCT_scRNAseq <- FindMarkers(object = C10SVEC_rawdata, only.pos = FALSE, test.use = "wilcox_limma", min.cells.feature = 5, ident.1 = "C10", ident.2 = "SVEC")
go_test_scT <- all.markers_SCT_scRNAseq %>% dplyr::filter(p_val_adj < 0.05)
go_test_scT <- go_test_scT %>% rownames_to_column(var = "gene") %>% select(gene) %>% as.list()
bg_scT <- normalized_scRNAseq_scT %>% rownames_to_column(var = "gene") %>% pull(gene)
bg_common <- read.csv("bg_GO_common.csv", header = TRUE, row.names = 1)
# Manually created a background file (bg_GO_common.csv) with all the genes/proteins identified by scProteomics and scRNA-seq in C10 and SVEC. I'll come back to this and write a code to do the same.
bg_common <- c(t(bg_common))
GO_scRNAseq_scT <- gost(query = go_test_scT,
                        organism = "mmusculus", ordered_query = F, 
                        multi_query = FALSE, significant = TRUE, 
                        exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, evcodes = TRUE, 
                        user_threshold = 0.001, correction_method = "fdr", 
                        domain_scope = "annotated", custom_bg = bg_common, 
                        numeric_ns = "", sources = c("GO", "KEGG"), 
                        as_short_link = FALSE)
flattened_GO_scRNASeq_scT <- as.data.frame(GO_scRNAseq_scT$result) %>%
  mutate(across(everything(), ~ map_chr(., toString)))
#write.csv(as.data.frame(flattened_GO_scRNASeq), file = "GO_scRNAseq_ScT_SVEC_vs_C10.csv")
GO_plot1 <- gostplot(GO_scRNAseq_scT, capped = FALSE, interactive = FALSE)
GO_plot1_1 <- publish_gostplot(GO_plot1, highlight_terms = c("GO:0001816", "GO:0001819", "GO:0002376", "GO:0002682", "GO:0005102", "GO:0005198", "GO:0006950", "GO:0006952", "GO:0006955", "GO:0007015", "GO:0007155", "GO:0008283", "GO:0009605", "GO:0009607", "GO:0009615", "GO:0009617", "GO:0010033", "GO:0010628", "GO:0015399", "GO:0016477", "GO:0022900", "GO:0022904", "GO:0030029", "GO:0030036", "GO:0030155", "GO:0030312", "GO:0030334", "GO:0031012", "GO:0032101", "GO:0032103", "GO:0034097", "GO:0035456", "GO:0040011", "GO:0040012", "GO:0040017", "GO:0042221", "GO:0043207", "GO:0043254", "GO:0044419", "GO:0045069", "GO:0045087", "GO:0045785", "GO:0048870", "GO:0051128", "GO:0051239", "GO:0051240", "GO:0051607", "GO:0051707", "GO:0062023", "GO:0065003", "GO:0070887", "GO:0071310", "GO:0071345", "GO:0097435", "GO:0098542", "GO:0098800", "GO:0140546", "GO:1903900", "GO:2000145", "GO:2000147", "KEGG:00190", "KEGG:04145", "KEGG:05020", "KEGG:05415"), width = NA, height = NA, filename = NULL) + theme_minimal(base_size = 18)
# manually found the overlapping "GO" and "KEGG" terms. I'll come back to this and write a code to find the overlapping terms for scRNA and scProt. 
#ggsave(filename = "GO_results_scR.png", plot = gostplot1_1, width = 10, height = 25, bg = "white" )

DefaultAssay(C10SVEC_rawdata) <- "Prot"
all.markers_SCT_scProt <- FindMarkers(object = C10SVEC_rawdata, assay = "Prot", only.pos = FALSE, test.use = "wilcox_limma", min.cells.feature = 5, ident.1 = "SVEC", ident.2 = "C10")
go_test_scp <- all.markers_SCT_scProt %>% dplyr::filter(p_val_adj < 0.05)
go_test_scp <- go_test_scp %>% rownames_to_column(var = "gene") %>% select(gene) %>% as.list()
bg_scp <- TMT_impute %>% rownames_to_column(var = "gene") %>% pull(gene)
GO_scProt<- gost(query = go_test_scp,
                 organism = "mmusculus", ordered_query = F, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.001, correction_method = "fdr", 
                 domain_scope = "annotated", custom_bg = bg_common, 
                 numeric_ns = "", sources = c("GO", "KEGG"), as_short_link = FALSE)
flattened_GO_scProt <- as.data.frame(GO_scProt$result) %>%
  mutate(across(everything(), ~ map_chr(., toString)))
#write.csv(as.data.frame(flattened_GO_scProt), file = "GO_scProt_SVEC_vs_C10.csv")
GOplot_scP <- gostplot(GO_scProt, capped = FALSE, interactive = FALSE)
GOplot_scP_1 <- publish_gostplot(GOplot_scP, highlight_terms = c("GO:0001816", "GO:0001819", "GO:0002376", "GO:0002682", "GO:0005102", "GO:0005198", "GO:0006950", "GO:0006952", "GO:0006955", "GO:0007015", "GO:0007155", "GO:0008283", "GO:0009605", "GO:0009607", "GO:0009615", "GO:0009617", "GO:0010033", "GO:0010628", "GO:0015399", "GO:0016477", "GO:0022900", "GO:0022904", "GO:0030029", "GO:0030036", "GO:0030155", "GO:0030312", "GO:0030334", "GO:0031012", "GO:0032101", "GO:0032103", "GO:0034097", "GO:0035456", "GO:0040011", "GO:0040012", "GO:0040017", "GO:0042221", "GO:0043207", "GO:0043254", "GO:0044419", "GO:0045069", "GO:0045087", "GO:0045785", "GO:0048870", "GO:0051128", "GO:0051239", "GO:0051240", "GO:0051607", "GO:0051707", "GO:0062023", "GO:0065003", "GO:0070887", "GO:0071310", "GO:0071345", "GO:0097435", "GO:0098542", "GO:0098800", "GO:0140546", "GO:1903900", "GO:2000145", "GO:2000147", "KEGG:00190", "KEGG:04145", "KEGG:05020", "KEGG:05415"), width = NA, height = NA, filename = NULL) + theme_minimal(base_size = 18)
#ggsave(filename = "GO_results_scP.png", plot = gostplot2_1, width = 10, height = 25, bg = "white")

# Finding correlations at library level as well as for every (same) mRNA-Protein pairs.
normalized_scRNAseq_SCT_long <- normalized_scRNAseq_scT %>% 
  rownames_to_column(var = "Gene") %>%  
  pivot_longer(!Gene, names_to = "SampleID", 
               values_to = "TPM") %>%  
  mutate(TPM = ifelse(TPM == 0, NA, TPM)) %>% 
  filter(!is.na(TPM))

normalized_scRNAseq_scT_RNA <- normalized_scRNAseq_scT %>% 
  setNames(paste0(names(.), "_RNA")) %>% 
  rownames_to_column(var = "gene") %>%  
  mutate_all(~replace(., . == 0, NA))
#%>%  pivot_longer(!gene, names_to = "SampleID", values_to = "TPM") %>% mutate(TPM = ifelse(TPM == 0, NA, TPM))

TMT_norm_prot <- TMT_norm %>% 
  setNames(paste0(names(.), "_Prot")) %>% 
  rownames_to_column(var = "gene") 
#%>% pivot_longer(!gene, names_to = "SampleID", values_to = "intensity") %>% mutate(intensity = ifelse(intensity == 0, NA, intensity))

combined_new <- full_join(TMT_norm_prot, normalized_scRNAseq_scT_RNA) %>%
  column_to_rownames(var = "gene")
correlation_matrix_new <- cor(as.data.frame(lapply(combined_new, as.numeric)),
                              method = "pearson", 
                              use = "pairwise.complete.obs")
correlation_matrix_C10_new <- correlation_matrix_new[
  grep("^C10", rownames(correlation_matrix_new)), 
  grep("^C10", colnames(correlation_matrix_new))]
hc_C10_new <- hclust(dist(correlation_matrix_C10_new))
correlation_matrix_C10_new_reordered <- correlation_matrix_C10_new[hc_C10_new$order, hc_C10_new$order]
pheatmap(correlation_matrix_C10_new_reordered,
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         color = colorRampPalette(c("white", "darkred"))(100),
         main = "C10 Correlation heatmap", 
         fontsize_row = 15, 
         fontsize_col = 15, 
         cellwidth = 17, 
         cellheight = 17, 
         fontsize = 15, 
         border_color = NA) 
#filename = "./pheatmap_C10_new.png"

correlation_matrix_SVEC_new <- correlation_matrix_new[
  grep("^SVEC", rownames(correlation_matrix_new)), 
  grep("^SVEC", colnames(correlation_matrix_new))]
hc_SVEC_new <- hclust(dist(correlation_matrix_SVEC_new))
reordered_correlation_matrix_SVEC_new <- correlation_matrix_SVEC_new[hc_SVEC_new$order, hc_SVEC_new$order]
pheatmap(reordered_correlation_matrix_SVEC_new,
         clustering_distance_rows = "correlation", 
         clustering_distance_cols = "correlation", 
         color = colorRampPalette(c("white", "darkred"))(100),
         main = "SVEC Correlation heatmap", 
         fontsize_row = 15, fontsize_col = 15, cellwidth = 17, cellheight = 17, fontsize = 15, border_color = NA, #filename = "./pheatmap_SVEC_new.png"
         )

TMT_norm_long <- TMT_norm %>% 
  rownames_to_column(var = "gene") %>% 
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "Intensity")

normalized_scRNAseq_SCT_long_new <- normalized_scRNAseq_scT %>%
  rownames_to_column(var = "gene") %>%  
  pivot_longer(!gene, names_to = "SampleID", values_to = "TPM") %>%
  mutate(TPM = ifelse(TPM == 0, NA, TPM))

combined_long_new <- full_join(TMT_norm_long, normalized_scRNAseq_SCT_long_new) %>% 
  filter(!is.na(SampleID)) %>% 
  filter(!is.na(TPM)) %>% 
  filter(!is.na(Intensity)) %>% 
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10", 
                           grepl("SVEC", SampleID) ~ "SVEC")) 
#%>% group_by(Group, Gene) %>% filter(sd(Intensity) != 0 & sd(TPM) != 0) %>%   ungroup()

combined_C10_new <- combined_long_new %>% filter(Group == "C10")
combined_SVEC_new <- combined_long_new %>% filter(Group == "SVEC")

test10 <- combined_C10_new %>%
  select(SampleID, gene, Intensity, TPM) %>%
  group_by(gene) %>%
  filter(sd(Intensity) != 0 & sd(TPM) != 0) %>%
  mutate(number = n()) %>%
  mutate(correlation = cor(Intensity, TPM, 
                           method = "pearson", 
                           use = "pairwise.complete.obs")) %>% 
  #distinct(gene, correlation) %>% view()
  #filter(n() > 4) %>% # Filter out groups with less than 2 observations
  mutate(pvalue = ifelse(n() > 2, 
                         cor.test(Intensity, TPM, 
                                  method = "pearson")$p.value, 
                         NA)) %>% 
  select(gene, correlation, pvalue) %>%  
  distinct() 
#mutate(color_group = case_when(pvalue < 0.05 ~ "red")) %>%
#mutate(both = gene %in% gene_list)

test10 %>%
  mutate(padj = p.adjust(pvalue, 
                         n = nrow(test10), method = "BH")) %>%
  mutate(group = ifelse(pvalue < 0.05, 
                        "significant", "non-significant")) %>% 
  ggplot(aes(x = correlation, 
             y = -log10(pvalue), 
             color = factor(group))) + 
  geom_point(na.rm = TRUE) + 
  geom_text_repel(data = . %>% 
                    filter(pvalue < 0.05), aes(label = gene)) + 
                  #vjust = -1, direction = "x", size = 4) +
  xlab("correlation coefficients") + 
  ylab("-log10(p-value)") + 
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.5)) +
  scale_y_continuous(limits = c(0, 7.5), 
                     breaks = seq(from = 0, to = 7.5, by = 2.5)) + 
  scale_color_manual(values = c("significant" = "red", 
                                "non-significant" = "black"), 
                     name = "padj") + 
  theme_minimal(base_size = 14)

#ggsave("correlation_RNA_protein_non-imputed_C10_pvalue.png", width = 8, height = 6, bg = "white")

test12 <- combined_SVEC_new %>%
  select(SampleID, gene, Intensity, TPM) %>%
  group_by(gene) %>%
  filter(sd(Intensity) != 0 & sd(TPM) != 0) %>%
  mutate(number = n()) %>% 
  mutate(correlation = cor(Intensity, TPM, 
                           method = "pearson", 
                           use = "pairwise.complete.obs")) %>% 
  #distinct(gene, correlation) %>% view()
  #filter(n() > 4) %>% # Filter out groups with less than 2 observations
  mutate(pvalue = ifelse(n() > 2, 
                         cor.test(Intensity,
                                  TPM, method = "pearson")$p.value, NA)) %>%
  select(gene, correlation, pvalue, number) %>%  
  distinct() 
#mutate(color_group = case_when(pvalue < 0.05 ~ "red")) %>%
#mutate(both = gene %in% gene_list)

test12 %>%
  mutate(padj = p.adjust(pvalue, n = nrow(test12), 
                         method = "BH")) %>% 
  mutate(group = ifelse(pvalue < 0.05, "significant", "non-significant")) %>%
  #mutate(color_group = case_when(padj < 0.05 ~ "red")) %>%
  #mutate(both = gene %in% gene_list) %>%
  ggplot(aes(x = correlation, 
             y = -log10(pvalue), 
             color = factor(group))) + 
  geom_point(na.rm = TRUE) + 
  #geom_hline(yintercept = -log10(0.05),linetype = "dashed", linewidth = 0.75, color = "red") +
  geom_text_repel(data = . %>% filter(pvalue < 0.05), 
                  aes(label = gene)) + #vjust = -5, direction = "x", size = 4, ) +
  xlab("correlation coefficients") + 
  ylab("-log10(p-value)") + 
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.5)) +
  scale_y_continuous(limits = c(0, 7.5), 
                     breaks = seq(from = 0, to = 7.5, by = 2.5)) + 
  scale_color_manual(values = c("significant" = "red", "non-significant" = "black"), name = "padj") + theme_minimal(base_size = 14)
#ggsave("correlation_RNA_protein_non-imputed_SVEC_pvalue.png", width = 8, height = 6, bg = "white")

test10 %>%
  select(gene, correlation, pvalue) %>%
  mutate(group = ifelse(pvalue < 0.05, "significant", "non-significant")) %>%
  group_by(group) %>%
  filter(!is.na(pvalue)) %>%
  ggplot()+
  aes(x = correlation, fill = group)+
  geom_histogram(bins = 100, alpha = 0.5)+
  xlab("Pearson correlation (r)")+
  ylab("mRNA-Protein pair (n)") + theme_minimal()
#ggsave("distribution of Pearson correlations.png", width = 8, height = 6, bg = "white")

#For venn diagram (used evenn - an online platform)
C10_proteins <- TMT_long %>% 
  filter(!is.na(Intensity)) %>%
  filter(!Group == "Blank") %>% 
  select(Gene, Group) %>% 
  distinct() %>%
  filter(Group == "C10")
#write.csv(C10_proteins, "C10_protein_list.csv")
SVEC_proteins <- TMT_long %>% 
  filter(!is.na(Intensity)) %>%
  filter(!Group == "Blank") %>% 
  select(Gene, Group) %>% 
  distinct() %>%
  filter(Group == "SVEC")
#write.csv(SVEC_proteins, "SVEC_protein_list.csv")
venn_diagram_input <- list(C10_proteins$Gene, SVEC_proteins$Gene)
SVEC_genes <- normalized_scRNAseq_SCT_long_new %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10", grepl("SVEC", SampleID) ~ "SVEC")) %>% 
  mutate(TPM = ifelse(TPM == 0, NA, TPM)) %>% 
  filter(!is.na(TPM)) %>%
  select(gene, Group) %>% 
  distinct() %>%
  filter(Group == "SVEC")
#write.csv(SVEC_genes, "SVEC_gene_list.csv")
C10_genes <- normalized_scRNAseq_SCT_long_new %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10", grepl("SVEC", SampleID) ~ "SVEC")) %>% 
  mutate(TPM = ifelse(TPM == 0, NA, TPM)) %>% 
  filter(!is.na(TPM)) %>%
  select(gene, Group) %>% 
  distinct() %>%
  filter(Group == "C10")
#write.csv(C10_genes, "C10_gene_list.csv")
