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
# scProteomics data analysis
filter_prot <- read.csv("./filter.csv")
peptide_input <- read.delim("../nanoSPINS_data/results_101624/tmt-report/abundance_peptide_GN.tsv") %>%
  select(contains("C10")| contains("SVEC")|
           contains("Control") | contains ("Index") | contains("Index_2") | contains("Peptide")) %>%
  filter(!Index_2 %in% filter_prot$Index) %>%
  dplyr::select(-c(Index, Index_2)) %>%
  pivot_longer(!Peptide, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC",
                           grepl("Control", SampleID) ~ "Blank"))
peptide_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  select(Peptide, SampleID, Batch) %>%
  unique() %>%
  group_by(SampleID, Batch) %>% ## Can add Batch for studying batch effects 
  add_count(name = "n") %>% 
  select(Batch, n) %>% 
  distinct(Batch, n) %>%
  ggplot()+
  aes(x = Batch, y = n)+
  geom_violin()+ 
  geom_jitter(width = 0.1) +
  geom_hline(yintercept = 4609.396, 
             linetype = "dashed", 
             color = "blue",
             linewidth = 1) +
  stat_summary(fun = "mean",
               geom = "point",
               colour = "red") +
  scale_y_continuous(limits = c(0,8000), 
                     breaks = seq(0,8000, by = 1000))+
  ylab("Number of peptides (n)")+
  xlab("TMT batch") +
  theme_classic(base_size = 25)
#ggtitle("scProteomics") +
#geom_text(data = . %>% filter(Group == "SVEC" & n < 750), aes(label = SampleID), hjust = -0.2, vjust = 0.2, size = 3, check_overlap = TRUE) +
#geom_text(data = . %>% filter(Group == "C10" & n < 750), aes(label = SampleID), hjust = -0.2, vjust = 0.2, size = 3, check_overlap = TRUE)
#ggsave("./number_of_peptides_per_TMT_batch.png", width = 15, height = 6, bg = "white")

mean_peptides_per_batch <- peptide_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  select(Peptide, SampleID, Batch) %>%
  unique() %>%
  group_by(SampleID, Batch) %>% ## Can add Batch for studying batch effects 
  add_count(name = "n") %>% 
  ungroup() %>%
  select(Batch, n) %>% 
  distinct(Batch, n) %>%
  ungroup() %>%
  group_by(Batch) %>%
  mutate(mean = mean(n)) %>%
  distinct(Batch, mean) %>%
  ungroup() %>%
  mutate(mean2 = mean(mean)) %>%
  select(mean2) %>% unique()

#unique_peptides_per_batch <- peptide_input %>%   filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>% mutate(Batch = Batch[,2]) %>% select(Peptide, Batch) %>% unique() %>% group_by(Batch) %>% ## Can add Batch for studying batch effects add_count(name = "n") %>% ungroup() %>% select(Batch, n) %>% distinct(Batch, n) %>% ungroup()
#write.csv(unique_peptides_per_batch, file = "unique_peptides_per_batch.csv")

protein_input <- read.delim("../nanoSPINS_data/results_101624/tmt-report/abundance_protein_GN.tsv") %>%
  select(contains("C10")| contains("SVEC")| contains("Control") | 
           contains ("Index") | contains("Gene")) %>%
  filter(!Index %in% filter_prot$Index) %>%
  dplyr::select(-c(Index)) %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC",
                           grepl("Control", SampleID) ~ "Blank"))

protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  select(Gene, SampleID, Batch) %>%
  unique() %>%
  group_by(SampleID, Batch) %>% 
  add_count(name = "n") %>% 
  select(Batch, n) %>% 
  distinct(Batch, n) %>%
  ggplot()+
  aes(x = Batch, y = n)+
  geom_violin()+ 
  geom_jitter(width = 0.2) +
  scale_y_continuous(limits = c(0,2000), 
                     breaks = seq(0,2000, by = 250)) +
  geom_hline(yintercept = 1105.926, 
             linetype = "dashed", 
             color = "blue", 
             linewidth = 1) +
  stat_summary(fun = "mean",
               geom = "point",
               colour = "red") +
  ylab("Number of proteins (n)") +
  xlab("TMT batch") +
  theme_classic(base_size = 25)
#ggsave("./number_of_proteins_per_TMT_batch.png", width = 15, height = 6, bg = "white")

mean_proteins_per_batch <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  select(Gene, SampleID, Batch) %>%
  unique() %>%
  group_by(SampleID, Batch) %>%
  add_count(name = "n") %>% 
  ungroup() %>%
  select(Batch, n) %>% 
  distinct(Batch, n) %>%
  ungroup() %>%
  group_by(Batch) %>%
  mutate(mean = mean(n)) %>%
  distinct(Batch, mean)  %>%
  ungroup() %>%
  mutate(mean2 = mean(mean)) %>%
  select(mean2) %>% unique()

#unique_proteins_per_batch <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>% mutate(Batch = Batch[,2]) %>% select(Gene, Batch) %>%   unique() %>% group_by(Batch) %>% add_count(name = "n") %>% ungroup() %>%  select(Batch, n) %>% distinct(Batch, n) %>% ungroup()
#write.csv(unique_proteins_per_batch, file = "unique_proteins_per_batch.csv")

protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  group_by(SampleID) %>%  
  add_count(name = "n") %>% 
  select(SampleID, Group, n) %>% 
  distinct(SampleID, Group, n) %>%
  ggplot()+
  aes(x = fct_reorder(Group, n), y = n)+
  geom_boxplot(aes(colour = Group), outlier.size = -1, width = 0.45, na.rm = TRUE) +
  geom_jitter(aes(colour = Group), width = 0.2, na.rm = TRUE) +
  stat_summary(fun.y="mean",color="black", shape=4) +
  scale_y_continuous(limits = c(0,1500), breaks = seq(0,1500, by = 250)) +
  ylab("Number of Proteins Identified (n)")+
  xlab("cell_type")+ ggtitle("scProteomics") +
  theme_minimal(base_size = 14)
#ggsave("./number_of_proteins_identified_C10_SVEC.png", width = 6, height = 4, bg = "white")

protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  group_by(SampleID) %>%  
  add_count(name = "n") %>% 
  select(SampleID, Group, n) %>% 
  distinct(SampleID, Group, n) %>%
  ungroup() %>% 
  select(Group, n) %>%
  group_by(Group) %>%
  summarise(mean = mean(n))

proteins_in_cell_lines <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  select(Gene, Group) %>% 
  distinct(Gene, Group)
#write.csv(proteins_in_cell_lines, file = "proteins_in_cell_lines.csv")

#peptide_input %>% 
#  filter(!is.na(Intensity)) %>%
#  filter(Group != "Blank") %>%
#  group_by(SampleID) %>%  
#  add_count(name = "n") %>% 
#  select(SampleID, Group, n) %>% 
#  distinct(SampleID, Group, n) %>%
#  ggplot()+
#  aes(x = fct_reorder(Group, n), y = n)+
#  geom_boxplot(aes(colour = Group), outlier.size = -1, width = 0.45) + 
#  geom_jitter(aes(colour = Group), width = 0.2) +
#  stat_summary(fun.y="mean",color="black", shape=4) +
#  scale_y_continuous(limits = c(0,8000), 
#                     breaks = seq(0,8000, by = 1000))+
#  ylab("Number of peptides Identified (n)")+
#  xlab("cell_type")+ ggtitle("scProteomics") +
#  theme_minimal(base_size = 14)

#peptide_input %>% 
#  filter(!is.na(Intensity)) %>%
#  filter(Group != "Blank") %>%
#  group_by(SampleID) %>%  
#  add_count(name = "n") %>% 
#  select(SampleID, Group, n) %>% 
#  distinct(SampleID, Group, n) %>%
#  ungroup() %>% 
#  select(Group, n) %>%
#  group_by(Group) %>%
#  summarise(mean = mean(n)) %>% view()

protein_input %>%
  filter(!is.na(Intensity)) %>%
  #filter(Group != "Blank") %>%
  group_by(SampleID) %>%  
  add_count(name = "n") %>% 
  select(SampleID, Group, n) %>% 
  distinct(SampleID, Group, n) %>%
  ggplot()+
  aes(x = fct_reorder(Group, n), y = n)+
  geom_boxplot(aes(colour = Group), outlier.size = -1, width = 0.45, na.rm = TRUE) +
  geom_jitter(aes(colour = Group), width = 0.2, na.rm = TRUE) +
  stat_summary(fun.y="mean",color="black", shape=4) +
  scale_y_continuous(limits = c(0,1500), breaks = seq(0,1500, by = 250)) +
  scale_colour_manual(values = c("Blank" = "black", "C10" = "#F8766D", "SVEC" = "#00BFC4"))+
  ylab("Number of Proteins Identified (n)")+
  xlab("cell_type")+ 
  ggtitle("scProteomics")+
  theme_minimal(base_size = 14)
#ggsave("./number_of_proteins_identified_Blank_C10_SVEC.png", width = 8, height = 4, bg = "white")

protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  ggplot() +
  aes(x = Group, y = Intensity, fill = Group) +
  geom_violin(aes(colour = Group), width = 0.40)+ 
  #geom_jitter(aes(colour = Group), width = 0.15)+
  stat_summary(fun = "median",
               geom = "crossbar",
               width = 0.15,
               colour = "black")+
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50, by = 10)) +
  xlab("")+
  ylab("log2(Protein Intensity)") +
  theme_minimal(base_size = 14)
#ggsave("median_sample_intensities_C10_SVEC.png", width = 5, height = 4, bg = "white")

#median_group_intensities <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>%   select(Group, Intensity) %>% group_by(Group) %>% summarise(median = median(Intensity))

protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  group_by(Group, Gene) %>%
  summarize(mean = mean(Intensity)) %>%
  ggplot()+
  aes(x = Group, y = mean, fill = Group)+
  geom_violin(aes(colour = Group), width = 0.65)+ 
  #geom_jitter(aes(colour = Group), width = 0.15)+
  stat_summary(fun = "median",
               geom = "crossbar",
               width = 0.15,
               colour = "black")+
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50, by = 10))+
  xlab("")+
  ylab("log2(Protein Intensity)")+
  theme_minimal(base_size = 14)
#ggsave("median_protein_intensity_distribution_C10_SVEC.png", width = 5, height = 4, bg = "white")

median_protein_intensities_pergroup <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  select(Group, Gene, Intensity) %>%
  group_by(Group, Gene) %>%
  mutate(mean = mean(Intensity)) %>%
  distinct(Group, Gene, mean) %>%
  ungroup() %>%
  group_by(Group) %>%
  summarize(median = median(mean))

length(unique(protein_input$SampleID)[grepl("C10", unique(protein_input$SampleID))])
length(unique(protein_input$SampleID)[grepl("SVEC", unique(protein_input$SampleID))])
length(unique(protein_input$SampleID)[grepl("Control", unique(protein_input$SampleID))])

#seventy_percent_cutoff <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% group_by(Gene) %>% add_count(name = "n") %>% select(Gene, Group, n) %>% distinct() %>% filter( n >= 57.6) %>% ungroup() %>% select(Gene) %>% distinct()
#write.csv(seventy_percent_cutoff, file = "seventy_percent_cutoff.csv")

#thirty_percent_cutoff <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% group_by(Gene) %>% add_count(name = "n") %>% filter( n >= 57.6) %>% ungroup() %>% select(Gene) %>% distinct()
#write.csv(thirty_percent_cutoff, file = "thirty_percent_cutoff.csv")
  
#seventy_percent_cutoff <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% group_by(Gene) %>% add_count(name = "n") %>% filter( n >= 134.4) %>% ungroup() %>% select(Gene) %>% distinct()
#write.csv(seventy_percent_cutoff, file = "seventy_percent_cutoff.csv")

#fifty_percent_cutoff <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% group_by(Gene) %>% add_count(name = "n") %>% filter( n >= 96) %>% ungroup() %>% select(Gene) %>% distinct()
#write.csv(fifty_percent_cutoff, file = "fifty_percent_cutoff.csv")

No_cutoff <- protein_input %>% filter(!is.na(Intensity)) %>% filter(Group != "Blank") %>% select(Gene) %>% distinct()
write.csv(No_cutoff, file = "no_cutoff.csv")

protein_input %>%
  filter(!is.na(Intensity)) %>%
  select(Gene) %>%
  unique() %>%
  add_count(name = "n") %>%
  select(n) %>%
  unique()

protein_input_750wide <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  mutate(Channel = paste(Group, str_sub(SampleID,-1,-1), sep = "_")) %>%
  group_by(SampleID) %>%
  add_count(name = "n") %>%
  filter(Channel != "C10_2") %>%
  filter(n >= 750) %>% 
  ungroup() %>%
  filter(Group != "Blank") %>%
  dplyr::distinct(Gene, SampleID, .keep_all = TRUE) %>% 
  dplyr::select(Gene, SampleID, Intensity) %>%
  pivot_wider(names_from = SampleID, 
              values_from = Intensity) %>%
  column_to_rownames(var= "Gene")

protein_input_750wide %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", 
               values_to = "Intensity") %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ 
                             "C10",
                           grepl("SVEC", SampleID) ~ 
                             "SVEC")) %>%
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  filter(!is.na(Intensity)) %>%
  group_by(Gene, Batch, Group) %>%
  mutate(Intensity = 2^Intensity) %>%
  mutate(Avg = mean(Intensity)) %>%
  mutate(CV = (sd(Intensity)/(Avg))) %>%
  filter(!is.na(CV)) %>%
  distinct(Gene, CV, Group) %>% 
  ggplot()+
  aes(x = reorder(Group, CV), y = CV, fill = Group)+
  ylab("Coefficient of Variation") +
  xlab("") +  
  scale_y_continuous(limits = c(0,3), 
                     breaks = seq(0,3, by = 0.5)) +
  geom_violin() + 
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  theme_minimal(base_size = 14) + ggtitle("scProteomics (normalized; non-imputed)")
#ggsave("CVs.png", bg = "white", height = 4, width = 5)

#median_CVs <- protein_input_750wide %>% rownames_to_column(var = "Gene") %>% pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>% mutate(Group = case_when(grepl("C10", SampleID) ~  "C10", grepl("SVEC", SampleID) ~ "SVEC")) %>% mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>% mutate(Batch = Batch[,2]) %>%  filter(!is.na(Intensity)) %>%  group_by(Gene, Batch, Group) %>%  mutate(Intensity = 2^Intensity) %>%  mutate(Avg = mean(Intensity)) %>%  mutate(CV = (sd(Intensity)/(Avg))) %>%  filter(!is.na(CV)) %>% distinct(Gene, CV, Group) %>% ungroup() %>%  group_by(Group) %>%  summarise(median = median(CV))

protein_input_750wide_norm <- as.data.frame(median_normalization(as.matrix(protein_input_750wide)))
sorted_colnames <- sort(colnames(protein_input_750wide_norm))
protein_input_750wide_norm <- protein_input_750wide_norm[, sorted_colnames]

meta <- data.frame(colnames(protein_input_750wide_norm)) %>%
  mutate(SampleID = `colnames.protein_input_750wide_norm.`) %>%
  mutate(Group = case_when(grepl("C10", SampleID) ~ "C10",
                           grepl("SVEC", SampleID) ~ "SVEC")) %>%
  #grepl("Blank", SampleID) ~ "Blank"))
  mutate(Batch = str_match(SampleID, "_\\s*(.*?)\\s*_")) %>%
  mutate(Batch = Batch[,2]) %>%
  mutate(Channel = paste(Group, str_sub(SampleID,-1,-1), sep = "_"))

protein_input_batch_corr <- ComBat.NA(protein_input_750wide_norm, 
                                      meta$Batch, 
                                      par.prior = TRUE,
                                      mean.only = FALSE,
                                      prior.plots = FALSE)

protein_input_batch_corr <- as.data.frame(protein_input_batch_corr$`corrected data`)

protein_input_batch_corr.1 <- protein_input_batch_corr %>% 
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

fifty_percent_cutoff <- protein_input_batch_corr %>% 
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
  dplyr::select(Gene) %>% 
  dplyr::distinct(Gene, .keep_all = TRUE)
#write.csv(fifty_percent_cutoff, file = "fifty_percent_cutoff.csv")

seventy_percent_cutoff <- protein_input_batch_corr %>% 
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
  mutate(Impute = case_when(n >= 64 & Group == "SVEC"  & n3 >= 48 ~ "Keep",
                            n >= 48 & Group == "C10" & n3 >= 64 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>% 
  ungroup() %>% 
  dplyr::select(Gene) %>% 
  dplyr::distinct(Gene, .keep_all = TRUE)
#write.csv(seventy_percent_cutoff, file = "seventy_percent_cutoff.csv")

thirty_percent_cutoff <- protein_input_batch_corr %>% 
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
  mutate(Impute = case_when(n >= 27 & Group == "SVEC"  & n3 >= 21 ~ "Keep",
                            n >= 21 & Group == "C10" & n3 >= 27 ~ "Keep", 
                            TRUE ~ "Discard")) %>%
  filter(Impute == "Keep") %>% 
  ungroup() %>% 
  dplyr::select(Gene) %>% 
  dplyr::distinct(Gene, .keep_all = TRUE)
#write.csv(thirty_percent_cutoff, file = "thirty_percent_cutoff.csv")

protein_input_impute <-  DreamAI(protein_input_batch_corr.1, 
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

protein_input_impute %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(!Gene, names_to = "SampleID", values_to = "Intensity") %>%
  ggplot()+
  aes(x = SampleID, y = Intensity)+
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8, face = "bold", colour = "black"))

protein_input_impute %>%
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
  ggplot()+
  aes(x = reorder(Group, CV), y = CV, fill = Group)+
  ylab("Coefficient of Variation")+
  xlab("")+  
  scale_y_continuous(limits = c(0,3), breaks = seq(0,3, by = 0.5)) +
  geom_violin() +
  stat_summary(fun = "median",
               geom = "point",
               color = "black") +
  theme_minimal(base_size = 14)
#ggsave(filename = "CVs.png", bg = "white", height = 4, width = 5)

median_CVs_imputedData <- protein_input_impute %>%
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
  ungroup() %>%
  group_by(Group) %>%
  summarise(median = median(CV))

protein_input_impute <- protein_input_impute[, sorted_colnames]

pca_test <- PCAtools::pca(protein_input_impute, scale = T, center = T)
pca_test$metadata <- meta
PCAtools::biplot(pca_test, x = "PC2", y =  "PC3", 
                 lab = 'Group',
                 showLoadings = FALSE,
                 labSize = 0, 
                 drawConnectors = FALSE,
                 legendPosition = 'right', 
                 colby = 'Group',
                 #encircle = TRUE, 
                 ellipse = TRUE,
                 ellipseLevel = 0.95)
#ggsave("scProteomics_PCA_plot.png", width = 8, height = 4, bg = "white")

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
  geom_point(size = 3, alpha = 0.7)+
  stat_ellipse()+
  theme_minimal(base_size = 22)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = "white"))+
  scale_color_manual(values = c("#F8766D", "#00BFC4"))
#ggsave("umap_scProteomics.png", width = 8, height = 5, bg = "white")
write.csv(protein_input_impute, "nanoSPINS_protein_imputed_results102924.csv")
