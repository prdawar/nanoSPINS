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

#Data integration, multimodal analysis, and marker analysis
RNA_input <- read.delim("SB22_12_PNNL_21.umicount.inex.all.tsv", 
                        row.names = "Gene") %>% 
  dplyr::select(one_of(RNA_BCs$Cell.Barcode))
#write.csv(RNA_input, file = "RNA_input.csv")
RNA_input <- read.csv("RNA_input.csv", row.names = 1) %>%
  rownames_to_column(var = "gene_id") %>%
  inner_join(., convert) %>% 
  mutate(Avg = rowSums(.[2:141]/140)) %>% 
  group_by(gene_name) %>% 
  slice_max(Avg, n = 1, with_ties = F ) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select( -gene_id, -Avg)
scRNAseq_rawdata <- RNA_input %>%
  as.sparse()
C10SVEC_rawdata <- CreateSeuratObject(counts = scRNAseq_rawdata, 
                                      min.cells = 5, 
                                      min.features = 1000)
scProt <- protein_input_impute
scProt <- scProt %>% 
  dplyr::select(one_of(RNA_BCs$Annotation)) %>% 
  as.sparse()
all.equal(colnames(scRNAseq_rawdata), colnames(scProt))
prot_assay <- CreateAssayObject(counts = scProt)
C10SVEC_rawdata[["Prot"]] <- prot_assay
Assays(C10SVEC_rawdata)
DefaultAssay(C10SVEC_rawdata)
C10SVEC_rawdata <- PercentageFeatureSet(C10SVEC_rawdata, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")
C10SVEC_rawdata$celltype <- RNA_BCs$Cell.type
VlnPlot(C10SVEC_rawdata, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "celltype", 
        pt.size = 1)

#as.data.frame(C10SVEC_rawdata$percent.mt) %>% rownames_to_column(var = "SampleID") %>% mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC")) %>% group_by(celltype) %>% summarise(median(`C10SVEC_rawdata$percent.mt`))

C10SVEC_rawdata <- SCTransform(C10SVEC_rawdata, 
                               vst.flavor = "v2", 
                               verbose = FALSE)
C10SVEC_rawdata <- RunPCA(C10SVEC_rawdata, 
                          verbose = FALSE)
VizDimLoadings(C10SVEC_rawdata, dims = c(1:10), 
               reduction = "pca")
DimPlot(C10SVEC_rawdata, reduction = "pca", 
        dims = c(1, 2), 
        pt.size = 2, 
        alpha = 0.9)
DimHeatmap(C10SVEC_rawdata, dims = 1, balanced = TRUE)
DimHeatmap(C10SVEC_rawdata, dims = 1:30, balanced = TRUE, assays = 'SCT')
ElbowPlot(C10SVEC_rawdata)
C10SVEC_rawdata <- FindNeighbors(C10SVEC_rawdata, 
                                 dims = 1:5)
C10SVEC_rawdata <- FindClusters(C10SVEC_rawdata, 
                                resolution = 1, 
                                verbose = FALSE)
C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           dims = 1:5, n.neighbors = 25)
DimPlot(C10SVEC_rawdata, 
        label = TRUE, 
        reduction = "umap",
        group.by = "celltype",
        shape.by = "seurat_clusters", 
        pt.size = 2, 
        alpha = 0.9)
#ggsave("C10SVEC_SCTRNA_UMAP.png", width = 6, height = 5, bg = "white")
C10SVEC_rawdata_C10markers_RNA <- FindMarkers(C10SVEC_rawdata, 
                                           assay = "SCT", 
                                           ident.1 = "C10", 
                                           ident.2 = "SVEC",
                                           group.by = "celltype",
                                           only.pos = FALSE,
                                           test.use = "wilcox_limma" ,
                                           # min.diff.pct = 0.2, 
                                           logfc.threshold = 0)
#write.csv(C10SVEC_rawdata_C10markers_RNA, file = "C10SVEC_rawdata_C10markers_RNA.csv")
EVP1_RNAseq <- EnhancedVolcano(C10SVEC_rawdata_C10markers_RNA, 
                               lab = rownames(C10SVEC_rawdata_C10markers_RNA), 
                               x = 'avg_log2FC', 
                               y = 'p_val', 
                               pCutoff = 3.541945e-06, 
                               FCcutoff = 1, 
                               pointSize = 3.0, 
                               labSize = 5, 
                               selectLab = c("Actb", "Actg1", "Bst2", "Col3a1", "Fdps", "Fn1", "H2-D1", "H2-K1", "Hmgb2", "Hsp90b1", "Hspb1", "Ifit1", "Ifit3", "Isg15", "Lgals3", "Rpl22", "Rplp1", "S100a4", "Stmn1", "Thbs1", "Top2a", "Vim"), 
                               drawConnectors = TRUE, 
                               xlim = c(-10,10), 
                               ylim = c(0, 25), 
                               cutoffLineWidth = 1, 
                               cutoffLineCol = "grey",
                               colConnectors = 'black',
                               arrowheads = FALSE, 
                               widthConnectors = 1,
                               parseLabels = TRUE, labCol = "black")
EVP1_RNAseq
DefaultAssay(C10SVEC_rawdata) <- "Prot"
VariableFeatures(C10SVEC_rawdata) <- rownames(C10SVEC_rawdata[["Prot"]])
C10SVEC_rawdata <-  ScaleData(C10SVEC_rawdata) %>% 
  RunPCA(reduction.name = 'apca')
C10SVEC_rawdata <- FindMultiModalNeighbors(C10SVEC_rawdata, 
                                           reduction.list = list("pca", "apca"),
                                           dims.list = list(1:5, 1:10),
                                           knn.range = 50, 
                                           smooth = FALSE,
                                           return.intermediate = TRUE,
                                           modality.weight.name = c("RNA.weight", "Prot.weight"))

C10SVEC_rawdata_C10markers_prot <- FindMarkers(C10SVEC_rawdata, 
                                              assay = "Prot", 
                                              ident.1 = "C10", 
                                              ident.2 = "SVEC",
                                              group.by = "celltype",
                                              only.pos = FALSE,
                                              test.use = "wilcox_limma" ,
                                              # min.diff.pct = 0.2, 
                                              logfc.threshold = 0)
#write.csv(C10SVEC_rawdata_C10markers_prot, file = "C10SVEC_rawdata_C10markers_prot.csv")
EVP1_Prot <- EnhancedVolcano(C10SVEC_rawdata_C10markers_prot, 
                             lab = rownames(C10SVEC_rawdata_C10markers_prot), 
                             x = 'avg_log2FC', 
                             y = 'p_val', 
                             pCutoff = 0.05, 
                             FCcutoff = 0.5, 
                             pointSize = 3.0, 
                             labSize = 5, 
                             selectLab = c("Actb", "Actg1", "Bst2", "Col3a1", "Fdps", "Fn1", "H2-D1", "H2-K1", "Hmgb2", "Hsp90b1", "Hspb1", "Ifit1", "Ifit3", "Isg15", "Lgals3", "Rpl22", "Rplp1", "S100a4", "Stmn1", "Thbs1", "Top2a", "Vim"), 
                             drawConnectors = TRUE, 
                             xlim = c(-5, 5), 
                             ylim = c(0, 25), 
                             cutoffLineWidth = 1, 
                             cutoffLineCol = "grey",
                             colConnectors = 'black',
                             arrowheads = FALSE, 
                             widthConnectors = 1, 
                             parseLabels = TRUE)
EVP1_Prot
EVP1_RNAseq + EVP1_Prot
ggsave("volcanoplots_scRNAseq_plus_scProt.png", width = 16, height = 8)

C10SVEC_rawdata <- FindClusters(C10SVEC_rawdata, 
                                graph.name = "wsnn", 
                                algorithm = 1, 
                                resolution = 1, 
                                verbose = FALSE
)

C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           nn.name = "weighted.nn", 
                           reduction.name = "wnn.umap", 
                           reduction.key = "wnnUMAP_"
)

p1 <- DimPlot(C10SVEC_rawdata, 
              group.by = "celltype",
              reduction = 'wnn.umap',
              label = TRUE, 
              repel = TRUE,
              label.size = 5) + NoLegend()
C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           reduction = 'pca', 
                           dims = 1:5, 
                           assay = 'RNA',
                           reduction.name = 'rna.umap', 
                           reduction.key = 'rnaUMAP_'
)
p2 <- DimPlot(C10SVEC_rawdata, reduction = 'rna.umap', 
              group.by = "celltype",
              label = TRUE, 
              repel = TRUE, 
              label.size = 5) + NoLegend()

C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata,
                           reduction = 'apca',
                           dims = 2:5,
                           assay = 'Prot',
                           reduction.name = 'prot.umap',
                           reduction.key = 'protUMAP_'
)

p3 <- DimPlot(C10SVEC_rawdata, group.by = "celltype",
              reduction = 'prot.umap', label = TRUE, 
              repel = TRUE, label.size = 5) + NoLegend()
p1+p2+p3
ggsave("integrative_umaps_all.png", width = 8, height = 5)

features <- c("Actb", "Actg1", "Bst2", "Col3a1", "Fdps", "Fn1", "H2-D1", "H2-K1", "Hmgb2", "Hsp90b1", "Hspb1", "Ifit1", "Ifit3", "Isg15", "Lgals3", "Rpl22", "Rplp1", "S100a4", "Stmn1", "Thbs1", "Top2a", "Vim")

features_rearrange <- c("Bst2", "H2-D1", "H2-K1", "Hmgb2", "Ifit3", "Isg15", "Lgals3", "Stmn1", "Actb", "Actg1", "Col3a1", "Fn1", "S100a4", "Vim")
  
#feature_test <- c("Actb", "Actg1")
#RidgePlot(C10SVEC_rawdata, features = feature_test, ncol = 2)
#FeaturePlot(C10SVEC_rawdata, features = feature_test)
#DoHeatmap(C10SVEC_rawdata, features = feature_test, size = 3, group.by = "celltype")
#DotPlot(C10SVEC_rawdata, features = features, split.by = "celltype") + RotatedAxis()

C10SVEC_rawdata_C10markers_RNA_features <- C10SVEC_rawdata_C10markers_RNA %>% 
  rownames_to_column("Gene") %>% 
  filter(Gene %in% features_rearrange)
C10SVEC_rawdata_C10markers_RNA_features <- C10SVEC_rawdata_C10markers_RNA_features[match(features_rearrange, C10SVEC_rawdata_C10markers_RNA_features$Gene, nomatch = 0), ]
colorpalette <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
heat1 <- DoHeatmap(C10SVEC_rawdata, 
                   group.by = "celltype", 
                   assay  ="SCT", 
                   features = as.character(C10SVEC_rawdata_C10markers_RNA_features$Gene), 
                   group.colors = c("#F8766D", "#00BFC4"),
                   label = TRUE, 
                   draw.lines = FALSE, size = 5) + 
  scale_fill_gradientn(colours = rev(colorpalette)) + 
  guides(color="none")

C10SVEC_rawdata_C10markers_Prot_features <- C10SVEC_rawdata_C10markers_prot %>%
  rownames_to_column("Gene") %>% 
  filter(Gene %in% features)
C10SVEC_rawdata_C10markers_Prot_features <- C10SVEC_rawdata_C10markers_Prot_features[match(features_rearrange, C10SVEC_rawdata_C10markers_Prot_features$Gene, nomatch = 0), ]
heat2 <- DoHeatmap(C10SVEC_rawdata, 
                   group.by = "celltype", 
                   assay  ="Prot", 
                   features = as.character(C10SVEC_rawdata_C10markers_Prot_features$Gene), 
                   group.colors = c("#F8766D", "#00BFC4"),
                   label = TRUE, 
                   draw.lines = FALSE, size = 5) + 
  scale_fill_gradientn(colours = rev(colorpalette)) + 
  guides(color="none")

heat1 + heat2
#ggsave("heatmap_integrative_analysis.png", width = 14, height = 5)

protein_list <- protein_input_impute %>%
  rownames_to_column("Gene") %>%
  select("Gene")
#write.csv(protein_list, file = "protein_list.csv")

nanoSPINS_RNA_long <- df_normalized_reads_nanoSPINS_new %>%
  pivot_longer(!gene, names_to = "SampleID", values_to = "cpm") %>%
  mutate(cpm = ifelse(cpm == 0, NA, cpm)) %>%
  mutate(SampleID = str_remove(SampleID, "_RNA"))

nanoSPINS_protein_long <- protein_input_750wide_norm_new %>%
  pivot_longer(!gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(SampleID = str_remove(SampleID, "_Prot"))
  
combined_corr_pair <- inner_join(nanoSPINS_protein_long, nanoSPINS_RNA_long) %>%
  filter(!is.na(SampleID) & !is.na(cpm) & !is.na(Intensity)) %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC",
    TRUE ~ "Other"
  ))

combined_corr_pair_C10 <- combined_corr_pair %>% 
  filter(Group == "C10") %>% 
  group_by(gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 3) %>%
  ungroup() %>%
  as_tibble()

combined_corr_pair_SVEC <- combined_corr_pair %>% 
  filter(Group == "SVEC") %>% 
  group_by(gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 4) %>%
  ungroup() %>%
  as_tibble()

combined_corr_pair_C10_nest <- group_by(combined_corr_pair_C10, gene) %>% nest()

combined_corr_pair_SVEC_nest <- group_by(combined_corr_pair_SVEC, gene) %>% nest()

cor_fun <- function(df) cor.test(df$Intensity, df$cpm, method = "pearson") %>% broom::tidy()
data_C10 <- combined_corr_pair_C10 %>% 
  select(gene) %>% 
  unique()
test_C10 <- mutate(combined_corr_pair_C10_nest, model = map(data, cor_fun))
corr_C10 <- dplyr::select(test_C10, -data) %>% 
  unnest(cols = c(model))
corr_C10 <- corr_C10 %>%
  ungroup() %>%
  mutate(FDR = p.adjust(`p.value`, method = "BH")) %>%
  mutate(Type2 = "mRNA-Protein_C10")

data_SVEC <- combined_corr_pair_SVEC %>% 
  select(gene) %>% 
  unique()
test_SVEC <- mutate(combined_corr_pair_SVEC_nest, model = map(data, cor_fun))
corr_SVEC <- dplyr::select(test_SVEC, -data) %>% 
  unnest(cols = c(model))
corr_SVEC <- corr_SVEC %>%
  ungroup() %>%
  mutate(FDR = p.adjust(`p.value`, method = "BH")) %>%
  mutate(Type2 = "mRNA-Protein_SVEC")

combined_pearson <- full_join(corr_C10, corr_SVEC)

combined_pearson %>%
  ggplot()+
  aes(x = estimate, fill = Type2)+
  geom_histogram(alpha = 0.6,
                 position = "identity",
                 bins = 75)+
  xlab("Pearson Correlation (r)")+
  ylab("Gene-Protein Pair (n)")+
  theme_minimal(base_size = 16)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())+
  scale_fill_brewer(name = "", palette = 6 ,type = "qual",
                    direction = -1)+
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal") +
  geom_vline(xintercept = 0.026, col = "black") +
  geom_vline(xintercept = 0.066, col = "purple")

write.csv(combined_pearson, file = "combined_pearson.csv")

#ggsave("mRNA-protein_correlations.png", width = 8, height = 6, bg = "white")

#correlation_results_C10 <- combined_corr_pair_C10 %>% select(gene, Intensity, cpm) %>% group_by(gene) %>% filter(sd(Intensity) != 0 & sd(cpm) != 0) %>% mutate(correlation = cor(Intensity, cpm, method = "pearson", use = "pairwise.complete.obs")) %>% select(gene, correlation) %>% distinct()

#combined_corr_pair_C10 %>% select(SampleID, gene, Intensity, cpm) %>% group_by(gene) %>% filter(sd(Intensity) != 0 & sd(cpm) != 0) %>% mutate(number = n()) %>% mutate(correlation = cor(Intensity, cpm, method = "pearson", use = "pairwise.complete.obs")) %>% mutate(pvalue = ifelse(n() > 2, cor.test(Intensity,  cpm, method = "pearson")$p.value, NA)) %>%  select(gene, correlation, pvalue) %>%  distinct() %>%  mutate(color_group = case_when(pvalue < 0.05 ~ "red")) %>% ggplot(aes(x = correlation,  y = -log(pvalue), color = factor(color_group))) + geom_point(na.rm = TRUE, ) +  geom_hline(yintercept = -log(0.05), linetype = "dashed",  color = "red") + xlab("correlation coefficients") + ylab("-log(p-value)") + scale_x_continuous(limits = c(-1, 1), breaks = seq(from = -1, to = 1, by = 0.2)) + scale_y_continuous(limits = c(0, 20),  breaks = seq(from = 0, to = 20, by = 5)) +  scale_color_manual(values = c("red" = "red")) + theme_minimal(base_size = 14)

genes_and_sampleID <- df_normalized_reads_nanoSPINS_new_0toNA %>% 
  pivot_longer(!gene, names_to = "SampleID", values_to = "cpm") %>% 
  filter(!is.na(cpm)) %>% 
  select(gene, SampleID) %>% 
  mutate(SampleID = str_remove(SampleID, "_RNA"))

proteins_and_sampleID <- protein_input_750wide_norm_new %>% 
  pivot_longer(!gene, names_to = "SampleID", values_to = "Intensity") %>% 
  filter(!is.na(Intensity)) %>% 
  select(gene, SampleID) %>% 
  mutate(SampleID = str_remove(SampleID, "_Prot"))

overlapping_genes_and_proteins <- genes_and_sampleID %>% 
  inner_join(proteins_and_sampleID, by = c("gene", "SampleID")) %>% 
  group_by(SampleID) %>% 
  mutate(n = n()) %>% 
  select(SampleID, n) %>% 
  unique() %>% 
  ungroup() %>% 
  mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC")) %>% 
  group_by(celltype) %>% 
  summarise(mean(n))

#percent mitochondrial reads scRNAseq datasets (cut-offs qualified)
RNA_input_all <- read.delim("SB22_12_PNNL_21.umicount.inex.all.tsv", 
                        row.names = "Gene") %>% 
  dplyr::select(one_of(RNA_selected_samples_new$Cell.Barcode))

RNA_input_all <- read.csv("RNA_input_all.csv", row.names = 1) %>%
  rownames_to_column(var = "gene_id") %>%
  inner_join(., convert) %>% 
  mutate(Avg = rowSums(.[2:303]/302)) %>% 
  group_by(gene_name) %>% 
  slice_max(Avg, n = 1, with_ties = F ) %>%
  ungroup() %>%
  column_to_rownames(var = "gene_name") %>% 
  dplyr::select( -gene_id, -Avg)
scRNAseq_rawdata_all <- RNA_input_all %>%
  as.sparse()
C10SVEC_rawdata_all <- CreateSeuratObject(counts = scRNAseq_rawdata_all, 
                                      min.cells = 5, 
                                      min.features = 1000)
C10SVEC_rawdata_all <- PercentageFeatureSet(C10SVEC_rawdata_all, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")
C10SVEC_rawdata_all$celltype <- RNA_selected_samples_new$Cell.type
VlnPlot(C10SVEC_rawdata_all, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        group.by = "celltype", 
        pt.size = 1)
percent_mt <- as.data.frame(C10SVEC_rawdata_all$percent.mt) %>% rownames_to_column(var = "SampleID") %>% mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC"))
#write.csv(percent_mt, file = "percent_mt.csv")
