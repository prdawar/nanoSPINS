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
library(gprofiler2)

getwd()
setwd("C:/Users/dawa726/OneDrive - PNNL/Desktop/nanoSPINS_Pranav")

load("nanoSPINS_proteome_and_transcriptome.RData")
#scProt-scRNA correlation analysis at sample level
protein_input_750wide_norm_new <- protein_input_750wide_norm %>% 
  setNames(paste0(names(.), "_Prot")) %>% 
  rownames_to_column(var = "gene") 

df_normalized_counts_1_nanoSPINS_protMatch <- df_normalized_counts_1 %>%
  select(one_of(colnames(protein_input_750wide_norm)))

protein_input_750wide_norm_RNAmatch <- protein_input_750wide_norm %>%
  select(one_of(colnames(df_normalized_counts_1_nanoSPINS_protMatch)))

df_normalized_counts_1_nanoSPINS_protMatch <- df_normalized_counts_1_nanoSPINS_protMatch %>%
  setNames(paste0(names(.), "_RNA")) %>%
  rownames_to_column(var = "gene")

protein_input_750wide_norm_RNAmatch <- protein_input_750wide_norm %>%
  setNames(paste0(names(.), "_Prot")) %>%
  rownames_to_column(var = "gene")

combined_new <- full_join(protein_input_750wide_norm_RNAmatch, df_normalized_counts_1_nanoSPINS_protMatch) %>% 
  column_to_rownames(var = "gene")

correlation_matrix_RNA_prot <- cor(as.data.frame(lapply(combined_new, as.numeric)), method = "pearson", use = "pairwise.complete.obs")
hc_RNA_prot <- hclust(dist(correlation_matrix_RNA_prot))
reordered_correlation_matrix_RNA_prot <- correlation_matrix_RNA_prot[hc_RNA_prot$labels, hc_RNA_prot$labels]
png("correlation_plot_crossmodality_RNA_prot.png", width = 3000, height = 3000, bg = "white", res = 300)
corrplot::corrplot(reordered_correlation_matrix_RNA_prot, 
                   method = 'shade', 
                   diag = FALSE, 
                   tl.cex = 0.2, 
                   col.lim = c(-1, 1)) %>% 
  corrplot::corrRect(c(1,70,161,210,282), lwd = 2, col = "black")
dev.off()

correlation_long <- melt(correlation_matrix_RNA_prot, 
                         varnames = c("Var1", "Var2"), 
                         value.name = "Correlation")
correlation_long %>%
  filter(grepl("_Prot$", Var1)) %>%
  filter(grepl("_RNA$", Var2)) %>%
  summarize(mean = mean(Correlation))

combined_overlapping <- inner_join(protein_input_750wide_norm_RNAmatch, df_normalized_counts_1_nanoSPINS_protMatch)

cpm_distribution <- df_normalized_counts_1_nanoSPINS_protMatch %>%
  pivot_longer(!gene, names_to = "SampleID", 
               values_to = "normalized_reads") %>%
  mutate(celltype = case_when(grepl("C10_[A-Z]", 
                                    SampleID) ~ "C10",
                              grepl("SVEC_[A-Z]", 
                                    SampleID) ~ "SVEC")) %>%
  filter(normalized_reads != "NA") %>%
  group_by(gene) %>%
  mutate(mean_cpm = mean(normalized_reads)) %>%
  select(gene, mean_cpm) %>% 
  unique() %>% 
  ungroup() %>%
  mutate(overlapping = ifelse(gene %in% combined_overlapping$gene, "Yes", "No"))

test4 <- nanoSPINS_proteins_list
test5 <- df_normalized_counts_1_nanoSPINS_protMatch %>% 
  select(gene) %>% 
  unique()
write.csv(test4, "protein_list.csv")
write.csv(test5, "gene_list.csv")

cpm_distribution %>% group_by(overlapping) %>% summarize(median = median(mean_cpm))
t.test(mean_cpm ~ overlapping, data = cpm_distribution)
cpm_distribution %>% 
  ggplot()+
  aes(x = overlapping, y = mean_cpm, fill = overlapping)+
  ylab("mean CPM") +
  xlab("") + 
  geom_violin() +
  scale_y_continuous(limits = c(0,20),
                     breaks = seq(0,20, by = 5)) +
  stat_summary(fun = "median_sd1",
               geom = "point",
               color = "black") +
  scale_fill_manual(values = c("Khaki", "RosyBrown")) +
  theme_minimal(base_size = 14)
#ggsave("median_cpm_values.png", width = 6, height = 4, bg = "white")


ggplot(cpm_distribution, 
       aes(x = mean_cpm, 
           fill = overlapping)) +
  geom_histogram(bins = 175) +
  #facet_wrap(~ celltype) +
  scale_y_continuous(limits = c(0, 1200), 
                     breaks = seq(0, 1200, by = 300)) +
  scale_x_continuous(limits = c(4, 16), 
                     breaks = seq(4,16, by = 2)) +
  labs(x = "mean log10(CPM)",
       y = "Genes (n)") +
  scale_fill_manual(values = c("Khaki", "RosyBrown")) +
#, fill = "Overlapping") +
  theme_minimal(base_size = 16)
#ggsave("histogram_overalpping_celltype.png", width = 6, height = 4, bg = "white")

nanoSPINS_RNA_long <- df_normalized_counts_1_nanoSPINS_protMatch %>%
  pivot_longer(!gene, names_to = "SampleID", values_to = "cpm") %>%
  mutate(cpm = ifelse(cpm == 0, NA, cpm)) %>%
  mutate(SampleID = str_remove(SampleID, "_RNA"))

nanoSPINS_protein_long <- protein_input_750wide_norm_RNAmatch %>%
  pivot_longer(!gene, names_to = "SampleID", values_to = "Intensity") %>%
  mutate(SampleID = str_remove(SampleID, "_Prot"))

combined_corr_pair <- inner_join(nanoSPINS_protein_long, nanoSPINS_RNA_long) %>%
  filter(!is.na(SampleID) & !is.na(cpm) & !is.na(Intensity)) %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC",
    TRUE ~ "Other"))

combined_corr_pair_C10 <- combined_corr_pair %>% filter(Group == "C10") %>% group_by(gene) %>% add_count(name = "Obs") %>% filter(Obs >= 4) %>% ungroup() %>% as_tibble()
combined_corr_pair_SVEC <- combined_corr_pair %>% filter(Group == "SVEC") %>% group_by(gene) %>% add_count(name = "Obs") %>% filter(Obs >= 4) %>% ungroup() %>% as_tibble()
combined_corr_pair_C10_nest <- group_by(combined_corr_pair_C10, gene) %>% nest()
combined_corr_pair_SVEC_nest <- group_by(combined_corr_pair_SVEC, gene) %>% nest()

cor_fun <- function(df) cor.test(df$Intensity, 
                                 df$cpm, 
                                 method = "pearson") %>% 
  broom::tidy()

data_C10 <- combined_corr_pair_C10 %>% select(gene) %>% unique()
test_C10 <- mutate(combined_corr_pair_C10_nest, model = map(data, cor_fun))
corr_C10 <- dplyr::select(test_C10, -data) %>% unnest(cols = c(model))
corr_C10 <- corr_C10 %>% ungroup() %>% mutate(FDR = p.adjust(`p.value`, method = "BH")) %>% mutate(Type2 = "mRNA-Protein_C10")
corr_C10 %>% summarize(mean(estimate))
data_SVEC <- combined_corr_pair_SVEC %>% select(gene) %>% unique()
test_SVEC <- mutate(combined_corr_pair_SVEC_nest, model = map(data, cor_fun))
corr_SVEC <- dplyr::select(test_SVEC, -data) %>% unnest(cols = c(model))
corr_SVEC <- corr_SVEC %>% ungroup() %>% mutate(FDR = p.adjust(`p.value`, method = "BH")) %>% mutate(Type2 = "mRNA-Protein_SVEC")
#filter(p.value < 0.05 & estimate < -0.5)
corr_SVEC %>% summarize(mean(estimate))
combined_pearson <- full_join(corr_C10, corr_SVEC)

corr_all %>% select(estimate) %>% summarise(mean = median(estimate))

#combined_pearson %>% ggplot() + aes(x = estimate, fill = Type2) + geom_histogram(alpha = 0.6, position = "identity", bins = 50) + xlab("Pearson Correlation (r)") + ylab("Gene-Protein Pair (n)") + theme_minimal(base_size = 16) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank(), axis.line = element_line()) + scale_fill_brewer(name = "", palette = 6, type = "qual", direction = -1) + geom_vline(xintercept = 0.022, col = "black") + theme(legend.position = "bottom", legend.direction = "vertical", legend.box = "horizontal")
#ggsave("correlation_plot_RNA_prot.png", width = 6, height = 5, bg = "white")

#combined_pearson %>% mutate(group = ifelse(p.value < 0.05, "significant", "non-significant")) %>% ggplot(aes(x = estimate, y = -log10(p.value), color = factor(group))) + geom_point(na.rm = TRUE, size = 3, alpha = 0.5) + #geom_text_repel(data = . %>% filter(p.value < 0.05), aes(label = gene), vjust = -1.5, direction = "x", size = 5) + xlab("correlation coefficients") + ylab("-log10(p-value)") + facet_wrap(~Type2) + scale_x_continuous(limits = c(-1, 1), breaks = seq(from = -1, to = 1, by = 0.5)) + scale_y_continuous(limits = c(0, 8), breaks = seq(from = 0, to = 8, by = 2)) + scale_color_manual(values = c("significant" = "red", "non-significant" = "darkgrey"), name = "p-value") + theme_minimal(base_size = 18)
#ggsave("combined_pearson.png", width = 10, height = 5, bg = "white")

#combined_pearson_sig_genes <- combined_pearson %>% select(gene, estimate, p.value, FDR, Type2) %>% filter(p.value < 0.05 & abs(estimate) > 0.5)
#combined_pearson_sig_genes_positive <- combined_pearson_sig_genes %>% select(gene, estimate, p.value, FDR, Type2) %>% filter(p.value < 0.05 & estimate >= 0.5) %>% select(gene) %>% unique()
#combined_pearson_sig_genes_positive <- combined_pearson_sig_genes_positive$gene
#combined_pearson_sig_genes_negative <- combined_pearson_sig_genes %>% select(gene, estimate, p.value, FDR, Type2) %>% filter(p.value < 0.05 & estimate <= -0.5) %>% select(gene) %>% unique()
#combined_pearson_sig_genes_negative <- combined_pearson_sig_genes_negative$gene
#nanoSPINS_protein_list <-  protein_input_750long %>% select(Gene) %>% unique() %>% as.list()
#nanoSPINS_protein_list <- nanoSPINS_protein_list$Gene
#nanoSPINS_RNA_list <- df_normalized_counts_1_nanoSPINS_protMatch %>% select(gene) %>% unique()
#nanoSPINS_RNA_list <- nanoSPINS_RNA_list$gene
#gene_list <- combined_new %>% rownames_to_column(var = "gene") %>% select(gene)
#gene_list <- gene_list$gene  
#GO_combined_pearson_bg_gene_list_negative <- gost(query = combined_pearson_sig_genes_negative, organism = "mmusculus", ordered_query = F, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated", custom_bg = gene_list, numeric_ns = "", sources = c("GO:CC"), as_short_link = FALSE)
#flattened_GO_combined_pearson_bg_gene_list_negative <- as.data.frame(GO_combined_pearson_bg_gene_list_negative$result) %>% mutate(across(everything(), ~ map_chr(., toString)))
#flattened_GO_combined_pearson_bg_gene_list_negative <- flattened_GO_combined_pearson_bg_gene_list_negative %>% select(-evidence_codes, -intersection)
#flattened_GO_combined_pearson_bg_gene_list_negative_subset <- flattened_GO_combined_pearson_bg_gene_list_negative %>% select(term_id, p_value, term_name)
#GO_combined_pearson_bg_gene_list_positive <- gost(query = combined_pearson_sig_genes_positive, organism = "mmusculus", ordered_query = F, multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, measure_underrepresentation = FALSE, evcodes = TRUE, user_threshold = 0.05, correction_method = "fdr", domain_scope = "annotated",custom_bg = gene_list, numeric_ns = "", sources = c("GO:CC"), as_short_link = FALSE)
#flattened_GO_combined_pearson_bg_gene_list_positive <- as.data.frame(GO_combined_pearson_bg_gene_list_positive$result) %>% mutate(across(everything(), ~ map_chr(., toString)))
#flattened_GO_combined_pearson_bg_gene_list_positive <- flattened_GO_combined_pearson_bg_gene_list_positive %>% select(-evidence_codes, -intersection)
#flattened_GO_combined_pearson_bg_gene_list_positive_subset <- flattened_GO_combined_pearson_bg_gene_list_positive %>% select(term_id, p_value, term_name)

combined_corr_pair_all <- inner_join(nanoSPINS_protein_long, nanoSPINS_RNA_long) %>%
  filter(!is.na(SampleID) & !is.na(cpm) & !is.na(Intensity)) %>%
  mutate(Group = case_when(
    grepl("C10", SampleID) ~ "C10",
    grepl("SVEC", SampleID) ~ "SVEC",
    TRUE ~ "Other")) %>%
  group_by(gene) %>%
  add_count(name = "Obs") %>%
  filter(Obs >= 4) %>%
  ungroup() %>%
  as_tibble()
combined_corr_pair_all_nest <- group_by(combined_corr_pair_all, gene) %>%
  nest()
data_all <- combined_corr_pair_all_nest %>% 
  select(gene) %>% 
  unique()
test_all <- mutate(combined_corr_pair_all_nest, 
                   model = map(data, cor_fun))
corr_all <- dplyr::select(test_all, -data) %>% 
  unnest(cols = c(model))
corr_all <- corr_all %>%
  ungroup() %>%
  mutate(FDR = p.adjust(`p.value`, method = "BH"))

corr_all %>%
  mutate(group = ifelse(p.value < 0.01, "significant", "non-significant")) %>% 
  ggplot(aes(x = estimate, 
             y = -log10(p.value), 
             color = factor(group))) + 
  geom_point(na.rm = TRUE, size = 3, alpha = 0.5) + 
  geom_text_repel(data = . %>% filter(FDR < 0.05 & gene %in% features_rearrange), aes(label = gene[gene %in% features_rearrange]), vjust = -1.5, direction = "x", size = 5) +
  xlab("correlation coefficients") + 
  ylab("-log10(p-value)") +
  scale_x_continuous(limits = c(-1, 1), 
                     breaks = seq(from = -1, to = 1, by = 0.5)) +
  scale_y_continuous(limits = c(0, 15), 
                     breaks = seq(from = 0, to = 15, by = 5)) + 
  scale_color_manual(values = c("significant" = "red", 
                                "non-significant" = "darkgrey"), 
                     name = "p-value") + 
  theme_minimal(base_size = 16)
#ggsave("correlation_plot_RNA_prot_all.png", width = 7, height = 4, bg = "white")

corr_all_sig_genes <- corr_all %>% 
  select(gene, estimate, p.value, FDR) %>%
  filter(p.value < 0.05 & abs(estimate) > 0.5)

nanoSPINS_proteins_list <- nanoSPINS_protein_long %>% 
  select(gene) %>% 
  unique()

GO_corr_all <- gost(query = corr_all_sig_genes_positive$gene,
                                                  organism = "mmusculus", 
                                                  ordered_query = F, 
                                                  multi_query = FALSE, 
                                                  significant = TRUE, 
                                                  exclude_iea = FALSE,
                                                  measure_underrepresentation = FALSE,
                                                  evcodes = TRUE,
                                                  user_threshold = 0.05,
                                                  correction_method = "fdr",
                                                  domain_scope = "annotated", custom_bg = nanoSPINS_proteins_list$gene,
                                                  numeric_ns = "",
                                                  sources = c("GO:CC"),
                                                  as_short_link = FALSE)
flattened_GO_corr_all <- as.data.frame(GO_corr_all$result) %>%
  mutate(across(everything(), ~ map_chr(., toString)))
flattened_GO_corr_all <- flattened_GO_corr_all %>%
  select(-evidence_codes, -intersection)
flattened_GO_corr_all <- flattened_GO_corr_all %>% select(term_id, p_value, term_name)

corr_all_sig_genes_positive <- corr_all_sig_genes %>% 
  filter(estimate >= 0.5)
corr_all_sig_genes_negative <- corr_all_sig_genes %>% 
  filter(estimate <= -0.5)

###
#Data integration, multimodal analysis, and marker analysis
scRNAseq_rawdata <- RNA_input_all_new2.2.1 %>%
  select(one_of(protein_samples)) %>% 
  as.sparse()
C10SVEC_rawdata <- CreateSeuratObject(counts = scRNAseq_rawdata, 
                                      min.cells = 0, 
                                      min.features = 0)
common_samples <- colnames(RNA_input_all_new2.2.1 %>% 
                             select(one_of(protein_samples)))
scProt <- protein_input_impute %>% 
  select(one_of(common_samples)) %>% as.sparse()

all.equal(colnames(scRNAseq_rawdata), colnames(scProt))
prot_assay <- CreateAssayObject(counts = scProt)
C10SVEC_rawdata[["Prot"]] <- prot_assay
Assays(C10SVEC_rawdata)
DefaultAssay(C10SVEC_rawdata)
C10SVEC_rawdata <- PercentageFeatureSet(C10SVEC_rawdata, 
                                        pattern = "^mt-", 
                                        col.name = "percent.mt")

RNA_BCs_test <- RNA_BCs_test %>% arrange(Annotation)
C10SVEC_rawdata$celltype <- RNA_BCs_test$Cell.type
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
        dims = c(1, 3), 
        pt.size = 2, 
        alpha = 0.9, 
        group.by = "celltype")
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
#write.csv(C10SVEC_rawdata_C10markers_RNA, file = "C10SVEC_rawdata_C10markers_RNA_032325.csv")
#EVP1_RNAseq <- EnhancedVolcano(C10SVEC_rawdata_C10markers_RNA, lab = rownames(C10SVEC_rawdata_C10markers_RNA), x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, FCcutoff = 1, pointSize = 3.0, labSize = 5, selectLab = c(""), drawConnectors = TRUE, xlim = c(-10,10), ylim = c(0, 25), cutoffLineWidth = 1, cutoffLineCol = "grey", colConnectors = 'black', arrowheads = FALSE, widthConnectors = 1, parseLabels = TRUE, labCol = "black")
#EVP1_RNAseq
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
#write.csv(C10SVEC_rawdata_C10markers_prot, file = "C10SVEC_rawdata_C10markers_prot_032325.csv")
#EVP1_Prot <- EnhancedVolcano(C10SVEC_rawdata_C10markers_prot, lab = rownames(C10SVEC_rawdata_C10markers_prot), x = 'avg_log2FC', y = 'p_val_adj', pCutoff = 0.05, FCcutoff = 1, pointSize = 3.0, labSize = 5, selectLab = c(""), drawConnectors = TRUE, xlim = c(-5, 5), ylim = c(0, 25), cutoffLineWidth = 1, cutoffLineCol = "grey", colConnectors = 'black', arrowheads = FALSE, widthConnectors = 1, parseLabels = TRUE)
#EVP1_Prot
#EVP1_RNAseq + EVP1_Prot
#ggsave("volcanoplots_scRNAseq_plus_scProt.png", width = 16, height = 8)

C10SVEC_rawdata <- FindClusters(C10SVEC_rawdata, 
                                graph.name = "wsnn", 
                                algorithm = 1, 
                                resolution = 1, 
                                verbose = FALSE)

C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           nn.name = "weighted.nn", 
                           reduction.name = "wnn.umap", 
                           reduction.key = "wnnUMAP_")

p1 <- DimPlot(C10SVEC_rawdata, 
              group.by = "celltype",
              reduction = 'wnn.umap',
              label = TRUE, 
              repel = TRUE,
              label.size = 0, 
              pt.size = 3, 
              alpha = 0.75) + NoLegend()
p1 + theme_cowplot(font_size = 16)
ggsave("wnn.png", width = 6, height = 4)

C10SVEC_rawdata <- RunUMAP(C10SVEC_rawdata, 
                           reduction = 'pca', 
                           dims = 1:5, 
                           assay = 'RNA',
                           reduction.name = 'rna.umap', 
                           reduction.key = 'rnaUMAP_')

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
              repel = TRUE, label.size = 5) + 
  NoLegend()

features_rearrange <- c("Bst2", "H2-D1", "H2-K1", "Hmgb2", "Ifit3", "Isg15", "Lgals3", "Stmn1", "Actg1", "Col3a1", "Fn1", "Fscn1", "Gsta4", "Manf", "S100a4", "Thbs1", "Vim")

#feature_test <- c("Actb", "Actg1")
#RidgePlot(C10SVEC_rawdata, features = feature_test, ncol = 2)
#FeaturePlot(C10SVEC_rawdata, features = feature_test)
#DoHeatmap(C10SVEC_rawdata, features = feature_test, size = 3, group.by = "celltype")
#DotPlot(C10SVEC_rawdata, features = features, split.by = "celltype") + RotatedAxis()

C10SVEC_rawdata_C10markers_RNA_features <- C10SVEC_rawdata_C10markers_RNA %>% 
  rownames_to_column("Gene") %>% filter(Gene %in% features_rearrange)
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
  rownames_to_column("Gene") %>% filter(Gene %in% features_rearrange)
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
#ggsave("heatmap_integrative_analysis_032425.png", width = 14, height = 5)

#percent mitochondrial reads scRNAseq datasets (cut-offs qualified)
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
                                          min.cells = 0, 
                                          min.features = 0)
C10SVEC_rawdata_all <- PercentageFeatureSet(C10SVEC_rawdata_all, 
                                            pattern = "^mt-", 
                                            col.name = "percent.mt")
percent_mt <- as.data.frame(C10SVEC_rawdata_all$percent.mt) %>% rownames_to_column(var = "SampleID") %>% mutate(celltype = case_when(str_detect(SampleID, "C10") ~ "C10", str_detect(SampleID, "SVEC") ~ "SVEC"))
#write.csv(percent_mt, file = "percent_mt_032425.csv")


