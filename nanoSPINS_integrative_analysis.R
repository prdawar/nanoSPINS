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
