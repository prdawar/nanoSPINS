#For GO analysis Figure 4A

protein_list <- protein_input %>%
  filter(!is.na(Intensity)) %>%
  filter(Group != "Blank") %>%
  select(Gene) %>% unique()
protein_list <- protein_list$Gene

gene_list <- unique_genes_scRNAseq_NS_1%>%
  select(gene) %>% 
  unique()
gene_list <- gene_list$gene

library(gprofiler2)
GO_protein_list_0.01 <- gost(query = protein_list,
                        organism = "mmusculus", 
                        ordered_query = F, 
                        multi_query = FALSE, 
                        significant = TRUE, 
                        exclude_iea = FALSE, 
                        measure_underrepresentation = FALSE, 
                        evcodes = TRUE, 
                        user_threshold = 0.01, 
                        correction_method = "fdr", 
                        domain_scope = "annotated", 
                        numeric_ns = "", 
                        sources = c("GO:CC"), 
                        as_short_link = FALSE)
flattened_GO_protein_list_0.01 <- as.data.frame(GO_protein_list_0.01$result) %>%
  mutate(across(everything(), ~ map_chr(., toString)))
flattened_GO_protein_list_0.01 <- flattened_GO_protein_list_0.01 %>%
  select(-evidence_codes, -intersection)
flattened_GO_protein_list_subset_0.01 <- flattened_GO_protein_list_0.01 %>% 
  select(term_id, p_value, term_name)

nuclear_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "nucleolus|nucleoplasm|nuclear|nucleus|splice|snRNP|nuclease|DNA|chromatin|chromosom|polymerase|histone|nucleo|transcription|RNA")) %>%
  distinct() %>%
  mutate(Component = "Nuclear")

cytoplasmic_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "cytoplasm|cytosol")) %>%
  distinct() %>%
  mutate(Component = "Cytoplasmic")

mitochondrial_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "mitochondria|mitochondrion")) %>%
  distinct() %>%
  mutate(Component = "Mitochondrial")

plasmamembrane_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "membrane")) %>%
  distinct() %>%
  mutate(Component = "Plasma membrane")

ER_golgi_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>% 
  filter(str_detect(term_name, "Golgi|endoplasmic")) %>%
  distinct() %>%
  mutate(Component = "ER and Golgi")

extraCel_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(grepl("extracellular", term_name)) %>%
  distinct() %>%
  mutate(Component = "Extracellular")

vacuole_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(grepl("vacuol", term_name)) %>%
  distinct() %>%
  mutate(Component = "Vacuole")

cytoskeleton_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "cytoskeleton|microtubule|actin")) %>%
  distinct() %>%
  mutate(Component = "Cytoskeleton")

synapse_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "synapse|synaptic|neuro|dendrite|dendritic")) %>%
  distinct() %>%
  mutate(Component = "synapse")

vesicle_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "vesicle|vesicular")) %>%
  distinct() %>%
  mutate(Component = "vesicle")

ribosome_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(str_detect(term_name, "ribosome|ribosomal")) %>%
  distinct() %>%
  mutate(Component = "ribosome")

misc_protein_0.01 <- flattened_GO_protein_list_subset_0.01 %>%
  filter(!term_id %in% nuclear_protein_0.01$term_id) %>%
  filter(!term_id %in% cytoplasmic_protein_0.01$term_id) %>%
  filter(!term_id %in% mitochondrial_protein_0.01$term_id) %>%
  filter(!term_id %in% plasmamembrane_protein_0.01$term_id) %>%
  filter(!term_id %in% ER_golgi_protein_0.01$term_id) %>%
  filter(!term_id %in% extraCel_protein_0.01$term_id) %>%
  filter(!term_id %in% cytoskeleton_protein_0.01$term_id) %>%
  filter(!term_id %in% vacuole_protein_0.01$term_id) %>%
  filter(!term_id %in% synapse_protein_0.01$term_id) %>%
  filter(!term_id %in% vesicle_protein_0.01$term_id)  %>%
  filter(!term_id %in% ribosome_protein_0.01$term_id)  %>%
  distinct() %>%
  mutate(Component = "Misc")

protein_component <- do.call("rbind",list(nuclear_protein_0.01,cytoplasmic_protein_0.01, mitochondrial_protein_0.01,plasmamembrane_protein_0.01, ER_golgi_protein_0.01,extraCel_protein_0.01, cytoskeleton_protein_0.01, vacuole_protein_0.01, synapse_protein_0.01, vesicle_protein_0.01, ribosome_protein_0.01, misc_protein_0.01)) %>%
  mutate(Component = as.factor(Component)) %>%
  mutate(Component = fct_relevel(Component, "Nuclear", "Cytoplasmic","Plasma membrane","Mitochondria","ER & Golgi Appratus", "Vacuole", "Ribosome", "Vesicle", "Synapse", "Extracellular"))
cols <- c("#fbb4ae", "#fddaec", "#b3cde3", "#ccebc5",
          "#decbe4", "#fed9a6", "#aef5fb", "#ffffcc", "#e5d8bd", "#aefbb4",  "#fbaecf", "#ddbf70" )

# Calculate the count for each component
component_counts_protein <- protein_component %>%
  group_by(Component) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

# Ensure the Component factor levels are ordered as desired
component_counts_protein <- component_counts_protein %>%
  mutate(Component = factor(Component, levels = c("Nuclear", "Cytoplasmic", "Plasma membrane", "Vacuole", "ER and Golgi", "Mitochondrial", "ribosome", "synapse", "vesicle", "Extracellular", "Cytoskeleton", "Misc")))

# Compose the pie chart
ggplot(component_counts_protein, aes(x = "", y = Count, fill = Component)) +
  geom_bar(stat = "identity", width = 1, show.legend = TRUE, color = NA) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cols) +
  theme_void() +
  theme(text = element_text(family = "Helvetica")) +
  labs(fill = "Component")
ggsave("nanoSPINS_protein_CCs.png", width = 8, height = 6, bg = "white")


GO_NS_0.01 <- gost(query = gene_list,
              organism = "mmusculus", 
              ordered_query = F, 
              multi_query = FALSE, 
              significant = TRUE, 
              exclude_iea = FALSE, 
              measure_underrepresentation = FALSE, 
              evcodes = TRUE, 
              user_threshold = 0.01, 
              correction_method = "fdr", 
              domain_scope = "annotated",
              numeric_ns = "", 
              sources = c("GO:CC"), 
              as_short_link = FALSE)
flattened_GO_NS_0.01 <- as.data.frame(GO_NS_0.01$result) %>%
  mutate(across(everything(), ~ map_chr(., toString)))
flattened_GO_NS_0.01 <- flattened_GO_NS_0.01 %>%
  select(-evidence_codes, -intersection)
flattened_GO_NS_subset_0.01 <- flattened_GO_NS_0.01 %>% select(term_id, p_value, term_name)

nuclear_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "nucleolus|nucleoplasm|nuclear|nucleus|splice|snRNP|nuclease|DNA|chromatin|chromosom|polymerase|histone|nucleo|transcription|RNA")) %>%
  distinct() %>%
  mutate(Component = "Nuclear")

cytoplasmic_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "cytoplasm|cytosol")) %>%
  distinct() %>%
  mutate(Component = "Cytoplasmic")

mitochondrial_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "mitochondria|mitochondrion")) %>%
  distinct() %>%
  mutate(Component = "Mitochondrial")

plasmamembrane_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "membrane")) %>%
  distinct() %>%
  mutate(Component = "Plasma membrane")

ER_golgi_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>% 
  filter(str_detect(term_name, "Golgi|endoplasmic")) %>%
  distinct() %>%
  mutate(Component = "ER and Golgi")

extraCel_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(grepl("extracellular", term_name)) %>%
  distinct() %>%
  mutate(Component = "Extracellular")

vacuole_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(grepl("vacuol", term_name)) %>%
  distinct() %>%
  mutate(Component = "Vacuole")

cytoskeleton_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "cytoskeleton|microtubule|actin")) %>%
  distinct() %>%
  mutate(Component = "Cytoskeleton")

synapse_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "synapse|synaptic|neuro|dendrite|dendritic")) %>%
  distinct() %>%
  mutate(Component = "synapse")

vesicle_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "vesicle|vesicular")) %>%
  distinct() %>%
  mutate(Component = "vesicle")

ribosome_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(str_detect(term_name, "ribosome|ribosomal")) %>%
  distinct() %>%
  mutate(Component = "ribosome")

misc_RNA_0.01 <- flattened_GO_NS_subset_0.01 %>%
  filter(!term_id %in% nuclear_RNA_0.01$term_id) %>%
  filter(!term_id %in% cytoplasmic_RNA_0.01$term_id) %>%
  filter(!term_id %in% mitochondrial_RNA_0.01$term_id) %>%
  filter(!term_id %in% plasmamembrane_RNA_0.01$term_id) %>%
  filter(!term_id %in% ER_golgi_RNA_0.01$term_id) %>%
  filter(!term_id %in% extraCel_RNA_0.01$term_id) %>%
  filter(!term_id %in% cytoskeleton_RNA_0.01$term_id) %>%
  filter(!term_id %in% vacuole_RNA_0.01$term_id) %>%
  filter(!term_id %in% synapse_RNA_0.01$term_id) %>%
  filter(!term_id %in% vesicle_RNA_0.01$term_id)  %>%
  filter(!term_id %in% ribosome_RNA_0.01$term_id)  %>%
  distinct() %>%
  mutate(Component = "Misc")

RNA_component <- do.call("rbind",list(nuclear_RNA_0.01,cytoplasmic_RNA_0.01, mitochondrial_RNA_0.01,plasmamembrane_RNA_0.01, ER_golgi_RNA_0.01,extraCel_RNA_0.01, cytoskeleton_RNA_0.01, vacuole_RNA_0.01, synapse_RNA_0.01, vesicle_RNA_0.01, ribosome_RNA_0.01, misc_RNA_0.01)) %>%
  mutate(Component = as.factor(Component)) %>%
  mutate(Component = fct_relevel(Component, "Nuclear", "Cytoplasmic","Plasma membrane","Mitochondria","ER & Golgi Appratus", "Vacuole", "Ribosome", "Vesicle", "Synapse", "Extracellular"))
cols <- c("#fbb4ae", "#fddaec", "#b3cde3", "#ccebc5",
          "#decbe4", "#fed9a6", "#aef5fb", "#ffffcc", "#e5d8bd", "#aefbb4",  "#fbaecf", "#ddbf70" )


png("RNA_component_pie_chart.png", width = 800, height = 800)
PieChart(Component, data = RNA_component,
         stat = "%",
         fill = cols,
         lwd = 2,
         lty = 1,
         values_size = 0.8,
         values_color = "black",
         color = "black",
         main = NULL,
         width = 10,
         heigth = 10)+
  theme(text=element_text(family="Helvetica"))
dev.off()

# Calculate the count for each component
component_counts <- RNA_component %>%
  group_by(Component) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = (Count / sum(Count)) * 100)

# Ensure the Component factor levels are ordered as desired
component_counts <- component_counts %>%
  mutate(Component = factor(Component, levels = c("Nuclear", "Cytoplasmic", "Plasma membrane", "Vacuole", "ER and Golgi", "Mitochondrial", "ribosome", "synapse", "vesicle", "Extracellular", "Cytoskeleton", "Misc")))

# Compose the pie chart
ggplot(component_counts, aes(x = "", y = Count, fill = Component)) +
  geom_bar(stat = "identity", width = 1, show.legend = TRUE, color = NA) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cols) +
  theme_void() +
  theme(text = element_text(family = "Helvetica")) +
  labs(fill = "Component")
ggsave("nanoSPINS_RNA_CCs.png", width = 8, height = 6, bg = "white")
