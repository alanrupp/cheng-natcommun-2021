# - Analysis of count data from Calcr NTS/AP TRAP-seq -------------------------

library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)

# - Read in count data and metadata -------------------------------------------
counts <- read_tsv("data/GSE176202_counts.tsv.gz")
metadata <- read_csv("data/GSE176202_metadata.csv.gz")
metadata$Cells <- factor(metadata$Cells, levels = c("Sup", "Bead"))

# save gene info for later and work with gene_id
genes <- select(counts, starts_with("gene"))

# convert count data into matrix with genes as rownames
counts <- select(counts, gene_id, starts_with("Sample")) %>%
  as.data.frame() %>%
  column_to_rownames("gene_id")

# convert metadata into matrix with samples as rownames
metadata <- as.data.frame(metadata) %>%
  column_to_rownames("Sample_ID")

# put metadata and counts in same order
counts <- counts[, rownames(metadata)]

# - Enrichment ----------------------------------------------------------------
# looking for enrichment in Bead / Sup
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts[, rownames(metadata)],
  colData = metadata,
  ~ Run_ID + Cells
) %>%
  DESeq2::DESeq()

# get DE between Bead and Sup and add back gene_name
enrichment <- DESeq2::results(dds) %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(genes, by = "gene_id")

# - Heatmap of enriched / de-enriched genes -----------------------------------
# get genes that are significantly different between bead and sup
enriched <- filter(enrichment, padj < 0.05) %>%
  filter(!gene_name %in% c("Cre", "Gfp_L10a")) %>%
  arrange(desc(log2FoldChange)) %>%
  pull(gene_id)

# get log2CPM z-scores for each sample
expression <- DESeq2::fpm(dds) %>% as.data.frame()
log2cpm <- apply(expression, c(1,2), function(x) log2(x+1))
zscores <- apply(log2cpm, 1, function(x) (x-mean(x))/sd(x)) %>% t()

# add metadata and organize for plotting
zscores <- as.data.frame(zscores) %>% rownames_to_column("gene_id") %>%
  filter(gene_id %in% enriched) %>%
  left_join(genes, by = "gene_id") %>%
  mutate(gene_id = factor(gene_id, levels = enriched)) %>%
  arrange(gene_id) %>%
  mutate("gene_name" = factor(gene_name, levels = unique(gene_name))) %>%
  pivot_longer(starts_with("Sample"), names_to = "Sample", values_to = "z") %>%
  left_join(rownames_to_column(metadata, "Sample"), by = "Sample") %>%
  arrange(Cells, Sample) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# plot heatmap
p <- ggplot(zscores, aes(x = Sample, y = as.numeric(gene_name), fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#798234", mid = "#fdfbe4", high = "#d46780",
                       name = "z-score") +
  annotate("text", x = 1.5, y = -0.2, label = "Sup", size = 3) +
  annotate("text", x = 3.5, y = -0.2, label = "Bead", size = 3) +
  theme_void() +
  scale_y_continuous(trans = "reverse",
                     breaks = length(enriched):1,
                     labels = levels(zscores$gene_name)[length(enriched):1]) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(hjust = 1, color = "black", size = 8,
                                   face = "italic"),
        legend.key.width = unit(0.02, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))

# - Saving for the manuscript -------------------------------------------------
# figure
if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/heatmap.pdf", p, width = 1.8, height = 4, units = "in", dpi = 400)

# table of raw data from plot
metadata$Pair <- recode(metadata$Run_ID, "Run_1905" = "A", "Run_2294" = "B")
metadata$Sample <- paste(metadata$Cells, metadata$Pair, sep = "_")
zscores$ID <- metadata$Sample[match(zscores$Sample, rownames(metadata))]
zscores$ID <- factor(zscores$ID, levels = c("Sup_A", "Sup_B", "Bead_A", "Bead_B"))
zscores %>% select(gene_name, z, ID) %>%
  pivot_wider(names_from = "ID", values_from = "z") %>%
  rename("Gene" = "gene_name") %>%
  mutate_if(is.numeric, ~ round(.x, 4)) %>%
  writexl::write_xlsx("figures/zscores.xlsx")
