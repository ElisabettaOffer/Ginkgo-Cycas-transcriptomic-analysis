# ==============================================================================
# SUPPLEMENTARY SCRIPT S2
# Targeted comparative transcriptomic analysis of Cycas circinalis male cone development without Mfuzz clustering
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

# Input file generated from the Novogene de novo transcriptome workflow.
# The table is expected to contain raw counts and annotation fields in a single file.
input_file <- "Cycas_circinalis_Novogene_counts_annotations.csv"

# Output directory.
output_dir <- "Cycas_Supplementary_Script_S2_Output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Statistical thresholds.
padj_cutoff <- 0.005
lfc_cutoff <- 1.5
min_count <- 10
min_samples <- 3

# Mfuzz clustering is intentionally not performed for Cycas because the sampled
# developmental stages are not exact one-to-one equivalents of the Ginkgo M1-M2-M3 series.

# Candidate heatmap settings.
max_genes_per_heatmap_category <- 4
include_leaf_in_candidate_heatmap <- TRUE

# If the expanded GO-category RDS file generated for the Ginkgo analysis is present,
# it will be reused to enforce direct comparability between Ginkgo and Cycas.
ginkgo_go_categories_rds <- "GO_categories_expanded_list.rds"

# Annotation-supported AGP/FLA-like candidate search.
run_annotation_supported_agp_screen <- TRUE
make_agp_fla_heatmap <- FALSE

# ------------------------------------------------------------------------------
# 1. Environment setup and package loading
# ------------------------------------------------------------------------------

cran_packages <- c(
  "ggplot2", "ggrepel", "pheatmap", "dplyr", "reshape2",
  "tibble", "stringr", "tidyr", "readr", "RColorBrewer", "data.table",
  "ontologyIndex", "openxlsx"
)

bioc_packages <- c("DESeq2", "Mfuzz", "Biobase")

install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (isTRUE(bioc)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

invisible(lapply(cran_packages, install_if_missing, bioc = FALSE))
invisible(lapply(bioc_packages, install_if_missing, bioc = TRUE))

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(reshape2)
library(tibble)
library(stringr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(data.table)
library(ontologyIndex)
library(openxlsx)

# ------------------------------------------------------------------------------
# 2. Data import and organization
# ------------------------------------------------------------------------------

raw_data <- data.table::fread(input_file, data.table = FALSE, check.names = FALSE)

# Standardize column names for stable downstream access.
names(raw_data) <- make.names(names(raw_data), unique = TRUE)

if (!"gene_id" %in% names(raw_data)) {
  stop("The input table must contain a column named 'gene_id'.")
}

sample_cols <- c(
  paste0("Male_S1_", 1:5),
  paste0("Male_S2_", 1:5),
  paste0("Male_S3_", 1:5),
  paste0("Leaf_", 1:5)
)

missing_sample_cols <- setdiff(sample_cols, names(raw_data))
if (length(missing_sample_cols) > 0) {
  stop("Missing expected sample columns: ", paste(missing_sample_cols, collapse = ", "))
}

count_matrix <- raw_data %>%
  dplyr::select(all_of(c("gene_id", sample_cols))) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

mode(count_matrix) <- "numeric"
count_matrix[is.na(count_matrix)] <- 0
count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

sample_metadata <- tibble(
  sample = colnames(count_matrix),
  condition = dplyr::case_when(
    str_detect(sample, "^Male_S1_") ~ "C1",
    str_detect(sample, "^Male_S2_") ~ "C2",
    str_detect(sample, "^Male_S3_") ~ "C3",
    str_detect(sample, "^Leaf_") ~ "Leaf",
    TRUE ~ NA_character_
  ),
  tissue = ifelse(condition == "Leaf", "Leaf", "Microsporophyll"),
  replicate = as.integer(str_extract(sample, "\\d+$"))
) %>%
  mutate(
    condition = factor(condition, levels = c("C1", "C2", "C3", "Leaf")),
    tissue = factor(tissue, levels = c("Microsporophyll", "Leaf"))
  ) %>%
  as.data.frame()

rownames(sample_metadata) <- sample_metadata$sample
write_csv(sample_metadata, file.path(output_dir, "S2_sample_metadata.csv"))

annotation_cols <- setdiff(names(raw_data), sample_cols)

annotation_table <- raw_data %>%
  dplyr::select(all_of(annotation_cols)) %>%
  distinct(gene_id, .keep_all = TRUE)

safe_column <- function(df, colname) {
  if (colname %in% names(df)) {
    as.character(df[[colname]])
  } else {
    rep(NA_character_, nrow(df))
  }
}

annotation_table <- annotation_table %>%
  mutate(
    best_description = dplyr::coalesce(
      safe_column(., "Swissprot.Description"),
      safe_column(., "NR.Description"),
      safe_column(., "PFAM.description"),
      safe_column(., "KO.Description"),
      safe_column(., "KOG.Description")
    ),
    annotation_search_text = paste(
      safe_column(., "Swissprot.ID"),
      safe_column(., "Swissprot.Description"),
      safe_column(., "NR.ID"),
      safe_column(., "NR.Description"),
      safe_column(., "PFAM.ID"),
      safe_column(., "PFAM.description"),
      safe_column(., "KO.ID"),
      safe_column(., "KO.Name"),
      safe_column(., "KO.Description"),
      safe_column(., "KOG.ID"),
      safe_column(., "KOG.Description"),
      sep = " | "
    ),
    annotation_search_text_upper = toupper(annotation_search_text)
  )

write_csv(annotation_table, file.path(output_dir, "S2_annotation_table_cleaned.csv"))


# ------------------------------------------------------------------------------
# 2B. UniProt-based homolog annotation and GO retrieval
# ------------------------------------------------------------------------------

# Novogene homology identifiers are used only as accessions to query UniProt.
# Common gene names, protein names, and GO identifiers are downloaded from UniProt
# to avoid relying on the GO annotations directly provided in the Novogene table.

uniprot_dir <- file.path(output_dir, "UniProt_Homolog_Annotation")
dir.create(uniprot_dir, showWarnings = FALSE, recursive = TRUE)

uniprot_source_columns <- intersect(c("Swissprot.ID", "NR.ID"), names(annotation_table))
if (length(uniprot_source_columns) == 0) {
  stop("No candidate UniProt accession source columns were found. Expected at least 'Swissprot ID' or 'NR ID' in the Novogene table.")
}

extract_uniprot_accessions <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  # UniProt accession regex supporting standard 6- and 10-character accessions.
  pattern <- "([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]([A-Z][A-Z0-9]{2}[0-9])?)"
  stringr::str_extract_all(x, pattern)
}

homolog_accession_table <- bind_rows(lapply(uniprot_source_columns, function(colname) {
  acc_list <- extract_uniprot_accessions(annotation_table[[colname]])
  tibble(
    gene_id = rep(annotation_table$gene_id, lengths(acc_list)),
    uniprot_accession = unlist(acc_list, use.names = FALSE),
    accession_source_column = colname
  )
})) %>%
  filter(!is.na(uniprot_accession), uniprot_accession != "") %>%
  mutate(uniprot_accession = toupper(uniprot_accession)) %>%
  distinct(gene_id, uniprot_accession, .keep_all = TRUE)

if (nrow(homolog_accession_table) == 0) {
  stop("No valid UniProt-like accessions could be extracted from the Novogene homology columns.")
}

write_csv(homolog_accession_table, file.path(uniprot_dir, "S2_gene_to_uniprot_accession_mapping.csv"))

fetch_uniprot_batch <- function(accessions, batch_size = 100, sleep_seconds = 0.3) {
  accessions <- unique(na.omit(accessions))
  accessions <- accessions[accessions != ""]
  if (length(accessions) == 0) return(tibble())

  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))

  out <- lapply(seq_along(batches), function(i) {
    batch <- batches[[i]]
    query <- paste0("(", paste0("accession:", batch, collapse = " OR "), ")")
    url <- paste0(
      "https://rest.uniprot.org/uniprotkb/search?query=",
      utils::URLencode(query, reserved = TRUE),
      "&fields=accession,id,gene_names,protein_name,organism_name,go_id&format=tsv&size=500"
    )

    message("Querying UniProt batch ", i, " of ", length(batches), " (", length(batch), " accessions)")
    Sys.sleep(sleep_seconds)

    tryCatch(
      readr::read_tsv(url, show_col_types = FALSE, progress = FALSE),
      error = function(e) {
        warning("UniProt query failed for batch ", i, ": ", conditionMessage(e))
        tibble()
      }
    )
  })

  bind_rows(out)
}

uniprot_cache_file <- file.path(uniprot_dir, "S2_UniProt_downloaded_annotations.tsv")

if (file.exists(uniprot_cache_file)) {
  uniprot_raw <- readr::read_tsv(uniprot_cache_file, show_col_types = FALSE)
} else {
  uniprot_raw <- fetch_uniprot_batch(unique(homolog_accession_table$uniprot_accession))
  if (nrow(uniprot_raw) == 0) {
    stop("No UniProt annotations were downloaded. Check internet access and the accession format in the Novogene table.")
  }
  readr::write_tsv(uniprot_raw, uniprot_cache_file)
}

names(uniprot_raw) <- make.names(names(uniprot_raw), unique = TRUE)

get_first_available_col <- function(df, candidates) {
  matched <- intersect(candidates, names(df))
  if (length(matched) == 0) return(rep(NA_character_, nrow(df)))
  as.character(df[[matched[1]]])
}

uniprot_clean <- tibble(
  uniprot_accession = get_first_available_col(uniprot_raw, c("Entry", "accession")),
  uniprot_entry_name = get_first_available_col(uniprot_raw, c("Entry.Name", "id")),
  uniprot_gene_names = get_first_available_col(uniprot_raw, c("Gene.Names", "gene_names")),
  uniprot_protein_name = get_first_available_col(uniprot_raw, c("Protein.names", "protein_name")),
  uniprot_organism = get_first_available_col(uniprot_raw, c("Organism", "organism_name")),
  uniprot_go_ids = get_first_available_col(uniprot_raw, c("Gene.Ontology.IDs", "go_id"))
) %>%
  mutate(
    uniprot_accession = toupper(uniprot_accession),
    primary_gene_symbol = stringr::str_split_fixed(dplyr::coalesce(uniprot_gene_names, ""), "\\s+", 2)[, 1],
    primary_gene_symbol = ifelse(primary_gene_symbol == "", NA_character_, primary_gene_symbol),
    entry_name_symbol = stringr::str_replace(dplyr::coalesce(uniprot_entry_name, ""), "_.*$", ""),
    entry_name_symbol = ifelse(entry_name_symbol == "", NA_character_, entry_name_symbol),
    common_gene_name = dplyr::coalesce(primary_gene_symbol, entry_name_symbol, uniprot_accession)
  ) %>%
  distinct(uniprot_accession, .keep_all = TRUE)

write_csv(uniprot_clean, file.path(uniprot_dir, "S2_UniProt_cleaned_annotations.csv"))

first_non_empty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

collapse_non_empty <- function(x) {
  x <- as.character(x)
  x <- sort(unique(x[!is.na(x) & x != ""]))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ";")
}

gene_uniprot_summary <- homolog_accession_table %>%
  left_join(uniprot_clean, by = "uniprot_accession") %>%
  group_by(gene_id) %>%
  summarise(
    homolog_uniprot_accessions = collapse_non_empty(uniprot_accession),
    homolog_entry_names = collapse_non_empty(uniprot_entry_name),
    common_gene_name = first_non_empty(common_gene_name),
    homolog_gene_names = collapse_non_empty(uniprot_gene_names),
    homolog_protein_name = first_non_empty(uniprot_protein_name),
    homolog_organism = first_non_empty(uniprot_organism),
    homolog_go_ids = collapse_non_empty(unlist(stringr::str_extract_all(paste(uniprot_go_ids, collapse = ";"), "GO:\\d{7}"))),
    .groups = "drop"
  )

write_csv(gene_uniprot_summary, file.path(uniprot_dir, "S2_gene_level_UniProt_annotation_summary.csv"))

annotation_table <- annotation_table %>%
  left_join(gene_uniprot_summary, by = "gene_id") %>%
  mutate(
    final_gene_label = dplyr::coalesce(common_gene_name, homolog_protein_name, best_description, gene_id),
    best_description = dplyr::coalesce(homolog_protein_name, best_description),
    annotation_search_text = paste(
      annotation_search_text,
      homolog_uniprot_accessions,
      homolog_entry_names,
      common_gene_name,
      homolog_gene_names,
      homolog_protein_name,
      homolog_organism,
      homolog_go_ids,
      sep = " | "
    ),
    annotation_search_text_upper = toupper(annotation_search_text)
  )

write_csv(annotation_table, file.path(uniprot_dir, "S2_annotation_table_with_UniProt_common_names_and_GO.csv"))

uniprot_go_mapping_tidy <- gene_uniprot_summary %>%
  filter(!is.na(homolog_go_ids), homolog_go_ids != "") %>%
  mutate(go_id = stringr::str_extract_all(homolog_go_ids, "GO:\\d{7}")) %>%
  tidyr::unnest(go_id) %>%
  distinct(gene_id, go_id)

if (nrow(uniprot_go_mapping_tidy) == 0) {
  stop("No GO identifiers were retrieved from UniProt homolog annotations. GO enrichment cannot proceed.")
}

write_csv(uniprot_go_mapping_tidy, file.path(uniprot_dir, "S2_UniProt_derived_GO_mapping_tidy.csv"))

# ------------------------------------------------------------------------------
# 3. DESeq2 normalization, low-count filtering, and VST transformation
# ------------------------------------------------------------------------------

dds_all <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ condition
)

keep_all <- rowSums(counts(dds_all) >= min_count) >= min_samples
dds_all <- dds_all[keep_all, ]
dds_all <- DESeq(dds_all)
vsd_all <- vst(dds_all, blind = FALSE)
vsd_all_mat <- assay(vsd_all)

male_samples <- rownames(sample_metadata)[sample_metadata$tissue == "Microsporophyll"]

dds_rep <- DESeqDataSetFromMatrix(
  countData = count_matrix[, male_samples],
  colData = droplevels(sample_metadata[male_samples, ]),
  design = ~ condition
)

keep_rep <- rowSums(counts(dds_rep) >= min_count) >= min_samples
dds_rep <- dds_rep[keep_rep, ]
dds_rep <- DESeq(dds_rep)
vsd_rep <- vst(dds_rep, blind = FALSE)
vsd_rep_mat <- assay(vsd_rep)

dataset_summary <- tibble(
  raw_transcripts = nrow(count_matrix),
  filtered_transcripts_all_samples = nrow(dds_all),
  filtered_transcripts_reproductive_stages = nrow(dds_rep),
  n_libraries = ncol(count_matrix),
  n_C1 = sum(sample_metadata$condition == "C1"),
  n_C2 = sum(sample_metadata$condition == "C2"),
  n_C3 = sum(sample_metadata$condition == "C3"),
  n_Leaf = sum(sample_metadata$condition == "Leaf")
)

write_csv(dataset_summary, file.path(output_dir, "S2_dataset_filtering_summary.csv"))

# ------------------------------------------------------------------------------
# 4. Quality control and PCA
# ------------------------------------------------------------------------------

qc_dir <- file.path(output_dir, "QC_PCA")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

svg(file.path(qc_dir, "S2_boxplot_vst_all_samples.svg"), width = 10, height = 6)
boxplot(vsd_all_mat, las = 2, main = "Cycas VST distribution of counts", col = rainbow(ncol(vsd_all_mat)))
dev.off()

melted_vst_all <- reshape2::melt(vsd_all_mat)

svg(file.path(qc_dir, "S2_densityplot_vst_all_samples.svg"), width = 10, height = 6)
print(
  ggplot(melted_vst_all, aes(x = value, color = Var2)) +
    geom_density(linewidth = 0.7) +
    theme_minimal(base_size = 14) +
    labs(title = "Cycas VST density plot", x = "VST expression", y = "Density")
)
dev.off()

plot_pca_custom <- function(vsd_object, intgroup, output_svg, title_label) {
  pca_data <- plotPCA(vsd_object, intgroup = intgroup, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  svg(output_svg, width = 8, height = 6)
  print(
    ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[intgroup]], label = name)) +
      geom_point(size = 4) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
      theme_minimal(base_size = 14) +
      labs(
        title = title_label,
        x = paste0("PC1: ", percent_var[1], "%"),
        y = paste0("PC2: ", percent_var[2], "%"),
        color = intgroup
      )
  )
  dev.off()

  list(pca_data = pca_data, percent_var = percent_var)
}

pca_all <- plot_pca_custom(
  vsd_object = vsd_all,
  intgroup = "condition",
  output_svg = file.path(qc_dir, "S2_PCA_all_samples_C1_C2_C3_Leaf.svg"),
  title_label = "Cycas PCA: microsporophyll stages and leaf"
)

pca_rep <- plot_pca_custom(
  vsd_object = vsd_rep,
  intgroup = "condition",
  output_svg = file.path(qc_dir, "S2_PCA_reproductive_stages_only.svg"),
  title_label = "Cycas PCA: microsporophyll stages only"
)

sample_dists_all <- dist(t(vsd_all_mat))
sample_dist_matrix_all <- as.matrix(sample_dists_all)
rownames(sample_dist_matrix_all) <- colnames(vsd_all_mat)
colnames(sample_dist_matrix_all) <- colnames(vsd_all_mat)

ann_col_all <- data.frame(
  Condition = sample_metadata[colnames(vsd_all_mat), "condition"],
  Tissue = sample_metadata[colnames(vsd_all_mat), "tissue"]
)
rownames(ann_col_all) <- colnames(vsd_all_mat)

svg(file.path(qc_dir, "S2_sample_distance_heatmap_all_samples.svg"), width = 8, height = 7)
pheatmap(
  sample_dist_matrix_all,
  annotation_col = ann_col_all,
  annotation_row = ann_col_all,
  clustering_distance_rows = sample_dists_all,
  clustering_distance_cols = sample_dists_all,
  main = "Cycas sample-to-sample distance"
)
dev.off()

# ------------------------------------------------------------------------------
# 5. Differential expression analysis
# ------------------------------------------------------------------------------

deg_dir <- file.path(output_dir, "DESeq2_DEGs")
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)

filter_deg <- function(res, contrast_name, padj_cutoff = 0.005, lfc_cutoff = 1.5) {
  as.data.frame(res) %>%
    mutate(gene_id = rownames(.), contrast = contrast_name) %>%
    filter(!is.na(padj)) %>%
    filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
}

summarise_deg_contrast <- function(res, contrast_name, padj_cutoff = 0.005, lfc_cutoff = 1.5) {
  df <- as.data.frame(res) %>% mutate(gene_id = rownames(.))
  sig <- df %>% filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  tibble(
    contrast = contrast_name,
    tested_transcripts_with_nonNA_padj = sum(!is.na(df$padj)),
    significant_DEGs_total = nrow(sig),
    upregulated_in_numerator_stage = sum(sig$log2FoldChange > lfc_cutoff),
    downregulated_in_numerator_stage = sum(sig$log2FoldChange < -lfc_cutoff),
    padj_cutoff = padj_cutoff,
    abs_log2FC_cutoff = lfc_cutoff
  )
}

# Developmental contrasts among reproductive stages.
res_C2_vs_C1 <- results(dds_rep, contrast = c("condition", "C2", "C1"))
res_C3_vs_C2 <- results(dds_rep, contrast = c("condition", "C3", "C2"))
res_C3_vs_C1 <- results(dds_rep, contrast = c("condition", "C3", "C1"))

deg_sig_C2_vs_C1 <- filter_deg(res_C2_vs_C1, "C2_vs_C1", padj_cutoff, lfc_cutoff)
deg_sig_C3_vs_C2 <- filter_deg(res_C3_vs_C2, "C3_vs_C2", padj_cutoff, lfc_cutoff)
deg_sig_C3_vs_C1 <- filter_deg(res_C3_vs_C1, "C3_vs_C1", padj_cutoff, lfc_cutoff)

write_csv(deg_sig_C2_vs_C1, file.path(deg_dir, "S2_DEGs_C2_vs_C1.csv"))
write_csv(deg_sig_C3_vs_C2, file.path(deg_dir, "S2_DEGs_C3_vs_C2.csv"))
write_csv(deg_sig_C3_vs_C1, file.path(deg_dir, "S2_DEGs_C3_vs_C1.csv"))

developmental_deg_ids <- unique(c(
  deg_sig_C2_vs_C1$gene_id,
  deg_sig_C3_vs_C2$gene_id,
  deg_sig_C3_vs_C1$gene_id
))

all_degs_table <- bind_rows(
  deg_sig_C2_vs_C1,
  deg_sig_C3_vs_C2,
  deg_sig_C3_vs_C1
) %>%
  filter(!is.na(gene_id), !is.na(padj)) %>%
  arrange(gene_id, padj, desc(abs(log2FoldChange))) %>%
  group_by(gene_id) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  left_join(annotation_table, by = "gene_id")

write_csv(all_degs_table, file.path(deg_dir, "S2_all_unique_developmental_DEGs_best_contrast.csv"))

# Reproductive-vs-leaf contrasts.
res_C1_vs_Leaf <- results(dds_all, contrast = c("condition", "C1", "Leaf"))
res_C2_vs_Leaf <- results(dds_all, contrast = c("condition", "C2", "Leaf"))
res_C3_vs_Leaf <- results(dds_all, contrast = c("condition", "C3", "Leaf"))

deg_sig_C1_vs_Leaf <- filter_deg(res_C1_vs_Leaf, "C1_vs_Leaf", padj_cutoff, lfc_cutoff)
deg_sig_C2_vs_Leaf <- filter_deg(res_C2_vs_Leaf, "C2_vs_Leaf", padj_cutoff, lfc_cutoff)
deg_sig_C3_vs_Leaf <- filter_deg(res_C3_vs_Leaf, "C3_vs_Leaf", padj_cutoff, lfc_cutoff)

write_csv(deg_sig_C1_vs_Leaf, file.path(deg_dir, "S2_DEGs_C1_vs_Leaf.csv"))
write_csv(deg_sig_C2_vs_Leaf, file.path(deg_dir, "S2_DEGs_C2_vs_Leaf.csv"))
write_csv(deg_sig_C3_vs_Leaf, file.path(deg_dir, "S2_DEGs_C3_vs_Leaf.csv"))

reproductive_enriched_ids <- unique(c(
  deg_sig_C1_vs_Leaf %>% filter(log2FoldChange > lfc_cutoff) %>% pull(gene_id),
  deg_sig_C2_vs_Leaf %>% filter(log2FoldChange > lfc_cutoff) %>% pull(gene_id),
  deg_sig_C3_vs_Leaf %>% filter(log2FoldChange > lfc_cutoff) %>% pull(gene_id)
))

deg_summary <- bind_rows(
  summarise_deg_contrast(res_C2_vs_C1, "C2_vs_C1", padj_cutoff, lfc_cutoff),
  summarise_deg_contrast(res_C3_vs_C2, "C3_vs_C2", padj_cutoff, lfc_cutoff),
  summarise_deg_contrast(res_C3_vs_C1, "C3_vs_C1", padj_cutoff, lfc_cutoff),
  summarise_deg_contrast(res_C1_vs_Leaf, "C1_vs_Leaf", padj_cutoff, lfc_cutoff),
  summarise_deg_contrast(res_C2_vs_Leaf, "C2_vs_Leaf", padj_cutoff, lfc_cutoff),
  summarise_deg_contrast(res_C3_vs_Leaf, "C3_vs_Leaf", padj_cutoff, lfc_cutoff),
  tibble(
    contrast = "Union_developmental_contrasts",
    tested_transcripts_with_nonNA_padj = NA_integer_,
    significant_DEGs_total = length(developmental_deg_ids),
    upregulated_in_numerator_stage = NA_integer_,
    downregulated_in_numerator_stage = NA_integer_,
    padj_cutoff = padj_cutoff,
    abs_log2FC_cutoff = lfc_cutoff
  ),
  tibble(
    contrast = "Union_reproductive_enriched_vs_leaf",
    tested_transcripts_with_nonNA_padj = NA_integer_,
    significant_DEGs_total = length(reproductive_enriched_ids),
    upregulated_in_numerator_stage = NA_integer_,
    downregulated_in_numerator_stage = NA_integer_,
    padj_cutoff = padj_cutoff,
    abs_log2FC_cutoff = lfc_cutoff
  )
)

write_csv(deg_summary, file.path(deg_dir, "S2_DEG_summary_by_contrast.csv"))

# ------------------------------------------------------------------------------
# 6. UniProt-derived GO mapping and Ginkgo-compatible GO categories
# ------------------------------------------------------------------------------

go_dir <- file.path(output_dir, "GO_Enrichment")
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

# GO enrichment is based on GO identifiers downloaded from UniProt using the
# Novogene homolog accessions. GO columns directly provided by Novogene are not
# used for enrichment.
go_mapping_tidy <- uniprot_go_mapping_tidy %>%
  filter(gene_id %in% rownames(vsd_all_mat)) %>%
  distinct(gene_id, go_id)

write_csv(go_mapping_tidy, file.path(go_dir, "S2_Cycas_GO_mapping_tidy_UniProt_derived.csv"))

# Seed terms match the Ginkgo analysis.
go_categories_seed <- list(
  Auxin = c("GO:0009850", "GO:0009851", "GO:0009733", "GO:0009734", "GO:0009926", "GO:0060918", "GO:0060919"),
  Gibberellin = c("GO:0009686", "GO:0009739", "GO:0009740"),
  Cytokinin = c("GO:0009691", "GO:0009735", "GO:0009736"),
  Ethylene = c("GO:0009692", "GO:0009693", "GO:0009723", "GO:0009873"),
  ABA = c("GO:0009687", "GO:0009688", "GO:0009737", "GO:0009738"),
  Brassinosteroids = c("GO:0016131", "GO:0016132", "GO:0009741", "GO:0009742"),
  Cell_Cycle = c("GO:0007049", "GO:0051301", "GO:0000278", "GO:0007067", "GO:0051321", "GO:0007126", "GO:0007059", "GO:0006260", "GO:0000082", "GO:0000086"),
  Sugar_Transport = c("GO:0005975", "GO:0005996", "GO:0005985", "GO:0005982", "GO:0006006", "GO:0006000", "GO:0008643", "GO:1901474", "GO:0015770", "GO:0055056"),
  Cell_Wall = c("GO:0071554", "GO:0009832", "GO:0009664", "GO:0042545", "GO:0009827", "GO:0044036", "GO:0030244", "GO:0010383", "GO:0010413", "GO:0010411", "GO:0009834", "GO:0071669"),
  Reproductive_Development = c("GO:0003006", "GO:0009908", "GO:0048443", "GO:0009555", "GO:0009556", "GO:0048235", "GO:0009860", "GO:0009553", "GO:0009554", "GO:0010154", "GO:0048316", "GO:0010431"),
  Dehiscence_Desiccation = c("GO:0009269", "GO:0009414", "GO:0009901", "GO:0080166", "GO:0120194", "GO:0048767", "GO:0010150", "GO:0010432")
)

build_go_categories_from_ontology <- function(go_categories_seed) {
  go_obo_file <- "go-basic.obo"
  go_obo_url <- "https://current.geneontology.org/ontology/go-basic.obo"

  if (!file.exists(go_obo_file)) {
    tryCatch(
      download.file(go_obo_url, destfile = go_obo_file, mode = "wb"),
      error = function(e) warning("GO ontology download failed. Seed-only GO categories will be used.")
    )
  }

  if (!file.exists(go_obo_file)) return(go_categories_seed)

  go <- ontologyIndex::get_ontology(go_obo_file, extract_tags = "everything")

  get_descendants_safe <- function(go, go_id) {
    if (!go_id %in% names(go$name)) return(character(0))
    descendants <- tryCatch(
      ontologyIndex::get_descendants(go, go_id),
      error = function(e) character(0)
    )
    unique(c(go_id, descendants))
  }

  clean_go_terms <- function(go, terms, namespace = "biological_process") {
    terms <- unique(terms)
    terms <- terms[terms %in% names(go$name)]
    terms <- terms[go$namespace[terms] %in% namespace]

    if (!is.null(go$obsolete)) {
      terms <- terms[!go$obsolete[terms] %in% TRUE]
    }

    too_generic_names <- c(
      "biological_process",
      "biological process",
      "cellular process",
      "metabolic process",
      "cellular metabolic process",
      "primary metabolic process",
      "developmental process",
      "response to stimulus"
    )

    terms <- terms[!tolower(go$name[terms]) %in% too_generic_names]
    unique(terms)
  }

  lapply(go_categories_seed, function(seeds) {
    expanded_terms <- unique(unlist(lapply(seeds, function(x) get_descendants_safe(go, x))))
    clean_go_terms(go, expanded_terms)
  })
}

if (file.exists(ginkgo_go_categories_rds)) {
  go_categories <- readRDS(ginkgo_go_categories_rds)
  message("Ginkgo expanded GO-category list loaded from: ", ginkgo_go_categories_rds)
} else {
  go_categories <- build_go_categories_from_ontology(go_categories_seed)
  saveRDS(go_categories, file.path(go_dir, "S2_GO_categories_expanded_or_seed_list.rds"))
  message("GO-category list generated from seed terms and current ontology when available.")
}

go_category_table <- bind_rows(lapply(names(go_categories), function(cat) {
  tibble(Category = cat, GO_ID = go_categories[[cat]])
}))
write_csv(go_category_table, file.path(go_dir, "S2_GO_categories_used_for_Cycas.csv"))

# ------------------------------------------------------------------------------
# 7. GO enrichment analysis
# ------------------------------------------------------------------------------

go_enrichment_analysis <- function(deg_list, mapping_df, go_categories, label) {
  deg_go <- mapping_df %>% filter(gene_id %in% deg_list)
  total_genes_with_go <- length(unique(deg_go$gene_id))
  total_background <- length(unique(mapping_df$gene_id))

  results <- tibble(
    Group = character(),
    Category = character(),
    n_DEG_hits = integer(),
    n_background_hits = integer(),
    p_value = numeric()
  )

  for (cat in names(go_categories)) {
    terms <- go_categories[[cat]]
    deg_hits <- length(unique(deg_go$gene_id[deg_go$go_id %in% terms]))
    bg_hits <- length(unique(mapping_df$gene_id[mapping_df$go_id %in% terms]))

    if (bg_hits > 0 && total_genes_with_go > 0) {
      mat <- matrix(
        c(deg_hits, total_genes_with_go - deg_hits, bg_hits, total_background - bg_hits),
        nrow = 2
      )
      pval <- fisher.test(mat, alternative = "greater")$p.value
    } else {
      pval <- 1
    }

    results <- bind_rows(
      results,
      tibble(
        Group = label,
        Category = cat,
        n_DEG_hits = deg_hits,
        n_background_hits = bg_hits,
        p_value = pval
      )
    )
  }

  results
}

pairwise_res_C2_vs_C1 <- go_enrichment_analysis(deg_sig_C2_vs_C1$gene_id, go_mapping_tidy, go_categories, "C2_vs_C1")
pairwise_res_C3_vs_C2 <- go_enrichment_analysis(deg_sig_C3_vs_C2$gene_id, go_mapping_tidy, go_categories, "C3_vs_C2")
pairwise_res_C3_vs_C1 <- go_enrichment_analysis(deg_sig_C3_vs_C1$gene_id, go_mapping_tidy, go_categories, "C3_vs_C1")

global_enrichment <- bind_rows(
  pairwise_res_C2_vs_C1,
  pairwise_res_C3_vs_C2,
  pairwise_res_C3_vs_C1
) %>%
  group_by(Group) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    log10_p = -log10(p_value),
    log10_padj = -log10(p_adj_BH),
    Group = factor(Group, levels = c("C2_vs_C1", "C3_vs_C2", "C3_vs_C1"))
  )

write_csv(global_enrichment, file.path(go_dir, "S2_global_GO_enrichment_pairwise_Cycas.csv"))

svg(file.path(go_dir, "S2_global_GO_enrichment_bubble_plot.svg"), width = 10, height = 6)
print(
  ggplot(global_enrichment, aes(x = Group, y = Category)) +
    geom_point(aes(size = n_DEG_hits, color = log10_p), alpha = 0.85) +
    scale_color_gradient(low = "yellow", high = "red", name = "-log10(p-value)") +
    scale_size_continuous(range = c(3, 12), name = "DEG hits") +
    theme_bw(base_size = 12) +
    labs(title = "Cycas targeted GO-category enrichment", x = "Developmental contrast", y = "GO category") +
    theme(
      axis.text.x = element_text(face = "bold", size = 11),
      axis.text.y = element_text(face = "bold", size = 11),
      panel.grid.major = element_line(color = "gray92"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    )
)
dev.off()

# ------------------------------------------------------------------------------
# 8. Mfuzz clustering omitted by design
# ------------------------------------------------------------------------------

# Mfuzz clustering is not performed in this Cycas-focused comparative workflow.
# The Cycas sampling series is developmentally informative but not strictly equivalent
# to the Ginkgo M1-M2-M3 framework: C2 corresponds to a more advanced reproductive
# phase than the Ginkgo M2 microsporogenesis window, and C3 includes cones in which
# some sporangia were already undergoing or approaching dehiscence.
# Therefore, the Cycas analysis is intentionally restricted to:
#   1. PCA-based quality control and stage/tissue separation;
#   2. pairwise differential expression;
#   3. targeted GO-category enrichment using Ginkgo-compatible categories;
#   4. a Ginkgo-informed representative reproductive candidate heatmap.

# ------------------------------------------------------------------------------
# 9. Focused reproductive candidate-gene workflow
# ------------------------------------------------------------------------------

candidate_dir <- file.path(output_dir, "Reproductive_Candidates")
dir.create(candidate_dir, showWarnings = FALSE, recursive = TRUE)

pollen_categories_seed <- list(
  Tapetum = c("GO:0010234", "GO:0048658", "GO:0048657", "GO:0048656", "GO:0048655"),
  Microsporogenesis = c("GO:0055046", "GO:0010480", "GO:0009556", "GO:0010152", "GO:0048235"),
  PollenWall = c("GO:0062075", "GO:0010584", "GO:0010208", "GO:0160030", "GO:0080110"),
  Dehiscence = c("GO:0120194", "GO:0080166", "GO:0009901")
)

# Load the GO ontology once. Re-loading go-basic.obo inside grouped operations is
# very slow for de novo transcriptome-scale datasets, so all GO depth and descendant
# calculations are precomputed in this section.
load_go_ontology_once <- function(go_obo_file = "go-basic.obo") {
  if (!file.exists(go_obo_file)) {
    warning("go-basic.obo was not found. Reproductive GO categories will use seed terms only, and GO depth will be set to NA.")
    return(NULL)
  }

  ontologyIndex::get_ontology(go_obo_file, extract_tags = "everything")
}

go_ontology <- load_go_ontology_once("go-basic.obo")

expand_reproductive_go_terms_fast <- function(seed_list, go) {
  if (is.null(go)) {
    return(seed_list)
  }

  lapply(seed_list, function(seeds) {
    expanded_terms <- unique(unlist(lapply(seeds, function(x) {
      if (!x %in% names(go$name)) return(x)

      descendants <- tryCatch(
        ontologyIndex::get_descendants(go, x),
        error = function(e) character(0)
      )

      unique(c(x, descendants))
    })))

    expanded_terms <- expanded_terms[expanded_terms %in% names(go$name)]

    # Reproductive candidate categories are biological-process categories.
    if (!is.null(go$namespace)) {
      expanded_terms <- expanded_terms[go$namespace[expanded_terms] %in% "biological_process"]
    }

    if (!is.null(go$obsolete)) {
      expanded_terms <- expanded_terms[!go$obsolete[expanded_terms] %in% TRUE]
    }

    unique(expanded_terms)
  })
}

pollen_categories_expanded <- expand_reproductive_go_terms_fast(pollen_categories_seed, go_ontology)

write_csv(
  bind_rows(lapply(names(pollen_categories_expanded), function(cat) {
    tibble(Category = cat, GO_ID = pollen_categories_expanded[[cat]])
  })),
  file.path(candidate_dir, "S2_reproductive_GO_terms_used_by_category.csv")
)

# Precompute GO depth once for all GO terms that may be used in the reproductive
# candidate workflow. This avoids repeated ontology traversal inside dplyr::summarise().
all_reproductive_go_ids <- unique(unlist(pollen_categories_expanded, use.names = FALSE))
all_mapped_go_ids <- unique(go_mapping_tidy$go_id)
go_ids_for_depth <- unique(c(all_reproductive_go_ids, all_mapped_go_ids))

compute_go_depth_table <- function(go_ids, go) {
  if (is.null(go)) {
    return(tibble(go_id = go_ids, go_depth = NA_real_))
  }

  go_ids <- unique(go_ids)
  go_ids <- go_ids[!is.na(go_ids) & go_ids != ""]

  depths <- vapply(go_ids, function(x) {
    if (!x %in% names(go$name)) return(NA_real_)

    ancestors <- tryCatch(
      ontologyIndex::get_ancestors(go, x),
      error = function(e) character(0)
    )

    as.numeric(length(unique(ancestors)))
  }, numeric(1))

  tibble(go_id = go_ids, go_depth = depths)
}

go_depth_table <- compute_go_depth_table(go_ids_for_depth, go_ontology)

write_csv(go_depth_table, file.path(candidate_dir, "S2_GO_depth_table_precomputed.csv"))

# Build a gene-category support table with precomputed GO-depth values.
category_go_long <- bind_rows(lapply(names(pollen_categories_expanded), function(cat) {
  tibble(
    Candidate_Category = cat,
    go_id = pollen_categories_expanded[[cat]]
  )
})) %>%
  distinct(Candidate_Category, go_id)

target_gene_scores_all_categories <- go_mapping_tidy %>%
  inner_join(category_go_long, by = "go_id") %>%
  distinct(gene_id, Candidate_Category, go_id) %>%
  left_join(go_depth_table, by = "go_id") %>%
  group_by(gene_id, Candidate_Category) %>%
  summarise(
    n_matching_go_terms = n_distinct(go_id),
    matching_go_ids = paste(sort(unique(go_id)), collapse = ";"),
    mean_go_depth = mean(go_depth, na.rm = TRUE),
    max_go_depth = max(go_depth, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_go_depth = ifelse(is.nan(mean_go_depth) | is.infinite(mean_go_depth), NA_real_, mean_go_depth),
    max_go_depth = ifelse(is.nan(max_go_depth) | is.infinite(max_go_depth), NA_real_, max_go_depth)
  )

write_csv(
  target_gene_scores_all_categories,
  file.path(candidate_dir, "S2_reproductive_GO_gene_support_scores_all_categories.csv")
)

candidate_df <- all_degs_table %>%
  filter(!is.na(gene_id), !is.na(padj)) %>%
  inner_join(target_gene_scores_all_categories, by = "gene_id") %>%
  mutate(
    abs_log2FC = abs(log2FoldChange),
    reproductive_enriched_vs_leaf = gene_id %in% reproductive_enriched_ids
  ) %>%
  group_by(Candidate_Category) %>%
  arrange(
    padj,
    desc(abs_log2FC),
    desc(reproductive_enriched_vs_leaf),
    desc(n_matching_go_terms),
    desc(mean_go_depth),
    desc(max_go_depth),
    .by_group = TRUE
  ) %>%
  slice_head(n = 300) %>%
  ungroup()

category_priority <- c("Tapetum", "PollenWall", "Dehiscence", "Microsporogenesis")

if (nrow(candidate_df) == 0) {
  warning("No reproductive candidates were recovered from UniProt-derived GO annotations. Empty candidate tables will be written.")
  candidate_df_unique <- candidate_df
  candidate_summary <- tibble(
    Candidate_Category = names(pollen_categories_expanded),
    n_unique_candidates = 0L
  )
} else {
  candidate_df_unique <- candidate_df %>%
    mutate(
      category_rank = match(Candidate_Category, category_priority),
      category_rank = ifelse(is.na(category_rank), 999, category_rank),
      padj_for_sort = ifelse(is.na(padj), Inf, padj),
      abs_log2FC_for_sort = ifelse(is.na(abs_log2FC), -Inf, abs_log2FC),
      n_go_for_sort = ifelse(is.na(n_matching_go_terms), -Inf, n_matching_go_terms),
      mean_depth_for_sort = ifelse(is.na(mean_go_depth), -Inf, mean_go_depth)
    ) %>%
    group_by(gene_id) %>%
    arrange(
      category_rank,
      desc(reproductive_enriched_vs_leaf),
      desc(n_go_for_sort),
      desc(mean_depth_for_sort),
      padj_for_sort,
      desc(abs_log2FC_for_sort),
      .by_group = TRUE
    ) %>%
    slice_head(n = 1) %>%
    ungroup()

  candidate_summary <- candidate_df_unique %>%
    group_by(Candidate_Category) %>%
    summarise(n_unique_candidates = n(), .groups = "drop") %>%
    right_join(
      tibble(Candidate_Category = names(pollen_categories_expanded)),
      by = "Candidate_Category"
    ) %>%
    mutate(n_unique_candidates = ifelse(is.na(n_unique_candidates), 0L, n_unique_candidates)) %>%
    arrange(Candidate_Category)
}

write_csv(candidate_df, file.path(candidate_dir, "S2_all_reproductive_candidates_before_duplicate_resolution.csv"))
write_csv(candidate_df_unique, file.path(candidate_dir, "S2_all_reproductive_candidates_unique_priority_based.csv"))
write_csv(candidate_summary, file.path(candidate_dir, "S2_reproductive_candidate_summary_by_category.csv"))

# ------------------------------------------------------------------------------
# 10. Stage-level expression summaries for candidate filtering
# ------------------------------------------------------------------------------

compute_stage_means_for_genes <- function(gene_ids, mat, metadata) {
  gene_ids <- intersect(gene_ids, rownames(mat))

  if (length(gene_ids) == 0) {
    return(tibble(
      gene_id = character(),
      C1_mean = numeric(),
      C2_mean = numeric(),
      C3_mean = numeric(),
      Leaf_mean = numeric(),
      peak_reproductive_stage = character(),
      max_reproductive_mean = numeric(),
      leaf_minus_max_reproductive = numeric()
    ))
  }

  cy1_samples <- rownames(metadata)[metadata$condition == "C1"]
  cy2_samples <- rownames(metadata)[metadata$condition == "C2"]
  cy3_samples <- rownames(metadata)[metadata$condition == "C3"]
  leaf_samples <- rownames(metadata)[metadata$condition == "Leaf"]

  out <- tibble(
    gene_id = gene_ids,
    C1_mean = rowMeans(mat[gene_ids, cy1_samples, drop = FALSE]),
    C2_mean = rowMeans(mat[gene_ids, cy2_samples, drop = FALSE]),
    C3_mean = rowMeans(mat[gene_ids, cy3_samples, drop = FALSE]),
    Leaf_mean = rowMeans(mat[gene_ids, leaf_samples, drop = FALSE])
  )

  rep_mean_matrix <- as.matrix(out[, c("C1_mean", "C2_mean", "C3_mean")])
  max_index <- max.col(rep_mean_matrix, ties.method = "first")

  out$max_reproductive_mean <- apply(rep_mean_matrix, 1, max, na.rm = TRUE)
  out$leaf_minus_max_reproductive <- out$Leaf_mean - out$max_reproductive_mean
  out$peak_reproductive_stage <- c("C1", "C2", "C3")[max_index]

  out
}

stage_means_all <- compute_stage_means_for_genes(rownames(vsd_all_mat), vsd_all_mat, sample_metadata)

# ------------------------------------------------------------------------------
# 11. Legacy annotation-supported marker search omitted
# ------------------------------------------------------------------------------

# The previous broad keyword-based marker search is omitted to avoid generating
# parallel candidate tables based on older criteria. Candidate selection for the
# main heatmap and for diagnostic outputs is now handled by the compact
# family-based workflow below.

# ------------------------------------------------------------------------------
# 11. Ginkgo-informed compact family-based Cycas heatmap
# ------------------------------------------------------------------------------

heatmap_dir <- file.path(output_dir, "Representative_Heatmap")
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

agp_fla_dir <- file.path(output_dir, "AGP_FLA_PFAM_Screen")
dir.create(agp_fla_dir, showWarnings = FALSE, recursive = TRUE)

row_zscore <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z
}

# Main logic:
#   - Use broad but biologically targeted annotation patterns.
#   - Apply the same reproductive filters to all families:
#       max_reproductive_mean > Leaf_mean
#       reproductive_enriched_vs_leaf == TRUE
#   - Do not filter by peak stage.
#   - Keep the heatmap compact by selecting a limited number of candidates per family.
#   - Preserve subfamily diversity so that key targets such as TPD1, AMS, and
#     sporopollenin-related genes are not lost simply because another subfamily
#     has more or stronger paralogous candidates.

deg_stats_best <- all_degs_table %>%
  dplyr::select(any_of(c(
    "gene_id", "contrast", "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj"
  ))) %>%
  distinct(gene_id, .keep_all = TRUE)

candidate_source_pool <- annotation_table %>%
  filter(gene_id %in% rownames(vsd_all_mat)) %>%
  left_join(deg_stats_best, by = "gene_id") %>%
  dplyr::select(-any_of(c(
    "C1_mean", "C2_mean", "C3_mean", "Leaf_mean",
    "peak_reproductive_stage", "max_reproductive_mean",
    "leaf_minus_max_reproductive"
  ))) %>%
  left_join(stage_means_all, by = "gene_id") %>%
  mutate(
    abs_log2FC = abs(log2FoldChange),
    reproductive_enriched_vs_leaf = gene_id %in% reproductive_enriched_ids,
    microsporophyll_enriched_relative_to_leaf = !is.na(max_reproductive_mean) & !is.na(Leaf_mean) & max_reproductive_mean > Leaf_mean,
    max_rep_minus_leaf = max_reproductive_mean - Leaf_mean,
    has_uniprot_name = !is.na(common_gene_name) & common_gene_name != "",
    has_swissprot_hit = !is.na(Swissprot.ID) & Swissprot.ID != "" & Swissprot.ID != "--",
    has_pfam = !is.na(PFAM.ID) & PFAM.ID != "" & PFAM.ID != "--",
    annotation_confidence_score = as.integer(has_uniprot_name) + as.integer(has_swissprot_hit) + as.integer(has_pfam),
    common_gene_name_upper = toupper(dplyr::coalesce(common_gene_name, "")),
    search_text_upper = toupper(dplyr::coalesce(annotation_search_text, "")),
    final_gene_label_clean = dplyr::case_when(
      !is.na(common_gene_name) & common_gene_name != "" ~ common_gene_name,
      !is.na(final_gene_label) & final_gene_label != "" ~ final_gene_label,
      !is.na(homolog_protein_name) & homolog_protein_name != "" ~ homolog_protein_name,
      !is.na(best_description) & best_description != "" ~ best_description,
      TRUE ~ gene_id
    )
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

write_csv(candidate_source_pool, file.path(heatmap_dir, "S2_candidate_source_pool_all_expressed_annotated_transcripts.csv"))

family_heatmap_rules <- tibble::tribble(
  ~group_order, ~Heatmap_Category, ~Family_Target, ~search_pattern, ~max_n, ~max_per_subfamily, ~curated_note,
  1, "SPL candidates", "SPL", "\\bSPL\\d*\\b|SQUAMOSA PROMOTER-BINDING-LIKE|SQUAMOSA", 5, 1,
     "Top reproductive-enriched SPL-like candidates.",
  2, "MADS/AGL candidates", "MADS_AGL", "MADS|AGAMOUS|\\bAGL\\d*\\b|SRF-TYPE TRANSCRIPTION FACTOR|GGM", 5, 1,
     "Top reproductive-enriched MADS-box/AGL/GGM-like candidates.",
  3, "EMS/TPD signaling candidates", "EMS_TPD", "\\bEMS1\\b|EXCESS MICROSPOROCYTES|\\bTPD1\\b|TAPETUM DETERMINANT", 5, 3,
     "Top reproductive-enriched EMS1- and TPD1-like candidates. The broad EXS acronym alone is intentionally not used.",
  4, "Tapetum regulation / PCD candidates", "AMS_MYB80_EAT_PTC", "\\bAMS\\b|ABORTED MICROSPORES|\\bMYB80\\b|\\bEAT1\\b|ETERNAL TAPETUM|\\bPTC1\\b|PERSISTENT TAPETAL", 5, 1,
     "Top reproductive-enriched AMS, MYB80, EAT1, and PTC1-like candidates.",
  5, "Advanced maturation candidates", "MYB101", "\\bMYB101\\b", 5, 2,
     "Top reproductive-enriched MYB101-like candidates.",
  6, "CALS5 candidates", "CALS5", "\\bCALS5\\b|CALLOSE SYNTHASE 5", 3, 2,
     "Top reproductive-enriched CALS5-like candidates.",
  7, "Sporopollenin candidates", "SPOROPOLLENIN", "CYP703|CYP704|\\bPKSB\\b|POLYKETIDE SYNTHASE|\\bTKPR1\\b|\\bTKPR2\\b|TETRAKETIDE|\\bABCG26\\b", 5, 1,
     "Top reproductive-enriched CYP703/CYP704/PKSB/TKPR/ABCG26-like candidates.",
  8, "LAC/EXPA wall remodeling candidates", "LAC_EXPA", "\\bLAC4\\b|\\bLAC17\\b|\\bLAC\\d*\\b|LACCASE|\\bEXPA1\\b|\\bEXPA\\d*\\b|EXPANSIN|\\bEXLA\\d*\\b", 5, 1,
     "Top reproductive-enriched LAC/EXPA/EXLA-like candidates."
)

write_csv(family_heatmap_rules, file.path(heatmap_dir, "S2_family_based_heatmap_rules.csv"))

# Manual exclusions after biological/annotation review.
# These transcripts remain in the diagnostic report but are not eligible for the
# main compact heatmap. This allows automatic replacement by the next best
# candidate passing the same filters.
manual_exclusion_after_review <- tibble::tribble(
  ~gene_id, ~manual_exclusion_reason,
  "Cluster-29530.4399", "Excluded from the main sporopollenin module after manual annotation review: CXP;2-3 label is ambiguous for canonical sporopollenin interpretation."
)

write_csv(
  manual_exclusion_after_review,
  file.path(heatmap_dir, "S2_manually_excluded_from_main_heatmap_after_review.csv")
)

assign_target_subfamily <- function(df, family_target) {
  df %>%
    mutate(
      target_subfamily = dplyr::case_when(
        family_target == "SPL" & str_detect(search_text_upper, "\\bSPL2\\b") ~ "SPL2",
        family_target == "SPL" & str_detect(search_text_upper, "\\bSPL8\\b") ~ "SPL8",
        family_target == "SPL" ~ dplyr::coalesce(common_gene_name, "SPL_other"),

        family_target == "MADS_AGL" & str_detect(search_text_upper, "\\bAGL104\\b") ~ "AGL104",
        family_target == "MADS_AGL" & str_detect(search_text_upper, "\\bAGL65\\b") ~ "AGL65",
        family_target == "MADS_AGL" & str_detect(search_text_upper, "\\bAGL66\\b") ~ "AGL66",
        family_target == "MADS_AGL" & str_detect(search_text_upper, "\\bAGL") ~ dplyr::coalesce(common_gene_name, "AGL_other"),
        family_target == "MADS_AGL" ~ dplyr::coalesce(common_gene_name, "MADS_other"),

        family_target == "EMS_TPD" & str_detect(search_text_upper, "\\bEMS1\\b|EXCESS MICROSPOROCYTES") ~ "EMS1",
        family_target == "EMS_TPD" & str_detect(search_text_upper, "\\bTPD1\\b|TAPETUM DETERMINANT") ~ "TPD1",

        family_target == "AMS_MYB80_EAT_PTC" & str_detect(search_text_upper, "\\bAMS\\b|ABORTED MICROSPORES") ~ "AMS",
        family_target == "AMS_MYB80_EAT_PTC" & str_detect(search_text_upper, "\\bMYB80\\b") ~ "MYB80",
        family_target == "AMS_MYB80_EAT_PTC" & str_detect(search_text_upper, "\\bEAT1\\b|ETERNAL TAPETUM") ~ "EAT1",
        family_target == "AMS_MYB80_EAT_PTC" & str_detect(search_text_upper, "\\bPTC1\\b|PERSISTENT TAPETAL") ~ "PTC1",

        family_target == "MYB101" ~ "MYB101",

        family_target == "CALS5" ~ "CALS5",

        family_target == "SPOROPOLLENIN" & str_detect(search_text_upper, "CYP703") ~ "CYP703",
        family_target == "SPOROPOLLENIN" & str_detect(search_text_upper, "CYP704") ~ "CYP704",
        family_target == "SPOROPOLLENIN" & str_detect(search_text_upper, "\\bPKSB\\b|POLYKETIDE SYNTHASE") ~ "PKSB",
        family_target == "SPOROPOLLENIN" & str_detect(search_text_upper, "\\bTKPR1\\b|\\bTKPR2\\b|TETRAKETIDE") ~ "TKPR",
        family_target == "SPOROPOLLENIN" & str_detect(search_text_upper, "\\bABCG26\\b") ~ "ABCG26",

        family_target == "LAC_EXPA" & str_detect(search_text_upper, "\\bLAC17\\b") ~ "LAC17",
        family_target == "LAC_EXPA" & str_detect(search_text_upper, "\\bLAC4\\b") ~ "LAC4",
        family_target == "LAC_EXPA" & str_detect(search_text_upper, "\\bLAC\\d*\\b|LACCASE") ~ "LAC_other",
        family_target == "LAC_EXPA" & str_detect(search_text_upper, "\\bEXPA1\\b") ~ "EXPA1",
        family_target == "LAC_EXPA" & str_detect(search_text_upper, "\\bEXPA\\d*\\b|EXPANSIN") ~ "EXPA_other",
        family_target == "LAC_EXPA" & str_detect(search_text_upper, "\\bEXLA\\d*\\b") ~ "EXLA_other",

        TRUE ~ dplyr::coalesce(common_gene_name, family_target, "unclassified")
      )
    )
}

select_family_candidates <- function(rule, source_pool) {
  candidates <- source_pool %>%
    filter(str_detect(search_text_upper, regex(rule$search_pattern, ignore_case = TRUE))) %>%
    mutate(
      group_order = rule$group_order,
      Heatmap_Category = rule$Heatmap_Category,
      Family_Target = rule$Family_Target,
      search_pattern = rule$search_pattern,
      curated_note = rule$curated_note,
      passes_microsporophyll_expression_filter = microsporophyll_enriched_relative_to_leaf,
      passes_reproductive_enrichment_filter = reproductive_enriched_vs_leaf,
      manually_excluded_after_review = gene_id %in% manual_exclusion_after_review$gene_id,
      manual_exclusion_reason = dplyr::case_when(
        manually_excluded_after_review ~ manual_exclusion_after_review$manual_exclusion_reason[
          match(gene_id, manual_exclusion_after_review$gene_id)
        ],
        TRUE ~ NA_character_
      ),
      passes_main_filters = passes_microsporophyll_expression_filter &
        passes_reproductive_enrichment_filter &
        !manually_excluded_after_review
    ) %>%
    assign_target_subfamily(rule$Family_Target)

  candidates_ranked <- candidates %>%
    arrange(
      desc(passes_main_filters),
      target_subfamily,
      padj,
      desc(abs_log2FC),
      desc(max_rep_minus_leaf),
      desc(annotation_confidence_score),
      common_gene_name,
      gene_id
    )

  # Keep all matched candidates in the diagnostic report, but select a compact,
  # diverse subset for the main heatmap. The number of retained candidates per
  # target_subfamily is rule-specific. This avoids over-compression of key
  # two-component groups such as EMS1/TPD1 while keeping broad families readable.
  selected <- candidates_ranked %>%
    filter(passes_main_filters) %>%
    group_by(target_subfamily) %>%
    arrange(
      padj,
      desc(abs_log2FC),
      desc(max_rep_minus_leaf),
      desc(annotation_confidence_score),
      .by_group = TRUE
    ) %>%
    slice_head(n = as.integer(rule$max_per_subfamily)) %>%
    ungroup() %>%
    arrange(
      padj,
      desc(abs_log2FC),
      desc(max_rep_minus_leaf),
      target_subfamily,
      gene_id
    )

  if (is.finite(rule$max_n)) {
    selected <- selected %>% slice_head(n = as.integer(rule$max_n))
  }

  list(candidates = candidates_ranked, selected = selected)
}

family_selection_results <- lapply(seq_len(nrow(family_heatmap_rules)), function(i) {
  select_family_candidates(family_heatmap_rules[i, ], candidate_source_pool)
})

all_family_candidate_report <- bind_rows(lapply(family_selection_results, function(x) x$candidates))
selected_heatmap_genes <- bind_rows(lapply(family_selection_results, function(x) x$selected)) %>%
  distinct(gene_id, Heatmap_Category, .keep_all = TRUE)

duplicated_family_assignments <- selected_heatmap_genes %>%
  group_by(gene_id) %>%
  filter(n() > 1) %>%
  ungroup() %>%
  arrange(gene_id, group_order)

selected_heatmap_genes <- selected_heatmap_genes %>%
  group_by(gene_id) %>%
  arrange(group_order, padj, desc(abs_log2FC), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

targets_not_recovered <- family_heatmap_rules %>%
  left_join(
    selected_heatmap_genes %>%
      group_by(Family_Target) %>%
      summarise(n_selected = n_distinct(gene_id), .groups = "drop"),
    by = "Family_Target"
  ) %>%
  mutate(n_selected = ifelse(is.na(n_selected), 0L, n_selected)) %>%
  filter(n_selected == 0L) %>%
  dplyr::select(group_order, Heatmap_Category, Family_Target, curated_note, n_selected)

heatmap_category_levels <- family_heatmap_rules$Heatmap_Category

selected_heatmap_genes <- selected_heatmap_genes %>%
  mutate(
    Heatmap_Category = factor(Heatmap_Category, levels = heatmap_category_levels),
    display_label = dplyr::case_when(
      !is.na(common_gene_name) & common_gene_name != "" ~ paste0(gene_id, " (", common_gene_name, ")"),
      TRUE ~ paste0(gene_id, " (", str_trunc(final_gene_label_clean, 30), ")")
    )
  ) %>%
  arrange(group_order, target_subfamily, padj, desc(abs_log2FC), gene_id)

write_csv(all_family_candidate_report, file.path(heatmap_dir, "S2_all_family_candidate_report_all_matches.csv"))
write_csv(selected_heatmap_genes, file.path(heatmap_dir, "S2_selected_family_based_heatmap_genes.csv"))
write_csv(targets_not_recovered, file.path(heatmap_dir, "S2_family_targets_not_recovered_after_filters.csv"))
write_csv(duplicated_family_assignments, file.path(heatmap_dir, "S2_selected_candidates_with_multiple_family_assignments.csv"))

ems_tpd_diagnostic_report <- all_family_candidate_report %>%
  filter(Family_Target == "EMS_TPD") %>%
  mutate(
    selected_in_main_heatmap = gene_id %in% selected_heatmap_genes$gene_id,
    exclusion_reason = dplyr::case_when(
      selected_in_main_heatmap ~ "selected_in_main_heatmap",
      manually_excluded_after_review ~ paste0("excluded after manual review: ", manual_exclusion_reason),
      !passes_microsporophyll_expression_filter & !passes_reproductive_enrichment_filter ~ "excluded: not microsporophyll-enriched and not reproductive-enriched vs leaf",
      !passes_microsporophyll_expression_filter ~ "excluded: max_reproductive_mean <= Leaf_mean",
      !passes_reproductive_enrichment_filter ~ "excluded: not reproductive_enriched_vs_leaf",
      TRUE ~ "excluded by compact top-N/subfamily ranking"
    )
  ) %>%
  arrange(
    target_subfamily,
    desc(selected_in_main_heatmap),
    padj,
    desc(abs_log2FC),
    desc(max_rep_minus_leaf),
    gene_id
  )

write_csv(ems_tpd_diagnostic_report, file.path(heatmap_dir, "S2_EMS_TPD_candidate_diagnostic_report.csv"))

candidate_inspection_dir <- file.path(heatmap_dir, "Candidate_Inspection_Heatmaps")
dir.create(candidate_inspection_dir, showWarnings = FALSE, recursive = TRUE)

plot_family_heatmap <- function(group_name, candidate_df, mat, metadata, out_file, title_label) {
  group_df <- candidate_df %>%
    filter(Heatmap_Category == group_name) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(
      desc(passes_main_filters),
      target_subfamily,
      padj,
      desc(abs_log2FC),
      desc(max_rep_minus_leaf),
      common_gene_name,
      gene_id
    )

  if (nrow(group_df) == 0) return(invisible(NULL))

  group_df <- group_df %>%
    filter(gene_id %in% rownames(mat)) %>%
    mutate(
      row_label = dplyr::case_when(
        !is.na(common_gene_name) & common_gene_name != "" ~ paste0(gene_id, " (", common_gene_name, ")"),
        TRUE ~ paste0(gene_id, " (", str_trunc(final_gene_label_clean, 30), ")")
      ),
      row_label = make.unique(row_label)
    )

  if (nrow(group_df) == 0) return(invisible(NULL))

  diagnostic_mat <- mat[group_df$gene_id, , drop = FALSE]
  diagnostic_mat_z <- row_zscore(diagnostic_mat)
  rownames(diagnostic_mat_z) <- group_df$row_label

  sample_order_diag <- rownames(metadata)[order(metadata$condition, rownames(metadata))]
  sample_order_diag <- intersect(sample_order_diag, colnames(diagnostic_mat_z))
  diagnostic_mat_z <- diagnostic_mat_z[, sample_order_diag, drop = FALSE]

  ann_col_diag <- data.frame(
    Condition = metadata[colnames(diagnostic_mat_z), "condition"],
    Tissue = metadata[colnames(diagnostic_mat_z), "tissue"]
  )
  rownames(ann_col_diag) <- colnames(diagnostic_mat_z)

  ann_row_diag <- data.frame(
    Passes_main_filters = ifelse(group_df$passes_main_filters, "yes", "no"),
    Manual_exclusion = ifelse(group_df$manually_excluded_after_review, "yes", "no"),
    Target_subfamily = group_df$target_subfamily,
    Peak_stage = group_df$peak_reproductive_stage,
    Reproductive_enriched = ifelse(group_df$reproductive_enriched_vs_leaf, "yes", "no"),
    Microsporophyll_enriched = ifelse(group_df$microsporophyll_enriched_relative_to_leaf, "yes", "no")
  )
  rownames(ann_row_diag) <- rownames(diagnostic_mat_z)

  svg(out_file, width = 14, height = max(6, 0.25 * nrow(diagnostic_mat_z) + 3))
  pheatmap(
    diagnostic_mat_z,
    annotation_col = ann_col_diag,
    annotation_row = ann_row_diag,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 9,
    angle_col = 45,
    main = title_label,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
    border_color = "white"
  )
  dev.off()

  invisible(NULL)
}

invisible(lapply(heatmap_category_levels, function(group_name) {
  safe_group_name <- gsub("[^A-Za-z0-9]+", "_", group_name)
  safe_group_name <- gsub("_+$", "", safe_group_name)

  plot_family_heatmap(
    group_name = group_name,
    candidate_df = all_family_candidate_report,
    mat = vsd_all_mat,
    metadata = sample_metadata,
    out_file = file.path(candidate_inspection_dir, paste0("S2_diagnostic_all_matches_", safe_group_name, ".svg")),
    title_label = paste0("Diagnostic all matches: ", group_name)
  )
}))

combined_diagnostic_candidates <- all_family_candidate_report %>%
  distinct(Heatmap_Category, gene_id, .keep_all = TRUE) %>%
  arrange(
    group_order,
    desc(passes_main_filters),
    target_subfamily,
    padj,
    desc(abs_log2FC),
    desc(max_rep_minus_leaf),
    common_gene_name,
    gene_id
  )

if (nrow(combined_diagnostic_candidates) > 0) {
  combined_diagnostic_candidates <- combined_diagnostic_candidates %>%
    filter(gene_id %in% rownames(vsd_all_mat)) %>%
    mutate(
      row_label = dplyr::case_when(
        !is.na(common_gene_name) & common_gene_name != "" ~ paste0(gene_id, " (", common_gene_name, ")"),
        TRUE ~ paste0(gene_id, " (", str_trunc(final_gene_label_clean, 30), ")")
      ),
      row_label = make.unique(paste0(Heatmap_Category, " | ", row_label))
    )

  combined_mat <- vsd_all_mat[combined_diagnostic_candidates$gene_id, , drop = FALSE]
  combined_mat_z <- row_zscore(combined_mat)
  rownames(combined_mat_z) <- combined_diagnostic_candidates$row_label

  sample_order_combined <- rownames(sample_metadata)[order(sample_metadata$condition, rownames(sample_metadata))]
  sample_order_combined <- intersect(sample_order_combined, colnames(combined_mat_z))
  combined_mat_z <- combined_mat_z[, sample_order_combined, drop = FALSE]

  ann_col_combined <- data.frame(
    Condition = sample_metadata[colnames(combined_mat_z), "condition"],
    Tissue = sample_metadata[colnames(combined_mat_z), "tissue"]
  )
  rownames(ann_col_combined) <- colnames(combined_mat_z)

  ann_row_combined <- data.frame(
    Group = combined_diagnostic_candidates$Heatmap_Category,
    Target_subfamily = combined_diagnostic_candidates$target_subfamily,
    Passes_main_filters = ifelse(combined_diagnostic_candidates$passes_main_filters, "yes", "no"),
    Manual_exclusion = ifelse(combined_diagnostic_candidates$manually_excluded_after_review, "yes", "no"),
    Peak_stage = combined_diagnostic_candidates$peak_reproductive_stage,
    Reproductive_enriched = ifelse(combined_diagnostic_candidates$reproductive_enriched_vs_leaf, "yes", "no"),
    Microsporophyll_enriched = ifelse(combined_diagnostic_candidates$microsporophyll_enriched_relative_to_leaf, "yes", "no")
  )
  rownames(ann_row_combined) <- rownames(combined_mat_z)

  group_runs_diag <- rle(as.character(combined_diagnostic_candidates$Heatmap_Category))
  gaps_row_diag <- cumsum(group_runs_diag$lengths)
  gaps_row_diag <- gaps_row_diag[-length(gaps_row_diag)]

  svg(
    file.path(candidate_inspection_dir, "S2_diagnostic_heatmap_all_family_matches_combined.svg"),
    width = 18,
    height = max(10, 0.22 * nrow(combined_mat_z) + 4)
  )

  pheatmap(
    combined_mat_z,
    annotation_col = ann_col_combined,
    annotation_row = ann_row_combined,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 6,
    fontsize_col = 8,
    angle_col = 45,
    main = "Diagnostic inspection: all matched candidates before/after filters",
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
    border_color = "white",
    gaps_row = gaps_row_diag
  )

  dev.off()
}

# Separate AGP/FLA PFAM/annotation-supported screen.
agp_fla_screen <- candidate_source_pool %>%
  mutate(
    pfam_agp_fla_support = str_detect(
      toupper(paste(dplyr::coalesce(PFAM.ID, ""), dplyr::coalesce(PFAM.description, ""), sep = " | ")),
      "PF02469|PF06376|FASCICLIN|ARABINOGALACTAN"
    ),
    annotation_agp_fla_support = str_detect(
      search_text_upper,
      "FASCICLIN-LIKE ARABINOGALACTAN|ARABINOGALACTAN|\\bFLA\\b|\\bAGP\\b"
    ),
    agp_fla_supported = pfam_agp_fla_support | annotation_agp_fla_support,
    keep_for_possible_supplementary_AGP_heatmap =
      agp_fla_supported & microsporophyll_enriched_relative_to_leaf & reproductive_enriched_vs_leaf
  ) %>%
  filter(agp_fla_supported) %>%
  arrange(
    desc(keep_for_possible_supplementary_AGP_heatmap),
    padj,
    desc(abs_log2FC),
    desc(max_rep_minus_leaf)
  )

agp_fla_screen_filtered <- agp_fla_screen %>%
  filter(keep_for_possible_supplementary_AGP_heatmap)

write_csv(agp_fla_screen, file.path(agp_fla_dir, "S2_AGP_FLA_PFAM_annotation_supported_screen_all_expressed.csv"))
write_csv(agp_fla_screen_filtered, file.path(agp_fla_dir, "S2_AGP_FLA_PFAM_annotation_supported_microsporophyll_enriched.csv"))

make_agp_fla_supplementary_heatmap <- TRUE

if (isTRUE(make_agp_fla_supplementary_heatmap) && nrow(agp_fla_screen_filtered) > 0) {
  agp_fla_heatmap_genes <- agp_fla_screen_filtered %>%
    mutate(agp_fla_label = dplyr::case_when(
      !is.na(common_gene_name) & common_gene_name != "" ~ common_gene_name,
      !is.na(final_gene_label) & final_gene_label != "" ~ str_trunc(final_gene_label, 30),
      TRUE ~ gene_id
    )) %>%
    group_by(agp_fla_label) %>%
    arrange(padj, desc(abs_log2FC), desc(max_rep_minus_leaf), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    slice_head(n = 20) %>%
    mutate(display_label = paste0(gene_id, " (", agp_fla_label, ")"))

  agp_mat <- vsd_all_mat[agp_fla_heatmap_genes$gene_id, , drop = FALSE]
  agp_mat_z <- row_zscore(agp_mat)
  rownames(agp_mat_z) <- agp_fla_heatmap_genes$display_label

  sample_order_agp <- rownames(sample_metadata)[order(sample_metadata$condition, rownames(sample_metadata))]
  sample_order_agp <- intersect(sample_order_agp, colnames(agp_mat_z))
  agp_mat_z <- agp_mat_z[, sample_order_agp, drop = FALSE]

  agp_ann_col <- data.frame(
    Condition = sample_metadata[colnames(agp_mat_z), "condition"],
    Tissue = sample_metadata[colnames(agp_mat_z), "tissue"]
  )
  rownames(agp_ann_col) <- colnames(agp_mat_z)

  svg(file.path(agp_fla_dir, "S2_optional_supplementary_AGP_FLA_microsporophyll_enriched_heatmap.svg"),
      width = 12, height = max(6, 0.25 * nrow(agp_mat_z) + 3))
  pheatmap(
    agp_mat_z,
    annotation_col = agp_ann_col,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 9,
    fontsize_col = 9,
    angle_col = 45,
    main = "Microsporophyll-enriched AGP/FLA-like candidates",
    color = colorRampPalette(RColorBrewer::brewer.pal(7, "RdPu"))(100),
    border_color = "white"
  )
  dev.off()
}

if (nrow(selected_heatmap_genes) == 0) {
  stop("No family-based candidates passed the main heatmap filters.")
}

if (isTRUE(include_leaf_in_candidate_heatmap)) {
  heatmap_mat_source <- vsd_all_mat
  heatmap_metadata <- sample_metadata
  output_heatmap_svg <- file.path(heatmap_dir, "S2_representative_Cycas_compact_family_based_reproductive_candidate_heatmap_with_leaf.svg")
} else {
  heatmap_mat_source <- vsd_rep_mat
  heatmap_metadata <- droplevels(sample_metadata[colnames(vsd_rep_mat), , drop = FALSE])
  output_heatmap_svg <- file.path(heatmap_dir, "S2_representative_Cycas_compact_family_based_reproductive_candidate_heatmap_reproductive_only.svg")
}

selected_heatmap_genes <- selected_heatmap_genes %>%
  filter(gene_id %in% rownames(heatmap_mat_source)) %>%
  mutate(display_label = make.unique(display_label)) %>%
  arrange(group_order, target_subfamily, padj, desc(abs_log2FC), gene_id)

heatmap_mat <- heatmap_mat_source[selected_heatmap_genes$gene_id, , drop = FALSE]
heatmap_mat_z <- row_zscore(heatmap_mat)
rownames(heatmap_mat_z) <- selected_heatmap_genes$display_label

sample_order <- rownames(heatmap_metadata)[order(heatmap_metadata$condition, rownames(heatmap_metadata))]
sample_order <- intersect(sample_order, colnames(heatmap_mat_z))
heatmap_mat_z <- heatmap_mat_z[, sample_order, drop = FALSE]

ann_col <- data.frame(
  Condition = heatmap_metadata[colnames(heatmap_mat_z), "condition"],
  Tissue = heatmap_metadata[colnames(heatmap_mat_z), "tissue"]
)
rownames(ann_col) <- colnames(heatmap_mat_z)

ann_row <- data.frame(
  Category = selected_heatmap_genes$Heatmap_Category,
  Target_subfamily = selected_heatmap_genes$target_subfamily,
  Peak_stage = selected_heatmap_genes$peak_reproductive_stage
)
rownames(ann_row) <- rownames(heatmap_mat_z)

category_runs <- rle(as.character(selected_heatmap_genes$Heatmap_Category))
gaps_row <- cumsum(category_runs$lengths)
gaps_row <- gaps_row[-length(gaps_row)]

condition_runs <- rle(as.character(ann_col$Condition))
gaps_col <- cumsum(condition_runs$lengths)
gaps_col <- gaps_col[-length(gaps_col)]

category_colors <- c(
  "SPL candidates" = "#6A3D9A",
  "MADS/AGL candidates" = "#984EA3",
  "EMS/TPD signaling candidates" = "#1F78B4",
  "Tapetum regulation / PCD candidates" = "#33A02C",
  "Advanced maturation candidates" = "#FF7F00",
  "CALS5 candidates" = "#66C2A5",
  "Sporopollenin candidates" = "#01665E",
  "LAC/EXPA wall remodeling candidates" = "#D95F02"
)

condition_colors <- c(
  "C1" = "#87CEEB",
  "C2" = "#90EE90",
  "C3" = "#FFD700",
  "Leaf" = "#BDBDBD"
)

tissue_colors <- c(
  "Microsporophyll" = "#756BB1",
  "Leaf" = "#74C476"
)

svg(output_heatmap_svg, width = 14, height = max(8, 0.30 * nrow(heatmap_mat_z) + 3))
pheatmap(
  heatmap_mat_z,
  annotation_col = ann_col,
  annotation_row = ann_row,
  annotation_colors = list(
    Condition = condition_colors,
    Tissue = tissue_colors,
    Category = category_colors
  ),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 9,
  angle_col = 45,
  main = "Compact family-based Cycas reproductive candidate module",
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
  border_color = "white",
  gaps_row = gaps_row,
  gaps_col = gaps_col
)
dev.off()

# ------------------------------------------------------------------------------
# 12. Manuscript-ready summary statistics
# ------------------------------------------------------------------------------

summary_dir <- file.path(output_dir, "Manuscript_Ready_Summary")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

pca_summary <- tibble(
  PCA = c("all_samples", "reproductive_stages_only"),
  PC1_percent_variance = c(pca_all$percent_var[1], pca_rep$percent_var[1]),
  PC2_percent_variance = c(pca_all$percent_var[2], pca_rep$percent_var[2])
)

write_csv(pca_summary, file.path(summary_dir, "S2_PCA_summary.csv"))

candidate_and_heatmap_summary <- bind_rows(
  candidate_summary %>%
    transmute(Summary_type = "GO_reproductive_candidate_catalogue", Group = Candidate_Category, n_genes = n_unique_candidates),
  selected_heatmap_genes %>%
    group_by(Heatmap_Category) %>%
    summarise(n_genes = n(), .groups = "drop") %>%
    transmute(Summary_type = "selected_representative_heatmap", Group = as.character(Heatmap_Category), n_genes = n_genes)
)

write_csv(candidate_and_heatmap_summary, file.path(summary_dir, "S2_candidate_and_heatmap_gene_summary.csv"))

if (exists("agp_fla_screen")) {
  agp_summary <- tibble(
    metric = c(
      "PFAM_or_annotation_supported_AGP_FLA_like_candidates_expressed",
      "PFAM_or_annotation_supported_AGP_FLA_like_candidates_microsporophyll_enriched",
      "PFAM_or_annotation_supported_AGP_FLA_like_candidates_reproductive_enriched_vs_leaf",
      "PFAM_or_annotation_supported_AGP_FLA_like_candidates_kept_for_supplementary_heatmap"
    ),
    value = c(
      nrow(agp_fla_screen),
      sum(agp_fla_screen$microsporophyll_enriched_relative_to_leaf, na.rm = TRUE),
      sum(agp_fla_screen$reproductive_enriched_vs_leaf, na.rm = TRUE),
      sum(agp_fla_screen$keep_for_possible_supplementary_AGP_heatmap, na.rm = TRUE)
    )
  )

  write_csv(agp_summary, file.path(summary_dir, "S2_AGP_FLA_like_summary.csv"))
} else {
  agp_summary <- tibble(
    metric = "PFAM_or_annotation_supported_AGP_FLA_like_candidates_screen_not_run",
    value = NA_real_
  )
  write_csv(agp_summary, file.path(summary_dir, "S2_AGP_FLA_like_summary.csv"))
}

short_report_file <- file.path(summary_dir, "S2_manuscript_numbers_short_report.txt")

# Use a small report-printing helper instead of tibble-specific print arguments.
# This avoids errors when an object is a data.frame rather than a tibble inside sink().
print_table_for_report <- function(x) {
  if (is.null(x)) {
    cat("[NULL]\n")
    return(invisible(NULL))
  }

  if (nrow(as.data.frame(x)) == 0) {
    cat("[empty table]\n")
    return(invisible(NULL))
  }

  print(as.data.frame(x), row.names = FALSE)
  invisible(NULL)
}

sink(short_report_file)
tryCatch({
  cat("MANUSCRIPT-READY CYCAS TRANSCRIPTOME SUMMARY STATISTICS\n")
  cat("=======================================================\n\n")
  cat("Libraries: ", ncol(count_matrix), " total; C1 = ", sum(sample_metadata$condition == "C1"),
      ", C2 = ", sum(sample_metadata$condition == "C2"),
      ", C3 = ", sum(sample_metadata$condition == "C3"),
      ", Leaf = ", sum(sample_metadata$condition == "Leaf"), ".\n", sep = "")
  cat("Transcripts retained after low-count filtering, all samples: ", nrow(dds_all), " out of ", nrow(count_matrix), ".\n", sep = "")
  cat("Transcripts retained after low-count filtering, reproductive stages only: ", nrow(dds_rep), " out of ", nrow(count_matrix), ".\n", sep = "")
  cat("PCA all samples: PC1 = ", pca_all$percent_var[1], "%; PC2 = ", pca_all$percent_var[2], "%.\n", sep = "")
  cat("PCA reproductive stages only: PC1 = ", pca_rep$percent_var[1], "%; PC2 = ", pca_rep$percent_var[2], "%.\n\n", sep = "")

  cat("DEG summary:\n")
  print_table_for_report(deg_summary)

  cat("\nMfuzz clustering: omitted by design because Cycas stages are not strict one-to-one equivalents of Ginkgo M1-M2-M3.\n")

  cat("\nCandidate and heatmap summary:\n")
  print_table_for_report(candidate_and_heatmap_summary)

  if (exists("targets_not_recovered") && nrow(targets_not_recovered) > 0) {
    cat("\nFamily targets with no candidates after reproductive/manual-review filters:\n")
    print_table_for_report(targets_not_recovered)
  }

  if (exists("manual_exclusion_after_review") && nrow(manual_exclusion_after_review) > 0) {
    cat("\nManual exclusions from main heatmap after annotation review:\n")
    print_table_for_report(manual_exclusion_after_review)
  }

  if (exists("all_family_candidate_report")) {
    cat("\nFamily-based candidate inspection summary:\n")
    diagnostic_summary_for_report <- all_family_candidate_report %>%
      group_by(Heatmap_Category) %>%
      summarise(
        n_all_matches = n_distinct(gene_id),
        n_passing_main_filters = n_distinct(gene_id[passes_main_filters]),
        .groups = "drop"
      )
    print_table_for_report(diagnostic_summary_for_report)
  }

  if (exists("agp_summary")) {
    cat("\nAGP/FLA-like annotation-supported summary:\n")
    print_table_for_report(agp_summary)
  }

  cat("\nNote: GO enrichment results, candidate tables, and heatmap gene selections are saved in the output directory.\n")
}, finally = {
  sink()
})

writeLines(capture.output(sessionInfo()), file.path(summary_dir, "S2_R_sessionInfo.txt"))

cat("\nSupplementary Script S2 completed successfully.\n")
cat("Output directory: ", output_dir, "\n", sep = "")
