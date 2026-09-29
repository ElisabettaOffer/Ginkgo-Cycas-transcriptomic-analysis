# ==============================================================================
# SUPPLEMENTARY SCRIPT S1
# Transcriptomic Data Analysis of Ginkgo biloba Male Cone Development
# Publication-ready minimal-difference version
#
# Relative to the validated original workflow, this version only removes
# unused/exploratory code and applies Benjamini-Hochberg correction to the
# existing GO-enrichment results. DESeq2, DEG ordering, Mfuzz, GO-category
# construction, Fisher tests, candidate selection and manuscript heatmap gene
# sets remain analytically unchanged.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup and Package Loading
# ------------------------------------------------------------------------------
# This supplementary script intentionally does not install or update packages.
# Run it in an R environment where the packages listed below are already
# available. Package versions used for the analysis are recorded at the end of
# the workflow through sessionInfo().

required_packages <- c(
  "DESeq2", "Mfuzz", "Biobase",
  "ggplot2", "dplyr", "readxl", "pheatmap", "reshape2",
  "readr", "stringr", "tidyr", "patchwork", "tibble",
  "RColorBrewer", "e1071", "ontologyIndex", "openxlsx", "curl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running Supplementary Script S1."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readxl)
  library(pheatmap)
  library(reshape2)
  library(readr)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(tibble)
  library(RColorBrewer)
  library(e1071)
  library(ontologyIndex)
})

required_input_files <- c(
  "count_RNAseq_MaleCones.xlsx",
  "annotazione_geni.xlsx"
)
missing_input_files <- required_input_files[!file.exists(required_input_files)]
if (length(missing_input_files) > 0L) {
  stop(
    "Missing input files: ",
    paste(missing_input_files, collapse = ", "),
    ". Place them in the working directory before running the script."
  )
}

# ------------------------------------------------------------------------------
# 2. Data Import & DESeq2 Initialization
# ------------------------------------------------------------------------------
# Load raw count matrix
data <- read_excel("count_RNAseq_MaleCones.xlsx")
count_matrix <- as.matrix(data[ , -1])
rownames(count_matrix) <- data[[1]]

# Create metadata table (coldata)
sample_names <- colnames(count_matrix)
condition <- gsub("_\\d+", "", sample_names)
coldata <- data.frame(
  row.names = sample_names,
  condition = factor(condition, levels = c("M1", "M2", "M3"))
)

# Initialize DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

# Filtering: keep genes with at least 10 counts in at least 3 samples
filter_low_counts <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[filter_low_counts, ]

# Run differential expression pipeline
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# ------------------------------------------------------------------------------
# 3. Quality Control (QC) Visualizations
# ------------------------------------------------------------------------------
# VST Boxplot
svg("boxplot_vst.svg", width = 8, height = 6)
boxplot(assay(vsd), las = 2, main = "VST Distribution of Counts", col = rainbow(15))
dev.off()

# VST Density Plot
svg("densityplot_vst.svg", width = 8, height = 6)
melted <- melt(assay(vsd))
ggplot(melted, aes(x = value, color = Var2)) +
  geom_density() +
  ggtitle("Density plot (VST)") +
  theme_minimal(base_size = 15)
dev.off()

# Principal Component Analysis (PCA)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

svg("PCAlot_vst.svg", width = 8, height = 6)
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 5, position = position_jitter(width = 1, height = 1)) +
  theme_minimal(base_size = 16) +
  labs(
    title = "PCA of 15 replicates",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )
dev.off()

# Sample-to-Sample Distance Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

svg("heatmap_replicates_vst.svg", width = 8, height = 6)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample-to-Sample Distance (VST)")
dev.off()

# Cook's Distance Outlier Detection
cook <- assays(dds)[["cooks"]]
cook_plots <- vector("list", ncol(cook))

for (i in seq_len(ncol(cook))) {
  df <- tibble(gene = seq_len(nrow(cook)), cooks = cook[, i]) %>%
    filter(cooks > 1 | gene <= 1000)
  
  cook_plots[[i]] <- ggplot(df, aes(x = gene, y = cooks)) +
    geom_point(size = 0.3, alpha = 0.6) +
    geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
    theme_minimal(base_size = 10) +
    ggtitle(colnames(cook)[i]) +
    theme(axis.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5))
}

cook_panel <- wrap_plots(cook_plots, ncol = 5, byrow = TRUE) +
  plot_annotation(title = "Cook’s distance per sample (threshold = 10)")

svg("CookDistance_panel.svg", width = 12, height = 7)
print(cook_panel)
dev.off()

# ------------------------------------------------------------------------------
# 4. Differential Expression Analysis (Pairwise Contrasts)
# ------------------------------------------------------------------------------
res_M2_vs_M1 <- results(dds, contrast = c("condition", "M2", "M1"))
res_M3_vs_M2 <- results(dds, contrast = c("condition", "M3", "M2"))
res_M3_vs_M1 <- results(dds, contrast = c("condition", "M3", "M1"))

# Filter significant DEGs (padj < 0.005 & |log2FC| > 1.5)
filter_deg <- function(res) {
  as.data.frame(res) %>%
    mutate(gene_id = rownames(.)) %>%
    filter(padj < 0.005, abs(log2FoldChange) > 1.5)
}

deg_sig_M2_vs_M1 <- filter_deg(res_M2_vs_M1)
deg_sig_M3_vs_M2 <- filter_deg(res_M3_vs_M2)
deg_sig_M3_vs_M1 <- filter_deg(res_M3_vs_M1)
all_degs_ids <- unique(c(deg_sig_M2_vs_M1$gene_id, deg_sig_M3_vs_M2$gene_id, deg_sig_M3_vs_M1$gene_id))

# Combined significant DEG table across all pairwise contrasts
all_degs_table <- bind_rows(
  deg_sig_M2_vs_M1 %>% mutate(contrast = "M2_vs_M1"),
  deg_sig_M3_vs_M2 %>% mutate(contrast = "M3_vs_M2"),
  deg_sig_M3_vs_M1 %>% mutate(contrast = "M3_vs_M1")
) %>%
  filter(!is.na(gene_id), !is.na(padj)) %>%
  arrange(
    gene_id,
    padj,
    desc(abs(log2FoldChange))
  ) %>%
  group_by(gene_id) %>%
  slice_head(n = 1) %>%
  ungroup()


# ------------------------------------------------------------------------------
# 5. Annotation & GO Mapping
# ------------------------------------------------------------------------------
gene_annotations <- read_excel("annotazione_geni.xlsx") %>%
  mutate(
    uniprot_id = str_extract(gene_description, "(?<=\\|)[A-Z0-9]{6}(?=\\|)"),
    gene_symbol = str_extract(gene_description, "(?<=\\|[^|]{0,10}\\|)[^_[:space:]]+")
  ) %>%
  dplyr::select(gene_id, uniprot_id, gene_symbol)

get_go_from_uniprot <- function(id_list, chunk_size = 200) {
  results <- list()
  chunks <- split(id_list, ceiling(seq_along(id_list) / chunk_size))
  
  for (i in seq_along(chunks)) {
    ids_string <- paste(chunks[[i]], collapse = ",")
    url <- paste0("https://rest.uniprot.org/uniprotkb/accessions?accessions=", ids_string, "&format=tsv&fields=accession,go_id")
    
   temp_file <- tempfile()
    curl::curl_download(url, temp_file)
    results[[i]] <- read_tsv(temp_file, show_col_types = FALSE)
    unlink(temp_file)
  }
  
  final_df <- bind_rows(results)
  colnames(final_df) <- c("uniprot_id", "go_id")
  return(final_df)
}

valid_uniprot_ids <- unique(gene_annotations$uniprot_id[!is.na(gene_annotations$uniprot_id)])
uniprot_go_map <- get_go_from_uniprot(valid_uniprot_ids)

go_mapping_tidy <- gene_annotations %>%
  left_join(uniprot_go_map, by = "uniprot_id") %>%
  filter(!is.na(go_id)) %>%
  separate_rows(go_id, sep = "; ") %>%
  mutate(go_id = trimws(go_id))

# ------------------------------------------------------------------------------
# 6. Mfuzz Soft Clustering & ggplot2 Object Creation
# ------------------------------------------------------------------------------
mat_degs <- assay(vsd)[rownames(assay(vsd)) %in% all_degs_ids, ]
mat_mean <- data.frame(
  M1 = rowMeans(mat_degs[, 1:5]),
  M2 = rowMeans(mat_degs[, 6:10]),
  M3 = rowMeans(mat_degs[, 11:15])
)

eset <- Biobase::ExpressionSet(assayData = as.matrix(mat_mean))
eset_s <- Mfuzz::standardise(eset)
m_est <- Mfuzz::mestimate(eset_s)
set.seed(123)
cl <- Mfuzz::mfuzz(eset_s, c = 6, m = m_est)

# Extract gene lists per cluster for downstream enrichment analysis
cluster_genes <- list()
for(i in 1:6){
  cluster_name <- paste0("Cluster_", i)
  cluster_genes[[cluster_name]] <- names(cl$cluster[cl$cluster == i & cl$membership[,i] > 0.5])
}

# Reproducibility guard: the validated analysis run used for the manuscript
# produced the following core-gene counts at membership > 0.5. This check does
# not alter the clustering; it prevents downstream outputs from being generated
# if the Mfuzz input order or execution path changes unexpectedly.
expected_core_counts <- c(
  Cluster_2 = 1387L,
  Cluster_3 = 901L,
  Cluster_4 = 1906L,
  Cluster_6 = 1327L
)
observed_core_counts <- lengths(cluster_genes)[names(expected_core_counts)]

if (!identical(unname(observed_core_counts), unname(expected_core_counts))) {
  stop(
    "Mfuzz reproducibility check failed. Expected core-gene counts: ",
    paste(names(expected_core_counts), expected_core_counts, sep = "=", collapse = ", "),
    "; observed: ",
    paste(names(observed_core_counts), observed_core_counts, sep = "=", collapse = ", "),
    ". Check that the DEG matrix retains the original VST row order."
  )
}
message("Mfuzz reproducibility check passed.")

# Tidy cluster data structure for ggplot2 compatibility
expression_matrix <- as.data.frame(Biobase::exprs(eset_s))
expression_matrix$gene_id <- rownames(expression_matrix)
expression_matrix$Cluster <- paste0("Cluster_", cl$cluster)
expression_matrix$Membership <- apply(cl$membership, 1, max)

mfuzz_long <- expression_matrix %>%
  tidyr::pivot_longer(cols = c("M1", "M2", "M3"), names_to = "Stage", values_to = "Expression") %>%
  mutate(Stage = factor(Stage, levels = c("M1", "M2", "M3")))

centers_df <- as.data.frame(cl$centers) %>%
  mutate(Cluster = paste0("Cluster_", row_number())) %>%
  tidyr::pivot_longer(cols = c("M1", "M2", "M3"), names_to = "Stage", values_to = "Center_Expression") %>%
  mutate(Stage = factor(Stage, levels = c("M1", "M2", "M3")))

# Generate Row 1: Mfuzz temporal clusters alignment
p_mfuzz <- ggplot() +
  geom_line(data = mfuzz_long, aes(x = Stage, y = Expression, group = gene_id, color = Membership), alpha = 0.25) +
  geom_line(data = centers_df, aes(x = Stage, y = Center_Expression, group = Cluster), color = "red", linewidth = 1.2) +
  scale_color_gradient(low = "lightyellow", high = "coral1", name = "Membership") +
  facet_wrap(~ Cluster, nrow = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "Mfuzz Temporal Clustering Profiles", x = "Developmental Stage", y = "Standardized Expression") +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

# ------------------------------------------------------------------------------
# 7. GO Enrichment Analysis & Function Definition
# ------------------------------------------------------------------------------
# Expanded custom GO categories generated by GO semantic expansion
# ------------------------------------------------------------------------------
# This block replaces the old static go_categories <- list(...).
# It builds go_categories dynamically from the current Gene Ontology DAG.
# Strategy:
#   1) curated GO-Slim-like seed terms for each biological category;
#   2) optional narrow keyword rescue only to find additional seed GO terms;
#   3) expansion to all descendants in the GO graph;
#   4) removal of obsolete and overly generic terms;
#   5) export of the final category table for inspection.
#
# The same go_categories object is used by the global enrichment,
# cluster-specific enrichment and reproductive candidate-gene workflows.

# ---- 7.0 GO ontology download/loading ----
go_obo_file <- "go-basic.obo"
go_obo_url  <- "https://current.geneontology.org/ontology/go-basic.obo"

if (!file.exists(go_obo_file)) {
  message("Downloading latest Gene Ontology: ", go_obo_url)
  download.file(go_obo_url, destfile = go_obo_file, mode = "wb")
} else {
  message("Using existing local Gene Ontology file: ", go_obo_file)
}

go <- ontologyIndex::get_ontology(go_obo_file, extract_tags = "everything")

# ---- 7.1 Curated GO seed terms ----
# These are intentionally NOT exhaustive: they are biologically meaningful anchors.
# All their descendants are added automatically below.
go_categories_seed <- list(
  Auxin = c(
    "GO:0009850", # auxin metabolic process
    "GO:0009851", # auxin biosynthetic process
    "GO:0009733", # response to auxin
    "GO:0009734", # auxin-activated signaling pathway
    "GO:0009926", # auxin polar transport
    "GO:0060918", # auxin efflux
    "GO:0060919"  # auxin influx
  ),

  Gibberellin = c(
    "GO:0009686", # gibberellin biosynthetic process
    "GO:0009739", # response to gibberellin
    "GO:0009740"  # gibberellin mediated signaling pathway
  ),

  Cytokinin = c(
    "GO:0009691", # cytokinin biosynthetic process
    "GO:0009735", # response to cytokinin
    "GO:0009736"  # cytokinin-activated signaling pathway
  ),

  Ethylene = c(
    "GO:0009692", # ethylene metabolic process
    "GO:0009693", # ethylene biosynthetic process
    "GO:0009723", # response to ethylene
    "GO:0009873"  # ethylene-activated signaling pathway
  ),

  ABA = c(
    "GO:0009687", # abscisic acid metabolic process
    "GO:0009688", # abscisic acid biosynthetic process
    "GO:0009737", # response to abscisic acid
    "GO:0009738"  # abscisic acid-activated signaling pathway
  ),

  Brassinosteroids = c(
    "GO:0016131", # brassinosteroid metabolic process
    "GO:0016132", # brassinosteroid biosynthetic process
    "GO:0009741", # response to brassinosteroid
    "GO:0009742"  # brassinosteroid mediated signaling pathway
  ),

  Cell_Cycle = c(
    "GO:0007049", # cell cycle
    "GO:0051301", # cell division
    "GO:0000278", # mitotic cell cycle
    "GO:0007067", # mitotic nuclear division
    "GO:0051321", # meiotic cell cycle
    "GO:0007126", # meiotic nuclear division
    "GO:0007059", # chromosome segregation
    "GO:0006260", # DNA replication
    "GO:0000082", # G1/S transition of mitotic cell cycle
    "GO:0000086"  # G2/M transition of mitotic cell cycle
  ),

  Sugar_Transport = c(
    "GO:0005975", # carbohydrate metabolic process
    "GO:0005996", # monosaccharide metabolic process
    "GO:0005985", # sucrose metabolic process
    "GO:0005982", # starch metabolic process
    "GO:0006006", # glucose metabolic process
    "GO:0006000", # fructose metabolic process
    "GO:0008643", # carbohydrate transport
    "GO:1901474", # carbohydrate transmembrane transport
    "GO:0015770", # sucrose transport
    "GO:0055056"  # D-glucose transmembrane transport
  ),

  Cell_Wall = c(
    "GO:0071554", # cell wall organization or biogenesis
    "GO:0009832", # plant-type cell wall biogenesis
    "GO:0009664", # plant-type cell wall organization
    "GO:0042545", # cell wall modification
    "GO:0009827", # plant-type cell wall modification
    "GO:0044036", # cell wall macromolecule metabolic process
    "GO:0030244", # cellulose biosynthetic process
    "GO:0010383", # lignin biosynthetic process
    "GO:0010413", # pectin metabolic process
    "GO:0010411", # xyloglucan metabolic process
    "GO:0009834", # plant-type secondary cell wall biogenesis
    "GO:0071669"  # plant-type cell wall loosening
  ),

  Reproductive_Development = c(
    "GO:0003006", # developmental process involved in reproduction
    "GO:0009908", # flower development
    "GO:0048443", # stamen development
    "GO:0009555", # pollen development
    "GO:0009556", # microsporogenesis
    "GO:0048235", # pollen sperm cell differentiation
    "GO:0009860", # pollen tube growth
    "GO:0009553", # embryo sac development
    "GO:0009554", # megasporogenesis
    "GO:0010154", # fruit development
    "GO:0048316", # seed development
    "GO:0010431"  # seed maturation
  ),

  Dehiscence_Dessication = c(
    "GO:0009269", # response to desiccation
    "GO:0009414", # response to water deprivation
    "GO:0009901", # anther dehiscence
    "GO:0080166", # regulation of anther dehiscence
    "GO:0120194", # positive regulation of anther dehiscence
    "GO:0048767", # abscission
    "GO:0010150", # leaf senescence
    "GO:0010432"  # seed maturation involved in desiccation tolerance
  )
)

# ---- 7.2 Optional narrow keyword rescue ----
# This is not used as a broad filter. It only finds extra seed terms whose GO names
# explicitly contain the biological concept. Descendants are then expanded via DAG.
use_keyword_rescue <- TRUE

keyword_seed_patterns <- list(
  Auxin = "auxin|indole-3-acetic acid|IAA",
  Gibberellin = "gibberellin",
  Cytokinin = "cytokinin|zeatin",
  Ethylene = "ethylene",
  ABA = "abscisic acid",
  Brassinosteroids = "brassinosteroid",
  Cell_Cycle = "cell cycle|mitotic|mitosis|meiotic|meiosis|cytokinesis|chromosome segregation",
  Sugar_Transport = "sucrose|glucose|fructose|starch|sugar|carbohydrate transport|carbohydrate transmembrane transport",
  Cell_Wall = "cell wall|cellulose|hemicellulose|pectin|lignin|xyloglucan",
  Reproductive_Development = "pollen|anther|tapetum|microsporogenesis|megasporogenesis|embryo sac|flower development|stamen development|seed development|pollen tube",
  Dehiscence_Dessication = "dehiscence|desiccation|water deprivation|dehydration|abscission|senescence"
)

find_keyword_seeds <- function(go, pattern, namespace = "biological_process") {
  ids <- names(go$name)
  hit <- ids[grepl(pattern, go$name[ids], ignore.case = TRUE)]
  hit <- hit[go$namespace[hit] %in% namespace]
  unique(hit)
}

# ---- 7.3 DAG expansion and cleaning functions ----
get_descendants_safe <- function(go, go_id) {
  if (!go_id %in% names(go$name)) return(character(0))
  descendants <- tryCatch(
    ontologyIndex::get_descendants(go, go_id),
    error = function(e) character(0)
  )
  unique(c(go_id, descendants))
}

expand_go_terms <- function(go, seeds) {
  seeds <- unique(seeds)
  seeds <- seeds[seeds %in% names(go$name)]
  unique(unlist(lapply(seeds, function(x) get_descendants_safe(go, x))))
}

clean_go_terms <- function(go, terms, namespace = "biological_process") {
  terms <- unique(terms)
  terms <- terms[terms %in% names(go$name)]

  # Keep only Biological Process terms for this biological-category enrichment.
  terms <- terms[go$namespace[terms] %in% namespace]

  # Remove obsolete terms when ontologyIndex provides an obsolete flag.
  if (!is.null(go$obsolete)) {
    terms <- terms[!go$obsolete[terms] %in% TRUE]
  }

  # Remove only very high-level generic terms.
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

as_character_go_field <- function(field, ids) {
  if (is.null(field)) return(rep(NA_character_, length(ids)))
  vals <- field[ids]
  vapply(vals, function(x) {
    if (length(x) == 0 || all(is.na(x))) {
      NA_character_
    } else {
      paste(as.character(x), collapse = "; ")
    }
  }, character(1))
}

make_go_term_table <- function(go, go_list) {
  ids <- unique(unlist(go_list, use.names = FALSE))
  ids <- ids[ids %in% names(go$name)]

  category_df <- bind_rows(lapply(names(go_list), function(cat) {
    tibble(
      Category = cat,
      GO_ID = go_list[[cat]]
    )
  }))

  definition_field <- if (!is.null(go$def)) go$def else go$definition

  term_df <- tibble(
    GO_ID = ids,
    GO_Name = as_character_go_field(go$name, ids),
    Namespace = as_character_go_field(go$namespace, ids),
    Definition = as_character_go_field(definition_field, ids)
  )

  category_df %>%
    left_join(term_df, by = "GO_ID") %>%
    arrange(Category, GO_ID)
}

# ---- 7.4 Build final go_categories object used by the rest of the script ----
go_categories <- list()

go_category_seed_report <- tibble(
  Category = character(),
  Seed_GO_ID = character(),
  Seed_GO_Name = character(),
  Source = character(),
  Status = character()
)

for (cat_name in names(go_categories_seed)) {
  manual_seeds <- unique(go_categories_seed[[cat_name]])
  valid_manual <- manual_seeds[manual_seeds %in% names(go$name)]
  missing_manual <- setdiff(manual_seeds, valid_manual)

  keyword_seeds <- character(0)
  if (isTRUE(use_keyword_rescue) && cat_name %in% names(keyword_seed_patterns)) {
    keyword_seeds <- find_keyword_seeds(go, keyword_seed_patterns[[cat_name]])
  }

  all_seeds <- unique(c(valid_manual, keyword_seeds))
  expanded_terms <- expand_go_terms(go, all_seeds)
  cleaned_terms <- clean_go_terms(go, expanded_terms)

  go_categories[[cat_name]] <- cleaned_terms

  go_category_seed_report <- bind_rows(
    go_category_seed_report,
    tibble(
      Category = cat_name,
      Seed_GO_ID = valid_manual,
      Seed_GO_Name = unname(go$name[valid_manual]),
      Source = "manual_seed",
      Status = "valid"
    ),
    tibble(
      Category = cat_name,
      Seed_GO_ID = keyword_seeds,
      Seed_GO_Name = unname(go$name[keyword_seeds]),
      Source = "keyword_rescue_seed",
      Status = "valid"
    ),
    tibble(
      Category = cat_name,
      Seed_GO_ID = missing_manual,
      Seed_GO_Name = NA_character_,
      Source = "manual_seed",
      Status = "not_found_in_current_GO"
    )
  )

  message("GO category built: ", cat_name, " -> ", length(cleaned_terms), " GO terms")
}

# Remove duplicate seed-report rows created when manual and keyword seeds overlap.
go_category_seed_report <- go_category_seed_report %>%
  distinct(Category, Seed_GO_ID, Source, .keep_all = TRUE)

# Export for manual inspection and reproducibility.
go_categories_table <- make_go_term_table(go, go_categories)
write_csv(go_categories_table, "GO_categories_expanded_table.csv")
write_csv(go_category_seed_report, "GO_categories_seed_report.csv")
saveRDS(go_categories, "GO_categories_expanded_list.rds")

# Optional Excel export. If openxlsx is available, this creates one workbook with
# one sheet per category plus summary sheets.
if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb_go <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb_go, "ALL_GO_CATEGORIES")
  openxlsx::writeData(wb_go, "ALL_GO_CATEGORIES", go_categories_table)

  openxlsx::addWorksheet(wb_go, "SEED_REPORT")
  openxlsx::writeData(wb_go, "SEED_REPORT", go_category_seed_report)

  for (cat_name in names(go_categories)) {
    sheet_name <- substr(cat_name, 1, 31)
    openxlsx::addWorksheet(wb_go, sheet_name)
    openxlsx::writeData(
      wb_go,
      sheet_name,
      go_categories_table %>% filter(Category == cat_name)
    )
  }

  openxlsx::saveWorkbook(wb_go, "GO_categories_expanded.xlsx", overwrite = TRUE)
}

# Sanity check: print the compact R list structure requested by the user.
cat("\nCompact go_categories object created. Example:\n")
print(lapply(go_categories, head, n = 10))
cat("\nFull expanded list saved as: GO_categories_expanded_list.rds\n")
cat("Full expanded table saved as: GO_categories_expanded_table.csv / GO_categories_expanded.xlsx\n\n")


# Core function: standardized to output 'Group' column for universal use
go_enrichment_analysis <- function(deg_list, mapping_df, go_categories, label) {
  deg_go <- mapping_df %>% filter(gene_id %in% deg_list)
  total_genes_with_go <- length(unique(deg_go$gene_id))
  total_background <- length(unique(mapping_df$gene_id))
  
  # Initialize tibble explicitly to prevent column missing errors
  results <- tibble(Group = character(), Category = character(), n_DEG_hits = integer(), p_value = numeric())
  
  for (cat in names(go_categories)) {
    terms <- go_categories[[cat]]
    deg_hits <- length(unique(deg_go$gene_id[deg_go$go_id %in% terms]))
    bg_hits  <- length(unique(mapping_df$gene_id[mapping_df$go_id %in% terms]))
    
    if (bg_hits > 0 && total_genes_with_go > 0) {
      mat <- matrix(c(deg_hits, total_genes_with_go - deg_hits, bg_hits, total_background - bg_hits), nrow = 2)
      pval <- fisher.test(mat, alternative = "greater")$p.value
    } else {
      pval <- 1
    }
    results <- bind_rows(results, tibble(Group = label, Category = cat, n_DEG_hits = deg_hits, p_value = pval))
  }
  return(results)
}

# ------------------------------------------------------------------------------
# 7A. Global Pairwise GO Enrichment Analysis (Bubble Plot)
# ------------------------------------------------------------------------------
pairwise_res_M2_vs_M1 <- go_enrichment_analysis(deg_sig_M2_vs_M1$gene_id, go_mapping_tidy, go_categories, "M2_vs_M1")
pairwise_res_M3_vs_M2 <- go_enrichment_analysis(deg_sig_M3_vs_M2$gene_id, go_mapping_tidy, go_categories, "M3_vs_M2")
pairwise_res_M3_vs_M1 <- go_enrichment_analysis(deg_sig_M3_vs_M1$gene_id, go_mapping_tidy, go_categories, "M3_vs_M1")

# Apply Benjamini-Hochberg correction independently within each pairwise
# contrast. Nominal p-values are retained in the exported table for traceability.
global_enrichment <- bind_rows(
  pairwise_res_M2_vs_M1,
  pairwise_res_M3_vs_M2,
  pairwise_res_M3_vs_M1
) %>%
  group_by(Group) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    log10_p = -log10(p_value),
    log10_p_adj_BH = -log10(p_adj_BH),
    significant_BH_0_05 = p_adj_BH < 0.05,
    Group = factor(Group, levels = c("M2_vs_M1", "M3_vs_M2", "M3_vs_M1"))
  )

readr::write_csv(
  global_enrichment,
  "GO_Enrichment_Global_results_with_BH.csv"
)

svg("GO_Enrichment_Global_Bubble_Plot_BH.svg", width = 10, height = 6)
ggplot(global_enrichment, aes(x = Group, y = Category)) +
  geom_point(aes(size = n_DEG_hits, color = log10_p_adj_BH), alpha = 0.85) +
  scale_color_gradient(
    low = "yellow",
    high = "red",
    name = "-log10(BH-adjusted p-value)"
  ) +
  scale_size_continuous(range = c(3, 12), name = "Number of Genes\n(DEG hits)") +
  theme_bw(base_size = 12) +
  labs(
    title = "Global Pairwise GO Enrichment Analysis",
    x = "Pairwise Comparison",
    y = "GO Category"
  ) +
  theme(
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 11),
    panel.grid.major = element_line(color = "gray92"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
dev.off()

# ------------------------------------------------------------------------------
# 7B. Cluster-specific GO Enrichment Analysis
# ------------------------------------------------------------------------------
all_cluster_enrichment <- tibble(Group = character(), Category = character(), n_DEG_hits = integer(), p_value = numeric())

for (i in 1:6) {
  cluster_name = paste0("Cluster_", i)
  cluster_deg_list <- cluster_genes[[cluster_name]]
  if (length(cluster_deg_list) > 0) {
    cluster_res <- go_enrichment_analysis(cluster_deg_list, go_mapping_tidy, go_categories, cluster_name)
    all_cluster_enrichment <- bind_rows(all_cluster_enrichment, cluster_res)
  }
}

# Apply Benjamini-Hochberg correction independently within each Mfuzz cluster.
# The underlying gene sets, GO categories, background, and Fisher tests remain
# unchanged relative to the validated original workflow.
all_cluster_enrichment <- all_cluster_enrichment %>%
  group_by(Group) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    log10_p = -log10(p_value),
    log10_p_adj_BH = -log10(p_adj_BH),
    significant_BH_0_05 = p_adj_BH < 0.05
  )

readr::write_csv(
  all_cluster_enrichment,
  "GO_Enrichment_Mfuzz_clusters_results_with_BH.csv"
)

mfuzz_BH_significant_category_summary <- all_cluster_enrichment %>%
  group_by(Group) %>%
  summarise(
    n_BH_significant_categories = sum(significant_BH_0_05),
    significant_categories = paste(
      Category[significant_BH_0_05],
      collapse = ";"
    ),
    .groups = "drop"
  )

readr::write_tsv(
  mfuzz_BH_significant_category_summary,
  "Mfuzz_BH_significant_category_summary.tsv"
)

plot_data <- all_cluster_enrichment %>%
  filter(p_adj_BH < 0.05)

# Safe conditional visualization layer
if (nrow(plot_data) > 0) {
  p_enrich <- ggplot(
    plot_data,
    aes(x = Category, y = n_DEG_hits, fill = log10_p_adj_BH)
  ) +
    geom_col(color = "black", width = 0.7, linewidth = 0.3) +
    scale_fill_gradient(
      low = "peachpuff",
      high = "coral1",
      name = "-log10(BH-adjusted p-value)"
    ) +
    facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 12) +
    labs(
      title = "BH-significant GO Category Enrichment Dynamics",
      x = "GO Category",
      y = "Number of Genes (DEG Hits)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.text.y = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      legend.position = "right"
    )
} else {
  p_enrich <- ggplot() +
    annotate(
      "text",
      x = 1,
      y = 1,
      label = "No GO categories remain significant after BH correction",
      size = 5,
      fontface = "bold"
    ) +
    theme_void() +
    labs(title = "BH-significant GO Category Enrichment Dynamics") +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
}

# ------------------------------------------------------------------------------
# 8. Unified Panel Layout Generation (Patchwork Assembly)
# ------------------------------------------------------------------------------
combined_manuscript_plot <- p_mfuzz / p_enrich + 
  plot_layout(heights = c(1, 1.3))

svg("Combined_Mfuzz_and_GO_Enrichment_Panel_BH.svg", width = 20, height = 12)
print(combined_manuscript_plot)
dev.off()

# ------------------------------------------------------------------------------
# 9. Targeted candidate-gene and heatmap analyses
# ------------------------------------------------------------------------------
vsd_mat <- assay(vsd)
ann_col <- data.frame(Stage = factor(condition, levels = c("M1", "M2", "M3")))
rownames(ann_col) <- colnames(vsd_mat)
ann_colors <- list(Stage = c(M1 = "#87ceeb", M2 = "#90ee90", M3 = "#ffff00"))

# Preserve VST row order while attaching gene annotations for the AGP/FLA panel.
expr_df_full <- tibble(gene_id = rownames(vsd_mat)) %>%
  left_join(gene_annotations, by = "gene_id")

# ------------------------------------------------------------------------------
# 9A. Arabinogalactans (AGP/FLA) Expression Dynamics
# ------------------------------------------------------------------------------
agp_fla_genes <- c(
  "Gb_00507",
  "Gb_01998",
  "Gb_02823",
  "Gb_04299",
  "Gb_05461",
  "Gb_06021",
  "Gb_06646",
  "Gb_06647",
  "Gb_09983",
  "Gb_12627",
  "Gb_12628",
  "Gb_23221",
  "Gb_23222",
  "Gb_27249",
  "Gb_27317",
  "Gb_31981",
  "Gb_34537",
  "Gb_36195",
  "Gb_37149",
  "Gb_40445",
  "Gb_40985",
  "novel.1842",
  "novel.1207"
)

agp_df <- expr_df_full %>% filter(gene_id %in% agp_fla_genes) %>%
  mutate(label = paste(gene_id, ifelse(grepl("FLA", gene_symbol), gene_symbol, 
                                       ifelse(grepl("AGP", gene_symbol), gene_symbol, "AGP/FLA")), sep = " - "))

mat_agp <- vsd_mat[agp_df$gene_id, ]
rownames(mat_agp) <- agp_df$label
mat_agp_z <- t(scale(t(mat_agp)))

svg("Heatmap_AGP_FLA_JIM.svg", width = 10, height = 6)
pheatmap(mat_agp_z, annotation_col = ann_col, annotation_colors = ann_colors,
         cluster_cols = FALSE, cluster_rows = TRUE, show_rownames = TRUE, fontsize_row = 12,
         main = "Expression Dynamics of AGPs and FLAs",
         color = colorRampPalette(brewer.pal(7, "RdPu"))(100), border_color = "white")
dev.off()

# ------------------------------------------------------------------------------
# 9B. Reproductive candidate-gene catalogue and category assignment
#     Complete catalogue with automatic duplicate resolution
# ------------------------------------------------------------------------------

# Maximum number of GO-supported candidates retained during the initial
# category-specific screening before duplicate resolution.
candidate_pool_per_category <- 200

# If TRUE, each GO term listed below is expanded to include all descendant GO terms.
# If you want to use only the exact GO terms listed in pollen_categories, set FALSE.
expand_pollen_go_terms <- TRUE

# Output folder
dir.create("Heatmaps_Pollen_Categories", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 9B.1. Make sure GO ontology is available for semantic expansion/depth scoring
# ------------------------------------------------------------------------------

if (expand_pollen_go_terms && !exists("go")) {
  if (!file.exists("go-basic.obo")) {
    download.file(
      "https://current.geneontology.org/ontology/go-basic.obo",
      destfile = "go-basic.obo",
      mode = "wb"
    )
  }
  
  go <- ontologyIndex::get_ontology(
    "go-basic.obo",
    extract_tags = "everything"
  )
}

# ------------------------------------------------------------------------------
# 9B.2. GO categories for reproductive / pollen-related processes
# ------------------------------------------------------------------------------

pollen_categories <- list(
  
  Tapetum = c(
    "GO:0010234",
    "GO:0048658",
    "GO:0048657",
    "GO:0048656",
    "GO:0048655"
  ),
  
  Microsporogenesis = c(
    "GO:0055046",
    "GO:0010480",
    "GO:0009556",
    "GO:0010152",
    "GO:0048235"
  ),
  
  PollenWall = c(
    "GO:0062075",
    "GO:0010584",
    "GO:0010208",
    "GO:0160030",
    "GO:0080110"
  ),
  
  Dehiscence = c(
    "GO:0120194",
    "GO:0080166",
    "GO:0009901"
  )
)

# ------------------------------------------------------------------------------
# 9B.3. Utility functions
# ------------------------------------------------------------------------------

expand_go_terms_9c <- function(go_terms, go, expand = TRUE) {
  
  if (!expand) {
    return(unique(go_terms))
  }
  
  expanded_terms <- unique(unlist(lapply(go_terms, function(x) {
    
    out <- x
    
    if (exists("go") && x %in% names(go$name)) {
      descendants <- tryCatch(
        ontologyIndex::get_descendants(go, x),
        error = function(e) character(0)
      )
      out <- c(out, descendants)
    }
    
    return(out)
  })))
  
  unique(expanded_terms)
}

get_go_depth_9c <- function(go, go_ids) {
  
  depths <- sapply(go_ids, function(x) {
    
    if (!exists("go") || !x %in% names(go$name)) {
      return(NA_real_)
    }
    
    ancestors <- tryCatch(
      ontologyIndex::get_ancestors(go, x),
      error = function(e) character(0)
    )
    
    length(unique(ancestors))
  })
  
  as.numeric(depths)
}

row_zscore_9c <- function(mat) {
  
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  
  return(mat_z)
}

# ------------------------------------------------------------------------------
# 9B.4. Function: get candidate genes for each category
# ------------------------------------------------------------------------------

get_candidate_genes_by_cat_9c <- function(go_terms,
                                          deg_table,
                                          mapping_df,
                                          category_name,
                                          go = NULL,
                                          max_candidates = 200) {
  
  df <- as.data.frame(deg_table)
  
  if (!"gene_id" %in% colnames(df)) {
    df$gene_id <- rownames(df)
  }
  
  # Keep only valid DEG rows
  df <- df %>%
    filter(!is.na(gene_id)) %>%
    filter(!is.na(padj))
  
  # Count how many category GO terms annotate each gene.
  # This is used to assign duplicated genes automatically.
  target_gene_scores <- mapping_df %>%
    filter(go_id %in% go_terms) %>%
    distinct(gene_id, go_id) %>%
    group_by(gene_id) %>%
    summarise(
      n_matching_go_terms = n_distinct(go_id),
      matching_go_ids = paste(sort(unique(go_id)), collapse = ";"),
      mean_go_depth = if (exists("go")) {
        mean(get_go_depth_9c(go, unique(go_id)), na.rm = TRUE)
      } else {
        NA_real_
      },
      max_go_depth = if (exists("go")) {
        max(get_go_depth_9c(go, unique(go_id)), na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      mean_go_depth = ifelse(is.infinite(mean_go_depth), NA_real_, mean_go_depth),
      max_go_depth = ifelse(is.infinite(max_go_depth), NA_real_, max_go_depth)
    )
  
  df %>%
    inner_join(target_gene_scores, by = "gene_id") %>%
    mutate(
      Category = category_name,
      abs_log2FC = abs(log2FoldChange)
    ) %>%
    arrange(
      padj,
      desc(abs_log2FC),
      desc(n_matching_go_terms),
      desc(mean_go_depth),
      desc(max_go_depth)
    ) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    slice_head(n = max_candidates)
}

# ------------------------------------------------------------------------------
# 9B.5. Expand GO terms for each category
# ------------------------------------------------------------------------------

pollen_categories_expanded <- lapply(
  pollen_categories,
  function(x) expand_go_terms_9c(
    go_terms = x,
    go = go,
    expand = expand_pollen_go_terms
  )
)

# Save expanded GO terms used for transparency
expanded_go_table_9c <- bind_rows(lapply(names(pollen_categories_expanded), function(cat_name) {
  data.frame(
    Category = cat_name,
    GO_ID = pollen_categories_expanded[[cat_name]],
    stringsAsFactors = FALSE
  )
}))

write.csv(
  expanded_go_table_9c,
  file = "Heatmaps_Pollen_Categories/GO_terms_used_for_each_category.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9B.6. Build candidate pool for all GO-based categories
# ------------------------------------------------------------------------------

candidate_df <- bind_rows(lapply(names(pollen_categories_expanded), function(cat_name) {
  
  get_candidate_genes_by_cat_9c(
    go_terms = pollen_categories_expanded[[cat_name]],
    deg_table = all_degs_table,
    mapping_df = go_mapping_tidy,
    category_name = cat_name,
    go = go,
    max_candidates = candidate_pool_per_category
  )
}))

# ------------------------------------------------------------------------------
# 9B.7. Add SDR genes as separate user-defined biological group
#      They are included only if not already assigned to a GO-based category.
# ------------------------------------------------------------------------------

sdr_genes <- tibble(
  gene_id = c("Gb_15883", "Gb_28589", "novel.855", "Gb_28585", "Gb_11536"),
  Category = "SDR",
  n_matching_go_terms = NA_integer_,
  matching_go_ids = NA_character_,
  mean_go_depth = NA_real_,
  max_go_depth = NA_real_,
  baseMean = NA_real_,
  log2FoldChange = NA_real_,
  lfcSE = NA_real_,
  stat = NA_real_,
  pvalue = NA_real_,
  padj = NA_real_,
  abs_log2FC = NA_real_
)

candidate_df <- bind_rows(candidate_df, sdr_genes)

# ------------------------------------------------------------------------------
# 9B.8. Export duplicated candidates before automatic assignment
# ------------------------------------------------------------------------------

duplicated_candidates_before <- candidate_df %>%
  group_by(gene_id) %>%
  filter(n_distinct(Category) > 1) %>%
  arrange(gene_id, Category) %>%
  ungroup()

write.csv(
  duplicated_candidates_before,
  file = "Heatmaps_Pollen_Categories/Duplicated_candidates_before_automatic_assignment.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9B.9. Automatic duplicate resolution with biological category priority
#      One gene is assigned to one category only.
#
#      Biological priority:
#      1. SDR, if manually defined
#      2. Tapetum
#      3. PollenWall
#      4. Dehiscence
#      5. Microsporogenesis
#
#      This avoids assigning ambiguous genes to the broad pollen/microsporogenesis
#      category when they are also annotated as tapetum- or pollen-wall-related.
# ------------------------------------------------------------------------------

category_priority_9c <- c(
  "SDR",
  "Tapetum",
  "PollenWall",
  "Dehiscence",
  "Microsporogenesis"
)

candidate_df_unique <- candidate_df %>%
  mutate(
    category_rank = match(Category, category_priority_9c),
    category_rank = ifelse(is.na(category_rank), 999, category_rank),
    
    padj_for_sort = ifelse(is.na(padj), Inf, padj),
    abs_log2FC_for_sort = ifelse(is.na(abs_log2FC), -Inf, abs_log2FC),
    n_go_for_sort = ifelse(is.na(n_matching_go_terms), -Inf, n_matching_go_terms),
    mean_depth_for_sort = ifelse(is.na(mean_go_depth), -Inf, mean_go_depth),
    max_depth_for_sort = ifelse(is.na(max_go_depth), -Inf, max_go_depth)
  ) %>%
  group_by(gene_id) %>%
  arrange(
    category_rank,
    desc(n_go_for_sort),
    desc(mean_depth_for_sort),
    desc(max_depth_for_sort),
    padj_for_sort,
    desc(abs_log2FC_for_sort),
    .by_group = TRUE
  ) %>%
  slice_head(n = 1) %>%
  ungroup()

write.csv(
  candidate_df_unique,
  file = "Heatmaps_Pollen_Categories/Automatic_unique_assignment_all_candidates_priority_based.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 9B.10. Annotate and export the complete non-redundant candidate catalogue
# ------------------------------------------------------------------------------
# All candidates are retained. No additional top-N category-specific subset is
# generated, because the complete catalogue is supplied as Supplementary Table S1
# and the two curated manuscript heatmaps are generated below.

gene_annotations_unique <- gene_annotations %>%
  distinct(gene_id, .keep_all = TRUE)

candidate_df_unique_annotated <- candidate_df_unique %>%
  left_join(gene_annotations_unique, by = "gene_id") %>%
  mutate(
    full_name = ifelse(
      !is.na(gene_symbol) & gene_symbol != "",
      paste0(gene_id, " (", gene_symbol, ")"),
      gene_id
    )
  )

write.csv(
  candidate_df_unique_annotated,
  file = "Heatmaps_Pollen_Categories/Supplementary_Table_All_Reproductive_Candidates_by_Category.csv",
  row.names = FALSE
)

candidate_category_summary <- candidate_df_unique_annotated %>%
  dplyr::count(Category, name = "n_unique_candidates") %>%
  dplyr::arrange(match(Category, category_priority_9c))

write.csv(
  candidate_category_summary,
  file = "Heatmaps_Pollen_Categories/Reproductive_candidate_category_summary.csv",
  row.names = FALSE
)

print(candidate_category_summary)
cat("\nCandidate-gene catalogue completed and exported without top-N exploratory subsets.\n")

# ------------------------------------------------------------------------------
# 9C. Curated manuscript heatmaps
# ------------------------------------------------------------------------------
# The expanded GO workflow above is used as an objective screening and annotation
# layer. For the main manuscript, we avoid plotting all category-specific heatmaps
# and instead generate two curated, biologically interpretable heatmaps:
#   Figure 5: structural genes for pollen wall, sporopollenin, endothecium, dehiscence
#   Figure 6: reproductive regulators, tapetum markers, microsporogenesis, SDR genes
#
# IMPORTANT MANUSCRIPT CHOICE:
# - The full non-redundant candidate catalogue is exported as a supplementary table.
# - The two main heatmaps are curated from GO-supported candidates plus a restricted
#   set of literature-supported pathway/regulatory markers.
# - The structural heatmap retains only the strongest M3-enriched members of the
#   cellulose/endothecium module. The full CESA/COBL4/CTL2 set is still searched
#   and reported, but only the top M3-expressed candidates are plotted to avoid
#   overloading the manuscript figure.
# - The plotting style uses pheatmap, with gene labels on the right and row
#   category annotation on the left, matching the graphical style used in the
#   earlier version of the script.

# ----------------------------------------------------------------------------
# 9C.1. Helper functions for curated manuscript heatmaps
# ----------------------------------------------------------------------------

find_genes_by_symbol_9d <- function(symbols,
                                    heatmap_category,
                                    source_label = "symbol_based_required_marker") {
  if (!exists("gene_annotations_unique")) {
    stop("gene_annotations_unique is not available. Run the annotation section before 9C.")
  }

  annot <- gene_annotations_unique
  if (!("gene_symbol" %in% names(annot))) annot$gene_symbol <- NA_character_
  if (!("full_name" %in% names(annot))) annot$full_name <- NA_character_

  annot <- annot %>%
    mutate(
      gene_symbol_clean = toupper(trimws(as.character(gene_symbol))),
      full_name_clean = toupper(trimws(as.character(full_name)))
    )

  symbols_clean <- toupper(trimws(symbols))

  hits <- annot %>%
    filter(
      gene_symbol_clean %in% symbols_clean |
        full_name_clean %in% symbols_clean
    ) %>%
    transmute(
      gene_id,
      preferred_symbol = dplyr::case_when(
        !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
        TRUE ~ symbols[match(gene_symbol_clean, symbols_clean)]
      ),
      Heatmap_Category = heatmap_category,
      source_for_manuscript_heatmap = source_label
    ) %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE)

  missing_symbols <- setdiff(
    symbols_clean,
    unique(annot$gene_symbol_clean[annot$gene_symbol_clean %in% symbols_clean])
  )

  if (length(missing_symbols) > 0) {
    warning(
      "No exact gene_symbol match found for: ",
      paste(missing_symbols, collapse = ", "),
      ". These markers will be reported in Curated_heatmap_missing_requested_markers.csv."
    )
  }

  attr(hits, "missing_symbols") <- missing_symbols
  hits
}

add_if_in_vsd_9d <- function(df) {
  df %>%
    filter(!is.na(gene_id), gene_id != "") %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    filter(gene_id %in% rownames(vsd_mat))
}

make_display_label_9d <- function(df) {
  df %>%
    left_join(gene_annotations_unique, by = "gene_id") %>%
    mutate(
      plot_symbol = dplyr::case_when(
        !is.na(preferred_symbol) & preferred_symbol != "" ~ preferred_symbol,
        !is.na(gene_symbol) & gene_symbol != "" ~ gene_symbol,
        TRUE ~ NA_character_
      ),
      display_label = ifelse(
        !is.na(plot_symbol) & plot_symbol != "",
        paste0(gene_id, " (", plot_symbol, ")"),
        gene_id
      )
    )
}

make_pheatmap_manuscript_9d <- function(gene_df,
                                     output_svg,
                                     plot_title,
                                     width = 14,
                                     height = 12,
                                     category_colors,
                                     category_levels,
                                     cluster_rows = FALSE,
                                     use_stage_means = FALSE,
                                     show_all_replicates = TRUE) {

  gene_df <- gene_df %>%
    mutate(Heatmap_Category = factor(Heatmap_Category, levels = category_levels)) %>%
    arrange(Heatmap_Category, gene_id) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    add_if_in_vsd_9d() %>%
    make_display_label_9d()

  missing_genes <- setdiff(gene_df$gene_id, rownames(vsd_mat))
  if (length(missing_genes) > 0) {
    warning("The following genes are not present in the VST matrix and will be skipped: ",
            paste(missing_genes, collapse = ", "))
  }

  if (nrow(gene_df) == 0) {
    stop("No genes available for heatmap: ", output_svg)
  }

  if (isTRUE(use_stage_means) && !isTRUE(show_all_replicates)) {
    # Compact manuscript view: one mean-expression column per developmental stage.
    # Kept as an optional mode, but not used for the final manuscript figures.
    mat <- data.frame(
      M1 = rowMeans(vsd_mat[gene_df$gene_id, ann_col$Stage == "M1", drop = FALSE]),
      M2 = rowMeans(vsd_mat[gene_df$gene_id, ann_col$Stage == "M2", drop = FALSE]),
      M3 = rowMeans(vsd_mat[gene_df$gene_id, ann_col$Stage == "M3", drop = FALSE])
    )
    rownames(mat) <- gene_df$gene_id

    ann_col_manuscript <- data.frame(
      Stage = factor(colnames(mat), levels = c("M1", "M2", "M3"))
    )
    rownames(ann_col_premium) <- colnames(mat)
    gaps_col <- c(1, 2)

  } else {
    # Final manuscript view: all individual RNA-seq libraries are displayed.
    # Columns are ordered by developmental stage, with visible gaps between M1, M2, and M3.
    sample_order <- rownames(ann_col)[order(ann_col$Stage, rownames(ann_col))]
    mat <- vsd_mat[gene_df$gene_id, sample_order, drop = FALSE]

    ann_col_premium <- ann_col[sample_order, , drop = FALSE]
    ann_col_premium$Stage <- factor(ann_col_premium$Stage, levels = c("M1", "M2", "M3"))

    stage_runs <- rle(as.character(ann_col_premium$Stage))
    gaps_col <- cumsum(stage_runs$lengths)
    gaps_col <- gaps_col[-length(gaps_col)]
  }

  # Row-wise Z-scores are calculated across the displayed columns.
  # In the final figure this means across all 15 biological replicates, not across stage means.
  mat_z <- row_zscore_9c(as.matrix(mat))
  rownames(mat_z) <- make.unique(gene_df$display_label)

  ann_row <- data.frame(
    Category = gene_df$Heatmap_Category
  )
  rownames(ann_row) <- rownames(mat_z)

  # Keep rows grouped by category and add visible gaps between categories.
  if (isFALSE(cluster_rows)) {
    category_runs <- rle(as.character(gene_df$Heatmap_Category))
    gaps_row <- cumsum(category_runs$lengths)
    gaps_row <- gaps_row[-length(gaps_row)]
  } else {
    gaps_row <- NULL
  }

  svg(output_svg, width = width, height = height)
  pheatmap::pheatmap(
    mat_z,
    annotation_col = ann_col_premium,
    annotation_row = ann_row,
    annotation_colors = list(
      Stage = ann_colors$Stage,
      Category = category_colors
    ),
    cluster_cols = FALSE,
    cluster_rows = cluster_rows,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 14,
    fontsize_col = 10,
    fontsize = 12,
    angle_col = 45,
    main = plot_title,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100),
    border_color = "white",
    gaps_row = gaps_row,
    gaps_col = gaps_col
  )
  dev.off()

  invisible(gene_df)
}

# ----------------------------------------------------------------------------
# 9C.2. Curated regulatory module for Figure 6
# ----------------------------------------------------------------------------

manuscript_regulatory_genes <- tibble::tribble(
  ~gene_id,     ~preferred_symbol, ~Heatmap_Category,                 ~source_for_manuscript_heatmap,
  "Gb_38764",   "EMS1",           "Early identity / tapetum",        "curated_GO_supported_or_literature",
  "Gb_21147",   "EMS1",           "Early identity / tapetum",        "curated_GO_supported_or_literature",
  "novel.1129", "TPD1",           "Early identity / tapetum",        "manual_literature_marker",
  "Gb_03922",   "SPL2",           "Early identity / tapetum",        "manual_literature_marker",
  "Gb_15605",   "SPL8",           "Early identity / tapetum",        "curated_GO_supported_or_literature",
  "novel.216",  "MADS6-like",     "Early identity / tapetum",        "manual_literature_marker",

  "novel.552",  "MMD1",           "Meiosis / microsporogenesis",     "curated_GO_supported_or_literature",
  "novel.943",  "MMD1",           "Meiosis / microsporogenesis",     "manual_literature_marker",
  "novel.431",  "PAIR1",          "Meiosis / microsporogenesis",     "manual_literature_marker",
  "novel.337",  "PAIR2",          "Meiosis / microsporogenesis",     "manual_literature_marker",
  "Gb_01637",   "AUG6",           "Meiosis / microsporogenesis",     "manual_literature_marker",
  "Gb_38139",   "JASON",          "Meiosis / microsporogenesis",     "curated_GO_supported_or_literature",
  "Gb_38936",   "JASON",          "Meiosis / microsporogenesis",     "curated_GO_supported_or_literature",

  "novel.1555", "AMS",            "Tapetum regulation / PCD",        "curated_GO_supported_or_literature",
  "novel.1876", "MYB80",          "Tapetum regulation / PCD",        "curated_GO_supported_or_literature",
  "novel.1801", "MYB80",          "Tapetum regulation / PCD",        "curated_GO_supported_or_literature",
  "novel.942",  "PTC1",           "Tapetum regulation / PCD",        "manual_literature_marker",
  "Gb_19452",   "EAT1",           "Tapetum regulation / PCD",        "curated_GO_supported_or_literature",
  "Gb_10444",   "CEP1",           "Tapetum regulation / PCD",        "curated_GO_supported_or_literature",
  "novel.295",  "API5",           "Tapetum regulation / PCD",        "manual_literature_marker",

  "Gb_23921",   "MYB101",         "Late pollen maturation / ABA",    "curated_GO_supported_or_literature",
  "Gb_31417",   "AGL104",         "Late pollen maturation / ABA",    "curated_GO_supported_or_literature",
  "novel.1622", "LBD27",          "Late pollen maturation / ABA",    "curated_GO_supported_or_literature",
  "Gb_16414",   "ZAT3",           "Late pollen maturation / ABA",    "curated_GO_supported_or_literature",
  "Gb_07367",   "NCED1",          "Late pollen maturation / ABA",    "manual_literature_marker",

  "Gb_15883",   "ORR24",          "SDR",                             "manual_SDR_marker",
  "Gb_28589",   "AB17B",          "SDR",                             "manual_SDR_marker",
  "novel.855",  "GGM13",          "SDR",                             "manual_SDR_marker",
  "Gb_28585",   "NFD4-like",      "SDR",                             "manual_SDR_marker",
  "Gb_11536",   "GAMYB",          "SDR",                             "manual_SDR_marker"
)

# ----------------------------------------------------------------------------
# 9C.3. Curated structural module for Figure 5
# ----------------------------------------------------------------------------

# Base structural list: pollen wall, sporopollenin, surface lipids, secondary wall,
# and dehiscence markers.
manuscript_structural_base <- tibble::tribble(
  ~gene_id,     ~preferred_symbol, ~Heatmap_Category,                              ~source_for_manuscript_heatmap,
  "Gb_38324",   "PME53",          "Early wall remodeling / callose",              "curated_GO_supported_or_literature",
  "Gb_12380",   "LRX3",           "Early wall remodeling / callose",              "curated_GO_supported_or_literature",
  "Gb_27977",   "B3GT3",          "Early wall remodeling / callose",              "curated_GO_supported_or_literature",
  "Gb_08787",   "CALS5",          "Callose / tetrad release",                     "curated_GO_supported_or_literature",
  "Gb_08788",   "CALS5",          "Callose / tetrad release",                     "curated_GO_supported_or_literature",
  "Gb_20660",   "QRT3",           "Callose / tetrad release",                     "curated_GO_supported_or_literature",
  "Gb_36473",   "A6-like",        "Callose / tetrad release",                     "curated_GO_supported_or_literature",
  "Gb_09137",   "INP1",           "Pollen aperture / wall patterning",            "manual_literature_marker",

  "Gb_15188",   "CYP703-like",    "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "Gb_34128",   "CYP704-like",    "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "Gb_02579",   "PKSB",           "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "Gb_29977",   "PKSB",           "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "Gb_15917",   "TKPR1",          "Sporopollenin synthesis / transport",          "manual_literature_marker",
  "Gb_31437",   "TKPR2",          "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "novel.2029", "4CLL1",          "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "Gb_31075",   "FACR2",          "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",
  "novel.949",  "SSL13",          "Sporopollenin synthesis / transport",          "manual_literature_marker",
  "Gb_00289",   "ABCG26-like",    "Sporopollenin synthesis / transport",          "curated_GO_supported_or_literature",

  "Gb_38453",   "HOTHEAD-like",   "Pollen coat / surface lipids",                 "manual_literature_marker",
  "novel.794",  "CER2-like",      "Pollen coat / surface lipids",                 "manual_literature_marker",

  "novel.1019", "LAC4",           "Endothecium secondary wall / lignin",           "curated_GO_supported_or_literature",
  "Gb_14613",   "LAC17",          "Endothecium secondary wall / lignin",           "curated_GO_supported_or_literature",

  "Gb_06127",   "EXPA1",          "Dehiscence / wall loosening",                  "curated_GO_supported_or_literature",
  "Gb_10500",   "XTH2",           "Dehiscence / wall loosening",                  "curated_GO_supported_or_literature",
  "Gb_12091",   "BGAL1",          "Dehiscence / wall loosening",                  "curated_GO_supported_or_literature",
  "novel.325",  "PGLR",           "Dehiscence / wall loosening",                  "curated_GO_supported_or_literature",
  "Gb_36896",   "ATL73",          "Dehiscence / wall loosening",                  "curated_GO_supported_or_literature",
  "Gb_24167",   "DTX35",          "Dehiscence / wall loosening",                  "curated_GO_supported_or_literature"
)

# Dynamically recover the cellulose/endothecium secondary-wall module by gene
# symbol, but retain in the manuscript heatmap only the strongest M3-enriched
# candidates. This keeps the CESA/COBL4/CTL2 narrative while avoiding an
# overloaded Figure 5. The full search/selection report is exported below.
required_cellulose_symbols <- c("CESA1", "CESA3", "CESA4", "CESA5", "CESA7", "CESA9", "COBL4", "CTL2")

required_cellulose_markers_all <- find_genes_by_symbol_9d(
  symbols = required_cellulose_symbols,
  heatmap_category = "Endothecium secondary wall / cellulose",
  source_label = "candidate_cellulose_secondary_wall_marker"
)

select_top_m3_cellulose_markers_9d <- function(df, max_markers = 4, min_m3_z = 0.5) {
  df <- df %>%
    add_if_in_vsd_9d() %>%
    distinct(gene_id, .keep_all = TRUE)

  if (nrow(df) == 0) return(df)

  stage_means <- data.frame(
    gene_id = df$gene_id,
    M1_mean = rowMeans(vsd_mat[df$gene_id, ann_col$Stage == "M1", drop = FALSE]),
    M2_mean = rowMeans(vsd_mat[df$gene_id, ann_col$Stage == "M2", drop = FALSE]),
    M3_mean = rowMeans(vsd_mat[df$gene_id, ann_col$Stage == "M3", drop = FALSE])
  )

  stage_z <- row_zscore_9c(as.matrix(stage_means[, c("M1_mean", "M2_mean", "M3_mean")]))
  colnames(stage_z) <- c("M1_z", "M2_z", "M3_z")

  stage_scores <- bind_cols(stage_means, as.data.frame(stage_z)) %>%
    rowwise() %>%
    mutate(
      peak_stage = c("M1", "M2", "M3")[which.max(c(M1_mean, M2_mean, M3_mean))]
    ) %>%
    ungroup()

  df_scored <- df %>%
    left_join(stage_scores, by = "gene_id")

  selected <- df_scored %>%
    filter(peak_stage == "M3", M3_z >= min_m3_z) %>%
    arrange(desc(M3_z), desc(M3_mean)) %>%
    slice_head(n = max_markers)

  # Fallback: if the z-score threshold is too strict for a particular dataset,
  # still keep the top M3-peaking candidates so that the cellulose module is shown.
  if (nrow(selected) == 0) {
    selected <- df_scored %>%
      filter(peak_stage == "M3") %>%
      arrange(desc(M3_z), desc(M3_mean)) %>%
      slice_head(n = max_markers)
  }

  selected %>%
    select(gene_id, preferred_symbol, Heatmap_Category, source_for_manuscript_heatmap)
}

# Adjust this number if the figure is still too dense or too sparse.
max_m3_cellulose_markers_to_plot <- 4
min_m3_z_for_cellulose_marker <- 0.5

required_cellulose_markers <- select_top_m3_cellulose_markers_9d(
  required_cellulose_markers_all,
  max_markers = max_m3_cellulose_markers_to_plot,
  min_m3_z = min_m3_z_for_cellulose_marker
)

manuscript_structural_genes <- bind_rows(
  manuscript_structural_base,
  required_cellulose_markers
) %>%
  distinct(gene_id, .keep_all = TRUE)

# Export a transparent report of which cellulose/endothecium markers were searched,
# which ones were selected for the main heatmap, and which symbols were not found.
cellulose_marker_selection_report <- required_cellulose_markers_all %>%
  add_if_in_vsd_9d() %>%
  left_join(
    data.frame(
      gene_id = .$gene_id,
      M1_mean = rowMeans(vsd_mat[.$gene_id, ann_col$Stage == "M1", drop = FALSE]),
      M2_mean = rowMeans(vsd_mat[.$gene_id, ann_col$Stage == "M2", drop = FALSE]),
      M3_mean = rowMeans(vsd_mat[.$gene_id, ann_col$Stage == "M3", drop = FALSE])
    ),
    by = "gene_id"
  ) %>%
  mutate(
    selected_for_manuscript_heatmap = gene_id %in% required_cellulose_markers$gene_id
  )

missing_required_cellulose_markers <- tibble::tibble(
  requested_symbol = attr(required_cellulose_markers_all, "missing_symbols")
)

write.csv(
  cellulose_marker_selection_report,
  file = "Curated_heatmap_cellulose_M3_marker_selection_report.csv",
  row.names = FALSE
)

write.csv(
  missing_required_cellulose_markers,
  file = "Curated_heatmap_missing_requested_markers.csv",
  row.names = FALSE
)

# ----------------------------------------------------------------------------
# 9C.4. Colors and plotting
# ----------------------------------------------------------------------------

regulatory_category_colors <- c(
  "Early identity / tapetum" = "#6A3D9A",
  "Meiosis / microsporogenesis" = "#1F78B4",
  "Tapetum regulation / PCD" = "#33A02C",
  "Late pollen maturation / ABA" = "#FF7F00",
  "SDR" = "#E31A1C"
)

structural_category_colors <- c(
  "Early wall remodeling / callose" = "#8C510A",
  "Callose / tetrad release" = "#BF812D",
  "Pollen aperture / wall patterning" = "#80CDC1",
  "Sporopollenin synthesis / transport" = "#01665E",
  "Pollen coat / surface lipids" = "#C51B7D",
  "Endothecium secondary wall / cellulose" = "#5E3C99",
  "Endothecium secondary wall / lignin" = "#8073AC",
  "Dehiscence / wall loosening" = "#D95F02"
)

manuscript_regulatory_annotated <- make_pheatmap_manuscript_9d(
  gene_df = manuscript_regulatory_genes,
  output_svg = "Figure6_reproductive_regulators_heatmap.svg",
  plot_title = "Expression Profile of Pollen Development Key Regulatory Genes",
  width = 14,
  height = 9,
  category_colors = regulatory_category_colors,
  category_levels = names(regulatory_category_colors),
  cluster_rows = FALSE,
  use_stage_means = FALSE,
  show_all_replicates = TRUE
)

manuscript_structural_annotated <- make_pheatmap_manuscript_9d(
  gene_df = manuscript_structural_genes,
  output_svg = "Figure5_structural_pollen_wall_dehiscence_heatmap.svg",
  plot_title = "Expression Profile of Cell Wall, Pollen Wall and Dehiscence Genes",
  width = 15,
  height = 10,
  category_colors = structural_category_colors,
  category_levels = names(structural_category_colors),
  cluster_rows = FALSE,
  use_stage_means = FALSE,
  show_all_replicates = TRUE
)

manuscript_gene_selection <- bind_rows(
  manuscript_regulatory_annotated %>% mutate(Figure = "Figure 6 - reproductive regulators"),
  manuscript_structural_annotated %>% mutate(Figure = "Figure 5 - structural pollen wall/dehiscence")
) %>%
  select(Figure, gene_id, preferred_symbol, gene_symbol, plot_symbol,
         display_label, Heatmap_Category, source_for_manuscript_heatmap, everything())

write.csv(
  manuscript_gene_selection,
  file = "Curated_heatmap_gene_selection.csv",
  row.names = FALSE
)

cat("\nCurated manuscript heatmaps generated with all 15 biological-replicate libraries shown.\n")
# ------------------------------------------------------------------------------
# 10. Manuscript-ready summary statistics
# ------------------------------------------------------------------------------
# This final block creates a compact set of internal summary files that can be used
# to report quantitative information in the Results section and figure legends.
# These files are intended as writing aids / quality-control summaries, not as a
# large mandatory set of supplementary tables. For the manuscript submission, the
# essential supplementary material can remain limited to:
#   1) this complete R script;
#   2) the raw count matrix;
#   3) the complete reproductive candidate-gene table used to curate Figures 5-6.

cat("\nStarting Section 10: manuscript-ready summary statistics...\n")

dir.create("Manuscript_Summary_Stats", showWarnings = FALSE)

# ---- 10.0 Robust helper functions ----
# Avoid ambiguous calls to count(), which can be masked by other packages, and
# safely convert possible list-columns to character vectors before summarising.
safe_as_character <- function(x) {
  if (is.list(x)) {
    vapply(x, function(z) paste(as.character(z), collapse = ";"), character(1))
  } else {
    as.character(x)
  }
}

safe_write_csv <- function(df, file) {
  readr::write_csv(as.data.frame(df), file)
  cat("Written:", file, "\n")
}

# ---- 10.1 Dataset, filtering and PCA summary ----
dataset_summary <- tibble::tibble(
  metric = c(
    "raw_genes_or_transcripts_in_count_matrix",
    "genes_or_transcripts_retained_after_low_count_filter",
    "genes_or_transcripts_removed_by_low_count_filter",
    "total_RNAseq_libraries",
    "M1_biological_replicates",
    "M2_biological_replicates",
    "M3_biological_replicates",
    "PCA_PC1_percent_variance",
    "PCA_PC2_percent_variance"
  ),
  value = c(
    nrow(count_matrix),
    nrow(dds),
    nrow(count_matrix) - nrow(dds),
    ncol(count_matrix),
    sum(coldata$condition == "M1"),
    sum(coldata$condition == "M2"),
    sum(coldata$condition == "M3"),
    percentVar[1],
    percentVar[2]
  )
)

safe_write_csv(
  dataset_summary,
  "Manuscript_Summary_Stats/01_dataset_filtering_PCA_summary.csv"
)

# ---- 10.2 Differential-expression summary by contrast ----
summarise_deg_contrast <- function(res,
                                   contrast_name,
                                   padj_cutoff = 0.005,
                                   lfc_cutoff = 1.5) {
  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id")

  sig <- df %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::filter(padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)

  tibble::tibble(
    contrast = contrast_name,
    tested_genes_with_nonNA_padj = sum(!is.na(df$padj)),
    significant_DEGs_total = nrow(sig),
    upregulated_in_numerator_stage = sum(sig$log2FoldChange > lfc_cutoff),
    downregulated_in_numerator_stage = sum(sig$log2FoldChange < -lfc_cutoff),
    padj_cutoff = padj_cutoff,
    abs_log2FC_cutoff = lfc_cutoff
  )
}

deg_summary <- dplyr::bind_rows(
  summarise_deg_contrast(res_M2_vs_M1, "M2_vs_M1"),
  summarise_deg_contrast(res_M3_vs_M2, "M3_vs_M2"),
  summarise_deg_contrast(res_M3_vs_M1, "M3_vs_M1"),
  tibble::tibble(
    contrast = "Union_all_pairwise_contrasts",
    tested_genes_with_nonNA_padj = NA_integer_,
    significant_DEGs_total = length(all_degs_ids),
    upregulated_in_numerator_stage = NA_integer_,
    downregulated_in_numerator_stage = NA_integer_,
    padj_cutoff = 0.005,
    abs_log2FC_cutoff = 1.5
  )
)

safe_write_csv(
  deg_summary,
  "Manuscript_Summary_Stats/02_DEG_summary_by_contrast.csv"
)

# ---- 10.3 Mfuzz cluster summary ----
cluster_assignment_table <- tibble::tibble(
  gene_id = names(cl$cluster),
  Cluster = paste0("Cluster_", as.integer(cl$cluster)),
  max_membership = apply(cl$membership, 1, max)
)

cluster_summary <- cluster_assignment_table %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(
    n_genes_assigned = dplyr::n(),
    n_genes_membership_gt_0_5 = sum(max_membership > 0.5),
    mean_max_membership = mean(max_membership),
    median_max_membership = median(max_membership),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Cluster)

safe_write_csv(
  cluster_summary,
  "Manuscript_Summary_Stats/03_Mfuzz_cluster_summary.csv"
)

# ---- 10.4 GO enrichment summary with Benjamini-Hochberg correction ----
# BH-adjusted values were calculated once, immediately after each enrichment
# analysis, and are reused here without recomputation.
format_enrichment_summary <- function(df, analysis_type) {
  df %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Analysis = analysis_type,
      Group = safe_as_character(Group),
      Category = safe_as_character(Category),
      significant_nominal_p_0_05 = p_value < 0.05
    ) %>%
    dplyr::select(
      Analysis, Group, Category, n_DEG_hits, p_value, p_adj_BH,
      log10_p, log10_p_adj_BH,
      significant_nominal_p_0_05, significant_BH_0_05
    ) %>%
    dplyr::arrange(Analysis, Group, p_adj_BH, p_value)
}

go_enrichment_summary <- dplyr::bind_rows(
  format_enrichment_summary(
    global_enrichment,
    "Pairwise_DEG_GO_category_enrichment"
  ),
  format_enrichment_summary(
    all_cluster_enrichment,
    "Mfuzz_cluster_GO_category_enrichment"
  )
)

safe_write_csv(
  go_enrichment_summary,
  "Manuscript_Summary_Stats/04_GO_enrichment_summary_with_BH.csv"
)

# ---- 10.5 Candidate-gene and curated heatmap summary ----
# The complete candidate-gene table is already exported in Section 9B as:
# Heatmaps_Pollen_Categories/Supplementary_Table_All_Reproductive_Candidates_by_Category.csv
# Here we only create a compact summary of candidate numbers and plotted genes.

candidate_summary <- candidate_df_unique_annotated %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Category = safe_as_character(Category)) %>%
  dplyr::group_by(Category) %>%
  dplyr::summarise(
    n_unique_candidates = dplyr::n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Category) %>%
  dplyr::mutate(Summary_type = "complete_candidate_catalogue") %>%
  dplyr::rename(Group = Category, n_genes = n_unique_candidates)

curated_heatmap_summary <- manuscript_gene_selection %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Figure = safe_as_character(Figure),
    Heatmap_Category = safe_as_character(Heatmap_Category)
  ) %>%
  dplyr::group_by(Figure, Heatmap_Category) %>%
  dplyr::summarise(
    n_genes = dplyr::n_distinct(gene_id),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Figure, Heatmap_Category) %>%
  dplyr::mutate(
    Summary_type = paste0("curated_heatmap_genes__", Figure),
    Group = Heatmap_Category
  ) %>%
  dplyr::select(Summary_type, Group, n_genes)

candidate_and_heatmap_summary <- dplyr::bind_rows(
  candidate_summary %>% dplyr::select(Summary_type, Group, n_genes),
  curated_heatmap_summary
)

safe_write_csv(
  candidate_and_heatmap_summary,
  "Manuscript_Summary_Stats/05_candidate_and_heatmap_gene_summary.csv"
)

# ---- 10.6 AGP/FLA candidate summary ----
agp_detected <- unique(agp_fla_genes)[unique(agp_fla_genes) %in% rownames(vsd_mat)]
agp_missing <- setdiff(unique(agp_fla_genes), rownames(vsd_mat))

agp_summary <- tibble::tibble(
  metric = c(
    "AGP_FLA_candidates_in_curated_domain_based_list",
    "AGP_FLA_candidates_detected_in_filtered_VST_matrix",
    "AGP_FLA_candidates_missing_from_filtered_VST_matrix",
    "detected_gene_ids",
    "missing_gene_ids"
  ),
  value = c(
    length(unique(agp_fla_genes)),
    length(agp_detected),
    length(agp_missing),
    paste(agp_detected, collapse = ";"),
    paste(agp_missing, collapse = ";")
  )
)

safe_write_csv(
  agp_summary,
  "Manuscript_Summary_Stats/06_AGP_FLA_candidate_summary.csv"
)

# ---- 10.7 Short text report: numbers to copy into Results / figure legends ----
get_metric <- function(summary_df, metric_name) {
  summary_df$value[summary_df$metric == metric_name][1]
}

report_lines <- c(
  "MANUSCRIPT-READY RNA-SEQ SUMMARY STATISTICS",
  "================================================",
  "",
  paste0("Libraries: ", get_metric(dataset_summary, "total_RNAseq_libraries"),
         " total; M1 = ", get_metric(dataset_summary, "M1_biological_replicates"),
         ", M2 = ", get_metric(dataset_summary, "M2_biological_replicates"),
         ", M3 = ", get_metric(dataset_summary, "M3_biological_replicates"), "."),
  paste0("Genes/transcripts retained after low-count filtering: ",
         get_metric(dataset_summary, "genes_or_transcripts_retained_after_low_count_filter"),
         " out of ", get_metric(dataset_summary, "raw_genes_or_transcripts_in_count_matrix"), "."),
  paste0("PCA: PC1 = ", get_metric(dataset_summary, "PCA_PC1_percent_variance"),
         "% variance; PC2 = ", get_metric(dataset_summary, "PCA_PC2_percent_variance"), "% variance."),
  "",
  "DEG summary:",
  paste(capture.output(print(deg_summary)), collapse = "\n"),
  "",
  "Mfuzz cluster summary:",
  paste(capture.output(print(cluster_summary)), collapse = "\n"),
  "",
  "Candidate and heatmap gene summary:",
  paste(capture.output(print(candidate_and_heatmap_summary)), collapse = "\n"),
  "",
  "AGP/FLA summary:",
  paste(capture.output(print(agp_summary)), collapse = "\n"),
  "",
  "Note: GO enrichment summary with BH correction is saved as 04_GO_enrichment_summary_with_BH.csv."
)

writeLines(
  report_lines,
  "Manuscript_Summary_Stats/00_manuscript_numbers_short_report.txt"
)
cat("Written: Manuscript_Summary_Stats/00_manuscript_numbers_short_report.txt\n")

# ---- 10.8 R session information ----
writeLines(
  capture.output(sessionInfo()),
  "Manuscript_Summary_Stats/R_sessionInfo.txt"
)
cat("Written: Manuscript_Summary_Stats/R_sessionInfo.txt\n")

cat("\nSection 10 completed: compact manuscript-ready summary statistics generated.\n")
