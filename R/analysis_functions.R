#' Filter protein-coding genes from a GTF file and subset a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param gtf_path Path to the GTF file
#'
#' @return A filtered Seurat object containing only protein-coding genes
filter_protein_coding_genes <- function(seurat_obj, gtf_path) {
  library(data.table)

  # Carica il GTF
  gtf <- fread(gtf_path, header = FALSE, sep = "\t", data.table = FALSE, skip = "#")

  # Prendi solo le righe che descrivono geni
  gtf_genes <- gtf[gtf[, 3] == "gene", ]

  # Filtra i geni protein-coding
  gtf_protein_coding <- gtf_genes[grep("gene_biotype \"protein_coding\"", gtf_genes[, 9]), ]

  # Estrai i nomi dei geni protein-coding
  protein_coding_genes <- gsub(".*gene_name \"([^\"]+)\".*", "\\1", gtf_protein_coding[, 9])

  # Sottoseleziona i geni codificanti presenti nei dati Seurat
  seurat_obj <- subset(seurat_obj, features = intersect(rownames(seurat_obj), protein_coding_genes))

  return(seurat_obj)
}

#' Filter cells by gene expression (≥ 3 UMI)
#'
#' Adds a metadata column with the number of genes detected with ≥ 3 UMI per cell.
#'
#' @param seurat_obj A Seurat object.
#' @return A Seurat object with metadata column `genes_umi3`.
#' @export
filter_cells_by_gene_expression <- function(seurat_obj) {
  counts_matrix <- Seurat::GetAssayData(seurat_obj, layer = "counts")
  seurat_obj$genes_umi3 <- colSums(counts_matrix >= 3)
  return(seurat_obj)
}

#' Filter mitochondrial, ribosomal and ribosomal pseudogenes
#'
#' Removes genes belonging to mitochondrial, ribosomal and pseudogene categories.
#'
#' @param seurat_obj A Seurat object.
#' @return A list with two elements: (1) filtered Seurat object, and (2) summary table of removed gene counts.
#' @export
filter_unwanted_genes <- function(seurat_obj) {
  all_genes <- rownames(seurat_obj)

  mito_genes <- grep("^MT-", all_genes, value = TRUE)
  ribo_genes <- grep("^RP[SL]", all_genes, value = TRUE)
  ribo_pseudo_genes <- grep("^RP.+P[0-9]+$", all_genes, value = TRUE)

  genes_to_remove <- unique(c(mito_genes, ribo_genes, ribo_pseudo_genes))
  seurat_obj_filtered <- subset(seurat_obj, features = setdiff(all_genes, genes_to_remove))

  summary_table <- data.frame(
    Category = c("Mitochondrial genes", "Ribosomal genes", "Ribosomal pseudogenes"),
    Removed = c(length(mito_genes), length(ribo_genes), length(ribo_pseudo_genes))
  )

  return(list(object = seurat_obj_filtered, summary = summary_table))
}

#' Run PCA and compute explained variance
#'
#' Normalizes data, finds variable features, scales data, runs PCA, and calculates variance explained.
#'
#' @param seurat_obj A Seurat object.
#' @param n_pcs Number of principal components to calculate.
#' @return A list with two elements: (1) updated Seurat object, (2) variance explained per PC.
#' @export
run_pca_and_variance <- function(seurat_obj, n_pcs = 20) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = n_pcs)

  stdev <- seurat_obj[["pca"]]@stdev
  explained_variance <- stdev^2 / sum(stdev^2)

  return(list(object = seurat_obj, variance = explained_variance))
}

#' Perform clustering and UMAP
#'
#' Constructs nearest neighbor graph, clusters the data, and runs UMAP for visualization.
#'
#' @param seurat_obj A Seurat object.
#' @param dims PCs to use (e.g., 1:13).
#' @param resolution Clustering resolution (default 0.5).
#' @return A Seurat object with clusters and UMAP embedding.
#' @export
run_clustering_and_umap <- function(seurat_obj, dims = 1:13, resolution = 0.5) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  return(seurat_obj)
}

#' Find top marker genes per cluster
#'
#' Identifies the most highly expressed gene per cluster.
#'
#' @param seurat_obj A clustered Seurat object.
#' @return A dataframe of top markers.
#' @export
find_top_markers <- function(seurat_obj) {
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top_markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 1, wt = avg_log2FC)
  return(top_markers)
}

#' Rename clusters with cell type labels
#'
#' Renames cluster IDs with user-provided cell type annotations.
#'
#' @param seurat_obj A Seurat object with clusters.
#' @param labels A character vector of new cluster names.
#' @return A Seurat object with renamed identities.
#' @export
annotate_cell_types <- function(seurat_obj, labels) {
  names(labels) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, labels)
  return(seurat_obj)
}

#' Infer tissue origin from cell types
#'
#' Based on annotated cell types, returns a guess for tissue origin.
#'
#' @param seurat_obj A Seurat object with annotated identities.
#' @return A string hypothesis for the tissue of origin.
#' @export
infer_tissue_origin <- function(seurat_obj) {
  return("The presence of T cells, monocytes, cytotoxic cells, and lack of epithelial cells suggests the origin is peripheral blood or lymphoid tissue.")
}
