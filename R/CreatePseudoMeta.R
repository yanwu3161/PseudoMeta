CreatePseudoMeta <- function(SeuratObj, CDSObj, Varnum = 1000, origindata = FALSE, PCA = FALSE) {
  required_packages <- c("Seurat", "dplyr", "ggplot2", "viridis")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package_", pkg, "_not installed or loaded. Please install it to use PseudoMeta."))
    } else {
      library(pkg, character.only = TRUE)
    }
  }
  if (!isClass("PseudoMeta")) {
    setClass(
      "PseudoMeta",
      slots = c(
        assays = "list",
        analysis = "list",
        reductions = "list",
        meta.data = "data.frame",
        ID = "character"
      )
    )
  }

  total_var_features <- length(SeuratObj@assays[["RNA"]]@var.features)
  if (Varnum > total_var_features) {
    stop("NUM ERROR: Not enough var_features in SeuratObj")
  }

  var_features <- SeuratObj@assays[["RNA"]]@var.features[1:Varnum]

  assays <- list(
    origin.data = if (origindata) SeuratObj@assays[["RNA"]]@data else NULL,
    var.matrix = as.data.frame(GetAssayData(SeuratObj, assay = "RNA")[var_features, , drop = FALSE]),
    var.features = var_features,
    var.feature.data = data.frame(gene = var_features, row.names = 1:Varnum)
  )
  analysis <- list(
    Spearman = list(),
    Kendall = list(),
    Mann_Kendall = list(),
    Page = list()
  )
  reductions <- list(
    PCA = if (PCA) as.data.frame(SeuratObj@reductions[["pca"]]@cell.embeddings) else NULL,
    UMAP = as.data.frame(SeuratObj@reductions[["umap"]]@cell.embeddings)
  )

  meta.data <-  mutate(as.data.frame(SeuratObj@meta.data),
                       setNames(
                         as.data.frame(
                           CDSObj@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]),"pseudotime"))

  PseudoMeta <- new(
    "PseudoMeta",
    assays = assays,
    analysis = analysis,
    reductions = reductions,
    meta.data = meta.data,
    ID = digest::digest(paste(var_features, collapse = ""))
  )
  return(PseudoMeta)
}
