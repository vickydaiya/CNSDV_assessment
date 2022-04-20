#install.packages("BiocManager")
#install.packages("devtools")

library(devtools)
#install_github("satijalab/azimuth")
#!/usr/bin/env Rscript

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
  stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
}

# Ensure glmGamPoi is installed
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}

# Ensure Azimuth is installed
if (packageVersion(pkg = "Azimuth") < package_version(x = "0.3.1")) {
  stop("Please install azimuth - remotes::install_github('satijalab/azimuth')", call. = FALSE)
}

library(Seurat)
library(Azimuth)
library(dplyr)
library(readr)
# Download the Azimuth reference and extract the archive
reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_motorcortex")

r <- function(file,outputfile)
{
  # Load the reference
  # Change the file path based on where the reference is located on your system.
  
  # Load the query object for mapping
  # Change the file path based on where the query file is located on your system.
  query <- LoadFileInput(path = file)
  
  # Calculate nCount_RNA and nFeature_RNA if the query does not
  # contain them already
  if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
    calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
    colnames(x = calcn) <- paste(
      colnames(x = calcn),
      "RNA",
      sep = '_'
    )
    query <- AddMetaData(
      object = query,
      metadata = calcn
    )
    rm(calcn)
  }
  
  # Calculate percent mitochondrial genes if the query contains genes
  # matching the regular expression "^MT-"
  if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
    query <- PercentageFeatureSet(
      object = query,
      pattern = '^MT-',
      col.name = 'percent.mt',
      assay = "RNA"
    )
  }
  
  # Preprocess with SCTransform
  query <- SCTransform(
    object = query,
    assay = "RNA",
    new.assay.name = "refAssay",
    residual.features = rownames(x = reference$map),
    reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
    method = 'glmGamPoi',
    ncells = 2000,
    n_genes = 2000,
    do.correct.umi = FALSE,
    do.scale = FALSE,
    do.center = TRUE
  )
  
  # Find anchors between query and reference
  anchors <- FindTransferAnchors(
    reference = reference$map,
    query = query,
    k.filter = NA,
    reference.neighbors = "refdr.annoy.neighbors",
    reference.assay = "refAssay",
    query.assay = "refAssay",
    reference.reduction = "refDR",
    normalization.method = "SCT",
    features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
    dims = 1:50,
    n.trees = 20,
    mapping.score.k = 100
  )
  
  # Transfer cell type labels and impute protein expression
  #
  # Transferred labels are in metadata columns named "predicted.*"
  # The maximum prediction score is in a metadata column named "predicted.*.score"
  # The prediction scores for each class are in an assay named "prediction.score.*"
  # The imputed assay is named "impADT" if computed
  
  refdata <- lapply(X = c("subclass", "class", "cluster", "cross_species_cluster"), function(x) {
    reference$map[[x, drop = TRUE]]
  })
  names(x = refdata) <- c("subclass", "class", "cluster", "cross_species_cluster")
  if (FALSE) {
    refdata[["impADT"]] <- GetAssayData(
      object = reference$map[['ADT']],
      slot = 'data'
    )
  }
  query <- TransferData(
    reference = reference$map,
    query = query,
    dims = 1:50,
    anchorset = anchors,
    refdata = refdata,
    n.trees = 20,
    store.weights = TRUE
  )
  
  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors,
    reference = reference$map,
    query = query,
    reductions = "pcaproject",
    reuse.weights.matrix = TRUE
  )
  
  # Calculate the query neighbors in the reference
  # with respect to the integrated embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference$map[["refDR"]]),
    query = Embeddings(query[["integrated_dr"]]),
    return.neighbor = TRUE,
    l2.norm = TRUE
  )
  
  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  query <- Azimuth:::NNTransform(
    object = query,
    meta.data = reference$map[[]]
  )
  
  # Project the query to the reference UMAP.
  query[["proj.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]],
    reduction.model = reference$map[["refUMAP"]],
    reduction.key = 'UMAP_'
  )
  
  
  # Calculate mapping score and add to metadata
  query <- AddMetaData(
    object = query,
    metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  
  data <- query@meta.data %>%
    tibble::rownames_to_column('cell') %>%
    select(-c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt','nCount_refAssay','nFeature_refAssay')) %>%
    select(cell,predicted.subclass,predicted.class,predicted.cluster,predicted.cross_species_cluster,predicted.subclass.score,predicted.class.score,predicted.cluster.score,predicted.cross_species_cluster.score,mapping.score)
  
  if (!file.exists("Azimuth local TSVs"))
  {
    dir.create("Azimuth local TSVs")
  }
  
  write_tsv(data,paste0("Azimuth local TSVs/",outputfile,".tsv"))
}
rds_files_dir <- paste0(getwd(),"/RDS_files/")
rds_files <- list.files(rds_files_dir)
for (rds_file in rds_files)
{
  azimuth_analysis(file = paste0(rds_files_dir,rds_file), outputfile = strsplit(rds_file,"\\.")[[1]][1])
}



