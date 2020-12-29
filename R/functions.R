#' Prepare the seurat object consisting of different datasets
#'
#' @param seurat_list a list of seurat objects
#' @param sample_names a vector of names for seurat objects or samples
#'
#' @return Returns a Seurat object that can be used for following analysis
#'
#' @import Seurat limma matrixStats stringr
#'
#' @export

Prepare_Seurat <- function(seurat_list, sample_names){
  common_genes <- rownames(seurat_list[[1]])
  for (i in 1:length(seurat_list)) {
    seurat_object <- seurat_list[[i]]
    common_genes <- intersect(rownames(seurat_object),common_genes)
    seurat_list[[i]]@meta.data$orig.ident = sample_names[i]
  }
  print(common_genes)
  data <- data.frame(matrix(ncol = 0, nrow = length(common_genes)))
  rownames(data) <- common_genes
  meta_data <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(meta_data) <- c("orig.ident","seurat_clusters","nCount_RNA","nFeature_RNA","percent.mt")

  for (i in 1:length(seurat_list)) {
    seurat_object <- seurat_list[[i]]
    raw_data = seurat_object@assays$RNA@counts[as.character(common_genes),]
    raw_meta_data = seurat_object@meta.data[,c("orig.ident","seurat_clusters","nCount_RNA","nFeature_RNA","percent.mt")]
    colnames(raw_data) <- str_c(sample_names[i],colnames(raw_data),sep = "_")
    rownames(raw_meta_data) <- str_c(sample_names[i],rownames(raw_meta_data),sep = "_")
    data = cbind(data,raw_data)
    meta_data = rbind(meta_data,raw_meta_data)
  }
  New_seurat <- CreateSeuratObject(counts = as.matrix(data), meta.data = meta_data)
  New_seurat <- NormalizeData(New_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  return(New_seurat)
}


#' Prepare metadata of the seurat object
#'
#' @param object seurat object from Prepare_seurat
#'
#' @return
#' @export

Prepare_Metadata <- function(object){
  cell.info <- object@meta.data[,c("orig.ident","seurat_clusters")]
  cell.info$seurat_clusters <- stringr::str_c(cell.info$orig.ident,cell.info$seurat_clusters,sep = "_")
  return(cell.info)
}

#' integration of different datasets
#'
#' @param object seurat object from Prepare_Seurat
#' @param metadata metadata from Prepare_Metadata
#' @param groups how many groups you expect to be in your final clusters
#'
#' @import stats
#'
#' @return
#' @export

clincluster <- function(object, metadata, groups){
  matrix <- expm1(object@assays$RNA@data)
  keep <- rowSums(as.matrix(matrix > 1)) > 5
  dge <- edgeR::DGEList(counts = matrix[keep,]) # make a edgeR object
  rm(matrix,keep)
  design <- model.matrix(~  0 + seurat_clusters, data = metadata)  # Use 0 because we do not need intercept for this linear model
  v <- voom(dge, design, plot = F)
  fit <- lmFit(v, design)
  initial.clusters <- colnames(design)
  nc <- ncol(design)
  ## Automating makeContrasts call in limma
  contrast_all <- gtools::permutations(v = initial.clusters, n = nc, r = 2)
  contrast_all <- apply(contrast_all, MARGIN = 1, function(x) return(paste(x[1],"-",x[2], sep = "")))
  cont.matrix <- makeContrasts(contrasts = contrast_all,
                               levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  n_deg <- matrix(0, ncol = nc, nrow = nc)  # number of DE genes
  colnames(n_deg) <- rownames(n_deg) <- gsub(x = colnames(design)[1:nc], pattern = "seurat_clusters",replacement = "")
  logcount <- object@assays$RNA@data[rownames(object@assays$RNA@data) %in% rownames(dge),]
  for(i in 1:(nc-1)) {
    for(j in (i+1):nc) {
      if(i == j) {
        n_deg[i,j] <- 0
      } else if (j < i) {
        coef_k = (i-1)*(nc-1)+j
      } else if (j > i) {
        coef_k = (i-1)*(nc-1)+j-1
      }

      if(i != j) {
        rls <- topTable(fit2, n = Inf, coef = coef_k, sort = "p", lfc = 0.6, p = 0.05 )
        if(nrow(rls) > 1) {
          v_expr <- logcount[match(rownames(rls),rownames(logcount)), metadata$seurat_clusters == rownames(n_deg)[i]]
          rls$ratio1 <- rowSums(as.matrix(v_expr > 0.5))/ncol(v_expr)
          v_expr <- logcount[match(rownames(rls),rownames(logcount)), metadata$seurat_clusters == colnames(n_deg)[j]]
          rls$ratio2 <- rowSums(as.matrix(v_expr > 0.5))/ncol(v_expr)
          rls$ratiomax <- rowMaxs(as.matrix(rls[,c("ratio1", "ratio2")]))
          rls$ratiomin <- rowMins(as.matrix(rls[,c("ratio1", "ratio2")]))
          rls <- rls[rls$ratiomax > 0.25, ]
          n_deg[i,j] <- sum(apply(rls, MARGIN = 1, function(x) return(abs(x[1]) * (x[9]+0.01)/(x[10]+0.01)))) ## 0.01 is used here to enhance the differences of on-off genes
        } else if (nrow(rls) == 1) {
          n_deg[i,j] <- sum(rls$logFC)
        }
        ## This eqaution take fold change and expression ratio into account
      }
    }
  }
  n_deg <- n_deg + t(n_deg)
  ## final cluster
  hc <- hclust(as.dist(n_deg))
  plot(hc); rect.hclust(hc, k = groups, border = "red")
  hc.cluster <- cutree(hc, k = groups)
  object@meta.data$clincluster <- hc.cluster[match(object@meta.data$seurat_clusters, names(hc.cluster))]
  object@meta.data$clincluster <- as.factor(object@meta.data$clincluster)
  return(object)
}
