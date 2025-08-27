Implement_VDJ <- function(seurat,
                          VDJ_data){
  
  ## initiate the data
  meta <- seurat@meta.data
  meta$clone_id <- NA
  meta$v_call <- NA
  
  ## implement the values
  for (i in 1:nrow(VDJ_data)) {
    cell <- VDJ_data$cell_id[i]
    if (cell %in% rownames(meta)) {
      meta[cell, "clone_id"] <- VDJ_data$clone_id[i]
      meta[cell, "v_call"] <- VDJ_data$v_call[i]
    }
  }
  
  ## update the object
  seurat@meta.data <- meta
  
  return(seurat)
}