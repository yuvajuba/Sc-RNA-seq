Implement_VDJ <- function(seurat,
                          VDJ_data){
  
  library(stringr)
  ## Manage VDJ data
  tmp <- VDJ_data %>% 
    group_by(cell_id) %>% 
    summarise(
      contig = n(),
      productif = case_when(
        n() == 1 ~ as.character(productive),
        n() == 2 ~ paste(unique(productive), sep = ""),
        TRUE ~ NA
      ),
      v_call = case_when(
        n() == 1 ~ as.character(v_call),
        n() == 2 ~ paste(unique(v_call), collapse = " | "),
        TRUE ~ NA
      ),
      c_call = case_when(
        n() == 1 ~ as.character(c_call),
        n() == 2 ~ paste(unique(c_call), collapse = " | "),
        TRUE ~ NA
      ),
      d_call = case_when(
        n() == 1 ~ as.character(d_call),
        n() == 2 ~ paste(unique(d_call), collapse = " | "),
        TRUE ~ NA
      ),
      j_call = case_when(
        n() == 1 ~ as.character(j_call),
        n() == 2 ~ paste(unique(j_call), collapse = " | "),
        TRUE ~ NA
      ),
      clone_id = case_when(
        n() == 1 ~ as.character(clone_id),
        n() == 2 ~ paste(unique(clone_id), collapse = " | "),
        TRUE ~ NA
      )
    ) %>% 
    unique()
  
  tmp <- tmp %>% 
    mutate(
      subunit = case_when(
        contig == 1 ~ substr(v_call,3,3),
        contig == 2 ~ paste0(substr(str_split_i(v_call, " | ", 1),3,3),
                             substr(str_split_i(v_call, " | ", 3),3,3)),
        TRUE ~ NA
      )
    )
  
  ## Update seuratobject metadata
  cell_metadata <- data.frame(cell_id = colnames(seurat), stringsAsFactors = FALSE)
  annotated_metadata <- dplyr::left_join(cell_metadata, tmp, by = "cell_id")
  
  seurat <- AddMetaData(seurat, metadata = annotated_metadata[, -1])
  
  
  
  return(seurat)
}
