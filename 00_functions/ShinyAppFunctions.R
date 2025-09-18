## 01 - Create a umap with centroids as labels for a group of cells

labeled_umap <- function(seuratobject,
                         metadata,
                         cluster,
                         palette){
  
  seurat = seuratobject
  grp = cluster
  
  if(all(c("umap_1","umap_2") %in% colnames(seurat@meta.data))){
    df = seurat@meta.data %>% 
      mutate(cond = seurat@meta.data[[metadata]] %in% grp,
             group = seurat@meta.data[[metadata]]) %>% 
      select(cond, group, umap_1, umap_2)
  } else {
    stop("UMAP coordinates are not available in your seura@metadata")
  }
  
  ## Initiate the centroids data and coordinates
  centroids <- df %>% 
    filter(group %in% grp) %>% 
    group_by(group) %>% 
    summarise(umap_1 = mean(umap_1),
              umap_2 = mean(umap_2))
  
  ## Setting the colors
  l = length(unique(seurat@meta.data[[metadata]]))
  cols <- palette[1:l] %>% 
    setNames(unique(seurat@meta.data[[metadata]]))
  
  lab_cols <- cols[names(cols) %in% grp]
  
  ## Creating the plot
  pl <- df %>% 
    ggplot(aes(x= umap_1,
               y= umap_2))+
    geom_point(aes(colour= group,
                   alpha= cond),
               size= 1.1)+
    geom_label_repel(data = centroids,
                     aes(label= group, 
                         fill = group),
                     size = 5,
                     fontface = "bold",
                     color = "black",
                     box.padding = 0.8,
                     label.padding = 0.4,
                     segment.color = NA)+
    labs(title = paste0(grp[1], "   vs   ", grp[2]))+
    scale_colour_manual(values = cols)+
    scale_fill_manual(values = lab_cols)+
    scale_alpha_manual(values = c("TRUE"=1, "FALSE"=0.02))+
    My_umap_shiny()+
    guides(colour = "none",
           alpha = "none",
           fill = "none")
  
  return(pl)
}




