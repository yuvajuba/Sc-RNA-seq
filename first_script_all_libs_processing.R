# Clean start :
# ----------- :
rm(list = ls())
graphics.off()
cat("\014")

# SETUP ENVIRONNEMENT :
# =================== :
source("00_setup.R")


# Processing : #### 
################# #
Libraries <- grep("Lib", list.files(path = "data/"), value = T)
List_SeuratObj <- list()
pct.mt.threshold <- 5  # Threshold for mitochondrial gene expression filtering
diagnosis <- c("M104", "M143", "M187") # sample labels at diagnosis
relapse <- c("M127", "M148", "M187r") # sample labels at relapse

t0 <- Sys.time()

for(Lib in Libraries){
  
  ##  Setting directories
  Raw_dir <- paste0("data/",Lib,"/")
  Seurat_outdir <- paste0(Lib,"/SeuratObject/")
  Out_dir_QC <- paste0(Lib,"/out_fig/QC/")
  Out_dir_fig <- paste0(Lib,"/out_fig/")
  
  ##  Importing 10x data      --------------------------------------
  Rawdata <- Read10X(paste0(Raw_dir,"filtered_feature_bc_matrix/"))
  VDJdata <- read.delim(paste0(Raw_dir,"airr_rearrangement.tsv"), header = T, sep = "\t") %>% 
    dplyr::select(-c(sequence, sequence_aa, sequence_alignment, germline_alignment))
  
  barcodes <- intersect(colnames(Rawdata$`Gene Expression`),
                        colnames(Rawdata$`Antibody Capture`))
  
  umi <- Rawdata$`Gene Expression`[,barcodes]
  hto <- Rawdata$`Antibody Capture`[,barcodes]
  
  ##  Creating seurat object & normalizing RNA data (log normalizing)   -------------
  SeuratObj <- CreateSeuratObject(counts = umi)
  SeuratObj <- NormalizeData(SeuratObj)
  
  ##  Find and scale variable features        ---------------------------------
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "mean.var.plot", nfeatures = 2000)
  SeuratObj <- ScaleData(SeuratObj, features = VariableFeatures(SeuratObj))
  
  ##  Add HTO data as an independent assay & normalize      --------------------------
  SeuratObj[["HTO"]] <- CreateAssayObject(counts = hto)
  SeuratObj <- NormalizeData(SeuratObj, assay = "HTO", normalization.method = "CLR")
  
  ##  Demultiplexing the HTO data     --------------------------------------------
  SeuratObj <- HTODemux(SeuratObj, 
                        assay = "HTO", positive.quantile = 0.99, nstarts = 100, nsamples = 100) # default params
  
  ##  QC HTODemux --> Import as PDF
  pdf(file = paste0(Out_dir_QC,Lib,"_QC_HTODemux.pdf"),
      width = 8,
      height = 6,
      title = "QC Demultiplexing")
  
  ##  Cell counts for the different HTOs    ---------------------------
  print(SeuratObj@meta.data %>% 
          dplyr::count(hash.ID) %>% 
          ggplot()+
          geom_bar(aes(x= hash.ID,
                       y= n,
                       fill= hash.ID),
                   width = 0.7,
                   stat = "identity",
                   position = "stack")+
          scale_fill_manual(values = MyPalette)+
          labs(x= "",
               y= "",
               fill= "HTO")+
          theme_bw()+
          theme(axis.text = element_text(size = 11, face = "bold"),
                legend.title = element_text(size = 12, face = "bold", colour = "darkred"),
                plot.title = element_text(size = 15, face = "bold", colour = "darkred")))
  
  ##  Compare number of UMI distributions across HTOs    --------------------
  print(VlnPlot(SeuratObj, features = "nCount_RNA", pt.size = 0.1, log = TRUE, group.by = "hash.ID")+
          labs(title = "UMI distributions",
               x= "")+
          scale_fill_manual(values = MyPalette))
  
  dev.off()
  
  ##  Subset only singlets      -------------------------------------------------------
  Idents(SeuratObj) <- SeuratObj$HTO_classification.global
  SeuratObj <- subset(SeuratObj, idents = "Singlet", invert = FALSE)
  
  ##  Mitochondiral gene expression : QC      ---------------------------------
  SeuratObj <- PercentageFeatureSet(SeuratObj,
                                    pattern = "^MT-",
                                    col.name = "percent.mt",
                                    assay = "RNA")
  
  ##  Visualization --> Importing as PDF    ---------------------------
  pdf(file = paste0(Out_dir_QC,Lib,"_preprocessing.pdf"), 
      width = 8, 
      height = 6,
      title = "Quality Control")
  
  print(VlnPlot(object = SeuratObj,
                features=c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size = 0.01, log = T))
  
  print(ggplot(SeuratObj@meta.data, 
               aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
          geom_point(size = 0.8) +
          scale_y_log10(breaks = c(90,100,300,1000,3000)) +
          scale_color_gradientn(colours = c("green","gold","orange","darkred"),
                                values = c(0,0.2,0.5,1),
                                limits = c(0,50)) +
          ggtitle("QC plot", "Number of detected genes in function of number of UMI") +
          labs(x = "Number of UMI by cell", y = "Number of detected genes by cell")+
          theme_bw())
  dev.off()
  
  ##  Filtering cells and visualization as PDF    -----------------------------
  SeuratObj <- subset(SeuratObj, subset = percent.mt < pct.mt.threshold & nFeature_RNA < 3500)
  
  pdf(file = paste0(Out_dir_QC,Lib,"_preprocessing_filtered.pdf"), 
      width = 8, 
      height = 6,
      title = "Quality Control - filtered")
  
  print(VlnPlot(object = SeuratObj,
                features=c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size = 0.01, log = T))
  
  print(ggplot(SeuratObj@meta.data, 
               aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
          geom_point(size = 0.8) +
          scale_y_log10(breaks = c(90,100,300,1000,3000)) +
          scale_color_gradientn(colours = c("green","gold","orange","darkred"),
                                values = c(0,0.2,0.5,1),
                                limits = c(0,50)) +
          ggtitle("QC plot", "Number of detected genes in function of number of UMI") +
          labs(x = "Number of UMI by cell", y = "Number of detected genes by cell")+
          theme_bw())
  dev.off()
  
  ##  Add vdj data          -----------------------------------------------------
  SeuratObj <- Implement_VDJ(SeuratObj, VDJdata)
  
  ##  SCtransforme the data     -------------------------------------------------
  SeuratObj <- SCTransform(SeuratObj, 
                           vars.to.regress = "percent.mt", 
                           verbose = T,
                           ncells = 5000,
                           variable.features.n = 3000)
  
  ##  Normalization         -----------------------------------------------
  SeuratObj <- NormalizeData(SeuratObj, 
                             normalization.method = "LogNormalize", 
                             scale.factor = median(SeuratObj$nCount_SCT))
  
  ##  Find variable features        -------------------------------------------
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
  SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj))
  
  ##   Clustering        ----------------------------------------------------------
  ###  Linear dim reduction      --------------------------------------------------
  SeuratObj <- RunPCA(SeuratObj, 
                      features = VariableFeatures(object = SeuratObj),
                      npcs = 50,
                      ndims.print = 1:2,
                      nfeatures.print = 10)
  
  ###  non-linear dim reduction      -----------------------------------------------
  SeuratObj <- RunUMAP(SeuratObj, dims = 1:30)
  
  ###  Clustering      -------------------------------------------------------------
  SeuratObj <- FindNeighbors(SeuratObj, dims = 1:30, verbose = FALSE)
  SeuratObj <- FindClusters(SeuratObj, verbose = FALSE, resolution = seq(0.1, 1.5, by= 0.2))
  
  ###  Add umap coordinates to metadata      ---------------------------------------
  SeuratObj@meta.data <- merge(SeuratObj@meta.data, 
                               as.data.frame(SeuratObj@reductions$umap@cell.embeddings),
                               by = 0)
  SeuratObj@meta.data <- SeuratObj@meta.data %>% 
    column_to_rownames(var = "Row.names")
  
  ###  Cell Cycle      --------------------------------------------------------
  SeuratObj <- CellCycleScoring(SeuratObj,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes)
  
  ###  Plotting the clustering as pdf    ------------------------------------------
  pdf(file = paste0(Out_dir_fig,Lib,"_clusters_before_filt.pdf"), 
      width = 8, 
      height = 6,
      title = "Clusters - no.filt")
  
  print(DimPlot(SeuratObj,
                group.by = "hash.ID",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8)+
          labs(title = "General conditions",
               colour = "Conditions"))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.1",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.5",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.9",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.1.5",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8))
  dev.off()
  
  
  ## Enrich metadata    ------------------------------------------------------
  
  ## Adding a proper label for the conditions and clonotypes
  SeuratObj@meta.data <- SeuratObj@meta.data %>% 
    mutate(Condition = case_when(hash.ID %in% diagnosis ~ "Diagnosis",
                                 hash.ID %in% relapse ~ "Relapse"),
           Clonotype = case_when(startsWith(clone_id,"clono") ~ sub("onotype","",clone_id),
                                 TRUE ~ clone_id))
  
  ## Adding some personalized clustering labels
  SeuratObj@meta.data <- SeuratObj@meta.data %>% 
    mutate(
      Clust_res0.1 = as.factor(paste0(SCT_snn_res.0.1, substr(SeuratObj$Condition,1,1))),
      Clust_res0.3 = as.factor(paste0(SCT_snn_res.0.3, substr(SeuratObj$Condition,1,1))),
      Clust_res0.5 = as.factor(paste0(SCT_snn_res.0.5, substr(SeuratObj$Condition,1,1))),
      Clust_res0.7 = as.factor(paste0(SCT_snn_res.0.7, substr(SeuratObj$Condition,1,1))),
      Clust_res0.9 = as.factor(paste0(SCT_snn_res.0.9, substr(SeuratObj$Condition,1,1))),
      Clust_res1.1 = as.factor(paste0(SCT_snn_res.1.1, substr(SeuratObj$Condition,1,1))),
      Clust_res1.3 = as.factor(paste0(SCT_snn_res.1.3, substr(SeuratObj$Condition,1,1))),
      Clust_res1.5 = as.factor(paste0(SCT_snn_res.1.5, substr(SeuratObj$Condition,1,1)))
    )
  
  ## Plotting some metadata
  pdf(file = paste0(Out_dir_fig,Lib,"_meta_before_filt.pdf"), 
      width = 8, 
      height = 6,
      title = "Metadata - no.filt")
  
  print(DimPlot(SeuratObj,
                group.by = "Condition",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8)+
          labs(title = "General conditions",
               colour = "Conditions"))
  
  print(DimPlot(SeuratObj,
                group.by = "Phase",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8)+
          labs(title = "Cell cycle",
               colour = "CC Phases"))
  
  print(DimPlot(SeuratObj,
                group.by = "Clonotype",
                pt.size = 0.7,
                alpha = 0.8)+
          labs(title = "Clonotypes")+
          theme(legend.position = "none"))
  
  print(DimPlot(SeuratObj, 
                group.by = "Clust_res0.1",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8))
  
  print(DimPlot(SeuratObj, 
                group.by = "Clust_res0.5",
                cols = MyPalette,
                pt.size = 0.7,
                alpha = 0.8))
  
  dev.off()
  
  
  ## Inspecting cell types   --------------------------------------------------------
  
  pdf(file = paste0(Out_dir_fig,Lib,"_cell_types_markers.pdf"), 
      width = 9, 
      height = 6,
      title = "Cell types")
  
  print(
    try(FeaturePlot(SeuratObj, 
                    features = c("CD3D","CD99","CD7","CD5"), 
                    cols = c("#EEEEEE","#221166"),
                    keep.scale = "all",
                    pt.size = 0.6,
                    alpha = 0.8)+
          labs(caption = "T-ALL blast marker genes")+
          theme(plot.caption = element_text(colour = "darkred", face = "bold", size = 12)))
  )
  
  print(
    try(FeaturePlot(SeuratObj, 
                    features = c("CD19","MS4A1"), 
                    cols = c("#EEEEEE","#221166"),
                    keep.scale = "all",
                    pt.size = 0.6,
                    alpha = 0.8)+
          labs(caption = "B cell markers")+
          theme(plot.caption = element_text(colour = "darkred", face = "bold", size = 12)))
  )
  
  print(
    try(FeaturePlot(SeuratObj, 
                    features = c("CD14", "CD68"), 
                    cols = c("#EEEEEE","#221166"),
                    keep.scale = "all",
                    pt.size = 0.6,
                    alpha = 0.8)+
          labs(caption = "Monocytes & macrophages markers")+
          theme(plot.caption = element_text(colour = "darkred", face = "bold", size = 12)))
  )
  
  print(
    try(FeaturePlot(SeuratObj, 
                    features = c("NKG7","GNLY"), 
                    cols = c("#EEEEEE","#221166"),
                    keep.scale = "all",
                    pt.size = 0.6,
                    alpha = 0.8)+
          labs(caption = "NK markers")+
          theme(plot.caption = element_text(colour = "darkred", face = "bold", size = 12)))
  )
  
  
  dev.off()
  
  ## Inspecting clonotypes    --------------------------------------------------
  
  SeuratObj@meta.data$TCR <- SeuratObj@meta.data$v_call
  for(c in c("Clonotype", "TCR")){
    SeuratObj@meta.data %>% 
      dplyr::count(.data[[c]]) %>% 
      dplyr::mutate(x = case_when(n == 1 ~ "unique",
                                  n >= 2 & n < 10 ~ "few",
                                  TRUE ~ .data[[c]])) -> tmp
    
    unique_cl <- tmp[[c]][which(tmp$x == "unique")]
    low_cl <- tmp[[c]][which(tmp$x == "few")]
    
    SeuratObj@meta.data <- SeuratObj@meta.data %>% 
      dplyr::mutate(!!sym(c) := case_when(is.na(.data[[c]]) ~ "missing",
                                          .data[[c]] %in% unique_cl ~ "unique",
                                          .data[[c]] %in% low_cl ~ "low_count",
                                          TRUE ~ .data[[c]]))
    
    # !!sym(c) converts the string "Clonotype" into a symbol (like Clonotype), 
    # and !! unquotes it so mutate() knows to use it dynamically.
    # However — R’s standard = doesn’t allow unquoting on the left-hand side 
    # of an assignment inside mutate()
    # That’s why tidyverse uses the special := operator (called the 
    # walrus operator 🦭 in this context) for dynamic assignment.
  }
  
  pdf(file = paste0(Out_dir_fig,Lib,"_clones&TCR_1_before_filt.pdf"), 
      width = 8, 
      height = 6,
      title = "Clonotypes & TCR")
  
  print(
    DimPlot(SeuratObj,
            group.by = "Clonotype",
            cols = c(MyPalette, colors(distinct = T)),
            pt.size = 1,
            alpha = .7,
            label = F,
            label.box = T,
            label.color = "#fff",
            repel = T)+
      My_umap_theme()+
      scale_fill_manual(values = MyPalette)
  )
  
  t <- 150
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(Condition, Clonotype) %>% 
      ggplot(aes(x= Condition, y= n, fill= Clonotype))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 4,
                fontface = "bold")+
      scale_fill_manual(values = MyPalette)+
      labs(title = "Clonotypes distribution",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()
  )
  
  print(
    DimPlot(SeuratObj,
            group.by = "TCR",
            cols = c(MyPalette, colors(distinct = T)),
            pt.size = 1,
            alpha = .7,
            label = F,
            label.box = T,
            label.color = "#fff",
            repel = T)+
      My_umap_theme()+
      scale_fill_manual(values = MyPalette)
  )
  
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(Condition, TCR) %>% 
      ggplot(aes(x= Condition, y= n, fill= TCR))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 4,
                fontface = "bold")+
      scale_fill_manual(values = MyPalette)+
      labs(title = "TCR distribution",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()
  )
  
  t <- 100
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(TCR, Clonotype, Condition) %>% 
      ggplot(aes(x= Clonotype, y= n, fill= TCR))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 4,
                fontface = "bold")+
      facet_grid(. ~ Condition)+
      scale_fill_manual(values = MyPalette)+
      labs(title = "TCR vs Clonotypes",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()+
      theme(axis.text.x = element_text(angle = 30),
            strip.background = element_rect(fill = "#051", colour = "#000", size = 0.7),
            strip.text = element_text(size = 12, face = "bold"))
  )
  
  print(
    DimPlot(SeuratObj,
            group.by = "SCT_snn_res.0.1",
            cols = c(MyPalette, colors(distinct = T)),
            pt.size = 1,
            alpha = .7,
            label = T,
            label.box = T,
            label.color = "#fff",
            repel = T)+
      My_umap_theme()+
      scale_fill_manual(values = MyPalette)
  )
  
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(SCT_snn_res.0.1, Clonotype) %>% 
      ggplot(aes(x= SCT_snn_res.0.1, y= n, fill= Clonotype))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 3,
                fontface = "bold")+
      scale_fill_manual(values = MyPalette)+
      labs(title = "Clonotypes distribution",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()
  )
  
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(Clust_res0.1, Clonotype) %>% 
      ggplot(aes(x= Clust_res0.1, y= n, fill= Clonotype))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 3,
                fontface = "bold")+
      scale_fill_manual(values = MyPalette)+
      labs(title = "Clonotypes distribution",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()
  )
  
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(SCT_snn_res.0.1, TCR) %>% 
      ggplot(aes(x= SCT_snn_res.0.1, y= n, fill= TCR))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 3,
                fontface = "bold")+
      scale_fill_manual(values = MyPalette)+
      labs(title = "TCR distribution",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()
  )
  
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(Clust_res0.1, TCR) %>% 
      ggplot(aes(x= Clust_res0.1, y= n, fill= TCR))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 3,
                fontface = "bold")+
      scale_fill_manual(values = MyPalette)+
      labs(title = "Clonotypes distribution",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()
  )
  
  dev.off()
  
  ## wave 2 
  
  pdf(file = paste0(Out_dir_fig,Lib,"_clones_clust_TCR_before_filt.pdf"), 
      width = 10, 
      height = 12,
      title = "Clonotypes & TCR & clust")
  
  t <- 100
  print(
    SeuratObj@meta.data %>% 
      dplyr::count(TCR, Clonotype, Clust_res0.1) %>% 
      ggplot(aes(x= Clonotype, y= n, fill= TCR))+
      geom_bar(position = "stack",
               stat = "identity",
               width = .8)+
      geom_text(aes(label= ifelse(n >= t, n, "")),
                position = position_stack(vjust = 0.5, reverse = F),
                colour = "white",
                size = 3.6,
                fontface = "bold")+
      facet_grid(Clust_res0.1 ~ ., 
                 scale = "free", 
                 space = "free")+
      coord_flip()+
      scale_fill_manual(values = MyPalette)+
      labs(title = "TCR - Clonotypes - Clusters",
           y = "Count",
           caption = paste0("Label count min threshold : ", t))+
      My_bplot_theme()+
      theme(axis.text.x = element_text(angle = 30),
            strip.background = element_rect(fill = "#051", colour = "#000", size = 0.7),
            strip.text = element_text(size = 11, face = "bold"),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank())
  )
  
  dev.off()
  
  ## wave 3
  
  pdf(file = paste0(Out_dir_fig,Lib,"_clones_cond_TCR_before_filt.pdf"), 
      width = 10, 
      height = 15,
      title = "Clonotypes & TCR & cond")
  
  print(
    SeuratObj@meta.data %>% 
      ggplot(aes(x= umap_1,
                 y= umap_2))+
      My_umap_theme()+
      geom_point(aes(colour= TCR),
                 size = .6,
                 alpha = .7)+
      facet_grid(Clonotype ~ Condition, 
                 margins = T, switch = "y")+
      theme(plot.margin = margin(t=0.05,b=0.05,l=0.05,r=0.05, unit = "in"),
            panel.background = element_rect(colour = "#555"))+
      guides(colour = guide_legend(override.aes = list(size = 3)))+
      scale_colour_manual(values = MyPalette)
  )
  
  dev.off()
  
  ## wave 4 -- cell cycle
  
  pdf(file = paste0(Out_dir_fig,Lib,"_cell_cycle_1_before_filt.pdf"), 
      width = 12, 
      height = 16,
      title = "CC 1")
  
  print(
    SeuratObj@meta.data %>% 
      ggplot(aes(x= umap_1,
                 y= umap_2))+
      My_umap_theme()+
      geom_point(aes(colour= TCR),
                 size = .6,
                 alpha = .7)+
      facet_grid(Clonotype ~ Phase, 
                 margins = T, switch = "y")+
      theme(plot.margin = margin(t=0.05,b=0.05,l=0.05,r=0.05, unit = "in"),
            panel.background = element_rect(colour = "#555"))+
      guides(colour = guide_legend(override.aes = list(size = 3)))+
      scale_colour_manual(values = MyPalette)
  )
  
  dev.off()
  
  
  pdf(file = paste0(Out_dir_fig,Lib,"_cell_cycle_2_before_filt.pdf"), 
      width = 12, 
      height = 18,
      title = "CC 2")
  
  
  print(
    SeuratObj@meta.data %>% 
      ggplot(aes(x= umap_1,
                 y= umap_2))+
      My_umap_theme()+
      geom_point(aes(colour= TCR),
                 size = .6,
                 alpha = .7)+
      facet_grid(Clust_res0.1 ~ Phase, 
                 margins = T, switch = "y")+
      theme(plot.margin = margin(t=0.05,b=0.05,l=0.05,r=0.05, unit = "in"),
            panel.background = element_rect(colour = "#555"))+
      guides(colour = guide_legend(override.aes = list(size = 3)))+
      scale_colour_manual(values = MyPalette)
  )
  
  dev.off()
  
  
  ## Additional metadata (Merging)    ------------------------------------------
  SeuratObj@meta.data <- SeuratObj@meta.data %>% 
    mutate(v2.cond_cc = paste0(substr(Condition,1,1),"_",Phase),
           v2.cond_cl = paste0(substr(Condition,1,1),"_",Clonotype),
           v2.cond_tcr = paste0(substr(Condition,1,1),"_",TCR),
           v2.clust_cc = paste0(Clust_res0.1,"_",Phase),
           v2.clust_cl = paste0(Clust_res0.1,"_",Clonotype),
           v2.clust_tcr = paste0(Clust_res0.1,"_",TCR),
           v2.cc_cl = paste0(Phase,"_",Clonotype),
           v2.cc_tcr = paste0(Phase,"_",TCR),
           v2.cl_tcr = paste0(Clonotype,"_",TCR),
           v3.cond_cc_cl = paste0(substr(Condition,1,1),"_",Phase,"_",Clonotype),
           v3.cond_cc_tcr = paste0(substr(Condition,1,1),"_",Phase,"_",TCR),
           v3.cond_cl_tcr = paste0(substr(Condition,1,1),"_",Clonotype,"_",TCR),
           v3.clust_cc_cl = paste0(Clust_res0.1,"_",Phase,"_",Clonotype),
           v3.clust_cc_tcr = paste0(Clust_res0.1,"_",Phase,"_",TCR),
           v3.clust_cl_tcr = paste0(Clust_res0.1,"_",Clonotype,"_",TCR),
           v3.cc_cl_tcr = paste0(Phase,"_",Clonotype,"_",TCR),
           v4.cond = paste0(substr(Condition,1,1),"_",Phase,"_",Clonotype,"_",TCR),
           v4.clust = paste0(Clust_res0.1,"_",Phase,"_",Clonotype,"_",TCR))
  
  
  
  ##  Save seuratobject       ------------------------------------------------
  saveRDS(SeuratObj, file = paste0(Seurat_outdir,"Seurat1.rds"))
  List_SeuratObj[[Lib]] <- SeuratObj
  
  sink(paste0(Seurat_outdir,"SeuratObjects_info"))
  cat(c("===========","Seurat1.rds","==========="), sep = "\n")
  print(SeuratObj)
  cat("=======================================================", sep = "\n")
  cat(c("Overview of the metadata :",
        "=========================="), sep = "\n")
  str(SeuratObj@meta.data)
  sink()
  
  
  t = Sys.time() - t0
  cat(c("===========","Duration  : ","==========="), sep = "\n")
  cat(as.numeric(t, units = "mins"), " minutes")
  
}

# Clean up environment :
# -------------------- :
rm(list = setdiff(ls(), c("List_SeuratObj", "MyPalette", 
                          "My_bplot_theme", "My_umap_theme")))








# # Saving the instance -----
# capture.output(devtools::session_info() ,file = "session_info.txt")






















