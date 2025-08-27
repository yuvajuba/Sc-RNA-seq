## Umap plot theme :  
## --------------- :
My_umap_theme <- function(){
  theme_void()+
    theme(# plot.
          plot.title = element_text(size = 17, face = "bold.italic", colour = "#024"),
          plot.subtitle = element_text(size = 14, face = "italic", colour = "#024", 
                                       margin = margin(b=0.3, unit = "in")),
          plot.caption = element_text(size = 11, face = "italic", colour = "#024",
                                      margin = margin(t=0.2, unit = "in")),
          plot.tag = element_text(size = 13, face = "bold", colour = "#000"),
          plot.tag.position = "topright",
          
          # legend.
          legend.margin = margin(l=0.15, unit = "in"),
          legend.title = element_text(size = 13, face = "bold", colour = "#024", 
                                      margin = margin(b=0.15, unit = "in")),
          legend.text = element_text(size = 11),
          legend.key = element_blank(),
          
          # facets.
          strip.text = element_text(size = 12, colour = "#fff", face = "bold"),
          strip.background = element_rect(fill = "#024"),
          strip.background.y = element_rect(fill = "#042"),
          strip.text.y = element_text(angle = 90))
}


My_umap_shiny <- function(){
  theme_void()+
    theme(# plot.
      plot.title = element_text(size = 20, face = "bold.italic", colour = "#024"),
      plot.subtitle = element_text(size = 16, face = "italic", colour = "#024", 
                                   margin = margin(b=0.3, unit = "in")),
      plot.caption = element_text(size = 13, face = "italic", colour = "#024",
                                  margin = margin(t=0.2, unit = "in")),
      plot.tag = element_text(size = 14, face = "bold", colour = "#000"),
      plot.tag.position = "topright",
      
      # legend.
      legend.margin = margin(l=0.15, unit = "in"),
      legend.title = element_text(size = 16, face = "bold", colour = "#024", 
                                  margin = margin(b=0.15, unit = "in")),
      legend.text = element_text(size = 13, face = "bold"),
      legend.key = element_blank(),
      
      # facets.
      strip.text = element_text(size = 15, colour = "#fff", face = "bold"),
      strip.background = element_rect(fill = "#024"),
      strip.background.y = element_rect(fill = "#042"),
      strip.text.y = element_text(angle = 90))
}




## Bar plot theme :  
## -------------- :
My_bplot_theme <- function(){
  theme_linedraw()+
    theme(
      plot.title = element_text(size = 16, face = "bold.italic", colour = "#024"),
      plot.subtitle = element_text(size = 13, face = "italic", colour = "#024", 
                                   margin = margin(b=0.1, unit = "in")),
      plot.caption = element_text(size = 11, face = "italic", colour = "#024",
                                  margin = margin(t=0.2, unit = "in")),
      plot.tag = element_text(size = 13, face = "bold", colour = "#000"),
      plot.tag.position = "topright",
      legend.margin = margin(l=0.15, unit = "in"),
      legend.title = element_text(size = 12, face = "bold", colour = "#024", 
                                  margin = margin(b=0.15, unit = "in")),
      legend.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, colour = "#024", face = "bold"),
      axis.title.x = element_text(margin = margin(t=0.15, unit = "in")),
      axis.text = element_text(size = 10),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(linetype = 2, colour = "#555"),
      panel.grid.minor.y = element_line(linetype = 2, colour = "#aaa")
    )
}
