#!/usr/bin/env Rscript
library(ggplot2)
custom_theme <- theme(plot.title = element_text(size = rel(1.3), face = "bold"),
                      plot.subtitle = element_text(size = rel(0.9), color = "grey40"),
                      plot.margin = margin(12, 12, 12, 12),
                      panel.grid.major.y = element_line(colour = 'gray'),
                      panel.grid.minor.y = element_line(colour = 'gray'),
                      panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      plot.background = element_rect(fill = NULL, colour = 'white'),
                      panel.background = element_rect(fill = 'white'),
                      
                      axis.line = element_line(colour = 'black', linewidth = 0.5),
                      axis.text = element_text(colour = "black", face = 'bold'),
                      axis.text.x = element_text(size = rel(1)),
                      axis.text.y = element_text(size = rel(1)),
                      axis.title = element_text(size = rel(1.2)),
                      axis.ticks = element_line(colour = 'black', linewidth = 0.8),
                           
                      legend.position = "bottom",
                      legend.title = element_text(face = 'bold'),
                      legend.background = element_blank(),
                      legend.key = element_blank())

plot_theme <- custom_theme
