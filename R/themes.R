library(ggplot2)
plot_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                    axis.title = element_text(face = "bold", hjust = 0.5, size = 12),
                    axis.text = element_text(hjust = 0.5, size = 11),
                    legend.position = "right",
                    legend.title = element_text(face = "bold", vjust = 0.85, size = 12),
                    legend.text = element_text(hjust = 0.5, vjust = 0.5, size = 10),
                    panel.background = element_rect(fill = "transparent"),
                    plot.background = element_rect(fill = "transparent"),
                    panel.grid = element_blank(),
                    panel.border = element_blank(),
                    axis.line = element_line(color = "#000000", size = 1))
