#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Boxplot for the figure F1D
### DATE: 13.04.2020
### AUTHOR: Artem Baranovskii

require(tidyverse)

set.seed(42)
p.diss <- file.path(args[1], "Data/PRJNA693557/")
p.out <- file.path(args[1], "Figures/")

res <- read_csv(file.path(p.diss, "ikarus/out/results_final.csv"))

pdf(file = paste0(p.out, "/", "F1_D.pdf"), onefile = T, width = 5, height = 5)
  ggplot(data = res %>% 
          pivot_longer(., c("Normal", "Tumor"), names_to = "gset", values_to = "aucell_score") %>% 
          mutate(gset = paste0(gset, " gene list")), 
        aes(x = tier_0, y = aucell_score)) + 
    geom_boxplot(aes(fill = tier_0)) + 
    facet_grid(~ gset) +
    scale_fill_manual(values = c("#AEC7E8", "#FF0000")) + 
    labs(fill = "Tier 0") +
    ylab("AUcell score") +
    theme_bw(base_size = 14) +
    theme(aspect.ratio = 5/4, 
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.background = element_blank(),
          legend.position = c(.15, .25))
dev.off()
