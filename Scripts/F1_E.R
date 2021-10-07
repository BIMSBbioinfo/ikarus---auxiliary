#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Boxplot for the figure F1D
### DATE: 13.04.2020
### AUTHOR: Artem Baranovskii

require(tidyverse)

set.seed(42)
p.diss <- file.path(args[1], "Data/ENCODE/")
p.out <- file.path(args[1], "Figures/")

### cancer cell lines selection
c.lines <- c("K562", "HepG2", "A549", "Hela-S3", "SK-N-SH", "MCF-7", "HT1080", "NCI-H460", 
             "SK-MEL-5", "SK-N-DZ", "BE2C", "Jurkat", "Clone E6-1", "Panc1")
            
res <- read_csv(file.path(p.diss, "ikarus/out/results_final.csv"))
res <- res %>% mutate(tier_0 = ifelse(raw %in% c.lines, "cancer cell line", tier_0),
                      tier_0 = factor(tier_0, levels = c("primary cell", "cell line", "cancer cell line")))

pdf(file = paste0(p.out, "/", "F1_E.pdf"), onefile = T, width = 5, height = 5)
ggplot(data = res %>% 
         pivot_longer(., c("Normal", "Tumor"), names_to = "gset", values_to = "aucell_score") %>% 
         mutate(gset = paste0(gset, " gene list")), aes(x = tier_0, y = aucell_score)) + 
  geom_boxplot(aes(fill = tier_0)) + 
  facet_grid(~ gset) +
  scale_fill_manual(values = c("#AEC7E8", "darkorange1", "#FF0000")) +
  labs(fill = "Tier 0") +
  ylab("AUcell score") +
  xlab("Tier 0") +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 5/4, 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        legend.position = c(.20, .21))
dev.off()