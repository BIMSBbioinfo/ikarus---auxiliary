#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Plot Ikarus's intersections with CN amplified genes from COSMIC
### DATE: 31.08.2021
### AUTHOR: Artem Baranovskii

require(tidyverse)

### setup
set.seed(42)
p.dat <- file.path(args[1], "Data/")
p.gs <- file.path(args[1], "Data/signatures/")
p.out <- file.path(args[1], "Figures/")


#### Read in the intersection between amplified genes and randomly generated genes
t.int <- purrr::map(list.files(file.path(p.dat, "CosmicCNA"),
                               full.names = T) %>% 
                      grep("_CNA", ., invert = T, value = T), 
                    ~ data.table::fread(.x)
                    ) %>% 
  purrr::reduce(., rbind) %>% 
  filter(tissue != "NS") %>% 
  mutate(gset = str_split_fixed(gset, "_sample", 2)[, 1],
         tissue = ifelse(str_detect(tissue, "central_"), "CNS",
                         ifelse(str_detect(tissue, "lymphoid"), "blood", 
                                ifelse(str_detect(tissue, "upper_aero"), "digestive_tract", tissue))),
         gset = ifelse(str_detect(gset, "Tumor"), "Tumor", gset))  
#### order by gene list int size
t.int[t.int$gset != "Random", ] %>% arrange(desc(int.size)) %>% pull(tissue) -> v.ord
t.int <- t.int %>% mutate(tissue = factor(tissue, levels = v.ord))


pdf(file = file.path(p.out, "F5_F.pdf"), onefile = T, width = 5, height = 5)
  print(
    ggplot(data = t.int[t.int$gset == "Random", ], aes(y = tissue, x = int.size / 162)) +
      #geom_vline(xintercept = 0) +
      geom_boxplot(outlier.size = 0.3) +
      geom_point(data = t.int[t.int$gset != "Random", ], aes(color = gset)) + 
      scale_color_manual(values = c("red")) +
      labs(x = "Intersection / gene set size", y = "Tumor tissue of origin", color = "Gene set") +
      theme_bw(base_size = 14) +
      theme(aspect.ratio = 5/3, axis.text.x = element_text(angle = 0, hjust = 1), 
            legend.position = "none",#c(0.65, 0.9), 
            legend.background = element_blank())
  )
dev.off()