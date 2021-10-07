#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Stoufer's Z across TCGA survival data (density plot)
### DATE: 12.09.2021
### AUTHOR: Artem Baranovskii

require(tidyverse)
require(patchwork)
require(broom)

### setup
set.seed(42)
p.dat <- file.path(args[1], "Data/")
p.gs <- file.path(args[1], "Data/signatures/")
p.out <- file.path(args[1], "Figures/")


#### Read Ikarus signatures into a gene list
l.gs <- map(read_lines(list.files(p.gs, pattern = "\\.gmt", full.names = T)),
            ~ str_split_fixed(.x, "\t", Inf)) %>%
  map(., ~ .x[which(str_length(.x) > 0)])
names(l.gs) <- map(l.gs, ~ .x[[1]])
l.gs <- map(l.gs, ~ .x[3:length(.x)])

#### Read in Stouffer's Z per gene
t.zsc <- readxl::read_xlsx(file.path(p.dat, "Gene expression.xlsx")) %>%
  filter(!grepl("\\?", Gene)) %>% 
  mutate(Gene = str_sub(Gene, 2, -1)) %>% 
  pivot_longer(., -contains("Gene"), names_to = "src", values_to = "z") %>% 
  mutate(gset = ifelse(Gene %in% l.gs$Tumor, "Tumor", "Else"))
t.zsc <- rbind(t.zsc %>% dplyr::filter(gset == "Tumor"), 
               t.zsc %>% mutate(gset = "Population"))

ordr <- t.zsc %>% filter(src != "Stouffer's Z") %>% 
  group_by(src, gset) %>% 
  summarize(z.mn = mean(z, na.rm = T), .groups = "drop") %>% 
  pivot_wider(., names_from = "gset", values_from = "z.mn") %>% 
  mutate(fc = Tumor - Population) %>% 
  arrange(fc) %>% 
  pull(src)

#### Render in pdf
pdf(file = file.path(p.out, "F5_C.pdf"), onefile = T, width = 6, height = 5)
  print(
    ggplot(data = t.zsc %>% filter(src == "Stouffer's Z"), aes(x = z, color = gset, fill = gset)) + 
    geom_density(alpha = 0.3) + 
    scale_fill_manual(values = c("#aaaaaa", "#FF0000")) + 
    scale_color_manual(values = c("#aaaaaa", "#FF0000")) + 
    geom_vline(xintercept = 0) +
    labs(y = "Kernel density", x = "Stouffer's Z", fill = "Gene set", color = "Gene set") +
    theme_bw(base_size = 14) +
    theme(aspect.ratio = 4/5, axis.text.x = element_text(angle = 0, vjust = 0.3, hjust = 0.5), legend.position = c(.8, .75))
  )
dev.off()