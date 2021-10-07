#!/opt/R/4.0/bin/Rscript
args = commandArgs(trailingOnly = TRUE)

### INFO: Upset plot; Intersection between CancerSEA signatures and Ikarus; Figure F4B
### DATE: 13.04.2020
### AUTHOR: Artem Baranovskii

require(tidyverse)
require(UpSetR)

### setup
set.seed(42)
p.cse <- file.path(args[1], "Data/CancerSEA/")
p.gs <- file.path(args[1], "Data/signatures/")
p.out <- file.path(args[1], "Figures/")


#### Read Ikarus signatures into a gene list
l.gs <- map(read_lines(list.files(p.gs, pattern = "\\.gmt", full.names = T)),
            ~ str_split_fixed(.x, "\t", Inf)) %>% 
  map(., ~ .x[which(str_length(.x) > 0)])
names(l.gs) <- map(l.gs, ~ .x[[1]])
l.gs <- map(l.gs, ~ .x[3:length(.x)])


### CacnerSEA
#### Quick download script if necessary
# l.pfx <- "http://biocc.hrbmu.edu.cn/CancerSEA/download/signature/"
n.sgs <- read_tsv(file.path(p.cse, "siglist.tsv"),
                  col_names = F, 
                  col_types = cols(.default = col_character())) %>% 
  pull(X1) %>% 
  str_replace_all(., " ", "_")
# dummy <- purrr::map(n.sgs,
#            ~ download.file(url = paste0(l.pfx, .x, ".txt"), 
#                            destfile = file.path("/local/abarano/Projects/ikarus/Data/auxiliary/CancerSEA/", paste0(.x, ".txt")), 
#                            method = "libcurl"))
#### Read in CancerSEA signatures
l.csea <- purrr::map(n.sgs, 
                     ~ read_tsv(file.path(p.cse, paste0(.x, ".txt")), 
                                col_names = T, col_types = cols(.default = col_character())) %>% 
                       pull(GeneName)) %>% 'names<-'(n.sgs)


#### Load in downloaded signatures
tmp.csea <- tibble(gene_name = reduce(c(list(l.gs$Tumor), l.csea), union))
t.csea <- purrr::map2(c(list(l.gs$Tumor), l.csea), 
                      c("Tumor_gene_list", names(l.csea)), 
                      ~ tmp.csea %>% mutate(!!sym(.y) := ifelse(gene_name %in% .x, 1, 0))) %>% 
  purrr::reduce(., inner_join, by = "gene_name")

#### Render in pdf
pdf(file = file.path(p.out, "F4_B.pdf"), onefile = T, width = 7.5, height = 5)
upset(as.data.frame(t.csea), 
            sets = colnames(t.csea)[-1], 
            mb.ratio = c(0.55, 0.45),
            sets.bar.color = "grey30", 
            mainbar.y.label = "Gene list intersection size", 
            sets.x.label = "Gene list size",
            order.by = "freq")
dev.off()