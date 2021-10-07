#!/opt/R/4.0/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)

### INFO: GO analysis of co-expressed genes identified by SEEK
### DATE: 13.04.2020
### AUTHOR: Artem Baranovskii

require(tidyverse)
require(gprofiler2)

### setup
set.seed(42)
p.see <- file.path(args[1], "Data/SEEK/")
p.out <- file.path(args[1], "Figures/")


#### Read in co-expression analysis results from SEEK
t.coexp <- read_tsv(file.path(p.see, "/coexpr.txt"),
                    skip = 77)
v.coexp <- t.coexp %>% slice(1:150) %>% pull(Gene)


gost_res <- gost(query = v.coexp,
                 organism = "hsapiens", ordered_query = FALSE,
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                 measure_underrepresentation = FALSE, evcodes = TRUE,
                 user_threshold = 2e-39, correction_method = "g_SCS",
                 domain_scope = "annotated", custom_bg = NULL,
                 numeric_ns = "", sources = c("GO"))

publish_gosttable(gost_res, 
                  highlight_terms = gost_res$result$term_id,
                  use_colors = TRUE, ggplot = F,
                  show_columns = c("source", "term_name", "term_size", "intersection_size", "p_value"),
                  filename = file.path(p.out, "F4_C.pdf"))