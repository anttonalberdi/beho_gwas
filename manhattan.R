# Libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help="Paath to the GWAS Rdata file produced by GridLMM")
parser$add_argument("-o", "--output", help="Path of the output manhatan plot")
args <- parser$parse_args()

# Load data
load(args$input)

# Prepare working table
gwas_result <- gwas$results %>%
  #slice_sample(n=100000) %>%
  separate(X_ID, into = c("chr", "bp"), sep = "_", remove = FALSE) %>% # separate chromosome and position
  mutate(bp=as.numeric(bp)) %>%
  mutate(chr=factor(chr, levels=c(c(1:19),"X"))) %>%
  rename(p=p_value_REML) %>%
  filter(!is.nan(p))

###
# Custom Manhattan
###
# https://r-graph-gallery.com/101_Manhattan_plot.html

# Append chromosome coordinates
gwas_result <- gwas_result %>%
    group_by(chr) %>%
    summarise(chr_len=max(bp)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas_result, ., by=c("chr"="chr")) %>%
    arrange(chr, bp) %>%
    mutate(bp_cum=bp+tot)

# Append significance level
gwas_result_unsig <- gwas_result %>%
    filter(p > 0.000001)

gwas_result_sig <- gwas_result %>%
    filter(p <= 0.000001)

# Generate chromosome number code
axisdf <- gwas_result %>%
    group_by(chr) %>%
    summarize(center=mean(bp_cum))

# Plot gwas
ggplot(gwas_result_unsig, aes(x=bp_cum, y=-log10(p))) +
        # Show all points
        geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
        scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
        # Show significant points
        geom_point(data=gwas_result_sig, aes(x=bp_cum, y=-log10(p)), color="#FFB91F", shape=17, alpha=0.9, size=1.6) +
        # custom X axis:
        scale_x_continuous(label = axisdf$chr, breaks= axisdf$center, expand = c(0.01, 0.01)) +
        scale_y_continuous(expand = c(0.01, 0.04)) +
        # Custom the theme:
        theme_bw() +
        labs(x="Chromosome") +
        theme(
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()
        )

# PRint GWAScat
ggsave(args$output,
    width = 20,
    height = 5,
    dpi = 300)
