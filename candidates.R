
# Libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help="Path to the base of the plink files .bed, .bim and .fam")
parser$add_argument("-o", "--output", help="Path of the output Rdata file")
args <- parser$parse_args()

# SNP data
load(args$input)

gwas$results %>%
  filter(p_value_REML < 0.000001) %>%
  write.table(., file=args$output, col.names=T, row.names=F, sep="\t", quote=FALSE)
