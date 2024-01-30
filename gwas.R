
# Libraries
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GridLMM))
suppressPackageStartupMessages(library(snpStats))
suppressPackageStartupMessages(library(qqman))

# Parse arguments
parser <- ArgumentParser()
parser$add_argument("-i", "--input", help="Path to the base of the plink files .bed, .bim and .fam")
parser$add_argument("-m", "--metadata", help="Path to the metadata tsv file")
parser$add_argument("-o", "--output", help="Path of the output Rdata file")
args <- parser$parse_args()

# SNP data
X = read.plink(args$input)
X = as(X$genotypes,'numeric')

# Metadata
data = read_table(args$metadata) %>%
      mutate(cage=factor(cage)) %>%
      mutate(mouse=factor(mouse)) %>%
      mutate(time=factor(time)) %>%
      as.data.frame()

# Run GWAS
gwas = GridLMM_GWAS(formula = dominance ~ 1 + (1|cage),
                                  test_formula = ~1,
                                  reduced_formula = ~0,
                                  data = data,
                                  X = X,
                                  X_ID = 'mouse',
                                  method = 'REML',
                                  fillNAX = TRUE,
                                  verbose = F,
                                  mc.cores = 8)

# Save GWAS object
save(gwas,file=args$output)
