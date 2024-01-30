# beho_gwas
GWAS analysis of BeHo genomic data

## Fetch VCF file from ERDA
```sh
wget https://sid.erda.dk/share_redirect/F0LbFHgFdk/variants_filtered.vcf.gz
```

## Subset VCF
Subset variant file for even representation of cages
```sh
#Create the batch file
cat <<EOF > 0_subset.sh
#!/bin/bash
#SBATCH --job-name=0_subset
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH --time=0:15:00

module load perl gsl/2.5 bcftools/1.16
bcftools view -s ^treatment_male_01_1,treatment_male_07_3,treatment_male_08_1 variants_filtered.vcf.gz > variants_filtered2.vcf
bcftools reheader -s rename.txt variants_filtered2.vcf > variants_filtered3.vcf
EOF

#Launch the batch file
sbatch 0_subset.sh
```


## Create plink bed files from VCF
```sh
#Create the batch file
cat <<EOF > 1_input.sh
#!/bin/bash
#SBATCH --job-name=1_input
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128gb
#SBATCH --time=1:00:00
module load openblas/0.3.24 plink/2.0.0
plink2 --vcf variants_filtered3.vcf --pheno beho_phenotype.tsv --update-sex beho_sex.tsv --max-alleles 2 --set-all-var-ids @_# --make-bed --out variants_filtered
# only desired chromosomes
#plink2 --vcf variants_filtered.vcf.gz --pheno beho_phenotype.tsv --update-sex beho_sex.tsv --max-alleles 2 --chr 1-2 --set-all-var-ids @_#  --make-bed --out variants_filtered_chr1-2
EOF

#Launch the batch file
sbatch 1_input.sh
```

## Run GWAS with GridLMM

### Against dominance

```sh
#Create the batch file
cat <<EOF > 2_gwas.sh
#!/bin/bash
#SBATCH --job-name=2_gwas
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=500gb
#SBATCH --time=24:00:00

# Load R (with all dependencies already installed)
module load gcc R/4.2.1

# Run R script
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_mean.tsv" -o "results/gwas_dominance_mean.Rdata"
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_cd.tsv" -o "results/gwas_dominance_cd.Rdata"
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_cr.tsv" -o "results/gwas_dominance_cr.Rdata"
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_dt.tsv" -o "results/gwas_dominance_dt.Rdata"
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_hr.tsv" -o "results/gwas_dominance_hr.Rdata"
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_ht.tsv" -o "results/gwas_dominance_ht.Rdata"
/maps/direct/software/R/4.2.1/bin/Rscript gwas.R -i "variants_filtered" -m "data/dominance_op.tsv" -o "results/gwas_dominance_op.Rdata"
EOF

#Launch the batch file
sbatch 2_gwas.sh
```

## Plot Manhattan chart

```sh
#Create the batch file
cat <<EOF > 3_manhattan.sh
#!/bin/bash
#SBATCH --job-name=3_manhattan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=1:00:00

# Load R (with all dependencies already installed)
module load gcc R/4.2.1

# Run R script
/maps/direct/software/R/4.2.1/bin/Rscript manhattan.R -i "gwas.Rdata" -o "manhattan.png"
EOF

#Launch the batch file
sbatch 3_manhattan.sh
```

## Plot Manhattan chart (weight)

```sh
#Create the batch file
cat <<EOF > 3_manhattan2.sh
#!/bin/bash
#SBATCH --job-name=3_manhattan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=1:00:00

# Load R (with all dependencies already installed)
module load gcc R/4.2.1

# Run R script
/maps/direct/software/R/4.2.1/bin/Rscript manhattan.R -i "gwas2.Rdata" -o "manhattan_weight.png"
EOF

#Launch the batch file
sbatch 3_manhattan2.sh
```

## Extract significant

```sh
#Create the batch file
cat <<EOF > 3_manhattan2.sh
#!/bin/bash
#SBATCH --job-name=3_manhattan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48gb
#SBATCH --time=1:00:00

# Load R (with all dependencies already installed)
module load gcc R/4.2.1

# Run R script
Rscript candidates.R -i "gwas2.Rdata" -o "gwas2_sig.tsv"
EOF

#Launch the batch file
sbatch 4_candidates.sh
```

## Significant variants

1_41965255  rs235926788 Gm29260
1_41987245  rs225058266 NA
1_42031168  rs33241262  NA
1_42031187  rs33242035  NA
1_42044011  rs33244413  NA
1_42076205  rs222796076 NA
1_42260929  rs264540779 NA
1_42691925  rs32194658  Pantr1
1_43054237  rs248553431 Gpr45
1_43293781  rs579943378 Gm29040
2_170316461 rs27603046  NA
2_170878051 rs241075516 NA
2_171104417 rs246720942 ENSMUSG00000084829
2_171109438 rs27620669  ENSMUSG00000084829
7_29506084  rs50258849  NA
9_54804919  rs248454421 Ireb2
9_69129503  NA  RORA
9_69129518  rs226488564 RORA
9_69129524  rs242479775 RORA
9_71703432  rs30426714  Cgnl1
15_3628694  rs49026775  NA
17_76537375 rs265589933 NA
17_76537379 NA  NA

RORA
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3797180/

## Variant exploration

```sh
plink2 --bfile variants_filtered --snp 9_69129524 --out 9_69129524 --export ped
```
