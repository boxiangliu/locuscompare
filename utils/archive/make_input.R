library(data.table)
library(stringr)

# PHACTR1
eqtl_fn='data/PHACTR1.eqtl.txt'
gwas_fn='data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'

eqtl=fread(eqtl_fn,select=c(2,4),col.names=c('SNP','pval'))
eqtl[,chr:=str_split_fixed(SNP,'_',5)[,1]]
eqtl[,pos:=str_split_fixed(SNP,'_',5)[,2]]

gwas=fread(paste0('gunzip -c ',gwas_fn),select=c(2:4,10),col.names=c('rsid','chr','pos','pval'))
