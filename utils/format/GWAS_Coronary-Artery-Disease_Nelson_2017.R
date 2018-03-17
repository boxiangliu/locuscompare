library(data.table)
source('config/config.R')

in_fn = '/mnt/data/shared/datasets/gwas_raw/GWAS_Coronary-Artery-Disease_Nelson_2017.txt'
gwas = fread(in_fn)
out = gwas[,list(rsid=snptestid,chr,pos=bp_hg19,pval=`p-value_gc`)]
fwrite(out,sprintf('%s/GWAS_Coronary-Artery-Disease_Nelson_2017.txt',gwas_dir),sep='\t')