library(data.table)
library(stringr)
source('config/config.R')
detach('package:locuscomparer',unload=TRUE)
library(locuscomparer)


in_fn = '/mnt/data/shared/datasets/gwas_raw/GWAS_Type-2-Diabetes_Zhao_2017.txt'
vcf_fn = '/mnt/data/shared/1KG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'

parse_chr_pos = function(chr_pos){
	split_chr_pos = str_split_fixed(chr_pos,':',2)
	chr = str_replace(split_chr_pos[,1],'chr','')
	pos = as.integer(split_chr_pos[,2])
	return(list(chr,pos))
}

gwas = fread(in_fn)
gwas[,c('chr','pos'):=parse_chr_pos(MarkerName)]
gwas2 = get_rsid(vcf_fn,gwas)
out = gwas2[,list(rsid,chr,pos,pval=PVAL)]
fwrite(out,sprintf('%s/GWAS_Type-2-Diabetes_Zhao_2017.txt',gwas_dir),sep='\t')