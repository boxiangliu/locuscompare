library(data.table)
library(stringr)
library(foreach)

setwd('/var/lib/mysql-files')

fn = list.files()
y=foreach(f=fn,.combine='rbind')%do%{
	x = fread(f,header=FALSE,col.names=c('gene_name','rsid','chr','pos','pval'))
	trait = str_replace_all(f,'eQTL_|_GTEx_2017.tsv','')
	x$trait = trait
	return(x)
}

fwrite(y[,list(trait,rsid,chr,pos,pval)],sprintf('PHACTR1.tsv'))