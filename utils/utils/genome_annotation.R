library(data.table)
library(stringr)

gencode_fn = '/srv/persistent/bliu2/rpe/data/reference/gencode.v19.annotation.gtf'
mean_expression_fn = '/srv/persistent/bliu2/rpe/data/rnaseq/mean_expression.txt'

read_gencode = function(fn=gencode_fn){
	x = fread(fn,select=c(1,3,4,5,7,9))
	setnames(x,c('chr','type','start','end','strand','annotation'))
	x = x[type=='gene',]
	x[,gene_id:=str_extract(annotation,'(?<=gene_id \\")(.+?)(?=\\";)')]
	x[,gene_name:=str_extract(annotation,'(?<=gene_name \\")(.+?)(?=\\";)')]
	x[,gene_type:=str_extract(annotation,'(?<=gene_type \\")(.+?)(?=\\";)')]
	x$annotation=NULL
	return(x)
}

read_mean_expression = function(){
	x = fread(mean_expression_fn)
	return(x)
}
