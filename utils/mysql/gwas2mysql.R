library(RMySQL)
library(pool)
library(data.table)
library(stringr)
library(DBI)
source('config/config.R')
source('utils/shared_functions.R')

read_gwas = function(gwas_fn){
	message('INFO - reading GWAS...')
	if (str_detect(gwas_fn,'gz')){
		command = sprintf('gunzip -c %s',paste0(gwas_dir,gwas_fn))
	} else {
		command = paste0(gwas_dir,gwas_fn)
	}

	message(command)

	gwas=tryCatch({
		fread(command,header=TRUE,fill=TRUE)
		}, 
		error = function(e) {NULL}
		)
	if (is.null(gwas)){
		return(gwas)
	}

	if ('snp_pos' %in% colnames(gwas)) {setnames(gwas,'snp_pos','pos')}
	if ('pvalue' %in% colnames(gwas)) {setnames(gwas,'pvalue','pval')}

	if (!('trait' %in% colnames(gwas))) {
		trait=str_replace(gwas_fn,'.txt','')
		trait=str_replace(gwas_fn,'.gz','')
		trait=str_split_fixed(trait,'_',4)[,2]
		gwas$trait=trait
	}
	return(gwas)
}

create_table = function(table_name){
	message('INFO - creating table...')
	if (dbExistsTable(locuscompare_pool,table_name)) {
		dbRemoveTable(locuscompare_pool,table_name)
	}

	dbExecute(
		conn = locuscompare_pool, 
		statement = sprintf("create table %s
		(%s_id int auto_increment primary key,
		trait varchar(100),
		rsid varchar(100),
		pval double);",table_name,table_name)
		)
}

upload_gwas = function(gwas,table_name){
	message('INFO - uploading GWAS')
	step = 1e6
	total_row = nrow(gwas)
	for (i in seq(1,total_row,step)) {
		start = i
		print(start)
		end = ifelse(i+step-1<=total_row,i+step-1,total_row)
		fwrite(gwas[start:end,list(trait, rsid, pval)],'gwas_tmp.csv')

		dbExecute(
			conn = locuscompare_pool,
			statement = sprintf("load data local infile 'gwas_tmp.csv' 
			into table %s 
			fields terminated by ','
			lines terminated by '\n'
			ignore 1 lines
			(trait, rsid, pval);",table_name)
			)

		unlink('gwas_tmp.csv')
	}
}

index_table = function(table_name){
	message('INFO - indexing table')
	dbExecute(
		conn = locuscompare_pool,
		statement = sprintf('alter table %s
			add index trait (trait),
			add index rsid (rsid)',table_name)
		)
}

upload_large_table = function(gwas_fn){
	out_fn = sprintf('processed_data/mysql/gwas2mysql/%s',str_replace(gwas_fn,'.gz',''))
	message(out_fn)
	table_name = str_replace_all(str_replace(gwas_fn,'.txt.gz',''),'-','_')
	command = sprintf("zcat %s/%s | cut -f1,2,5 | sed '1 s/pvalue/pval/' > %s",gwas_dir,gwas_fn,out_fn)
	system(command)

	create_table(table_name)

	dbExecute(
		conn = locuscompare_pool,
		statement = sprintf("load data local infile '%s' 
		into table %s 
		fields terminated by '\t'
		lines terminated by '\n'
		ignore 1 lines
		(trait, rsid, pval);",
		out_fn,
		table_name)
		)

	file.remove(out_fn)
}

locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

##
## The following batch 1 to 5 are on py-gy1: 
##

# gwas_dir = '/mnt/data/shared/datasets/gwas/batch1/'
# gwas_fn_list = list.files(gwas_dir)
# for (gwas_fn in gwas_fn_list){
# 	print(gwas_fn)
# 	gwas = read_gwas(gwas_fn)
# 	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
# 	create_table(table_name)
# 	upload_gwas(gwas,table_name)
# 	index_table(table_name)
# }

# gwas_dir = '/mnt/data/shared/datasets/gwas/batch2/'
# gwas_fn_list = list.files(gwas_dir)
# for (gwas_fn in gwas_fn_list){
# 	print(gwas_fn)
# 	gwas = read_gwas(gwas_fn)
# 	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
# 	create_table(table_name)
# 	upload_gwas(gwas,table_name)
# 	index_table(table_name)
# }


# gwas_dir = '/mnt/data/shared/datasets/gwas/batch3_japanese/'
# gwas_fn_list = list.files(gwas_dir)
# for (gwas_fn in gwas_fn_list){
# 	print(gwas_fn)
# 	gwas = read_gwas(gwas_fn)
# 	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
# 	create_table(table_name)
# 	upload_gwas(gwas,table_name)
# 	index_table(table_name)
# }

# gwas_dir = '/mnt/data/shared/datasets/gwas/batch4_ukbb/'
# gwas_fn_list = list.files(gwas_dir)
# for (gwas_fn in gwas_fn_list){
# 	print(gwas_fn)
# 	gwas = read_gwas(gwas_fn)
# 	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
# 	create_table(table_name)
# 	upload_gwas(gwas,table_name)
# 	index_table(table_name)
# }


# gwas_dir = '/mnt/data/shared/datasets/gwas/batch5/'
# gwas_fn_list = list.files(gwas_dir,pattern='txt.gz$')
# for (gwas_fn in gwas_fn_list){
# 	print(gwas_fn)
# 	gwas = read_gwas(gwas_fn)
# 	table_name = str_replace_all(str_replace(gwas_fn,'.txt.gz',''),'-','_')
# 	create_table(table_name)
# 	upload_gwas(gwas,table_name)
# 	index_table(table_name)
# }

##
## Mike has created a new set of GWAS to replace the old 
## batch 1 to 5. Note that the japanese ones and UKBB stayed the 
## same. 
## 
gwas_dir = '/users/mgloud/projects/gwas/data/munged/'
gwas_fn_list = list.files(gwas_dir,pattern='txt.gz$')
gwas_fn_list = gwas_fn_list[!gwas_fn_list %in% c('GWAS_Waist-Format-2_Shungin_2015.txt.gz','GWAS_Waist-Format-1_Shungin_2015.txt.gz','GWAS_Blood-Cell-Traits_Astle_2016.txt.gz','GWAS_Circulating-Metabolites_Kettunen_2016.txt.gz')]

for (gwas_fn in gwas_fn_list[c(134:length(gwas_fn_list))]){
	print(gwas_fn)
	gwas = read_gwas(gwas_fn)
	if (is.null(gwas)){
		write.table(gwas_fn,'utils/mysql/failed_gwas.txt',append=TRUE,quote=FALSE,col.names=FALSE,row.names=FALSE)
		next
	}
	table_name = str_replace_all(str_replace(gwas_fn,'.txt.gz',''),'-','_')
	create_table(table_name)
	upload_gwas(gwas,table_name)
	index_table(table_name)
}

# Two large tables:
upload_large_table('GWAS_Blood-Cell-Traits_Astle_2016.txt.gz')
upload_large_table('GWAS_Circulating-Metabolites_Kettunen_2016.txt.gz')