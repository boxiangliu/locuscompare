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
		fread(command)
		}, 
		error = function(e) {NULL}
		)
	if (is.null(gwas)){
		next
	}

	if ('snp_pos' %in% colnames(gwas)) {setnames(gwas,'snp_pos','pos')}
	if ('pvalue' %in% colnames(gwas)) {setnames(gwas,'pvalue','pval')}

	if (!('trait' %in% colnames(gwas))) {
		trait=str_replace(gwas_fn,'.txt','')
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


locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

gwas_dir = '/mnt/data/shared/datasets/gwas/batch1/'
gwas_fn_list = list.files(gwas_dir)
for (gwas_fn in gwas_fn_list){
	print(gwas_fn)
	gwas = read_gwas(gwas_fn)
	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
	create_table(table_name)
	upload_gwas(gwas,table_name)
	index_table(table_name)
}

gwas_dir = '/mnt/data/shared/datasets/gwas/batch2/'
gwas_fn_list = list.files(gwas_dir)
for (gwas_fn in gwas_fn_list){
	print(gwas_fn)
	gwas = read_gwas(gwas_fn)
	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
	create_table(table_name)
	upload_gwas(gwas,table_name)
	index_table(table_name)
}


gwas_dir = '/mnt/data/shared/datasets/gwas/batch3_japanese/'
gwas_fn_list = list.files(gwas_dir)
for (gwas_fn in gwas_fn_list){
	print(gwas_fn)
	gwas = read_gwas(gwas_fn)
	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
	create_table(table_name)
	upload_gwas(gwas,table_name)
	index_table(table_name)
}

gwas_dir = '/mnt/data/shared/datasets/gwas/batch4_ukbb/'
gwas_fn_list = list.files(gwas_dir)
for (gwas_fn in gwas_fn_list){
	print(gwas_fn)
	gwas = read_gwas(gwas_fn)
	table_name = str_replace_all(str_replace(gwas_fn,'.txt',''),'-','_')
	create_table(table_name)
	upload_gwas(gwas,table_name)
	index_table(table_name)
}

gwas_dir = '/mnt/data/shared/datasets/gwas/batch5/'
gwas_fn_list = list.files(gwas_dir,pattern='txt.gz$')
for (gwas_fn in gwas_fn_list){
	print(gwas_fn)
	gwas = read_gwas(gwas_fn)
	table_name = str_replace_all(str_replace(gwas_fn,'.txt.gz',''),'-','_')
	create_table(table_name)
	upload_gwas(gwas,table_name)
	index_table(table_name)
}
