library(data.table)
library(stringr)
library(RMySQL)
library(pool)
library(DBI)
source('config/config.R')

in_fn = 'processed_data/mysql/coloc2mysql/all_finemap_results_sorted.txt'


read_colocalization = function(fn){
    x = fread(fn)
    return(x)
}

munge_colocalization = function(x){
	split_snp = str_split_fixed(x$ref_snp,'_',2)
	x$chr = split_snp[,1]
	x$pos = as.integer(split_snp[,2])
	x$ref_snp = NULL

	setnames(x, 
		c('base_gwas_file','gwas_trait','eqtl_file','feature','-log_gwas_pval','-log_eqtl_pval'),
		c('gwas','trait','eqtl','gene_id','logp_gwas','logp_eqtl'))

	base_trait = basename(x$trait)
	split_trait = str_split_fixed(base_trait,'_',4)
	x$trait = ifelse(split_trait[,2]=='',split_trait[,1],split_trait[,2])

	x$gwas = str_replace_all(x$gwas,'-','_')
	x$gwas = str_replace(x$gwas,'_txt_gz','')

	x$eqtl = str_replace(x$eqtl,'_allpairs_txt_gz','')
	x$eqtl = paste0('eQTL_',x$eqtl,'_GTEx_v7')
	
	x = x[!is.na(clpp),]

	return(x[,list(gwas,trait,eqtl,gene_id,logp_gwas,logp_eqtl,clpp)])
}

get_gencode = function(){
	gencode = dbGetQuery(
        conn = locuscompare_pool,
        statement = "select gene_id, gene_name, chr, start 
			from gencode_v19_gtex_v6p;"
    )
    return(gencode)
}

merge_with_gencode = function(colocalization){
	colocalization = merge(colocalization,gencode,by='gene_id')
	setnames(colocalization,'start','pos')
	setcolorder(colocalization,c('gwas','trait','eqtl','gene_id','gene_name','chr','pos','logp_gwas','logp_eqtl','clpp'))
	return(colocalization)
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
		gwas varchar(100),
		trait varchar(64),
		eqtl varchar(64),
		gene_id varchar(20),
		gene_name varchar(20),
		chr varchar(5),
		pos integer,
		logp_gwas double,
		logp_eqtl double,
		clpp double);",table_name,table_name)
		)
}

upload_table = function(colocalization,table_name){
	message('INFO - uploading GWAS')
	step = 1e6
	total_row = nrow(colocalization)
	for (i in seq(1,total_row,step)) {
		start = i
		print(start)
		end = ifelse(i+step-1<=total_row,i+step-1,total_row)
		fwrite(colocalization[start:end,],'coloc_tmp.csv')

		dbExecute(
			conn = locuscompare_pool,
			statement = sprintf("load data local infile 'coloc_tmp.csv' 
			into table %s
			fields terminated by ','
			lines terminated by '\n'
			ignore 1 lines
			(gwas, trait, eqtl, gene_id, gene_name, chr, pos, logp_gwas, logp_eqtl, clpp);",table_name)
			)

		unlink('coloc_tmp.csv')
	}
}

index_table = function(table_name){
	message('INFO - indexing table')
	dbExecute(
		conn = locuscompare_pool,
		statement = sprintf('alter table %s
			add index gwas (gwas),
			add index trait (trait),
			add index eqtl (eqtl),
			add index gene_id (gene_id)',table_name)
		)
}

locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

colocalization = read_colocalization(in_fn)
colocalization = munge_colocalization(colocalization)
colocalization = colocalization[gwas != 'GWAS_Circulating_Metabolites_Kettunen_2016'] # temporarily remove before this table is uploaded.

gencode = get_gencode()
colocalization = merge_with_gencode(colocalization)

table_name = 'eCAVIAR'
create_table(table_name)
upload_table(colocalization,table_name)
index_table(table_name)