library(data.table)
library(stringr)
library(RMySQL)
library(pool)
library(DBI)
source('config/config.R')

in_fn = '/srv/persistent/bliu2/rpe/processed_data/finemap/manhattan/2018-06-23_21-47-08_rpe_amd_gtex/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'


read_colocalization = function(fn){
    x = fread(fn)
    return(x)
}

munge_colocalization = function(x){
	y = x[,list(gwas = base_gwas_file,
				trait = gwas_trait,
				eqtl = eqtl_file,
				gene_id = feature,
				clpp = clpp)]
	# update once Mike provides new format:
	y[,gwas := 'GWAS_Age_Related_Macular_Degeneration_Fritsche_2013']
	y[,trait := 'Advanced-vs-Controls']
	y[,eqtl := sprintf('eQTL_%s_GTEx_v6p',str_replace(eqtl,'_allpairs_txt_gz',''))]
	return(y)
}

get_gencode = function(){
	gencode = dbGetQuery(
        conn = locuscompare_pool,
        statement = "select gene_id, chr, start 
			from gencode_v19_gtex_v6p;"
    )
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
			(gwas, trait, eqtl, gene_id, clpp);",table_name)
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
gencode = get_gencode()
colocalization = merge(colocalization,gencode,by='gene_id')
setnames(colocalization,'start','pos')

table_name = 'eCAVIAR'
create_table(table_name)
upload_table(colocalization,table_name)
index_table(table_name)
