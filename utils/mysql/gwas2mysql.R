library(RMySQL)
library(pool)
library(data.table)
library(stringr)
library(DBI)
source('config/config.R')
source('utils/shared_functions.R')


locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

for (gwas_fn in c("GWAS_Asthma_Moffatt_2010.txt","GWAS_BMI-Europeans_Locke_2015.txt",
	"GWAS_BMI-Mixed_Locke_2015.txt","GWAS_Coronary-Heart-Disease_Nikpay_2015.txt",
	"GWAS_Ischemic-Stroke_Malik_2016.txt","GWAS_Rheumatoid-Arthritis-Asian_Okada_2014.txt",
	"GWAS_Rheumatoid-Arthritis-European_Okada_2014.txt",
	"GWAS_Rheumatoid-Arthritis-Mixed_Okada_2014.txt","GWAS_Type-2-Diabetes_Scott_2017.txt",
	"GWAS_Waist-Hip-Ratio-Europeans_Shungin_2015.txt",
	"GWAS_Waist-Hip-Ratio-Mixed_Shungin_2015.txt","GWAS_Coronary-Artery-Disease_Nelson_2017.txt",
	"GWAS_Type-2-Diabetes_Zhao_2017.txt")){

	print(gwas_fn)
	
	gwas=tryCatch({
		fread(paste0(gwas_dir,gwas_fn))
		}, 
		error = function(e) {NULL}
		)
	if (is.null(gwas)){
		next
	}

	trait=str_replace(gwas_fn,'.txt','')
	gwas$trait=trait

	table_name = str_replace_all(trait,'-','_')

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

	dbWriteTable(
		conn = locuscompare_pool,
		name = table_name,
		value = gwas[,list(trait,rsid,pval)],
		row.names=FALSE,
		append=TRUE
		)

	dbExecute(
		conn = locuscompare_pool,
		statement = sprintf('alter table %s
			add index trait (trait),
			add index rsid (rsid)',table_name)
		)
}
