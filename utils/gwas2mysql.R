library(RMySQL)
library(pool)
library(data.table)
library(stringr)
library(DBI)

# Variables:
gwas_dir='/mnt/data/shared/datasets/gwas/'

# Functions:
connect_database=function(dbname){
	locuscompare_db <- dbPool(
		RMySQL::MySQL(), 
		dbname = dbname,
		host = "localhost",
		username = "root",
		password = "admin"
	)
}



locuscompare_db=connect_database('locuscompare')

for (gwas_fn in c("GWAS_Asthma_Moffatt_2010.txt","GWAS_BMI-Europeans_Locke_2015.txt",
	"GWAS_BMI-Mixed_Locke_2015.txt","GWAS_Coronary-Heart-Disease_Nikpay_2015.txt",
	"GWAS_Ischemic-Stroke_Malik_2016.txt","GWAS_Rheumatoid-Arthritis-Asian_Okada_2014.txt",
	"GWAS_Rheumatoid-Arthritis-European_Okada_2014.txt",
	"GWAS_Rheumatoid-Arthritis-Mixed_Okada_2014.txt","GWAS_Type-2-Diabetes_Scott_2017.txt",
	"GWAS_Waist-Hip-Ratio-Europeans_Shungin_2015.txt",
	"GWAS_Waist-Hip-Ratio-Mixed_Shungin_2015.txt")){
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

	
	dbWriteTable(conn = locuscompare_db,
				 name = paste0(str_replace_all(trait,'-','_'),'_traits'),
				 value = data.frame(trait=trait),
				 field.types = c(trait='varchar(100)'),
				 row.names=FALSE,
				 overwrite=TRUE)
	
	dbWriteTable(conn = locuscompare_db,
				 name = str_replace_all(trait,'-','_'),
				 value = gwas,
				 field.types=c(trait='varchar(100)',rsid='varchar(100)',chr='varchar(4)',pos='int',pval='double'),
				 row.names=FALSE,
				 overwrite=TRUE)
	
	dbExecute(conn = locuscompare_db,
			  statement = sprintf('alter table %s add index chr_pos (chr,pos), add index rsid (rsid)',str_replace_all(trait,'-','_')))
}
