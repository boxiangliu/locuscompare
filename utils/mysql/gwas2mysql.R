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

# for (gwas_fn in c("GWAS_Asthma_Moffatt_2010.txt","GWAS_BMI-Europeans_Locke_2015.txt",
# 	"GWAS_BMI-Mixed_Locke_2015.txt","GWAS_Coronary-Heart-Disease_Nikpay_2015.txt",
# 	"GWAS_Ischemic-Stroke_Malik_2016.txt","GWAS_Rheumatoid-Arthritis-Asian_Okada_2014.txt",
# 	"GWAS_Rheumatoid-Arthritis-European_Okada_2014.txt",
# 	"GWAS_Rheumatoid-Arthritis-Mixed_Okada_2014.txt","GWAS_Type-2-Diabetes_Scott_2017.txt",
# 	"GWAS_Waist-Hip-Ratio-Europeans_Shungin_2015.txt",
# 	"GWAS_Waist-Hip-Ratio-Mixed_Shungin_2015.txt","GWAS_Coronary-Artery-Disease_Nelson_2017.txt",
# 	"GWAS_Type-2-Diabetes_Zhao_2017.txt")){

gwas_dir = '/mnt/data/shared/datasets/gwas/batch2/'
for (gwas_fn in c("GWAS_Allergies_Ferreira_2017.txt","GWAS_Asthma_Demenais_2017.txt",
	"GWAS_Atopic-Dermatitis_Hirota_2012.txt","GWAS_Atrial-Fibrillation_Low_2017.txt",
	"GWAS_BMI_Akiyama_2017.txt","GWAS_BMI-Exome_Turcot_2017.txt",
	"GWAS_Breast-Cancer-ER-Negative_Milne_2017.txt",
	"GWAS_Breast-Cancer-ER-Positive_Milne_2017.txt","GWAS_Breast-Cancer_Milne_2017.txt",
	"GWAS_Coronary-Artery-Disease_Howson_2017.txt",
	"GWAS_Estimated-Bone-Mineral-Density_Liu_2017.txt",
	"GWAS_Exfoliation-Syndrome_Aung_2017.txt",
	"GWAS_Insomnia_Hammerschlag_2017.txt","GWAS_Intelligence_Sniekers_2017.txt",
	"GWAS_Myocardial-Infarction_Hirokawa_2015.txt","GWAS_Neuroticism_Luciano_2017.txt",
	"GWAS_Plasma-Lipids-HDL-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-HDL-Euro-Exome_Liu_2017.txt",
	"GWAS_Plasma-Lipids-HDL-European-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-LDL-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-LDL-Euro-Exome_Liu_2017.txt",
	"GWAS_Plasma-Lipids-LDL-European-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-TC-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-TC-Euro-Exome_Liu_2017.txt",
	"GWAS_Plasma-Lipids-TC-European-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-TG-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Plasma-Lipids-TG-Euro-Exome_Liu_2017.txt",
	"GWAS_Plasma-Lipids-TG-European-EastAsian-Exome_Lu_2017.txt",
	"GWAS_Schizophrenia_Li_2017.txt","GWAS_Type-2-Diabetes-Study1_Imamura_2015.txt",
	"GWAS_Type-2-Diabetes-Study2_Imamura_2015.txt","GWAS_Facial-Shape_Claes_2018.txt")){

	print(gwas_fn)
	
	gwas=tryCatch({
		fread(paste0(gwas_dir,gwas_fn))
		}, 
		error = function(e) {NULL}
		)
	if (is.null(gwas)){
		next
	}

	if (!('trait' %in% colnames(gwas))) {
		trait=str_replace(gwas_fn,'.txt','')
		gwas$trait=trait
	}


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
