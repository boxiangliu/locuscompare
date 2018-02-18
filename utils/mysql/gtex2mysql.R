library(RMySQL)
library(pool)
library(data.table)
library(stringr)
source('utils/shared_functions.R')
source('config/config.R')

# Functions:
read_gtex=function(tissue_id){
	gtex_dir='/mnt/data/shared/datasets/gtex/GTEx_Analysis_2015-01-12/eqtl_updated_annotation/v6p_fastQTL_allpairs_FOR_QC_ONLY/'
	suffix='_Analysis.v6p.FOR_QC_ONLY.allpairs.txt.gz'
	gtex=fread(paste0('zcat ',gtex_dir,tissue_id,suffix),select=c(1,2,4),col.names=c('gene_id','snp','pval'))
	split_snp=str_split_fixed(gtex$snp,'_',5)
	gtex[,chr:=split_snp[,1]]
	gtex[,pos:=as.integer(split_snp[,2])]
	gtex[,snp:=NULL]
	return(gtex)
}


add_rsid=function(gtex,tkg){
	merged=merge(gtex,tkg,by=c('chr','pos'))
	return(merged)
}

add_gene_name=function(gtex,anno){
	merged=merge(gtex,anno,by='gene_id')
	return(merged)
}


# Main:
locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

tkg = read_1kg()
tkg = tkg[rsid!='.']

for (tissue_id in c(
	"Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta",
	"Artery_Coronary","Artery_Tibial","Brain_Anterior_cingulate_cortex_BA24",
	"Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere","Brain_Cerebellum",
	"Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus",
	"Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia",
	"Breast_Mammary_Tissue","Cells_EBV-transformed_lymphocytes",
	"Cells_Transformed_fibroblasts","Colon_Sigmoid","Colon_Transverse",
	"Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis",
	"Heart_Atrial_Appendage","Heart_Left_Ventricle","Liver","Lung","Muscle_Skeletal",
	"Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate",
	"Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg",
	"Small_Intestine_Terminal_Ileum","Spleen","Stomach","Testis","Thyroid","Uterus",
	"Vagina","Whole_Blood")){
	print(tissue_id)
	table_name = paste0('eQTL_',str_replace(tissue_id,'-','_'),'_GTEx_v6p')

	gtex=read_gtex(tissue_id)
	gtex=add_rsid(gtex,tkg)

	if (dbExistsTable(locuscompare_pool,table_name)) {
		dbRemoveTable(locuscompare_pool,table_name)
	}

	dbExecute(
		conn = locuscompare_pool, 
		statement = sprintf("create table %s
		(%s_id int auto_increment primary key,
		trait varchar(24),
		rsid varchar(100),
		pval double);",table_name,table_name)
		)

	fwrite(gtex[,list(trait = gene_id, rsid, pval)],'tmp.csv')
	
	dbExecute(
		conn = locuscompare_conn,
		statement = sprintf("load data local infile 'tmp.csv' 
		into table %s 
		fields terminated by ','
		lines terminated by '\n'
		ignore 1 lines
		(trait, rsid, pval);",table_name)
		)

	dbExecute(
		conn = locuscompare_pool,
		statement = sprintf('alter table %s
			add index trait (trait),
			add index rsid (rsid)',table_name)
		)
}