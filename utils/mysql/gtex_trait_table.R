library(RMySQL)
library(pool)
library(data.table)
library(stringr)
source('utils/shared_functions.R')
source('config/config.R')

locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

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


	table_name = paste0('eQTL_',str_replace(tissue_id,'-','_'),'_GTEx_v6p')

	dbExecute(
	conn = locuscompare_pool, 
	statement = sprintf("
		create table %s_trait as 
		(select distinct t2.gene_id,t2.gene_name 
		from %s as t1 
		join gencode_v19_gtex_v6p as t2 
		on t1.trait = t2.gene_id 
		order by t2.gene_name);",table_name,table_name)
	)
}