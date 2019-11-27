library(RMySQL)
library(pool)
library(data.table)
library(stringr)
library(dplyr)
source('utils/shared_functions.R')
source('config/config.R')

# Functions:
read_gtex_sex_specific=function(tissue_id, sex){

	# gtex_sex_dir should be specified in config/config.R

	suffix=paste0('.all_associations', sex, '.txt.gz')
	gtex=fread(
		input = paste0("zcat ",gtex_sex_dir,tissue_id,suffix," | sed 's/_/\t/g'"),
		sep = '\t',
		select = c(1,2,3,11),
		col.names = c('gene_id','chr','pos','pval'))
	return(gtex)
}


add_rsid=function(gtex,tkg){
	merged=inner_join(gtex,tkg,by=c('chr','pos'))
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

# Map unversioned gene IDs to the associated GENCODE v19 gene in a
# fast data structure
gencode_v19 = read_gene_anno()
gencode_v19$gene_base = sapply(gencode_v19$gene_id, function(x) {strsplit(x, "\\.")[[1]][1]})

tkg = read_1kg_hg38_chunked()
tkg = tkg[rsid!='.'&!str_detect(rsid,';')]
tkg$chr = paste0("chr", tkg$chr)

for (tissue_id in c("Breast_Mammary_Tissue", "Liver", "Muscle_Skeletal", "Adipose_Subcutaneous",
	"Skin_Not_Sun_Exposed_Suprapubic", "Artery_Aorta", "Brain_Cerebellum",
	"Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta",
	"Artery_Coronary","Artery_Tibial","Brain_Anterior_cingulate_cortex_BA24",
	"Brain_Caudate_basal_ganglia","Brain_Cerebellar_Hemisphere",
	"Brain_Cortex","Brain_Frontal_Cortex_BA9","Brain_Hippocampus","Brain_Hypothalamus",
	"Brain_Nucleus_accumbens_basal_ganglia","Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1",
	"Brain_Substantia_nigra","Cells_EBV-transformed_lymphocytes",
	"Cells_Cultured_fibroblasts","Colon_Sigmoid","Colon_Transverse",
	"Esophagus_Gastroesophageal_Junction","Esophagus_Mucosa","Esophagus_Muscularis",
	"Heart_Atrial_Appendage","Heart_Left_Ventricle","Lung", "Minor_Salivary_Gland",
	"Nerve_Tibial","Ovary","Pancreas","Pituitary","Prostate",
	"Skin_Sun_Exposed_Lower_leg",
	"Small_Intestine_Terminal_Ileum","Spleen","Stomach","Thyroid","Whole_Blood")){

	for (sex in c(".females", ".males_downsample"))
	{
		short_sex = gsub("\\.", "", sex)
		short_sex = gsub("_downsample", "", short_sex)

		print(tissue_id)
		print(sex)
		table_name = paste0('eQTL_',str_replace(tissue_id,'-','_'),'_', short_sex, '_GTEx_v8')

		# Read GTEx data and add relevant metadata
		gtex=read_gtex_sex_specific(tissue_id, sex)
		gtex=add_rsid(gtex,tkg)

		# Trim version numbers from genes
		gtex$new_gencode_gene_id = gtex$gene_id
		gtex$gene_id = sapply(gtex$new_gencode_gene_id, function(x) {strsplit(x, "\\.")[[1]][1]})
		
		# Select the corresponding GENCODE v19 gene names
		# Then upload these IDs instead, unless none exists
		gtex=inner_join(gencode_v19, gtex, by=c("gene_base" = "gene_id"))

		# Now upload the table to MySQL database in chunks

		# Delete old table from database
		if (dbExistsTable(locuscompare_pool,table_name)) {
			dbRemoveTable(locuscompare_pool,table_name)
		}
		
		# Create new table in database
		dbExecute(
			conn = locuscompare_pool, 
			statement = sprintf("create table %s
			(%s_id int auto_increment primary key,
			trait varchar(24),
			rsid varchar(100),
			pval double);",table_name,table_name)
			)

		# Upload data in chunks
		step = 1e6
		total_row = nrow(gtex)
		for (i in seq(1,total_row,step)) {
			start = i
			print(start)
			end = ifelse(i+step-1<=total_row,i+step-1,total_row)
			fwrite(gtex[start:end,c("gene_id", "rsid", "pval")],'gtex_sex_tmp.csv')

			dbExecute(
				conn = locuscompare_pool,
				statement = sprintf("load data local infile 'gtex_sex_tmp.csv' 
				into table %s 
				fields terminated by ','
				lines terminated by '\n'
				ignore 1 lines
				(trait, rsid, pval);",table_name)
				)

			unlink('gtex_sex_tmp.csv')
		}

		# Index the newly added table
		dbExecute(
			conn = locuscompare_pool,
			statement = sprintf('alter table %s
				add index trait (trait),
				add index rsid (rsid)',table_name)
			)

		# Create or update the table mapping Ensembl gene IDs
		# to gene names
		if (dbExistsTable(locuscompare_pool,sprintf("%s_trait", table_name))) {
			dbRemoveTable(locuscompare_pool,sprintf("%s_trait", table_name))
		}

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
		rm(gtex)
	}
}

