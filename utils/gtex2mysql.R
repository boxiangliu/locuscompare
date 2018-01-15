library(RMySQL)
library(pool)
library(data.table)
library(stringr)

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

read_1kg=function(){
  vcf=fread('zcat /mnt/data/shared/1KG/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz',
            select=1:3,
            col.names=c('chr','pos','rsid'),
            skip = 253)
  return(vcf)
}

read_gene_anno=function(){
  anno=fread('/mnt/data/shared/datasets/gtex/GTEx_Analysis_2015-01-12/extra/gencode.v19.genes.v6p.hg19.bed',
             select = 5:7,
             col.names = c('gene_id','gene_name','type'))
  return(anno)
}

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
locuscompare_db=connect_database('locuscompare')
tkg=read_1kg()
anno=read_gene_anno()

for (tissue_id in c("Adipose_Subcutaneous","Adipose_Visceral_Omentum","Adrenal_Gland","Artery_Aorta",
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
  
  gtex=read_gtex(tissue_id)
  gtex=add_rsid(gtex,tkg)
  gtex=add_gene_name(gtex,anno)
  
  dbWriteTable(conn = locuscompare_db,
               name = paste0('traits_',str_replace(tissue_id,'-','_'),'_GTEx_2017'),
               value = unique(gtex[,list(gene_name,gene_id,type)]),
               field.types = c(gene_name='varchar(19)',gene_id='varchar(18)',type='varchar(24)'),
               row.names=FALSE)
  
  dbWriteTable(conn = locuscompare_db,
               name = paste0('eQTL_',str_replace(tissue_id,'-','_'),'_GTEx_2017'),
               value = gtex,
               field.types=c(gene_name='varchar(19)',gene_id='varchar(24)',type='varchar(24)',rsid='varchar(100)',chr='varchar(4)',pos='int',pval='double'),
               row.names=FALSE)
  
}



             


