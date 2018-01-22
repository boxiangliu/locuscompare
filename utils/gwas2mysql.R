library(RMySQL)
library(pool)
library(data.table)
library(stringr)

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

# for (gwas_fn in c("GWAS_2hrGlu_unknown_2000.txt","GWAS_ADHD_unknown_2000.txt",
#                "GWAS_Adiponectin_unknown_2000.txt","GWAS_AgeAtMenopause_unknown_2000.txt",
#                "GWAS_Aggression_unknown_2000.txt","GWAS_BMI-Europeans_unknown_2000.txt",
#                "GWAS_BMI-Mixed_unknown_2000.txt","GWAS_IS_unknown_2000.txt",
#                "GWAS_T2D_unknown_2000.txt","GWAS_WHR-adjBMI-Europeans_unknown_2000.txt",
#                "GWAS_WHR-adjBMI-Mixed_unknown_2000.txt")){
    
for (gwas_fn in c("GWAS_WHR-adjBMI-Europeans_unknown_2000.txt")){
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
    
}
