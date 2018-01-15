library(RMySQL)
library(pool)

gwas_fn='data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz'
eqtl_fn='data/PHACTR1.eqtl.txt'

gwas=fread(paste0('gunzip -c ',gwas_fn),select=c(2:4,10),col.names=c('rsid','chr','pos','pval'))
eqtl=fread(eqtl_fn,select=c(2,4),col.names=c('SNP','pval'))

eqtl[,chr:=as.integer(str_split_fixed(SNP,'_',5)[,1])]
eqtl[,pos:=as.integer(str_split_fixed(SNP,'_',5)[,2])]

eqtl=merge(eqtl,gwas[,list(chr,pos,rsid)],by=c('chr','pos'))
eqtl=eqtl[,list(rsid,chr,pos,pval)]

locuscompare_db <- dbPool(
    RMySQL::MySQL(), 
    dbname = "locuscompare",
    host = "rds-mysql-locuscompare.cbhpzvkzr3rc.us-west-1.rds.amazonaws.com",
    username = "admin",
    password = "12345678"
)

locuscompare_db <- dbPool(
  RMySQL::MySQL(), 
  dbname = "locuscompare",
  host = "localhost",
  username = "root",
  password = "admin"
)


dbWriteTable(conn = locuscompare_db,
             name = 'GWAS_coronaryArteryDisease_nelson_2017',
             value = gwas[chr==6],
             field.types=c(rsid='varchar(100)',chr='varchar(4)',pos='int',pval='double(20,10)'),
             row.names=FALSE,
             overwrite=TRUE)

dbWriteTable(conn = locuscompare_db,
             name = 'eQTL_AdiposeSubcutaneous_GTEx_2017',
             value = eqtl,
             field.types=c(rsid='varchar(100)',chr='varchar(4)',pos='int',pval='double(20,10)'),
             row.names=FALSE,
             overwrite=TRUE)
