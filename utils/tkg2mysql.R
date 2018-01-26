# Upload 1000 genomes data onto mysql
# Boxiang Liu
# 01/25/2018

library(data.table)
library(RMySQL)
library(pool)
library(foreach)
library(doMC)
registerDoMC(20)
source('utils/shared_functions.R')

tkg=foreach(i=c(1:22,'X','Y'),.combine='rbind')%dopar%{
	fread(
		input = sprintf('/srv/persistent/bliu2/locuscompare/data/1KG/chr%s.txt',i),
		skip = 1,
		col.names = c('chr','pos','rsid','ref','alt','AF','EAS_AF','AMR_AF','AFR_AF','EUR_AF','SAS_AF')
	)
}


reference_db = connect_database('reference')


dbWriteTable(conn = reference_db,
           name = 'tkg_p3v5a',
           value = tkg,
           field.types = c(rsid='varchar(100)',chr='varchar(4)',pos='int',ref='varchar(130)',alt='varchar(662)',AF='double',EAS_AF='double',AMR_AF='double',AFR_AF='double',EUR_AF='double',SAS_AF='double'),
           row.names=FALSE)

dbExecute(
	conn = reference_db,
	statement = 'alter table reference.tkg_p3v5a add index rsid (rsid), add index chr (chr), add index chr_pos (chr,pos)')
