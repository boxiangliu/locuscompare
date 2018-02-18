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
source('config/config.R')

locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)


tkg=foreach(i=c(1:22,'X','Y'),.combine='rbind')%dopar%{
	fread(
		input = sprintf('%s/1KG/chr%s.txt',data_dir,i),
		skip = 1,
		col.names = c('chr','pos','rsid','ref','alt','AF','EAS_AF','AMR_AF','AFR_AF','EUR_AF','SAS_AF')
	)
}

tkg = tkg[rsid!='.']


dbExecute(
	conn = locuscompare_pool, 
	statement = "create table tkg_p3v5a
	(tkg_p3v5a_id int auto_increment primary key,
	rsid varchar(100),
	chr varchar(4),
	pos int,
	ref varchar(130),
	alt varchar(662),
	AF double,
	EAS_AF double,
	AMR_AF double,
	AFR_AF double,
	EUR_AF double,
	SAS_AF double);")

dbWriteTable(
	conn = locuscompare_pool,
	name = 'tkg_p3v5a',
	value = tkg,
	row.names=FALSE,
	append = TRUE
	)

dbExecute(
	conn = locuscompare_db,
	statement = 'alter table locuscompare.tkg_p3v5a 
		add index rsid (rsid), 
		add index chr (chr), 
		add index chr_pos (chr,pos);'
	)