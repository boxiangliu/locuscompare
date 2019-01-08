# Upload 1000 genomes data onto mysql
# Boxiang Liu
# 01/25/2018

library(data.table)
library(RMySQL)
library(pool)
library(foreach)
source('utils/shared_functions.R')
source('config/config.R')

locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)


dbExecute(
	conn = locuscompare_pool, 
	statement = "create table tkg_p3v5a_hg19
	(tkg_p3v5a_id int auto_increment primary key,
	chr varchar(4),
	pos int,
	rsid varchar(100),
	ref varchar(130),
	alt varchar(662),
	AF double,
	EAS_AF double,
	AMR_AF double,
	AFR_AF double,
	EUR_AF double,
	SAS_AF double);")

for (i in c(1:22,'X','Y')){
	dbExecute(
		conn = locuscompare_pool,
		statement = sprintf("load data local infile '%s/1KG/GRCh37/chr%s.txt' 
		into table tkg_p3v5a_hg19
		fields terminated by '\t'
		lines terminated by '\n'
		ignore 1 lines
		(chr, pos, rsid, ref, alt, AF, EAS_AF, AMR_AF, AFR_AF, EUR_AF, SAS_AF);",data_dir,i)
		)
}


dbExecute(
	conn = locuscompare_pool,
	statement = 'alter table locuscompare.tkg_p3v5a_hg19 
		add index rsid (rsid), 
		add index chr (chr), 
		add index chr_pos (chr,pos);'
	)