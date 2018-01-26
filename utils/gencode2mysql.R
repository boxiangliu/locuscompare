# Upload gencode v19 data onto mysql
# Boxiang Liu
# 01/25/2018

library(data.table)
library(RMySQL)
library(pool)
library(stringr)
source('utils/shared_functions.R')


reference_db = connect_database('reference')
gencode = read_gene_anno(select = 1:7, col.names = c('chr','start','end','strand','gene_id','gene_name','type'))
gencode[,chr:=str_replace(chr,'chr','')]


dbWriteTable(conn = reference_db,
           name = 'gencode_v19_gtex_v6p',
           value = gencode,
           field.types = c(
           				chr='varchar(4)',start='int',end='int',
           				strand='varchar(1)',gene_name='varchar(19)',
           				gene_id='varchar(18)',type='varchar(24)'),
           row.names=FALSE)


dbExecute(
	conn = reference_db,
	statement = 'alter table reference.gencode_v19_gtex_v6p add index gene_name (gene_name), add index gene_id (gene_id), add index type (type)')