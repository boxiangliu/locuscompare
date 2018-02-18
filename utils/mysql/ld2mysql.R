library(data.table)
library(stringr)
library(RMySQL)
library(pool)
source('utils/shared_functions.R')
source('config/config.R')

in_dir = 'data/ld/'

locuscompare_pool = dbPool(
	drv = RMySQL::MySQL(), 
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password
	)

dbExecute(
	conn = locuscompare_pool, 
	statement = "create table tkg_p3v5a_ld
	(tkg_p3v5a_ld_id int auto_increment primary key,
	SNP_A varchar(100),
	SNP_B varchar(100),
	R2_AFR double,
	R2_AMR double,
	R2_EAS double,
	R2_EUR double,
	R2_SAS double);"
	)


for (chrom in c(1:22,'X')){
	for (POP in c('AFR','AMR','EAS','EUR','SAS')){
		dt = fread(
			input = sprintf('%s/chr%s_%s.ld',in_dir,chrom,POP),
			select = c(3,6,7),
			col.names = c('SNP_A','SNP_B',sprintf('R2_%s',POP))
			)[!str_detect(SNP_A,';')&!str_detect(SNP_B,';')]
		assign(POP, dt)
	}


	ld = AFR
	for (POP in list(AMR,EAS,EUR,SAS)){
		ld = merge(
			x = ld,
			y = POP,
			by = c('SNP_A','SNP_B'),
			sort = FALSE,
			all = TRUE
			)
	}

	fwrite(ld,'ld_tmp.csv')

	dbExecute(
		conn = locuscompare_pool,
		statement = "load data local infile 'ld_tmp.csv' 
		into table tkg_p3v5a_ld
		fields terminated by ','
		lines terminated by '\n'
		ignore 1 lines
		(SNP_A, SNP_B, R2_AFR, R2_AMR, R2_EAS, R2_EUR, R2_SAS);"
		)

	unlink('ld_tmp.csv')
	rm(ld)
}


dbExecute(
	conn = locuscompare_pool,
	statement = 'alter table tkg_p3v5a_ld add index SNP_A (SNP_A), add index SNP_B (SNP_B);'
	)

