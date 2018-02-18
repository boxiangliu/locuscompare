library(data.table)
library(stringr)
library(RMySQL)
library(pool)
source('utils/shared_functions.R')
source('config/config.R')

in_dir = 'data/ld2/'
locuscompare_db=connect_database(
	dbname = 'locuscompare',
	host = aws_host,
	username = aws_username,
	password = aws_password)


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

	dbWriteTable(conn = locuscompare_db,
		name = 'ld',
		value = ld,
		field.types = c(
			SNP_A = 'varchar(100)',
			SNP_B = 'varchar(100)',
			R2_AFR = 'double',
			R2_AMR = 'double',
			R2_EAS = 'double',
			R2_EUR = 'double',
			R2_SAS = 'double'),
		row.names = FALSE,
		append = TRUE)
}

dbExecute(
	conn = locuscompare_db,
	statement = 'ALTER TABLE locuscompare.ld ADD COLUMN ld_id INT NOT NULL AUTO_INCREMENT FIRST, ADD PRIMARY KEY (ld_id);'
	)

dbExecute(
	conn = locuscompare_db,
	statement = 'alter table locuscompare.ld add index SNP_A (SNP_A), add index SNP_B (SNP_B);'
	)

