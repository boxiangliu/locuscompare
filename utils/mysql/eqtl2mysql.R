library(RMySQL)
library(pool)
library(data.table)
library(stringr)
source('utils/shared_functions.R')
source('config/config.R')

hcasmc_eqtl_cfg = list(
	fn = '/srv/persistent/bliu2/HCASMC_eQTL/processed_data/rasqual/output_merged/expressed_genes.pval.txt',
	table_name = 'eQTL_HCASMC_Liu_2018',
	rsid_type = 'rsid',
	pval_type = 'pval',
	trait_type = 'name',
	trait_col = 1, 
	rsid_col = 2,
	pval_col = 26)

hcasmc_sqtl_cfg = list(
	fn = '/srv/persistent/bliu2/HCASMC_eQTL/processed_data/sqtl/fastQTL/nominal/all.nominal.txt.gz',
	table_name = 'sQTL_HCASMC_Liu_2018',
	rsid_type = 'rsid',
	pval_type = 'pval',
	trait_type = 'name',
	trait_col = 1,
	rsid_col = 2,
	pval_col = 4)

rpe_eqtl_glucose_cfg = list(
	fn = '/srv/persistent/bliu2/rpe/processed_data/rasqual/output/glucose/joint/all_associations_filt.txt.gz',
	table_name = 'eQTL_RPE-glucose_Liu_2018',
	rsid_type = 'chr+pos',
	pval_type = 'chisq',
	trait_type = 'id_name',
	trait_col = 1,
	rsid_col = c(3,4),
	pval_col = 11)

rpe_eqtl_galactose_cfg = list(
	fn = '/srv/persistent/bliu2/rpe/processed_data/rasqual/output/galactose/joint/all_associations_filt.txt.gz',
	table_name = 'eQTL_RPE-galactose_Liu_2018',
	rsid_type = 'chr+pos',
	pval_type = 'chisq',
	trait_type = 'id_name',
	trait_col = 1,
	rsid_col = c(3,4),
	pval_col = 11)

rpe_sqtl_glucose_cfg = list(
	fn = '/srv/persistent/bliu2/rpe/processed_data/sqtl/fastQTL/nominal/glucose/all.nominal.txt.gz',
	table_name = 'sQTL_RPE-glucose_Liu_2018',
	rsid_type = 'chr_pos',
	pval_type = 'pval',
	trait_type = 'name',
	trait_col = 1,
	rsid_col = 2,
	pval_col = 4)

rpe_sqtl_galactose_cfg = list(
	fn = '/srv/persistent/bliu2/rpe/processed_data/sqtl/fastQTL/nominal/galactose/all.nominal.txt.gz',
	table_name = 'sQTL_RPE-galactose_Liu_2018',
	rsid_type = 'chr_pos',
	pval_type = 'pval',
	trait_type = 'name',
	trait_col = 1,
	rsid_col = 2,
	pval_col = 4)

geuvadis_eqtl_cfg = list(
	fn = '/srv/persistent/bliu2/locuscompare/data/qtl/eQTL_Blood_Geuvadis_sumstats.txt.gz',
	table_name = 'eQTL_Blood_Geuvadis_2013',
	rsid_type = 'chr+pos',
	pval_type = 'pval',
	trait_type = 'id',
	trait_col = 1,
	rsid_col = c(2,3),
	pval_col = 4
	)

eQTLgen_eqtl_cfg = list(
	fn = '/srv/persistent/bliu2/locuscompare/data/qtl/eQTL_Blood_eQTLGen_sumstats.txt.gz',
	table_name = 'eQTL_Blood_eQTLGen_2018',
	rsid_type = 'rsid',
	pval_type = 'pval',
	trait_type = 'name', # the file also has an ENSEMBL ID column
	trait_col = 9,
	rsid_col = 5,
	pval_col = 4
	)

BrainMeta_eqtl_cfg = list(
	fn = '/srv/persistent/bliu2/locuscompare/data/qtl/eQTL_BrainMeta_Qi2018_sumstats.txt.gz',
	table_name = 'eQTL_Brain_Qi_2018',
	rsid_type = 'chr+pos',
	pval_type = 'pval',
	trait_type = 'id',
	trait_col = 1,
	rsid_col = c(2,3),
	pval_col = 4)

BrainMeta_mqtl_cfg = list(
	fn = '/srv/persistent/bliu2/locuscompare/data/qtl/mQTL_BrainMeta_Qi2018_sumstats.txt.gz',
	table_name = 'mQTL_Brain_Qi_2018',
	rsid_type = 'chr+pos',
	pval_type = 'pval',
	trait_type = 'name',
	trait_col = 1,
	rsid_col = c(2,3),
	pval_col = 4)

FetalBrain_mqtl_cfg = list(
	fn = '/srv/persistent/bliu2/locuscompare/data/qtl/mQTL_FetalBrain_Hannon2015_sumstats.txt.gz',
	table_name = 'mQTL_Fetal-Brain_Hannon_2015',
	rsid_type = 'chr+pos',
	pval_type = 'pval',
	trait_type = 'name',
	trait_col = 1,
	rsid_col = c(2,3),
	pval_col = 4)


# Functions:
add_rsid=function(x,tkg){
	merged=merge(x,tkg,by=c('chr','pos'))
	return(merged)
}

add_gene_name=function(x,anno){
	merged=merge(x,anno,by='gene_id')
	return(merged)
}


# Main:
# locuscompare_pool = dbPool(
# 	drv = RMySQL::MySQL(), 
# 	dbname = 'locuscompare',
# 	host = aws_host,
# 	username = aws_username,
# 	password = aws_password
# 	)

conn = DBI::dbConnect(RMariaDB::MariaDB(),group = 'locuscompare')

tkg = read_1kg()
tkg = tkg[rsid!='.'&!str_detect(rsid,';')]

for (cfg in list(FetalBrain_mqtl_cfg))
	# ,hcasmc_sqtl_cfg,
	# rpe_eqtl_glucose_cfg,rpe_eqtl_galactose_cfg,
	# rpe_sqtl_glucose_cfg,rpe_sqtl_galactose_cfg))
{


	select_ = c(cfg$trait_col, cfg$rsid_col, cfg$pval_col)
	if (length(select_) == 3) {
		colnames_ = c('trait','rsid','pval')
	} else {
		colnames_ = c('trait','chr','pos','pval')
	}

	message('Reading ', cfg$fn)
	x = fread(cfg$fn, select = select_, col.names = colnames_)

	if (cfg$rsid_type == 'chr_pos') {
		message('rsID type is chr_pos.')
		chr_pos = str_split_fixed(x$rsid,'_',3)[,c(1,2)]
		x$chr = chr_pos[,1]
		x$pos = chr_pos[,2]
		x$rsid = NULL
	}

	if (cfg['rsid_type'] %in% c('chr+pos','chr_pos')) {
		message('rsID type is chr+pos')
		if (!is.character(x$chr)) mode(x$chr) = 'character'
		x=add_rsid(x,tkg)
		x$chr = NULL
		x$pos = NULL
	}

	if (cfg$trait_type == 'id_name'){
		message('trait type is id_name.')
		name = str_split_fixed(x$trait, '_', 2)[,2]
		x$trait = name
	}

	if (cfg$pval_type == 'chisq'){
		message('pval type is chisq.')
		pval = pchisq(x$pval,1,lower.tail=FALSE)
		x$pval = pval
	}

	table_name = cfg$table_name

	if (dbExistsTable(conn,table_name)) {
		dbRemoveTable(conn,table_name)
	}

	message('Uploading...')
	dbExecute(
		conn = conn, 
		statement = sprintf("create table %s
		(%s_id int auto_increment primary key,
		trait varchar(24),
		rsid varchar(100),
		pval double);",table_name,table_name)
		)


	step = 1e6
	# total_row = nrow(x)
	total_row = 2e6

	for (i in seq(1,total_row,step)) {
		start = i
		print(start)
		end = ifelse(i+step-1<=total_row,i+step-1,total_row)
		fwrite(x[start:end,list(trait = gene_id, rsid, pval)],'x_tmp.csv')

		dbExecute(
			conn = conn,
			statement = sprintf("load data local infile 'x_tmp.csv' 
			into table %s 
			fields terminated by ','
			lines terminated by '\n'
			ignore 1 lines
			(trait, rsid, pval);",table_name)
			)

		unlink('x_tmp.csv')
	}

	message('Indexing...')
	dbExecute(
		conn = conn,
		statement = sprintf('alter table %s
			add index trait (trait),
			add index rsid (rsid)',table_name)
		)

	rm(x)
}