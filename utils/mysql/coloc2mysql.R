library(data.table)
library(stringr)

in_fn = '/srv/persistent/bliu2/rpe/processed_data/finemap/manhattan/2018-06-23_21-47-08_rpe_amd_gtex/Fritsche_sorted_txt_gz_finemap_clpp_status.txt'

read_colocalization = function(fn){
    x = fread(fn)
    return(x)
}

munge_colocalization = function(x){
	y = x[,list(gwas = base_gwas_file,
				trait = gwas_trait,
				eqtl = eqtl_file,
				gene_id = feature,
				clpp = clpp)]
	# update once Mike provides new format:
	y[,gwas := 'GWAS_Age_Related_Macular_Degeneration_Fritsche_2013']
	y[,trait := 'Advanced-vs-Controls']
	y[,eqtl := sprintf('eQTL_%s_GTEx_v6p',str_replace(eqtl,'_allpairs_txt_gz',''))]
	return(y)
}

create_table = function(){
	message('INFO - creating table...')
	if (dbExistsTable(locuscompare_pool,table_name)) {
		dbRemoveTable(locuscompare_pool,table_name)
	}

	dbExecute(
		conn = locuscompare_pool, 
		statement = sprintf("create table %s
		(%s_id int auto_increment primary key,
		gwas varchar(100),
		rsid varchar(100),
		pval double);",table_name,table_name)
		)
}

upload_table = function(){
    
}


colocalization = read_colocalization(in_fn)
colocalization = munge_colocalization(colocalization)

table_name = 'eCAVIAR'
create_table(table_name)
upload_table(table_name)

