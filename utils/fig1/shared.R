get_gwas = function(coordinate,mysql_table='GWAS_Type_2_Diabetes_Zhao_2017',mysql_trait='Type-2-Diabetes'){
	gwas=dbGetQuery(
		conn = locuscompare_pool,
		statement = sprintf(
			"select t1.rsid, t1.pval 
			from %s as t1 
			join tkg_p3v5a as t2 
			on t1.rsid = t2.rsid 
			where t1.trait = '%s' 
			and t2.chr = '%s' 
			and t2.pos >= %s 
			and t2.pos <= %s;",
			mysql_table,
			mysql_trait,
			coordinate$chr,
			coordinate$start,
			coordinate$end
			)
		)
	return(gwas)
}

get_gene_id = function(gene_name){
	gene_id = unlist(dbGetQuery(
	conn = locuscompare_pool,
	statement = sprintf(
		"select gene_id 
		from gencode_v19_gtex_v6p
		where gene_name = '%s'",
		gene_name
		)
	))
	return(gene_id)
}

get_eqtl = function(tissue,gene_id,coordinate){
	eqtl = dbGetQuery(
		conn = locuscompare_pool,
		statement = sprintf(
			"select t1.rsid, t1.pval 
			from eQTL_%s_GTEx_v6p as t1 
			join tkg_p3v5a as t2 
			on t1.rsid = t2.rsid 
			where t1.trait = '%s' 
			and t2.chr = '%s' 
			and t2.pos >= %s 
			and t2.pos <= %s;",
			tissue,
			gene_id,
			coordinate$chr,
			coordinate$start,
			coordinate$end
			)
		)
	return(eqtl)
}
