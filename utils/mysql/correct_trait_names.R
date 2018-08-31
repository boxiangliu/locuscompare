# Correct trait names check if the trait column have the right trait names
# and correct them if not. 
source('global.R')

# read all GWAS tables 
get_gwas_tables = function(conn) {

	table_names = dbListTables(conn)
	idx = which(str_detect(table_names, 'GWAS') & !str_detect(table_names,'_trait'))
	gwas_tables = table_names[idx]
	
	return(gwas_tables)

}



# check trait table exist
check_trait_table = function(conn,gwas_table){
	
	trait_table = paste0(gwas_table,'_trait')

	exists = dbExistsTable(
		conn = locuscompare_pool,
		name = trait_table
	)

	return(exists)

}


# get distinct trait name
get_trait_name = function(conn,table){
	
	statement = sprintf('select distinct trait from %s;',table)
	trait = dbGetQuery(
		conn = conn,
		statement = statement
	)

	return(trait$trait)

}


# check whether trait name is part of table name
trait_in_table = function(trait,table){
	
	trait = str_replace_all(trait,'-','_')
	detected = str_detect(table,trait)
	return(detected)

}



########
# Main #
########
gwas_table = get_gwas_tables(locuscompare_pool)
counter = 1
for (table in gwas_table){

	counter = counter + 1
	print(counter)
	
	exists = check_trait_table(locuscompare_pool,table)
	if (!exists){
	
		print(paste0(table,'_trait does not exist!'))
	
	} else {
	
		trait = get_trait_name(locuscompare_pool,table)

		if (length(trait) == 1){

			detected = trait_in_table(trait,table)
			
			if (!detected){

				print(table)
				print(trait)

			}
		}
	}
}