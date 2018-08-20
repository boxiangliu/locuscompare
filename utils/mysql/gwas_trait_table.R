library(RMySQL)
library(pool)
library(data.table)
library(stringr)
source('utils/shared_functions.R')
source('config/config.R')

locuscompare_pool = dbPool(
    drv = RMySQL::MySQL(), 
    dbname = 'locuscompare',
    host = aws_host,
    username = aws_username,
    password = aws_password
)

# Read UKBB mapping: 
UKBB_mapping = fread('data/ukbb/UKBB_ICD10_mapping.tsv',header=FALSE)
setnames(UKBB_mapping,c('trait','display_trait'))
UKBB_mapping[,display_trait:=str_replace(display_trait,'Diagnoses - main ICD10: ','')]

# Get all GWAS tables:
table_names = dbGetQuery(
    conn = locuscompare_pool, 
    statement = "select table_name 
        from information_schema.tables
        where table_name like 'GWAS_%'
        and table_name not like '%_trait'
        and table_name not like '%UKBB%'
        and table_name not like '%Kanai%'"
    )

# Get trait:
for (i in 1:nrow(table_names)){
    tn = table_names[i,]
    trait = dbGetQuery(
        conn = locuscompare_pool,
        statement = sprintf("select distinct trait 
                            from %s",tn)
        )


    # Create GWAS trait table:
    if (str_detect(tn,'UKBB')){
        trait = merge(trait,UKBB_mapping)
    } else {
        result = str_split_fixed(trait$trait,'_',4)
        if (result[1,1]=='GWAS'){
            trait$display_trait = result[,2]
        } else {
            trait$display_trait = result[,1]
        }
    }
    
    print(i)
    print(trait)
    
    # Upload to MySQL:
    trait_table = sprintf('%s_trait',tn)
    if (dbExistsTable(locuscompare_pool,trait_table)) {
        dbRemoveTable(locuscompare_pool,trait_table)
    }
    dbWriteTable(
        conn = locuscompare_pool, 
        name = trait_table,
        value = trait,
        overwrite = TRUE
    )
}