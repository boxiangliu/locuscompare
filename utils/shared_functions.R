source('config/config.R')

read_1kg=function(select=1:3,col.names=c('chr','pos','rsid')){
  vcf=fread(cmd = sprintf('zcat %s/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz',tkg_dir),
            select=select,
            col.names=col.names,
            skip = 253)
  return(vcf)
}

read_gene_anno=function(select=5:7,col.names=c('gene_id','gene_name','type')){
  anno=fread('/mnt/data/shared/datasets/gtex/GTEx_Analysis_2015-01-12/extra/gencode.v19.genes.v6p.hg19.bed',
             select = select,
             col.names = col.names)
  return(anno)
}

connect_database=function(dbname,host='localhost',username='root',password='admin'){
  locuscompare_db <- dbPool(
    RMySQL::MySQL(), 
    dbname = dbname,
    host = host,
    username = username,
    password = password
  )
}