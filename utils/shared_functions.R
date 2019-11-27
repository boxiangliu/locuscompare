source('config/config.R')

read_1kg=function(select=1:3,col.names=c('chr','pos','rsid')){
  vcf=fread(cmd = sprintf('zcat %s/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz',tkg_dir),
            select=select,
            col.names=col.names,
            skip = 253)
  return(vcf)
}

read_1kg_hg38_chunked=function(select=1:3,col.names=c('chr','pos','rsid'))
{
  kg_list = list()
  for (chrom in c(as.character(1:22), "X", "Y"))
  {
	  print(chrom)
	  # tkg_hg38_dir should be specified in config/config.R
	  kg_list[[chrom]]=fread(sprintf('%s/ALL.chr%s_GRCh38.genotypes.20170504.vcf.gz',tkg_hg38_dir, chrom),
		    select=select,
		    col.names=col.names,
		    skip = 131)
  }
  vcf = do.call(rbind, kg_list)
  return(vcf)
}

read_gene_anno=function(select=5:7,col.names=c('gene_id','gene_name','type'), fn = 'data/gencode/gencode.v19.genes.v6p.hg19.bed'){
  anno=fread(fn,
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


