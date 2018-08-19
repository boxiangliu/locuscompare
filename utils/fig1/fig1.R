library(data.table)
# detach('package:locuscomparer',unload=TRUE)
# devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)
source('global.R')


# Reading data:
PHACTR1_coordinate=list(chr='6',start=11717037,end=13717037)
PCSK9_coordinate=list(chr='1',start=54505149,end=56505149)

PHACTR1_gwas=dbGetQuery(
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
		'GWAS_Coronary_Heart_Disease_Nikpay_2015',
		'GWAS_Coronary-Heart-Disease_Nikpay_2015',
		PHACTR1_coordinate$chr,
		PHACTR1_coordinate$start,
		PHACTR1_coordinate$end
		)
	)

PCSK9_gwas=dbGetQuery(
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
		'GWAS_Coronary_Heart_Disease_Nikpay_2015',
		'GWAS_Coronary-Heart-Disease_Nikpay_2015',
		PCSK9_coordinate$chr,
		PCSK9_coordinate$start,
		PCSK9_coordinate$end
		)
	)

PHACTR1 = unlist(dbGetQuery(
	conn = locuscompare_pool,
	statement = sprintf(
		"select gene_id 
		from gencode_v19_gtex_v6p
		where gene_name = '%s'",
		'PHACTR1'
		)
	))

PCSK9 = unlist(dbGetQuery(
	conn = locuscompare_pool,
	statement = sprintf(
		"select gene_id 
		from gencode_v19_gtex_v6p
		where gene_name = '%s'",
		'PCSK9'
		)
	))

Coronary_Artery_PHACTR1_eqtl = dbGetQuery(
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
		'eQTL_Artery_Coronary_GTEx_v6p',
		PHACTR1,
		PHACTR1_coordinate$chr,
		PHACTR1_coordinate$start,
		PHACTR1_coordinate$end
		)
	)

Lung_PHACTR1_eqtl = dbGetQuery(
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
		'eQTL_Lung_GTEx_v6p',
		PHACTR1,
		PHACTR1_coordinate$chr,
		PHACTR1_coordinate$start,
		PHACTR1_coordinate$end
		)
	)

Adipose_PCSK9_eqtl = dbGetQuery(
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
		'eQTL_Adipose_Visceral_Omentum_GTEx_v6p',
		PCSK9,
		PCSK9_coordinate$chr,
		PCSK9_coordinate$start,
		PCSK9_coordinate$end
		)
	)



p1 = locuscomparer::main(
	in_fn1 = PHACTR1_gwas,
	in_fn2 = Coronary_Artery_PHACTR1_eqtl,
	title1 = 'CAD GWAS',
	title2 = 'Coronary Artery eQTL',
	vcf_fn = sprintf('%s/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE
	)

p2 = locuscomparer::main(
	in_fn1 = PHACTR1_gwas,
	in_fn2 = Lung_PHACTR1_eqtl,
	title1 = 'CAD GWAS',
	title2 = 'Lung eQTL',
	vcf_fn = sprintf('%s/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE,
	snp = 'rs9349379'
	)

PCSK9_gwas_lead_SNP=PCSK9_gwas[which.min(PCSK9_gwas$pval),rsid]
PCSK9_eqtl_lead_SNP=Adipose_PCSK9_eqtl[which.min(Adipose_PCSK9_eqtl$pval),rsid]

p3 = locuscomparer::main(
	in_fn1 = PCSK9_gwas,
	in_fn2 = Adipose_PCSK9_eqtl,
	title1 = 'CAD GWAS',
	title2 = 'Adipose eQTL',
	vcf_fn = sprintf('%s/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE,
	snp = PCSK9_eqtl_lead_SNP
	)

p4 = locuscomparer::main(
	in_fn1 = PCSK9_gwas,
	in_fn2 = Adipose_PCSK9_eqtl,
	title1 = 'CAD GWAS',
	title2 = 'Adipose eQTL',
	vcf_fn = sprintf('%s/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE,
	snp = PCSK9_gwas_lead_SNP
	)


fig1c_i = p3$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
fig1c=plot_grid(fig1c_i,p4$locuszoom1,ncol=1,rel_heights=c(0.8,1))
p=plot_grid(p1$locuscompare,p2$locuscompare,fig1c,p3$locuscompare,labels=c('A','B','C','D'))
save_plot('figure/fig1.pdf',p,base_width=8,base_height=8)
