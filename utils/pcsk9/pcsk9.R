# detach('package:locuscomparer',unload=TRUE)
# devtools::install_github("boxiangliu/locuscomparer")
library(data.table)
library(locuscomparer)
source('global.R')
library(ggrepel)
library(cowplot)

fig_dir = 'figure/pcsk9/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

PCSK9_coordinate=list(chr='1',start=54505149,end=56505149)

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
setDT(PCSK9_gwas)

PCSK9 = unlist(dbGetQuery(
	conn = locuscompare_pool,
	statement = sprintf(
		"select gene_id 
		from gencode_v19_gtex_v6p
		where gene_name = '%s'",
		'PCSK9'
		)
	))


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
setDT(Adipose_PCSK9_eqtl)


PCSK9_gwas_lead_SNP=PCSK9_gwas[which.min(PCSK9_gwas$pval),rsid]
PCSK9_eqtl_lead_SNP=Adipose_PCSK9_eqtl[which.min(Adipose_PCSK9_eqtl$pval),rsid]

p3 = locuscomparer::main(
	in_fn1 = PCSK9_gwas,
	in_fn2 = Adipose_PCSK9_eqtl,
	title1 = 'CAD GWAS',
	title2 = 'Adipose eQTL',
	vcf_fn = sprintf('%s/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE,
	snp = PCSK9_eqtl_lead_SNP,
	legend_position = 'topright'
	)


p4 = locuscomparer::main(
	in_fn1 = PCSK9_gwas,
	in_fn2 = Adipose_PCSK9_eqtl,
	title1 = 'CAD GWAS',
	title2 = 'Adipose eQTL',
	vcf_fn = sprintf('%s/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE,
	snp = PCSK9_gwas_lead_SNP,
	legend_position = 'topright'
	)



p3$locuszoom1 = p3$locuszoom1 + ylab(expression(atop('CAD GWAS',-log[10]*'(P)')))
p3$locuszoom2 = p3$locuszoom2 + ylab(expression(atop('Adipose eQTL',-log[10]*'(P)')))

p4$locuszoom1 = p4$locuszoom1 + ylab(expression(atop('CAD GWAS',-log[10]*'(P)')))
p4$locuszoom2 = p4$locuszoom2 + ylab(expression(atop('Adipose eQTL',-log[10]*'(P)')))


fig1a_i = p3$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
fig1a=plot_grid(fig1a_i,p4$locuszoom1,ncol=1,rel_heights=c(0.8,1))
plot_blank = ggplot() + geom_blank()
bottom = plot_grid(p3$locuscompare,p4$locuscompare,labels=c('b','c'))
p=plot_grid(fig1a,bottom,nrow=2,labels=c('a',''))
save_plot(sprintf('%s/pcsk9.pdf',fig_dir),p,base_width=8,base_height=8)
