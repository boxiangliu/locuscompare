# detach('package:locuscomparer',unload=TRUE)
# devtools::install_github("boxiangliu/locuscomparer",force=TRUE)

library(data.table)
library(locuscomparer)
source('global.R')
library(ggrepel)
library(cowplot)

fig_dir = 'figure/pleotropy/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

SORT1_coordinate=list(chr='1',start=108852187,end=110852187)

SORT1_gwas=dbGetQuery(
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
		'GWAS_Coronary_Artery_Disease_Nelson_2017',
		'Coronary-Artery-Disease',
		SORT1_coordinate$chr,
		SORT1_coordinate$start,
		SORT1_coordinate$end
		)
	)
setDT(SORT1_gwas)


LDL_gwas=dbGetQuery(
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
		'GWAS_Metabolic_QTLs_Kanai_2017',
		'Low-Density-Lipoprotein-Cholesterol',
		SORT1_coordinate$chr,
		SORT1_coordinate$start,
		SORT1_coordinate$end
		)
	)
setDT(LDL_gwas)




lead_SNP='rs12740374'

p3 = locuscomparer::main(
	in_fn1 = SORT1_gwas,
	in_fn2 = LDL_gwas,
	title1 = 'CAD GWAS',
	title2 = 'LDL GWAS',
	vcf_fn = sprintf('%s/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
	combine = FALSE,
	snp = lead_SNP,
	legend_position = 'bottomright'
	)

p3$locuszoom1 = p3$locuszoom1 + ylab(expression(atop('CAD GWAS',-log[10]*'(P)')))
p3$locuszoom2 = p3$locuszoom2 + ylab(expression(atop('LDL eQTL',-log[10]*'(P)')))

fig1a_i = p3$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
fig1a=plot_grid(fig1a_i,p3$locuszoom1,ncol=1,rel_heights=c(0.8,1),labels=c('a','b'))
p=plot_grid(fig1a,p3$locuscompare,ncol=2,labels=c('','c'))
save_plot(sprintf('%s/sort1_and_ldl.pdf',fig_dir),p,base_width=8,base_height=4)
