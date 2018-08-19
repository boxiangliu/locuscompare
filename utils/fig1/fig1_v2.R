library(data.table)
library(locuscomparer)
source('global.R')
source('utils/fig1/shared.R')
library(ggrepel)
library(cowplot)



out_dir = 'processed_data/fig1_v2/'
fig_dir = 'figure/fig1_v2/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}


combine_plots = function(p,labels=letters[1:3]){
	a = p$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
	ab =plot_grid(a,p$locuszoom1,ncol=1,rel_heights=c(0.8,1),labels=labels[1:2])
	abc = plot_grid(ab,p$locuscompare,ncol=2,rel_widths=c(1,0.9),labels=c('',labels[3]))
	return(abc)
}

ARMS2_coordinate=list(chr='10',start=123214169,end=125214169)

ARMS2_gwas = get_gwas(ARMS2_coordinate)
ARMS2_ensg = get_gene_id('ARMS2')
ARMS2_coloc_tissues = c(
	Testis='Testis', 
	Esophagus_Muscularis='Esophagus (Musc.)', 
	Skin_Sun_Exposed_Lower_leg = 'Lower leg skin', 
	Artery_Tibial = 'Tibial artery', 
	Esophagus_Gastroesophageal_Junction = 'Esophagus (GEJ)'
	)

ARMS2_plots = list()
for (i in seq_along(ARMS2_coloc_tissues)){
	tissue_mysql = names(ARMS2_coloc_tissues)[i]
	tissue_plot = ARMS2_coloc_tissues[i]
	eqtl = get_eqtl(tissue_mysql,ARMS2_ensg,ARMS2_coordinate)

	p = locuscomparer::main(
		in_fn1 = ARMS2_gwas,
		in_fn2 = eqtl,
		title1 = 'T2D GWAS',
		title2 = paste(tissue_plot,'eQTL'),
		vcf_fn = sprintf('%s/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir),
		snp = 'rs2421016',
		combine = FALSE,
		legend_position = ifelse(tissue_mysql=='Testis','topright','bottomright')
	)

	p$locuszoom1 = p$locuszoom1 + ylab(bquote(atop('T2D GWAS',-log[10]*'(P)')))
	p$locuszoom2 = p$locuszoom2 + ylab(bquote(atop(paste(.(tissue_plot),' eQTL'),-log[10]*'(P)')))
	ARMS2_plots[[tissue_mysql]] = p
}

out_fn = sprintf('%s/ARMS2.rds',out_dir)
saveRDS(ARMS2_plots,out_fn)
# ARMS2_plots = readRDS(out_fn)


# Figure 1:
tissue = 'Artery_Tibial'
p = ARMS2_plots[[tissue]]
entire = combine_plots(p,c('b','c','d'))
save_plot(sprintf('%s/%s.pdf',fig_dir,tissue),entire,base_width=8.5,base_height=4)

# Figure S1:
plot_ls = list()
i=0
for (tissue in names(ARMS2_coloc_tissues)[c(5,2,3,1)]){
	i = i + 1
	plot_ls[[i]] = ARMS2_plots[[tissue]]$locuscompare
}

entire = plot_grid(plotlist = plot_ls,nrow=2,labels=letters[1:4])
save_plot(sprintf('%s/fig_s1.pdf',fig_dir),entire,base_width=8,base_height=8)
