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

plot_locuscompare = function(coordinate, gwas, gene_id, tissues){

	plots = list()

	for (i in seq_along(tissues)){

		tissue_mysql = names(tissues)[i]
		tissue_plot = tissues[i]
		print(tissue_plot)

		eqtl = get_eqtl(tissue_mysql,gene_id,coordinate)

		vcf_tmp_fn = sprintf('%s_tmp.vcf',format(Sys.time(), "%Y-%b-%d-%H:%M:%S"))
		vcf_fn = sprintf('%s/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir,coordinate$chr)

		command = sprintf('bcftools view %s %s:%s-%s > %s', 
			vcf_fn, 
			coordinate$chr, 
			coordinate$start, 
			coordinate$end,
			vcf_tmp_fn)

		system(command)

		p = locuscomparer::main(
			in_fn1 = gwas,
			in_fn2 = eqtl,
			title1 = 'T2D GWAS',
			title2 = paste(tissue_plot,'eQTL'),
			vcf_fn = vcf_tmp_fn,
			snp = 'rs2421016',
			combine = FALSE,
			legend_position = ifelse(tissue_mysql=='Testis','topright','bottomright')
		)

		p$locuszoom1 = p$locuszoom1 + ylab(bquote(atop('T2D GWAS',-log[10]*'(P)')))
		p$locuszoom2 = p$locuszoom2 + ylab(bquote(atop(paste(.(tissue_plot),' eQTL'),-log[10]*'(P)')))
		plots[[tissue_mysql]] = p

		file.remove(vcf_tmp_fn)
	}

	return(plots)

}

ARMS2_coordinate=list(chr='10',start=123214169,end=125214169)

ARMS2_gwas = get_gwas(coordinate = ARMS2_coordinate,
	mysql_table='GWAS_Type_2_Diabetes_Zhao_2017',
	mysql_trait='Type-2-Diabetes')

ARMS2_ensg = get_gene_id('ARMS2')

ARMS2_coloc_tissues = c(
	Testis='Testis', 
	Esophagus_Muscularis='Esophagus (Musc.)', 
	Skin_Sun_Exposed_Lower_leg = 'Lower leg skin', 
	Artery_Tibial = 'Tibial artery', 
	Esophagus_Gastroesophageal_Junction = 'Esophagus (GEJ)'
	)

ARMS2_plots = plot_locuscompare(ARMS2_coordinate,ARMS2_gwas,ARMS2_ensg,ARMS2_coloc_tissues)


PLEKHA1_coordinate=list(chr='10',start=123134094,end=125134094)

PLEKHA1_gwas = get_gwas(coordinate = ARMS2_coordinate,
	mysql_table='GWAS_Type_2_Diabetes_Zhao_2017',
	mysql_trait='Type-2-Diabetes')

PLEKHA1_ensg = get_gene_id('PLEKHA1')

PLEKHA1_coloc_tissues = c(
	Nerve_Tibial='Tibial nerve', 
	Adipose_Subcutaneous='SubQ adipose', 
	Lung = 'Lung', 
	Ovary = 'Ovary'
	)

PLEKHA1_plots = plot_locuscompare(PLEKHA1_coordinate,PLEKHA1_gwas,PLEKHA1_ensg,PLEKHA1_coloc_tissues)

# Figure 1:
p1 = ARMS2_plots[['Testis']]$locuscompare + annotate(geom = 'text', x = 0.5, y = 0.6, label = 'italic(ARMS2)', parse=TRUE)
p2 = PLEKHA1_plots[['Adipose_Subcutaneous']]$locuscompare + annotate(geom = 'text', x = 0.6, y = 0.6, label = 'italic(PLEKHA1)', parse=TRUE)
blank = ggplot() + geom_blank()

bottom = plot_grid(p1, p2, nrow = 1, labels = c('b','c'))
entire = plot_grid(blank, bottom, nrow = 2, labels = c('a',''))
fig_fn = sprintf('%s/missing_a.pdf',fig_dir)
save_plot(fig_fn,entire,base_width=8,base_height=8)

# Figure S1:
plot_ls = list()
i=0
for (tissue in names(ARMS2_coloc_tissues)[c(4,5,2,3)]){
	i = i + 1
	plot_ls[[i]] = ARMS2_plots[[tissue]]$locuscompare
}

entire = plot_grid(plotlist = plot_ls,nrow=2,labels=letters[1:4])
save_plot(sprintf('%s/fig_s1.pdf',fig_dir),entire,base_width=8,base_height=8)

# Figure S2:
plot_ls = list()
i=0
for (tissue in names(PLEKHA1_coloc_tissues)[c(3,1,4)]){
	i = i + 1
	plot_ls[[i]] = PLEKHA1_plots[[tissue]]$locuscompare
}

entire = plot_grid(plotlist = plot_ls,nrow=2,labels=letters[1:4])
save_plot(sprintf('%s/fig_s2.pdf',fig_dir),entire,base_width=8,base_height=8)
