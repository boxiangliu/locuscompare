library(data.table)
# detach('package:locuscomparer', unload = TRUE)
# devtools::install_github('boxiangliu/locuscomparer')
library(locuscomparer)
source('global.R')
source('utils/fig1/shared.R')

# Reading data:
PLEKHA1_coordinate=list(chr='10',start=123134094,end=125134094)
ARMS2_coordinate=list(chr='10',start=123214169,end=125214169)
ITFG3_coordinate=list(chr='16',start=0,end=1284545)


PLEKHA1_gwas=get_gwas(PLEKHA1_coordinate)
PLEKHA1_ensg = get_gene_id('PLEKHA1')
PLEKHA1_coloc_tissues = c(
	Nerve_Tibial = 'Nerve - Tibial',
	Adipose_Subcutaneous = 'Adipose - SubQ',
	Lung = 'Lung',
	Ovary = 'Ovary'
	)

PLEKHA1_plots = list()
for (i in seq_along(PLEKHA1_coloc_tissues)){
	tissue_mysql = names(PLEKHA1_coloc_tissues)[i]
	tissue_plot = PLEKHA1_coloc_tissues[i]
	eqtl = get_eqtl(tissue_mysql,PLEKHA1_ensg,PLEKHA1_coordinate)

	p = locuscomparer::main(
		in_fn1 = PLEKHA1_gwas,
		in_fn2 = eqtl,
		title1 = 'T2D GWAS',
		title2 = paste(tissue_plot,'eQTL'),
		snp = 'rs2421016',
		combine = FALSE
	)

	PLEKHA1_plots[[tissue_mysql]] = p

	# Locuscompare:
	save_plot(sprintf('figure/fig1/gwas_t2d_zhao_2017/PLEKHA1_%s_locuscompare.pdf',tissue_mysql),p$locuscompare)
	
	# Locuszoom:
	p1 = p$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
	p2 = plot_grid(p1,p$locuszoom1,ncol=1,align='v',rel_heights=c(0.8,1))
	save_plot(sprintf('figure/fig1/gwas_t2d_zhao_2017/PLEKHA1_%s_locuszoom.pdf',tissue_mysql),p2)
}


ARMS2_gwas=get_gwas(ARMS2_coordinate)
ARMS2_ensg = get_gene_id('ARMS2')
ARMS2_coloc_tissues = c(
	Testis='Testis', 
	Esophagus_Muscularis='Esophagus - Muscularis', 
	Skin_Sun_Exposed_Lower_leg = 'Skin - Lower Leg', 
	Artery_Tibial = 'Artery - Tibial', 
	Esophagus_Gastroesophageal_Junction = 'Esophagus - GEJ',
	Colon_Sigmoid = 'Colon - Sigmoid'
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
		snp = 'rs2421016',
		combine = FALSE
	)

	ARMS2_plots[[tissue_mysql]] = p

	# Locuscompare:
	save_plot(sprintf('figure/fig1/gwas_t2d_zhao_2017/ARMS2_%s_locuscompare.pdf',tissue_mysql),p$locuscompare)

	# Locuszoom:
	p1 = p$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
	p2 = plot_grid(p1,p$locuszoom1,ncol=1,align='v',rel_heights=c(0.8,1))
	save_plot(sprintf('figure/fig1/gwas_t2d_zhao_2017/ARMS2_%s_locuszoom.pdf',tissue_mysql),p2)
}


ITFG3_gwas=get_gwas(ITFG3_coordinate)
ITFG3_ensg = get_gene_id('ITFG3')
ITFG3_coloc_tissues = c(
	Nerve_Tibial='Nerve - Tibial', 
	Adipose_Subcutaneous='Adipose - SubQ', 
	Artery_Tibial = 'Artery - Tibial', 
	Muscle_Skeletal = 'Muscle - Skeletal',
	Artery_Aorta = 'Artery - Aorta'
	)

ITFG3_plots = list()
for (i in seq_along(ITFG3_coloc_tissues)){
	tissue_mysql = names(ITFG3_coloc_tissues)[i]
	tissue_plot = ITFG3_coloc_tissues[i]
	eqtl = get_eqtl(tissue_mysql,ITFG3_ensg,ITFG3_coordinate)

	p = locuscomparer::main(
		in_fn1 = ITFG3_gwas,
		in_fn2 = eqtl,
		title1 = 'T2D GWAS',
		title2 = paste(tissue_plot,'eQTL'),
		snp = 'rs9940149',
		combine = FALSE
	)

	ITFG3_plots[[tissue_mysql]] = p

	# Locuscompare:
	save_plot(sprintf('figure/fig1/gwas_t2d_zhao_2017/ITFG3_%s_locuscompare.pdf',tissue_mysql),p$locuscompare)

	# Locuszoom:
	p1 = p$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
	p2 = plot_grid(p1,p$locuszoom1,ncol=1,align='v',rel_heights=c(0.8,1))
	save_plot(sprintf('figure/fig1/gwas_t2d_zhao_2017/ITFG3_%s_locuszoom.pdf',tissue_mysql),p2)
}
save(PLEKHA1_plots,ARMS2_plots,ITFG3_plots,file='processed_data/fig1/gwas_t2d_zhao_2017.rda')