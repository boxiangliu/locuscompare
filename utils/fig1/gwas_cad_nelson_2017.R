library(data.table)
library(locuscomparer)
source('global.R')
source('utils/fig1/shared.R')

# Reading data:
TGFB1_coordinate=list(chr='19',start=40807492,end=42807492)


TGFB1_gwas=get_gwas(TGFB1_coordinate,mysql_table='GWAS_Coronary_Artery_Disease_Nelson_2017',mysql_trait='GWAS_Coronary-Artery-Disease_Nelson_2017')
TGFB1_ensg = get_gene_id('TGFB1')
TGFB1_coloc_tissues = c(
	Liver = 'Liver',
	Adipose_Visceral_Omentum = 'Adipose - Visceral',
	Adipose_Subcutaneous = 'Adipose_SubQ',
	Muscle_Skeletal = 'Muscle - Skeletal',
	Thyroid = 'Thyroid',
	Heart_Left_Ventricle = 'Heart - LV'
	)

TGFB1_plots = list()
for (i in seq_along(TGFB1_coloc_tissues)){
	tissue_mysql = names(TGFB1_coloc_tissues)[i]
	tissue_plot = TGFB1_coloc_tissues[i]
	eqtl = get_eqtl(tissue_mysql,TGFB1_ensg,TGFB1_coordinate)

	p = locuscomparer::main(
		in_fn1 = TGFB1_gwas,
		in_fn2 = eqtl,
		title1 = 'CAD GWAS',
		title2 = paste(tissue_plot,'eQTL'),
		vcf_fn = '/mnt/data/shared/1KG/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
		snp = 'rs56254331',
		combine = FALSE
	)

	TGFB1_plots[[tissue_mysql]] = p

	# Locuscompare:
	save_plot(sprintf('figure/fig1/gwas_cad_nelson_2017/TGFB1_%s_locuscompare.pdf',tissue_mysql),p$locuscompare)
	
	# Locuszoom:
	p1 = p$locuszoom2+theme(axis.text.x = element_blank(), axis.title.x = element_blank())
	p2 = plot_grid(p1,p$locuszoom1,ncol=1,align='v',rel_heights=c(0.8,1))
	save_plot(sprintf('figure/fig1/gwas_cad_nelson_2017/TGFB1_%s_locuszoom.pdf',tissue_mysql),p2)
}

save(TGFB1_plots,file='processed_data/fig1/gwas_cad_nelson_2017.rda')
