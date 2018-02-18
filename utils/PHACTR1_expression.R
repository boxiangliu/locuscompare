library(data.table)
library(cowplot)


gtex = fread('data/gtex/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct',sep='\t')
tmp = unlist(gtex[Description == 'PHACTR1', !colnames(gtex) %in% c('Name','Description'), with = F])
color = fread('data/gtex/gtex_tissue_colors.txt')
color_map = paste0('#',color$tissue_color_hex)
names(color_map) = color$tissue_site_detail

PHACTR1 = data.table(
	tissue = names(tmp),
	rpkm = tmp)

setorder(PHACTR1,rpkm)
PHACTR1[,tissue := factor(tissue,tissue)]


p=ggplot(PHACTR1,aes(x = tissue, y = rpkm, fill = tissue))+
	geom_bar(stat='identity')+
	xlab('') + ylab('Median RPKM') +
	coord_flip()+ 
	scale_fill_manual(values = color_map, guide = 'none')
save_plot('figure/PHACTR1_expression.png',p,base_width=8,base_height=9)