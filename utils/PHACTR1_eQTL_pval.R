library(data.table)
library(cowplot)
library(stringr)

phactr1 = fread('data/example/PHACTR1.tsv')
color = fread('data/gtex/gtex_tissue_colors.txt')
color_map = paste0('#',color$tissue_color_hex)
names(color_map) = color$tissue_site_detail

name_map = color$tissue_site_detail
names(name_map) = str_replace(color$tissue_site_detail_id,'-','_')

# rs9349379 is the lead SNP for PHACTR1 GWAS loci,
# rs66974016 is the only SNP in LD (r2>=0.5) with rs9349379.
eqtl_pval = phactr1[rsid %in% c('rs66974016','rs9349379'),list(min_p = min(pval)),by = 'trait']
setorder(eqtl_pval,-min_p)
eqtl_pval[,logp:=-log10(min_p)]
eqtl_pval[,trait:=name_map[trait]]
eqtl_pval[,trait:=factor(trait,trait)]


x_breaks=c(0,round(-log10(0.05/44),2),5,10,15)
p = ggplot(eqtl_pval, aes(x = trait, y = logp, fill = trait))+
	geom_bar(stat='identity')+
	xlab('') + ylab('-log10(P)') +
	geom_hline(yintercept = -log10(0.05/44), color = 'red', linetype = 2) +  
	coord_flip() +
	scale_fill_manual(values = color_map, guide = 'none') + 
	scale_y_continuous(breaks = x_breaks, labels = x_breaks)



save_plot('figure/PHACTR1_eQTL.png',p,base_width=8,base_height=9)