library(data.table)
library(stringr)
library(cowplot)

clpp_fn = 'processed_data/mysql/coloc2mysql/all_finemap_results_sorted.txt'
sample_size_fn = '/srv/persistent/bliu2/locuscompare/processed_data/colocalization/gtex_sample_size/sample_size.txt'
fig_dir = 'figure/colocalization/sample_size_vs_clpp/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

clpp = fread(clpp_fn)
clpp[, eqtl_file := str_split_fixed(eqtl_file,'_all',2)[,1]]
arms2 = clpp[str_detect(gwas_trait,'Zhao_2017') & str_detect(feature,'ENSG00000254636')]

sample_size = fread(sample_size_fn)
arms2 = merge(arms2,sample_size,by.x = 'eqtl_file', by.y = 'tissue')
arms2$eqtl_file = paste0('  ',arms2$eqtl_file)
p1 = ggplot(arms2, aes(x = size, y = clpp, label = eqtl_file)) + 
	geom_point() + 
	geom_text(hjust = 0,angle = 90) + 
	xlab('Size') + ylab('CLPP') +
	ylim(0,0.03)
fig_fn = sprintf('%s/arms2.pdf',fig_dir)
save_plot(fig_fn, p1)

plekha1 = clpp[str_detect(gwas_trait,'Zhao_2017') & str_detect(feature,'ENSG00000107679')]
plekha1 = merge(plekha1,sample_size,by.x = 'eqtl_file', by.y = 'tissue')
plekha1$eqtl_file = paste0('  ',plekha1$eqtl_file)
p2 = ggplot(plekha1, aes(x = size, y = clpp, label = eqtl_file)) + 
	geom_point() + 
	geom_text(hjust = 0,angle = 90) + 
	xlab('Size') + ylab('CLPP')
fig_fn = sprintf('%s/plekha1.pdf',fig_dir)
save_plot(fig_fn, p2)
