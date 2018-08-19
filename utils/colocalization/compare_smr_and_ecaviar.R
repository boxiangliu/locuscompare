library(data.table)
library(stringr)
source('utils/utils/genome_annotation.R')
library(manhattan)
library(cowplot)
library(ggrepel)

smr_fn = 'processed_data/colocalization/Nelson_gtexv7_coronaryartery_pvalthreshold5e-5.out.smr'
ecaviar_fn = 'processed_data/colocalization/UKBB_GWAS1KG_EXOME_CAD_SOFT_META_PublicRelease_300517_txt_gz_finemap_clpp_status.txt'
fig_dir = 'figure/colocalization/compare_smr_and_ecaviar/'
if (!dir.exists(fig_dir)) {dir.create(fig_dir,recursive=TRUE)}

munge_smr = function(x){
	y = x[,list(chrom=paste0('chr',ProbeChr),gene_name=Gene,pos=Probe_bp,y=-log10(p_SMR),gwas_logp=-log10(p_GWAS),eqtl_logp=-log10(p_eQTL),method='SMR')]
	y$gene_type = gencode$gene_type[match(y$gene_name,gencode$gene_name)]
	y = y[gene_type == 'protein_coding',]
	y$gene_type = NULL
	return(y)
}

munge_ecaviar = function(x){
	colnames(x) = c('snp','eqtl_file','gwas_trait','gene_name','n_tested_snps','clpp_score','gwas_logp','eqtl_logp','gwas_file','clpp_mod_score')
	split_snp = str_split_fixed(x$snp,'_',2)
	x$chrom = split_snp[,1]
	x$pos = as.integer(split_snp[,2])
	y = x[,list(chrom=paste0('chr',chrom),pos,gene_name,y=clpp_score,gwas_logp,eqtl_logp,method='eCAVIAR')]
	y = y[eqtl_logp>-log10(5e-5)&gwas_logp>-log10(5e-5),]
	y$gene_type = gencode$gene_type[match(y$gene_name,gencode$gene_id)]
	y = y[gene_type == 'protein_coding',]
	y$gene_type = NULL
	y$gene_name = gencode$gene_name[match(y$gene_name,gencode$gene_id)]
	return(y)
}

g_legend <- function(a.gplot){ 
	tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
	leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
	legend <- tmp$grobs[[leg]] 
	return(legend)
} 

gencode = read_gencode()

color_red='#F8766D'
color_green='#00BA38'
color_blue='#619CFF'

smr_cutoff = 3.4e-5
ecaviar_cutoff = 0.02

smr = fread(smr_fn)
ecaviar = fread(ecaviar_fn)

smr = munge_smr(smr)
ecaviar = munge_ecaviar(ecaviar)


smr_hits = smr[y > -log10(smr_cutoff),gene_name]
ecaviar_hits = ecaviar[y > ecaviar_cutoff,gene_name]
union = union(smr_hits,ecaviar_hits)
intersection = intersect(smr_hits,ecaviar_hits)
smr_only = setdiff(smr_hits,intersection)
ecaviar_only = setdiff(ecaviar_hits,intersection)

smr[,label:=ifelse(gene_name %in% union, gene_name, '')]
ecaviar[,label:=ifelse(gene_name %in% union, gene_name, '')]

smr[,color:=ifelse(gene_name %in% intersection,color_red,NA)]
smr[,color:=ifelse(gene_name %in% ecaviar_only,color_green,color)]
smr[,color:=ifelse(gene_name %in% smr_only,color_blue,color)]

ecaviar[,color:=ifelse(gene_name %in% intersection,color_red,NA)]
ecaviar[,color:=ifelse(gene_name %in% ecaviar_only,color_green,color)]
ecaviar[,color:=ifelse(gene_name %in% smr_only,color_blue,color)]

set.seed(10)
p1=manhattan(smr,build='hg19')+
	geom_hline(yintercept=-log10(smr_cutoff),color='pink',linetype=2)+
	scale_y_continuous(labels=function(x){abs(x)})+
	geom_text_repel(aes(label=label),ylim=c(-log10(smr_cutoff),NA),color='black')+
	annotate(geom='text',x=2.5e9,y=6.5,label='SMR',size=6) + 
	geom_point(data=function(x){x[label!='',]}) + 
	ylab(bquote(-log[10]*'(P)'))


set.seed(42)
p2=manhattan(ecaviar,build='hg19')+
	geom_hline(yintercept=ecaviar_cutoff,color='pink',linetype=2)+
	scale_y_sqrt(labels=function(x){abs(x)})+
	geom_text_repel(aes(label=label),ylim=c(ecaviar_cutoff+0.15,NA),color='black')+
	annotate(geom='text',x=2.5e9,y=0.35,label='eCAVIAR',size=6) + 
	geom_point(data=function(x){x[label!='',]}) + 
	ylab('CLPP')


dummy = data.table(x=1:6,y=1:6,Hits = c('SMR','eCAVIAR','Both'))
p3 = ggplot(data = dummy, aes(x,y,color=Hits)) + geom_point()
legend = g_legend(p3) 

left = plot_grid(p1,p2,ncol=1,align='v',labels=c('a','b'))
p = plot_grid(left,legend,ncol=2,rel_widths=c(0.85,0.15))
fig_fn = sprintf('%s/ecaviar_and_smr.pdf',fig_dir)
save_plot(fig_fn,p,base_height=6,base_width=8)