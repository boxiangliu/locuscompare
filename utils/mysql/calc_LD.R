source('config/config.R')

out_dir='/srv/persistent/bliu2/locuscompare/data/ld/'
if (!dir.exists(out_dir)) {dir.create(out_dir,recursive=TRUE)}

calc_LD=function(chr,pop,out_dir,vcf_fn,panel=NULL){
	pop_fn=sprintf('%s/population/%s.txt',data_dir,pop)

	command=sprintf('%s --vcf %s --keep-allele-order --maf 0.01 --keep %s --r2 --ld-window 9999999 --ld-window-kb 10000 --out %s/%s_%s',plink,vcf_fn,pop_fn,out_dir,chr,pop)
	print(command)
	system(command)
}


for (pop in c('EUR','EAS','SAS','AFR','AMR')){
	message('Population: ', pop)
	for(i in 1:22){
		chrom = paste0('chr',i)
		message('Chromosome: ', chrom)
		calc_LD(
			chr = chrom,
			pop = pop,
			out_dir = out_dir,
			vcf_fn = sprintf('%s/ALL.%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',tkg_dir,chrom)
		)
	}
	message('Chromosome: X')
	calc_LD(
		chr = 'chrX',
		pop = pop,
		out_dir = out_dir,
		vcf_fn = sprintf('%s/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz',tkg_dir)
	)
}