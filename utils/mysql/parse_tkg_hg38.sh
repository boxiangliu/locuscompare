for i in {1..22} X Y; do 
bcftools query -f \
"%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/EAS_AF\t%INFO/AMR_AF\t%INFO/AFR_AF\t%INFO/EUR_AF\t%INFO/SAS_AF\n" \
-H \
/srv/persistent/bliu2/shared/1000genomes/phase3v5a/GRCh38/ALL.chr${i}_GRCh38_sites.20170504.vcf.gz \
> /srv/persistent/bliu2/locuscompare/data/1KG/GRCh38/chr$i.txt &
done 
