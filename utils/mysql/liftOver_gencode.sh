awk '{print $1"\t"$2"\t"$3"\t"$4","$5","$6","$7}' data/gencode/gencode.v19.genes.v6p.hg19.bed > data/gencode/tmp.bed
liftOver data/gencode/tmp.bed /srv/persistent/bliu2/tools/ucsc_tools/hg19ToHg38.over.chain.gz data/gencode/tmp.hg38.bed data/gencode/tmp.bed.unMapped
sed 's/,/\t/g' data/gencode/tmp.hg38.bed > data/gencode/gencode.v19.genes.v6p.hg19.bed
rm data/gencode/tmp*