in_dir=data/ld/

create_table(){
POP=$1
mysql -h locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com \
-u root \
-pmontgomerylab \
-D locuscompare \
-e "create table tkg_p3v5a_ld_${POP} 
(tkg_p3v5a_ld_${POP}_id int auto_increment primary key, 
SNP_A varchar(100), 
SNP_B varchar(100),
R2 double);"
}
export -f create_table

drop_table(){
POP=$1
mysql -h locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com \
-u root \
-pmontgomerylab \
-D locuscompare \
-e "drop table tkg_p3v5a_ld_${POP};"
}
export -f drop_table

process_ld(){
in_dir=$1
i=$2
POP=$3
cat $in_dir/chr${i}_${POP}.ld | awk 'BEGIN{OFS="\t"}{if ($3 !~ /;|\./ && $6 !~ /;|\./) {print $3,$6,$7}}' > $in_dir/chr${i}_${POP}.mysql.ld
}
export -f process_ld

ld2mysql(){
in_dir=$1
i=$2
POP=$3
mysql -h locuscompare-us-west-2a.cvocub39nnri.us-west-2.rds.amazonaws.com \
-u root \
-pmontgomerylab \
-D locuscompare \
--local-infile \
--compress \
-e "load data local infile '$in_dir/chr${i}_${POP}.mysql.ld'
into table tkg_p3v5a_ld_${POP}
fields terminated by '\t'
lines terminated by '\n'
ignore 1 lines
(SNP_A, SNP_B, R2);"
}
export -f ld2mysql


echo 'creating tables..'
parallel -j5 create_table {} ::: AFR AMR EAS EUR SAS

for i in {1..22} X; do
echo chr$i
echo 'processing ld..'
parallel -j5 process_ld $in_dir $i {} ::: AFR AMR EAS EUR SAS
echo 'uploding to mysql..'
parallel -j5 ld2mysql $in_dir $i {} ::: AFR AMR EAS EUR SAS
done


