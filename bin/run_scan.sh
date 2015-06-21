# a bash script for quickly running all the commands in the 1000genomes_Selection_Exercise

CHR=$1
START=$2
END=$3
CLST1=$4
CLST2=$5
CLST3=$6
IHS_CLST=$7

echo -e "\n* Preparing clst and iid files"
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panels
awk -v c1=$CLST1 -v c2=$CLST2 -v c3=$CLST3 'NR>1($2 == c1 || $2 == c2 || $2 == c3) {print $1, $1, $2}' integrated_call_samples_v3.20130502.ALL.panel > ${CLST1}_${CLST2}_${CLST3}_clst.txt
awk '{print $1}' ${CLST1}_${CLST2}_${CLST3}_clst.txt > ${CLST1}_${CLST2}_${CLST3}_iids.txt
awk -v c1=$CLST1 -v c2=$CLST2 '($3 == c1 || $3 == c2) {print $0}' ${CLST1}_${CLST2}_${CLST3}_clst.txt > ${CLST1}_${CLST2}_clst.txt
awk -v c1=$CLST1 -v c3=$CLST3 '($3 == c1 || $3 == c3) {print $0}' ${CLST1}_${CLST2}_${CLST3}_clst.txt > ${CLST1}_${CLST3}_clst.txt
awk -v c2=$CLST2 -v c3=$CLST3  '($3 == c2 || $3 == c3) {print $0}' ${CLST1}_${CLST2}_${CLST3}_clst.txt > ${CLST2}_${CLST3}_clst.txt

echo -e "\n* Downloading data"
#tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR}:${START}-${END} | bgzip -c > 1000genomes_phase3_chr${CHR}_${START}_${END}.vcf.gz

echo -e "\n* Filtering data"
bcftools view -v snps \
         -m 2 -M 2 \
         -S ${CLST1}_${CLST2}_${CLST3}_iids.txt \
         -O z \
         1000genomes_phase3_chr${CHR}_${START}_${END}.vcf.gz > ${CLST1}_${CLST2}_${CLST3}_chr${CHR}_${START}_${END}.vcf.gz

echo -e "\nFst"      
plink-1.9 --vcf ${CLST1}_${CLST2}_${CLST3}_chr${CHR}_${START}_${END}.vcf.gz \
          --fst \
          --within ${CLST1}_${CLST2}_clst.txt \
          --double-id \
          --id-delim . \
          --out ${CLST1}_${CLST2}_chr${CHR}_${START}_${END}_filtered \
          --set-missing-var-ids @:#
 
echo -e "\n* PBS"      
plink-1.9 --vcf ${CLST1}_${CLST2}_${CLST3}_chr${CHR}_${START}_${END}.vcf.gz \
          --fst \
          --within ${CLST1}_${CLST3}_clst.txt \
          --double-id \
          --id-delim . \
          --out ${CLST1}_${CLST3}_chr${CHR}_${START}_${END}_filtered \
          --set-missing-var-ids @:#
          
plink-1.9 --vcf ${CLST1}_${CLST2}_${CLST3}_chr${CHR}_${START}_${END}.vcf.gz \
          --fst \
          --within ${CLST2}_${CLST3}_clst.txt \
          --double-id \
          --id-delim . \
          --out ${CLST2}_${CLST3}_chr${CHR}_${START}_${END}_filtered \
          --set-missing-var-ids @:#

echo -e "* iHS"
awk -v c=$IHS_CLST  '($3 == c) {print $1}' ${CLST1}_${CLST2}_${CLST3}_clst.txt > ${IHS_CLST}_iids.txt

bcftools view -S ${IHS_CLST}_iids.txt \
         -q '.05 minor' \
         -O z \
         ${CLST1}_${CLST2}_${CLST3}_chr${CHR}_${START}_${END}.vcf.gz > ${IHS_CLST}_chr${CHR}_${START}_${END}_filtered.vcf.gz

bcftools view ${IHS_CLST}_chr${CHR}_${START}_${END}_filtered.vcf.gz \
         -H | awk '{print $1, $3, "0", $2, $4, $5}' > chr${CHR}_${START}_${END}_filtered.bim

plink-1.9 --bim chr${CHR}_${START}_${END}_filtered.bim \
          --cm-map genetic_map_chr2_combined_b37.txt 2 \
          --make-just-bim \
          --out chr${CHR}_${START}_${END}_filtered_gm
  
awk '{print $1, $2, $3, $4}' chr${CHR}_${START}_${END}_filtered_gm.bim > chr${CHR}_${START}_${END}_filtered_gm.map

selscan --ihs \
        --vcf ${IHS_CLST}_chr${CHR}_${START}_${END}_filtered.vcf.gz \
        --map chr${CHR}_${START}_${END}_filtered_gm.map \
        --out ${IHS_CLST}_chr${CHR}_${START}_${END}_filtered   

norm --ihs --files ${IHS_CLST}_chr${CHR}_${START}_${END}_filtered.ihs.out

echo -e "\n* Plotting results"
Rscript ../bin/plot.R 2 108013601 111013601 CHB FIN YRI CHB
