# a bash script for quickly running all the commands in the 1000genomes_Selection_Exercise
# Example of a run (call from within data)
# ../bin/run_scan.sh 2 108513601 110513601 CHB FIN YRI CHB
CHR=$1
START=$2
STOP=$3
CLST1=$4
CLST2=$5
CLST3=$6
IHS_CLST=$7


tag=chr${CHR}_${START}_${STOP}

echo -e "\n* Preparing clst and iid files"
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panels
awk 'NR>1{print $1, $1, $2}' integrated_call_samples_v3.20130502.ALL.panel > 1kgv3_clst.txt

echo -e "\n* Downloading data"
#tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz ${CHR}:${START}-${STOP} | bgzip -c > 1000genomes_phase3_${tag}.vcf.gz

echo -e "\n* Filtering to obtain SNVs (30secs)"
#bcftools view -v snps -m 2 -M 2 -O z 1000genomes_phase3_${tag}.vcf.gz > snv_1000genomes_${tag}.vcf.gz

plink-1.9 --vcf snv_1000genomes_$tag.vcf.gz \
          --within 1kgv3_clst.txt \
          --hardy midp \
          --keep-cluster-names $CLST1 \
          --out tmp_$CLST1
          
plink-1.9 --vcf snv_1000genomes_$tag.vcf.gz \
          --within 1kgv3_clst.txt \
          --hardy midp \
          --keep-cluster-names $CLST2 \
          --out tmp_$CLST2
          
plink-1.9 --vcf snv_1000genomes_$tag.vcf.gz \
          --within 1kgv3_clst.txt \
          --hardy midp \
          --keep-cluster-names $CLST3 \
          --out tmp_$CLST3
          
# Pull out all SNPs that fail HWE (with a loose criterion) in at least one population test
awk '$9<1e-5{print $2}' tmp_$CLST1.hwe tmp_$CLST2.hwe tmp_$CLST3.hwe | sort | uniq > tmp.exclude.txt

# Filter out HW failures 
echo -e "\n* Excluding SNPs failing HW test... (30 sec)"
#bcftools view snv_1000genomes_$tag.vcf.gz --exclude 'ID=@tmp.exclude.txt' -O z > snv_1000genomes_${tag}_HWEfiltered.vcf.gz

echo -e "\n* Computing Fst's..."
plink-1.9 --vcf snv_1000genomes_${tag}_HWEfiltered.vcf.gz \
          --fst \
          --within 1kgv3_clst.txt \
          --out ${CLST1}_${CLST2}_${tag}_HWEfiltered \
          --keep-cluster-names $CLST1 $CLST2
          
plink-1.9 --vcf snv_1000genomes_${tag}_HWEfiltered.vcf.gz \
          --fst \
          --within 1kgv3_clst.txt \
          --out ${CLST1}_${CLST3}_${tag}_HWEfiltered \
          --keep-cluster-names $CLST1 $CLST3
          
plink-1.9 --vcf snv_1000genomes_${tag}_HWEfiltered.vcf.gz \
          --fst \
          --within 1kgv3_clst.txt \
          --out ${CLST2}_${CLST3}_${tag}_HWEfiltered \
          --keep-cluster-names $CLST2 $CLST3
  
echo -e "\n* Preparing selscan input..."
  
# write a file of just IIDs from IHS_CLST 
awk -v c=${IHS_CLST} '$3==c{print $1}' 1kgv3_clst.txt > ${IHS_CLST}_iids.txt

# keep IHS_CLST individuals from the filtered vcf and filter out rare variants (at maf < 5%)
# this will be used for haplotype based selection analysis.
bcftools view -S ${IHS_CLST}_iids.txt\
         -q '.05 minor' \
         -O z \
         snv_1000genomes_${tag}_HWEfiltered.vcf.gz > ${IHS_CLST}_${tag}_HWEfiltered.vcf.gz
         
# write bim file from vcf (-H supress header)
bcftools view ${IHS_CLST}_${tag}_HWEfiltered.vcf.gz \
          -H | awk '{print $1, $3, "0", $2, $4, $5}' > ${tag}_HWEfiltered.bim

# get genetic positions and interpolate genetic positions for variants not present in the genetic_map
plink-1.9 --bim CHB_${tag}_HWEfiltered.bim \
          --cm-map genetic_map_chr2_combined_b37.txt 2 \
          --make-just-bim \
          --out ${tag}_HWEfiltered_gm
          
# format as map file
awk '{print $1, $2, $3, $4}' ${tag}_HWEfiltered_gm.bim > ${tag}_HWEfiltered_gm.map

#selscan --ihs \
#        --vcf ${IHS_CLST}_${tag}_HWEfiltered.vcf.gz \
#        --map ${tag}_HWEfiltered_gm.map \
#        --out ${IHS_CLST}_${tag}_HWEfiltered

norm --ihs --bins 20 --files ${IHS_CLST}_${tag}_HWEfiltered.ihs.out 


Rscript ../bin/plot.R $CHR $START $STOP $CLST1 $CLST2 $CLST3 ${IHS_CLST}
