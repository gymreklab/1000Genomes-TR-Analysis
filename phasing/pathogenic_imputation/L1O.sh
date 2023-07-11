while getopts s:c:b:p: option
do
        case "${option}"
        in
                s) sss=${OPTARG};;
                c) chr=$OPTARG;;
                b) beagle=$OPTARG;;
                p) pop=$OPTARG;;
        esac
done

VCF=/projects/ps-gymreklab/helia/ensembl/ensemble_out/ensemble_chr"$chr"_filtered.vcf.gz ## To get initial TR call
snp=/projects/ps-gymreklab/yal084_storage/ensemble_phasing/SNP_chr"$chr".vcf.gz  ## To get SNPs for target sample
#procStr=/projects/ps-gymreklab/yal084_storage/ensemble_phasing/phased_merged_with_SNP/chr"$chr"_SNP_merged.vcf.gz ## To get nearby SNPs and TRs
procStr=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/phase_extra_repeats/concat_to_previous_files/chr"$chr"_final_SNP_merged_additional_TRs.vcf.gz
SAMPLEID=$sss
REGIONS=/projects/ps-gymreklab/helia/ensembl/1000Genomes-TR-Analysis/phasing/pathogenic_imputation/pathogenic_loci/loci_chr"$chr"_extra.txt


while read region; do
  start=$(echo $region | cut -d " " -f 2)
  end=$(echo $region | cut -d " " -f 3)
  window_start=$(($start - 50000))
  window_end=$(($end + 50000))
  REGION="chr$chr:$start-$end"
  WINDOW="chr$chr:$window_start-$window_end"
  echo $SAMPLEID
  echo $chr
  echo $WINDOW
  #extract snp IDs
  bcftools query -r $WINDOW -f '%ID\n' $snp > ID.${SAMPLEID}.${REGION}.txt

  #remove our target sample from snpstrs
  bcftools view -r $WINDOW $procStr -v snps,indels --samples ^$SAMPLEID --no-update --output-type z --output-file ref.${SAMPLEID}.${chr}.vcf.gz --force-samples

  #keep only our target sample from snps
  bcftools view $snp -r $WINDOW -v snps,indels --samples $SAMPLEID --no-update --output-type z --output-file exclude.${SAMPLEID}.${chr}.vcf.gz --force-samples
  bcftools index -f ref.${SAMPLEID}.${chr}.vcf.gz

  #impute strs for target sample
  java -Xmx4g -jar $beagle gt=exclude.${SAMPLEID}.${chr}.vcf.gz ref=ref.${SAMPLEID}.${chr}.vcf.gz out=imputed.${SAMPLEID}.${chr}
  bcftools index imputed.${SAMPLEID}.${chr}.vcf.gz

  #remove snps and extra STRs from imputed file
  bcftools view imputed.${SAMPLEID}.${chr}.vcf.gz -r $REGION --exclude ID=@ID.${SAMPLEID}.${REGION}.txt -O z -o imputed_${pop}/${chr}/imputed.str.${SAMPLEID}.${REGION}.vcf.gz


  # Get unphased STRs
  bcftools view -s $SAMPLEID -h $VCF --no-update > ensemble.${REGION}.${SAMPLEID}.vcf
  bcftools view -s $SAMPLEID -H $VCF --no-update -r $REGION | grep $start >> ensemble.${REGION}.${SAMPLEID}.vcf
  bgzip -f ensemble.${REGION}.${SAMPLEID}.vcf
  bcftools index -f ensemble.${REGION}.${SAMPLEID}.vcf.gz

  #write gts for both imputed and original strs
  python write_bases.py imputed_${pop}/${chr}/imputed.str.${SAMPLEID}.${REGION}.vcf.gz $SAMPLEID $SAMPLEID.${REGION}.imputeResult.txt
  python write_bases.py ensemble.${REGION}.${SAMPLEID}.vcf.gz $SAMPLEID $SAMPLEID.${REGION}.groundTruth.txt

  #sort and join
  sort -f ${SAMPLEID}.${REGION}.groundTruth.txt > ${SAMPLEID}.${REGION}.ground.sorted.txt
  sort -f ${SAMPLEID}.${REGION}.imputeResult.txt > ${SAMPLEID}.${REGION}.imputed.sorted.txt
  join -1 1 -2 1 ${SAMPLEID}.${REGION}.ground.sorted.txt ${SAMPLEID}.${REGION}.imputed.sorted.txt | awk 'NF==5{print}' > diff_${pop}/${chr}/${SAMPLEID}.${REGION}.diff.txt


  rm $SAMPLEID.${REGION}.groundTruth.txt $SAMPLEID.${REGION}.imputeResult.txt ${SAMPLEID}.${REGION}.ground.sorted.txt ${SAMPLEID}.${REGION}.imputed.sorted.txt ensemble.${REGION}.${SAMPLEID}.vcf.gz*

  echo "done"

  rm exclude.${SAMPLEID}.${chr}.vcf.gz* ref.${SAMPLEID}.${chr}.vcf.gz* imputed.${SAMPLEID}.${chr}.vcf.gz*
  rm ID.${SAMPLEID}.${REGION}.txt imputed.${SAMPLEID}.${chr}.log

done <$REGIONS
