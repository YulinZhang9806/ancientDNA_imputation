hg19ref=/mnt/data/Genomes/hg19_1000g/whole_genome.fa
GATK=/public/software/adna/GATK-3.5.0/GenomeAnalysisTK.jar
picard=/home/yulin_zhang/picard.jar
samtools=/public/software/adna/samtools-1.5/samtools
vcftools=/public/software/vcftools/src/cpp/vcftools
bcftools=/public/software/adna/bcftools/bcftools
GP=/mnt/data/public_data/aDNA_Imp/1kg_filtered
sample=/public/adna/yulin_zhang
beagle=/home/yulin_zhang/beagle.r1399.jar

-------------preparing_sorted_bam_files/1000GP_files/1240k_VCFfile----------------------
###sort bam
$samtools sort -@ 4 Loschbour.hg19_1000g.bam -o Loschbour_sorted.bam
###quality control--PCR duplicate removed, mapping quality >30
$samtools rmdup Loschbour_sorted.bam Loschbour_sorted_dupremoved.bam
$samtools view -@ 4 -q 30 -b Loschbour_sorted_dupremoved.bam > Loschbour_dupremoved_filtered.bam
###random subsample bam file
$samtools view -@ 4 -s 0.09 -b Loschbour_dupremoved_filtered.bam > Loschbour_dupremoved_filtered_2x.bam
###build bam index
$samtools index -b Loschbour_dupremoved_filtered_2x.bam > Loschbour_dupremoved_filtered_2x.bai


###filter 1000GP vcf, no indel, no multiallelic
$bcftools view --max-alleles 2 --exclude-types indels chr$chromo.1kg.phase3.v5a.vcf.gz --output-type v -o chr$chromo.filtered.vcf
#bgzip -c input.vcf > output.vcf.gz
#bcftools index input.vcf.gz
#bcftools merge input1.vcf.gz input2.vcf.gz input3.vcf.gz > merge.vcf.gz

###making 1240k vcf
python make_1240k_vcf.py

--------------UnifiedGenotyper--------------

###call genotype with 1000GP; 1240k by changing the --alleles to 1240k.vcf
java -jar $GATK \
-T UnifiedGenotyper \
-nct 4 \
-R $hg19ref \
-I $sample/LBK_dupremoved_filtered_1.9x.bam \
--min_base_quality_score 20 \
--output_mode EMIT_ALL_SITES \
--allSitePLs \
--alleles $GP/chr$chromo.filtered.vcf \
--genotyping_mode GENOTYPE_GIVEN_ALLELES \
-o $sample/LBK_1.9x.chr$chromo.vcf \

--------------Beagle_imputation-------------
###filter deamination/PL SNPs in vcf
python filter_deamination_vcf.py $chromosome_number

###split vcf by chromosome
$vcftools --vcf Loschbour_2x.1240k.filtered.vcf --chr Y --recode --recode-INFO-all --out Loschbour_2x.1240k.chrY

###split vcf for beagle
cat /public/adna/yulin_zhang/Loschbour_2x.chr$chromo.filtered.vcf | java -jar /home/yulin_zhang/splitvcf.jar $chromo 50000 25000 Loschbour_beagle_chr$chromo

###beagle
java -jar $beagle gprobs=true impute=true gl=$sample/LBK_beagle/$input ref=$GP/chr1.filtered.vcf map=$sample/plink_GRCh37_map/plink.chr1.GRCh37.map

###merge beagle output files by chromosome
java -jar /home/yulin_zhang/mergevcf.jar $chromo Loschbour_impute_0.5x.1240k.chr$chromo.*  > /public/adna/yulin_zhang/Loschbour_beagle/Loschbour_impute_0.5x_merge.1240k.chr$chromo.vcf

###merge beagle imputed chromosomes to one file
python merge_BeagleVCF.py

-----------evaluation-------------------

###filter genotype vcf of high coverage sequencing samples
cd /mnt/data/public_data/aDNA_Imp/LBK
$vcftools --vcf /mnt/data/public_data/aDNA_Imp/LBK/LBK.chr$chromo.vcf --minQ 20 --min-meanDP 15 --recode --recode-INFO-all --out LBK.chr$chromo.filtered

###filter vcf files with GP thresholds
python threshold_imputation_vcf.py

###compare high coverage file with imputed files
java -jar /home/yulin_zhang/SnpSift.jar concordance -v /mnt/data/public_data/aDNA_Imp/LBK/LBK.chr1.vcf LBK_impute_merge_chr1.vcf

###making figures for accuracy evaluation
python 