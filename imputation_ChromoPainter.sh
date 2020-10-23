hg19ref=/mnt/data/Genomes/hg19_1000g/whole_genome.fa
GATK=/public/software/adna/GATK-3.5.0/GenomeAnalysisTK.jar
picard=/home/yulin_zhang/picard.jar
samtools=/public/software/adna/samtools-1.5/samtools
vcftools=/public/software/vcftools/src/cpp/vcftools
bcftools=/public/software/adna/bcftools/bcftools
GP=/mnt/data/public_data/aDNA_Imp/1kg_filtered
sample=/public/adna/yulin_zhang
beagle=/home/yulin_zhang/beagle.r1399.jar
dbsnp=/public/adna/yulin_zhang/dbsnp_138.hg19.vcf
Mills=/public/adna/yulin_zhang/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
GP_phase1=/public/adna/yulin_zhang/1000G_phase1.indels.hg19.sites.vcf
plink=/public/software/adna/plink_1.9-linux_x86_64/plink
smartpca=/public/software/adna/EIG/bin/smartpca
root=NE1
root_1240k=v42.4.1240K

-----------------------preparing_sorted_bam_files---------------------------
###sort bam
$samtools sort -@ 4 $root.bam -o $root.sorted.bam
###quality control--PCR duplicate removed, mapping quality >30
$samtools view -@ 4 -q 30 -b $root.sorted.bam > $root.sorted.filteredq.bam
$samtools rmdup $root.sorted.filteredq.bam > $root.sorted.filteredq.rmdup.bam
###realign around indels
java -jar $GATK \
-nct 4 \
-R $hg19ref \
-T RealignerTargetCreator \
-I $root.sorted.filteredq.rmdup.bam \
-o $root.sorted.filteredq.rmdup.realn.intervals \
-known $Mills \
-known $GP_phase1 
java -jar $GATK \
-R $hg19ref \
-nct 4 \
-T IndelRealigner \
-targetIntervals $root.sorted.filteredq.rmdup.realn.intervals \
-I $root.sorted.filteredq.rmdup.bam \
-o $root.sorted.filteredq.rmdup.realn.bam \
-known $Mills \
-known $GP_phase1 
###BQSR
java -jar GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R $hg19ref \
-I $root.sorted.filteredq.rmdup.realn.bam \
-knownSites $dbsnp \
-knownSites $Mills \
-knownSites $GP_phase1 \
-o $root.sorted.filteredq.rmdup.realn.grp
java -jar GenomeAnalysisTK.jar \
-T PrintReads \
-R $hg19ref \
--disable_indel_quals \
--preserve_qscores_less_than 6 \
-I $root.sorted.filteredq.rmdup.realn.bam \
-BQSR $root.sorted.filteredq.rmdup.realn.grp \
-o $root.sorted.filteredq.rmdup.realn.rebq.bam \

#(note here: sort & index should use the same program, e.g. samtools or Picardtools)
###random subsample bam file (22x)
$samtools view -@ 4 -s 0.04545 -b $root.sorted.filteredq.rmdup.realn.rebq.bam > ${root}_1x.bam
###build bam index
$samtools index -b ${root}_1x.bam

-----------------------preparing_1000GP_reference&1240k_vcf---------------------------
###filter 1000GP vcf, no indel, no multiallelic, maf>0.01 (use vcftools, not bcftools, since bcftools can`t remove indel/multiallelic completely)
$vcftools --gzvcf chr$chrom.1kg.phase3.v5a.vcf.gz --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --recode --recode-INFO-all --out chr$chrom.filtered.vcf
#$bcftools view --max-alleles 2 --min-alleles 2 --exclude-types indels chr$chromo.1kg.phase3.v5a.vcf.gz | $bcftools norm -d none --output-type v -o chr$chromo.filtered.vcf
###making 1240k vcf
python eigenstrat2vcf.py -r $root_1240k

-------------------------------UnifiedGenotyper----------------------------
###call genotype with 1000GP; 1240k by changing the --alleles to 1240k.vcf
java -jar $GATK \
-T UnifiedGenotyper \
-nct 4 \
-R $hg19ref \
-I ${root}_1x.bam \
-L chr$chromo.filtered.vcf \
--output_mode EMIT_ALL_SITES \
--allSitePLs \
--alleles chr$chromo.filtered.vcf \
--genotyping_mode GENOTYPE_GIVEN_ALLELES \
-o ${root}_1x.vcf \

------------------------------Beagle_imputation-----------------------------
###keep only non-missing PL SNPs, filter deamination transversion to missing in vcf
python filter_deamination_vcf.py $chromosome_number
###split input vcf by chromosome
$vcftools --vcf ${root}_1x.1240k.vcf --chr $chrom --recode --recode-INFO-all --out ${root}_1x.1240k.chr$chrom.vcf
###split input vcf for quicker imputation
cat ${root}_1x.1240k.chr$chrom.vcf | java -jar splitvcf.jar $chrom 5000 3000 ${root}_1x.1240k.chr$chrom.beagle
###beagle
start=$(bcftools query -f '%CHROM\t%POS\n' ${root}_1x.1240k.chr$chrom.beagle.$input.vcf.gz |head -1|cut -f2)
end=$(bcftools query -f '%CHROM\t%POS\n' ${root}_1x.1240k.chr$chrom.beagle.$input.vcf.gz |tail -1|cut -f2)
#note here if imputing whole chromosome, the first input file should have chrom=[chr#chrom:-end] and the last should be chrom=[chr#chrom:start-]
java -Xmx120g -jar $beagle gprobs=true impute=true gl=${root}_1x.1240k.chr$chrom.beagle.$input.vcf.gz ref=chr$chrom.filtered.vcf map=plink.chr$chrom.GRCh37.map excludesamples=non.eur.excl excludemarkers=non1240k.excl chrom=chr$chrom:$start-$end out=${root}_1x.1240k.chr$chrom.impute.$input
###merge beagle output files by chromosome
java -jar /home/yulin_zhang/mergevcf.jar $chrom ${root}_1x.1240k.chr$chrom.impute.$input.*  > ${root}_1x.1240k.chr$chrom.impute.vcf.gz
###merge beagle imputed chromosomes to one file
#python merge_BeagleVCF.py
$bcftools index ${root}_1x.1240k.chr$chrom.impute.vcf.gz
$bcftools concat ${root}_1x.1240k.chr$chrom.impute.vcf.gz ${root}_1x.1240k.chr$chrom.impute.vcf.gz -Oz -o ${root}_1x.1240k.impute.vcf.gz

---------------------------GLIMPSE_imputation--------------------------
#chunking reference chromosome into pieces
bgzip chr$chrom.filtered.vcf
$bcftools index chr$chrom.filtered.vcf.gz
/home/yulin_zhang/GLIMPSE_chunk --input chr$chrom.filtered.vcf.gz --region $chrom --window-size 10000000 --buffer-size 1000000 --output chunks.$chrom.txt
#imputation
bgzip ${root}_1x.1240k.chr$chrom.vcf
$bcftools index ${root}_1x.1240k.chr$chrom.vcf.gz
start=$(less chunks.$chrom.txt|cut -f 1)
end=$(less chunks.$chrom.txt|cut -f 2)
/home/yulin_zhang/GLIMPSE_phase --input ${root}_1x.1240k.chr$chrom.vcf.gz --reference chr$chrom.filtered.vcf.gz --map /GLIMPSE/maps/genetic_maps.b37/$chrom.b37.gmap.gz --input-region $chrom:$start-$end --output-region chr1:10642-249240539 --output ${root}_1x.1240k.impute.$input.chr$chrom.vcf --thread 20
bgzip ${root}_1x.1240k.impute.$input.chr$chrom.vcf
$bcftools index ${root}_1x.1240k.impute.$input.chr$chrom.vcf.gz
#ligate
ls ${root}_1x.1240k.impute.$input.chr$chrom.vcf > ligate.$chrom.txt
/home/yulin_zhang/GLIMPSE_ligate --input ligate.$chrom.txt --output ${root}_1x.1240k.impute.chr$chrom.vcf
#phasing
/home/yulin_zhang/GLIMPSE_sample --input ${root}_1x.1240k.impute.$input.chr$chrom.vcf.gz --solve --output ${root}_1x.1240k.impute.$input.chr$chrom.phased.vcf

-----------------------------ChromoPainter-----------------------------
###filter imputed files by MAF>0.01, GP>0.99
python filter_GP&AF.py 
###merge vcf of different samples (for flist, one indexed .vcf.gz file per line)
$bcftools merge --file-list flist -Ov -o Worldpop_capture.vcf
###filter snps with missing genotypes,indel,multiallelic,maf<0.01
$vcftools --vcf Worldpop_capture.vcf --max-missing-count 0 --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01 --recode --recode-INFO-all --out Worldpop_capture.filsnp.vcf
###split by chromosome
$vcftools --vcf Worldpop_capture.filsnp.vcf --chr $chrom --recode --recode-INFO-all --out Worldpop_capture.filsnp.$chrom.vcf
###convert to .ped .map of plink
$plink --vcf Worldpop_capture.filsnp.$chrom.vcf --recode12 --out Worldpop_capture.filsnp.$chrom
###convert to .haps of chromopainter
./plink2chromopainter.pl -p Worldpop_capture.filsnp.$chrom.ped -m  Worldpop_capture.filsnp.$chrom.map -o  Worldpop_capture.filsnp.$chrom.phase
###make recombination file
./convertrecfile.pl -M hap Worldpop_capture.filsnp.$chrom.phase genetic_map_GRCh37_$chrom.txt Worldpop_capture.filsnp.$chrom.recombfile
###run fineStructure
./fs Worldpop_capture_$chrom.cp -phasefiles Worldpop_capture.filsnp.$chrom.phase -recombfiles Worldpop_capture.filsnp.$chrom.recombfile -idfile Worldpop_capture.id -s1minsnps 5000 -s3iters 10000 -s4iters 10000 -go
