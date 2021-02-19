#!/bin/sh

#  DNA_Mutation_Pipeline.sh
#  
#
#  Created by xyl on 2020/2/18.
#


############################################################## Set Directory and Data Path ##########################################################################

tumorfastq1=/home/xueyinglyu/project/renbo/rawdata/UM_45_FKDN202373527-1A_H5JFJDSXY_L4_1.fq.gz
tumorfastq2=/home/xueyinglyu/project/renbo/rawdata/UM_45_FKDN202373527-1A_H5JFJDSXY_L4_2.fq.gz
bloodfastq1=/home/xueyinglyu/project/renbo/rawdata/UM_44_FKDN202373526-1A_H5JFJDSXY_L2_1.fq.gz
bloodfastq2=/home/xueyinglyu/project/renbo/rawdata/UM_44_FKDN202373526-1A_H5JFJDSXY_L2_2.fq.gz
output_dir=/home/xueyinglyu/project/renbo/DNA/P15
sample=P15


DBSNP=/home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/dbsnp_146.hg38.vcf.gz
GENOME=/home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta
kgINDEL=/home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
kgSNP=/home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
bed=/home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/hg38_targets.bed


############################################################## Process Data ##########################################################################

fastp -i $tumorfastq1 -I $tumorfastq2 -o $outputdir/$sample_tumor_R1.fastq.gz -O $outputdir/$sample_tumor_R2.fastq.gz -q 20 -l 50 -g -w 8 --json $sample_tumor.json --html $sample_tumor.html
fastp -i $bloodfastq1 -I $bloodfastq2 -o $outputdir/$sample_blood_R1.fastq.gz -O $outputdir/$sample_blood_R2.fastq.gz -q 20 -l 50 -g -w 8 --json $sample_blood.json --html $sample_blood.html


bwa mem -t 16 -M -R '@RG\tID:Tumor\tSM:Tumor\tLB:WES\tPL:Illumina' /home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta $outputdir/$sample_tumor_R1.fastq.gz $outputdir/$sample_tumor_R2.fastq.gz 1 > $outputdir/$sample_tumor.sam
bwa mem -t 16 -M -R '@RG\tID:Tumor\tSM:Tumor\tLB:WES\tPL:Illumina' /home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta $outputdir/$sample_blood_R1.fastq.gz $outputdir/$sample_blood_R2.fastq.gz 1 > $outputdir/$sample_blood.sam


gatk --java-options "-Xmx25G -Djava.io.tmpdir=./" SortSam -SO coordinate -I $outputdir/$sample_tumor.sam -O $outputdir/$sample_tumor.bam
gatk --java-options "-Xmx25G -Djava.io.tmpdir=./" SortSam -SO coordinate -I $outputdir/$sample_blood.sam -O $outputdir/$sample_blood.bam

samtools sort -@ 16 $outputdir/$sample_tumor.bam -o $outputdir/$sample_tumor.sort.bam
samtools sort -@ 16 $outputdir/$sample_blood.bam -o $outputdir/$sample_blood.sort.bam

sambamba markdup --overflow-list-size 600000 --tmpdir='./' -r $outputdir/$sample_tumor.sort.bam $outputdir/$sample_tumor.dedup.bam
sambamba markdup --overflow-list-size 600000 --tmpdir='./' -r $outputdir/$sample_blood.sort.bam $outputdir/$sample_blood.dedup.bam

gatk --java-options "-Xmx25G -Djava.io.tmpdir=./" FixMateInformation -I $outputdir/$sample_tumor.dedup.bam -O $outputdir/$sample_tumor.fixed.bam -SO coordinate
gatk --java-options "-Xmx25G -Djava.io.tmpdir=./" FixMateInformation -I $outputdir/$sample_blood.dedup.bam -O $outputdir/$sample_blood.fixed.bam -SO coordinate

gatk --java-options -Xmx25G BaseRecalibrator -I $outputdir/$sample_tumor.fixed.bam -R $GENOME --output $outputdir/$sample_tumor.table --known-sites $kgSNP --known-sites $kgINDEL
gatk --java-options -Xmx25G BaseRecalibrator -I $outputdir/$sample_blood.fixed.bam -R $GENOME --output $outputdir/$sample_blood.table --known-sites $kgSNP --known-sites $kgINDEL

gatk --java-options -Xmx25G ApplyBQSR -I $outputdir/$sample_tumor.fixed.bam -R $GENOME --output $outputdir/$sample_tumor.recal.bam -bqsr $outputdir/$sample_tumor.table
gatk --java-options -Xmx25G ApplyBQSR -I $outputdir/$sample_blood.fixed.bam -R $GENOME --output $outputdir/$sample_blood.recal.bam -bqsr $outputdir/$sample_blood.table

############################################################## Three Tools are applied to call mutations ##########################################################################

#################################### GATK Mutect2 ##################################
gatk Mutect2 -R $GENOME -I $outputdir/$sample_tumor.recal.bam -I $outputdir/$sample_blood.recal.bam -tumor $sample_tumor -O $outputdir/$sample.mutect2.vcf.gz
gatk FilterMutectCalls -R $GENOME -V $outputdir/$sample.mutect2.vcf.gz -O $outputdir/$sample.mutect2.filter.vcf.gz --filtering-stats $outputdir/$sample.mutect2.stats
zcat $outputdir/$sample.mutect2.filter.vcf.gz | grep -E '^#|PASS' > $outputdir/$sample.mutect2.filtered.vcf

#################################### Strelka2 ####################################
configureStrelkaSomaticWorkflow.py --normalBam $outputdir/$sample_blood.recal.bam --tumorBam $outputdir/$sample_tumor.recal.bam --referenceFasta $GENOME --runDir $outputdir --exome --reportEVSFeatures --outputCallableRegions --callRegions $bed
runWorkflow.py -m local -j 8
cat <(zcat $outputdir/results/variants/somatic.snvs.vcf.gz | grep '#') <(zcat $outputdir/results/variants/somatic.snvs.vcf.gz | grep 'PASS') <(zcat $outputdir/results/variants/somatic.indels.vcf.gz | grep 'PASS')  > $outputdir/$sample.strelka2.filtered.vcf.gz

#################################### Lancet ####################################
lancet --tumor $outputdir/$sample_tumor.recal.bam --normal $outputdir/$sample_blood.recal.bam --ref $GENOME --bed $bed --num-threads 32


############################################################## Mutation Annotation ##########################################################################

perl vcf2maf.pl --input-vcf $outputdir/$sample.mutect2.filtered.vcf --output-maf $outputdir/$sample.mutect2.maf --tumor-id $sample_tumor --normal-id $sample_blood --ref-fasta /home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta --vep-forks 16 --vep-data /home/xueyinglyu/database/ --vep-path /home/xueyinglyu/anaconda2/envs/vep/bin/ --ncbi-build GRCh38 --filter-vcf /home/xueyinglyu/database/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
perl vcf2maf.pl --input-vcf $outputdir/$sample.strelka2.filtered.vcf --output-maf $outputdir/$sample.strelka2.maf --tumor-id TUMOR --normal-id NORMAL --ref-fasta /home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta --vep-forks 16 --vep-data /home/xueyinglyu/database/ --vep-path /home/xueyinglyu/anaconda2/envs/vep/bin/ --ncbi-build GRCh38 --filter-vcf /home/xueyinglyu/database/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
perl vcf2maf.pl --input-vcf $outputdir/$sample.lancet.filtered.vcf --output-maf $outputdir/$sample.lancet.maf --tumor-id $sample_tumor --normal-id $sample_normal --ref-fasta /home/xueyinglyu/biosoft/GATK/resources/bundle/hg38/Homo_sapiens_assembly38.fasta --vep-forks 16 --vep-data /home/xueyinglyu/database/ --vep-path /home/xueyinglyu/anaconda2/envs/vep/bin/ --ncbi-build GRCh38 --filter-vcf /home/xueyinglyu/database/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz


############################################################## Mutation Filter (https://github.com/mskcc/ngs-filters) ##########################################################################

/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_common_variants.R $outputdir/$sample.mutect2.maf $outputdir/$sample.mutect2.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_low_conf.R $outputdir/$sample.mutect2.tmp.maf $outputdir/$sample.mutect2.tmp2.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_blacklist_regions.R $outputdir/$sample.mutect2.tmp2.maf $outputdir/$sample.mutect2.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/tag_filters.py -itxt /home/xueyinglyu/biosoft/ngs-filters-1.4/data/hotspot-list-union-v1-v2.txt -m $outputdir/$sample.mutect2.tmp.maf -o $outputdir/$sample.mutect2.tmp2.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_vaf.R $outputdir/$sample.mutect2.tmp2.maf $outputdir/$sample.mutect2.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_ffpe.R $outputdir/$sample.mutect2.tmp.maf $outputdir/$sample.mutect2.tmp2.maf


/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_common_variants.R $outputdir/$sample.strelka2.maf $outputdir/$sample.strelka2.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_low_conf.R $outputdir/$sample.strelka2.tmp.maf $outputdir/$sample.strelka2.tmp2.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_blacklist_regions.R $outputdir/$sample.strelka2.tmp2.maf $outputdir/$sample.strelka2.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/tag_filters.py -itxt /home/xueyinglyu/biosoft/ngs-filters-1.4/data/hotspot-list-union-v1-v2.txt -m $outputdir/$sample.strelka2.tmp.maf -o $outputdir/$sample.strelka2.tmp2.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_vaf.R $outputdir/$sample.strelka2.tmp2.maf $outputdir/$sample.strelka2.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_ffpe.R $outputdir/$sample.strelka2.tmp.maf $outputdir/$sample.strelka2.tmp2.maf


/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_common_variants.R $outputdir/$sample.lancet.maf $outputdir/$sample.lancet.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_low_conf.R $outputdir/$sample.lancet.tmp.maf $outputdir/$sample.lancet.tmp2.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_blacklist_regions.R $outputdir/$sample.lancet.tmp2.maf $outputdir/$sample.lancet.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/tag_filters.py -itxt /home/xueyinglyu/biosoft/ngs-filters-1.4/data/hotspot-list-union-v1-v2.txt -m $outputdir/$sample.lancet.tmp.maf -o $outputdir/$sample.lancet.tmp2.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_vaf.R $outputdir/$sample.lancet.tmp2.maf $outputdir/$sample.lancet.tmp.maf
/home/xueyinglyu/biosoft/ngs-filters-1.4/applyFilter.sh /home/xueyinglyu/biosoft/ngs-filters-1.4/filter_ffpe.R $outputdir/$sample.lancet.tmp.maf $outputdir/$sample.lancet.tmp2.maf


