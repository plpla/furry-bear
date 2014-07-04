#!/bin/bash

#We need java and samtools
module load nne-790-ab/samtools/0.1.18-1
module load java/jdk1.7.0

PATH_TO_PICARD=/rap/nne-790-ab/projects/pplante/Metagenome_SNP_detection_pipeline/Software/picard-tools-1.114/
PATH_TO_GATK=/rap/nne-790-ab/projects/pplante/Metagenome_SNP_detection_pipeline/Software

#Use absolute paths
referenceSequence=$1
#reference must be one one line or all line must have the same length otherwise segfault!
bamFilesDirectory=$2
#Outpout dir name. It will be created and it will be used in files name
outputDir=$3



#Create directory for input Data
mkdir $outputDir
mkdir $outputDir/Reference
mkdir $outputDir/BamFiles


#Place symbolic links to data in directory
ln -s $referenceSequence $outputDir/Reference/.
ln -s $bamFilesDirectory/*.bam $outputDir/BamFiles/.

#Change CWD
cd $outputDir

#Prepare reference dictionary and index
referenceFile=$(basename $referenceSequence)
cd Reference/
java -jar $PATH_TO_PICARD/CreateSequenceDictionary.jar R=$(basename $referenceSequence) O=$(echo $referenceFile|cut -f1 -d'.').dict
samtools faidx $referenceFile 2>&1 | tee -a ${outputDir}.log

#Prepare bam files
cd ../BamFiles

#SoftClip to make sure eveything is OK (added 2013-11-12 because of an error... should not effect results).
java -jar -Xmx4g $PATH_TO_PICARD/CleanSam.jar  INPUT=$(ls *R1*.bam) OUTPUT=${outputDir}_1.bam 2>&1 | tee -a ${outputDir}.log
java -jar -Xmx4g $PATH_TO_PICARD/CleanSam.jar  INPUT=$(ls *R2*.bam) OUTPUT=${outputDir}_2.bam 2>&1 | tee -a ${outputDir}.log

#Merge
java -jar -Xmx4g $PATH_TO_PICARD/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT INPUT=${outputDir}_1.bam INPUT=${outputDir}_2.bam OUTPUT=$(echo $outputDir)_merged.bam 2>&1 | tee -a ${outputDir}.log

#Add header group
java -jar $PATH_TO_PICARD/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT INPUT=${outputDir}_merged.bam OUTPUT=${outputDir}_merged_Grouped.bam  \
 LB=$outputDir PL=illumina RGPU=$outputDir SM=$outputDir 2>&1 | tee -a ${outputDir}.log

#Index
samtools index ${outputDir}_merged_Grouped.bam 2>&1 | tee -a ${outputDir}.log

cd ..

#Run GATK!
java -jar -Xmx5g  \
 $PATH_TO_GATK -T DepthOfCoverage -baseCounts -R Reference/$referenceFile -nt 8 -rf BadCigar\
 -I BamFiles/${outputDir}_merged_Grouped.bam -o ${outputDir}.gatk 2>&1 | tee ${outputDir}.log

#Run the conversion from GATK output to my usual tab format
python /rap/nne-790-ab/projects/pplante/GATK_PICARD_PIPELINE/GATKToTab.py Reference/$referenceFile ${outputDir}.gatk > ${outputDir}.tab

#Run the SNPs research script:
python /rap/nne-790-ab/projects/pplante/GATK_PICARD_PIPELINE/GATKSearchSNPS.py Reference/$referenceFile ${outputDir}.gatk > ${outputDir}-LowVarSNPs.tab

#Run SNPs calling: is there a syn or non-syn mutation?
python /rap/nne-790-ab/projects/pplante/GATK_PICARD_PIPELINE/RemoveSynonymousMutation.py Reference/$referenceFile ${outputDir}-SNPS.tab >${outputDir}-SNPs_Type
