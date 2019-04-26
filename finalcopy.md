#A pipeline for identifying potential polymorphisms associated with a change in wing shape in following weak or strong selection. 
Katie Pelletier 
Biol722 Project 
April 26, 2019 

##Introduction

***

Many traits of interest to biologists, including many human diseases, are highly polygenic with alleles of small effect size contributing. However, many studies have focused on traits with a monogenic basis and alleles of large effect sizes (Flint and MacKay, 2009). Understanding how changes in genetic architecture relate to phenotypic changes for complex traits is important for understanding not only disease but also the evolution of complex traits. By studying the genetic architecture underlying phenotypic change in natural populations and changes to this architecture following artificial selection, we can broaden our understanding of the relationship between genotype and phenotype for complex traits. 

A common model to study the interaction between genetic architecture and complex traits is Drosophila wing shape. However, the observed variation across many Drosophild species is would predict a much faster rate of evolution than what is observed for these species (Houle *et al* 2017). This indicates that there is a much more complex relationship between
genetic architecture and wing shape phenotypes and a deeper understanding of the number and types of variants contributing to wing shape change is needed. Pitches *et al* performed a genome-wide association study to identify a number of variants contributing to changes in wing shape in D. melanogaster (2019). This study identified over 2000 polymorphisms implicated in shape of the wing using the Drosophila genetic resource panel (DGRP) (MacKay *et al* 2012). These lines were isolated from wild populations but have undergone many generations of rearing in the lab and thus, have accumulated lab adaptations and are not truly representative of natural populations. Because many of the genes associated with a change in wing shape are involved in developmental pathways, variant selection could be associated with this lab adaptation and not applicable in a natural context. By sampling natural populations, we are able to ask if the same variants, or variants in the same genes, are segregating in these populations, indicating that these variants are not just a by-product of lab adaptation. Because wing shape is a highly dimensional and complex trait, it seems unlikely that variants in a single gene can result in a single unique phenotype. This is further supported by the complexity of the wing developmental regulatory network, which incorporates many signaling pathways and crosstalk. By identifying the polymorphisms associated with selection of variant phenotypes, we ask if polymorphisms in many genes can be associated with a particular phenotype or if this is restricted to polymorphisms in a single gene. 

Of the identified polymorphisms, three variant phenotypes were selected for further analysis: *dachous* (*ds*), *neutralized* (*neur*) and *extra macrochaetae* (*emc*). All three genes are known members of the *Drosophila* wing developmental gene network (Pitchers *et al* 2019). To obtain a reference phenotype, an RNAi knockdown was performed for gene and this phenotype was used for comparisons for the selection experiments. In other words, we are defining the gene of interest phenotype as that of a weak knockdown of that gene. 

To ask about the association between polymorphisms in a natural context and the identified DGRP polymorphisms, flies were collected from three populations in northern Michigan: Fenn Valley winery (FVW, sampled in 2012 and 2013), Phillip Hills Orchard (PHO, 2013) and Country Mills Orchard (CMO, 2013). A shape score was calculated for each individual that described the similarity (or lack thereof) in wing shape to either that of the *ds* or *neur* knockdown phenotype. These scores were used to select the 75 most extreme individuals, that is either most similar to or least similar to, relative to the phenotype of interest. Least similar pools are referred to as the left-selected while most similar individuals are right-selected, because of their distribution when all shape scores are plotted. These pools of individuals were sequenced to identify natural polymorphisms. This design is equivalent to a single generation of selection, and thus we do not expect to see large changes in allele frequency between the right and left selected pools. Thus, in this data we will have to detect a low signal of ‘true’ variants relative to the noise from other variants in the population.  

To ask about the alleles under selection, a pool of genetic variants that could be selected upon were created by mixing 10 DGRP lines and selecting the animals that were most similar to (up selection) and least similar to (down selection) the knockdown phenotype of the gene of
interest for ~15 generations. A third population was also maintained without selection (control) to control for lab adaptations. Approximately 75 individuals from each pool were randomly selected for sequencing. Polymorphisms in this population will be analyzed to understand the number and types of alleles under selection for a specific phenotype. Because of the selection on these populations, we expect a high amount of noise in the data due to alleles captured by
drift. However, the signal of the alleles we are interested in will also be higher, compared to the wild-caught population experiment. In this case, we will have to identify polymorphisms with a high signal/high noise ratio as opposed to the wild-caught experiment where we will identify polymorphisms with a low signal/low noise ratio. 

Pool sequencing (pool-seq) is a whole genome sequencing approach that involves creating pools of individuals within a single library preparation (Schlotter *et al* 2014). This provides a low-cost option to get genomic information of sequence diversity within a population. However, there is a cost in that you lose information about individual genotypes and haplotype
information. Not all software designed for genome analysis performs well with pool-seq data because of the number of genomes present in each sample. One notable example is gatk, a commonly used variant caller for other sequencing methods, which have a high false positive rate for samples with pools less than 16 diploid genomes and fails to run for samples with greater than 16 diploid genomes (Huang *et al* 2015). Because of this, a number of tools specific to pool-seq data have been developed. Three specific tools I will be using in this project are PoPoolation, PoPoolation2 and CRISP (Reviewed in Schlotter *et al* 2014). PoPoolation is a tool designed to make inferences about genomic data within a single pooled sample while PoPoolation2 is designed to compare genomic information between two different pooled samples (Kofler *et al* 2011a and b). The final tool, CRISP, is a variant caller specifically designed to identify variants present in multiple pool-seq samples (Bansal et al 2011). 

Using the genome analysis pipeline used by other students in the Dworkin lab (modified from Paul Knoops and Sarah Marzec) with some modifications, I intend to identify regions of the genome associated with a change in wing phenotype in either wild-caught or artificially selected populations and begin to identify possible polymorphisms for future analysis. With this
project, I will also investigate the inclusion of an indel realignment step in the analysis, using the gatk indel realigner to make inferences about the impact such steps will have on our results. 


***

###1. Sequencing and Populations.

All natural population capture and selection experiments were performed by in the Dworkin lab and artificial selection experiments were performed in the Houle lab. 
Sequencing of all samples was performed at the All samples were prepared with the Rubicon Genomics ThruPLEX DNA Library Preparation kit which does not have a automatic size selection procedure. This may result in some shorter sequences in these samples.

###2. Initial Quality Control 

The files were transferred to the McMaster/Golding lab server from Ian’s storage space at MSU. The first step following transfer is to ensure that this happened correctly using an md5sum check.

``` 
md5sum - c md5.txt 
``` 
Flags: "-c" = report if checksums match contents of files.

All files were successfully transferred, so I can begin my analysis.

I next want to check the overall quality of these samples with FastQC. This program will identify problems with the quality of the reads that could have occurred during sequencing (Babraham Bioinformatics, 2015). FastQC is run with the following script.

I ran this within each directory (each directory corresponds to sequences from a single run). Example for a single directory shown below. 
``` 
fastqc -o /home/katie/dachsProject/WingShapeBSA/20141205_A_DNASeq_PE /home/katie/dachsProject/WingShapeBSA/beforeTrim/20141205_A_DNASeq_PE
``` 
Flags: "-o" = specify an output directory.

fastqc.html flies can then be pulled down to my local machine to view. The quality of these reads is generally ok, although there are some low quality bases at the ends of reads and an over representation of adapter sequences, as we would expect because we have not yet trimmed the adapters out of these reads.

``` 
scp -r katie@info.mcmaster.ca:/home/katie/dachsProject/WingShapeBSA/beforeTrim/ ~/Documents/phd/WingShapeBSA/fastQC/ 
``` 
Flags: "-r" = copies directories recursively

###3. Trimming of reads 

Trimmomatic is used to remove adapter sequences from reads and is called with the script below (Bolger *et al.* 2014).

```
#!/bin/bash

#create a variable for the run. This needs to be modified each time. 
runvar=20141205_A_DNASeq_PE

#raw reads
raw_dir=/home/katie/dachsProject/WingShapeBSA/${run_var}

#Trimmomatic
trim=/usr/local/trimmomatic/trimmomatic-0.36.jar

#Adapter file, used adapters specified in README file with sequencing 
adapt=/home/katie/adapterSequences/TruSeq3-PE.fa


trim_dir=/2/scratch/Katie/WingshapeBSA/trimmomatic/{runvar}

files=(${raw_dir}/*_R1_001.fastq.gz) 

for file in ${files[@]} 
do
name=${file} 
base=`basename ${name} _R1_001.fastq.gz` 
java -jar ${trim} PE -threads 5 -phred33 \ 
		-trimlog ${trim_dir}/trimlog.txt \
		${raw_dir}/${base}_R1_001.fastq.gz \ 
		${raw_dir}/${base}_R2_001.fastq.gz \ 
		${trim_dir}/${base}_R1_PE.fastq.gz \
		${trim_dir}/${base}_R1_SE.fastq.gz \ 
		${trim_dir}/${base}_R2_PE.fastq.gz \ 
		${trim_dir}/${base}_R2_SE.fastq.gz \ 
		ILLUMINACLIP:${adapt} LEADING:3 TRAILING:3 MAXINFO:40:0.5 MINLEN:36
done 
``` 
Flags: "-threads 5" = uses 5 threads. "-phred33" = Illumina sequencing. "trimlog" = path and name of log file. "ILLUMINACLIP" = adapter sequence file "LEADING" = Removal at beining of read, for low quality. "TRAILING" = Removal at end of read, for low quality. "MAXINFO" = balances read length and base quality to maximize information provided, intermediate values chosen. "MINLEN" = Smallest possible window, smaller reads after trimming are not kept.

Following trimming, FastQC is again run on these files and pulled down to my local machine to ensure that this step was successful. There are no longer adapter sequences detected in these files and the low-quality bases have been removed. In some files there is an over representation of some highly repetitive sequences that is detected. This may result in issues later in analysis.

These are modified from the scripts above.

``` 
fastqc -o /2/scratch/Katie/WingshapeBSA/trimmomatic/20141205_A_DNASeq_PE /home/katie/dachsProject/WingShapeBSA/afterTrim/20141205_A_DNASeq_PE 
```
``` 
scp -r katie@info.mcmaster.ca:/home/katie/dachsProject/WingShapeBSA/afterTrim/ ~/Documents/phd/WingShapeBSA/fastQC/
```

###4. Mapping 

Before I can map my reads, I need to download the genome. I chose to use the most recent Drosophila genome (v6.23). When using this genome, it is important to note that all contigs which are not assigned to a particular chromosome are listed separately, instead of being concatenated together into a single ‘super contig’ as was done in previous releases of the genome.
```
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.23_FB2018_04/fasta/dmel-all-chromosome-r6.23.fasta.gz
```
I will be using BWA (Burrows-Wheeler Alignment) to map reads for this pipeline (Li and Durbin 2010). 
Before I can map my reads, I have to index the genome so it is accessible by bwa.
```
bwa index dmel-all-chromosome-r6.23.fasta.gz
```
I am now able to begin mapping of my reads.

```
#!/bin/bash
#
#project directory (to keep each run separate)
run=20141205_A_DNASeq_PE

# directory of processed sequences with trimmomatic
trim_dir=/2/scratch/Katie/WingShapeBSA/trimmomatic/${run}

# variable for the reference genome
refGenome=/home/katie/flyGenome/dmel_r6.23/bwa/dmel-all-chromosome-r6.23
.fasta

# make output directory from mapping outputs
output=/2/scratch/Katie/WingShapeBSA/mapped/bwa/SAM_files/${run}

# make BWA directory path
bwa_dir=/usr/local/bwa/0.7.8

cd ${bwa_dir}

#list all files to be read (this selects the left end from each PE pair)
files=(${trim_dir}/*_L001_R1_PE.fastq.gz)

for 
file in ${files[@]} 
do 
name=${file} 
base=`basename ${name} _L001_R1_PE.fastq.gz`

bwa mem -t 8 -M ${refGenome} \ 
		${trim_dir}/${base}_L001_R1_PE.fastq.gz \
		${trim_dir}/${base}_L001_R2_PE.fastq.gz > ${output}/${base}_bwa_PE.SAM
done
 ``` 
Flags: "-t 8" = uses 8 threads "-M" = mark shorter, spit reads

###5. Merge reads 

In these experiments, each sample was run twice to increase read depth. I now want to merge these two groups. This is done with Samtools and the script below.

```
#!/bin/bash Only works if all lanes are L001/L002

#Each sample was run on 2 seperate runs, there is a mathing sequencing
#file in each directory (identical names) Within each directory, each
#library was run on 2 lanes marked as 001 and 002
run1=20141205_B_DNASeq_PE 
run2=20150107_B_DNASeq_PE

sam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/SAM_files
merged_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/SAM_files/merge

files=(${bam_dir}/${run1}/*_L001_bwa_PE)

for 
file in ${files[@]} 
do 
name=${file} 
base=`basename ${name} _L001_bwa_PE` 
samtools merge ${merged_dir}/${base}_merged_aligned_PE.bam \ 
	${bam_dir}/${run1}/${base}_L001_bwa_PE \
	${bam_dir}/${run1}/${base}_L002_$_L002_bwa_PE \
	${bam_dir}/${run2}/${base}_L001_bwa_PE \
	${bam_dir}/${run2}/${base}_L002_bwa_PE 
done 
```

Although I did not include it in this pipeline, in the future I would add read group information to the files before they are merged.

###6. Convert SAM to BAM files 

I can now convert the SAM files into the smaller BAM format, also using Samtools (Li et al. 2009) . When merging files, I also flag low quality reads and sort the files by position.

```
#!/bin/bash

#Specify input directory
sam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/SAM_files/merge

#Specify output directory
bam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge

files=(${sam_dir}/*.SAM)

for 
file in ${files[@]} 
do 
name=${file} base=`basename ${name} .SAM`

samtools view -b -q -@5 ${sam_dir}/${base}.SAM | samtools sort -o ${bam_dir}/${base}.bam 
done 
``` 
Flags: "-b" = BAM file output. "-q" = filters for low quality reads, BQ<20 id the default. "-o" = output file name. 

At this point, I can check the mapping of the reads. First, I can look at the statistics and flags associated with each read using the following scripts.

First, the coverage needs to be calculated at each base 
```
#!/bin/bash
#
#Specify output directory
bam_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge

#temp storage 
temp=/2/scratch/Katie/WingShapeBSA/temp

files=(${bam_dir}/* _PE.bam)

for file in ${files[@]} do name=${file} base=`basename ${name} .bam`
samtools depth ${bam_dir}/{${base}_PE.bam > ${temp}/${base}_PE.coverage
done 
``` 
Then this script can be used to call the R script below 
```
#!/bin/bash

#temp storage 
temp=/2/scratch/Katie/WingShapeBSA/temp

#Rscripts dir
Rscripts=/home/katie/bin/Rscripts

files=(${temp}/* _PE.bam)
for file in ${files[@]} 
do 
name=${file} 
base=`basename ${name} .bam`
Rscript ${Rscripts}/depth_histogram.R ${base}.bam ${base}
done 
```

```
##coverage_histogram.R
## need next line to call arguments:
args <- commandArgs(trailingOnly = TRUE)

#read in the first argument which should be the file
dat <- read.table(args[1])
#the title should be the second argument (the base name)
title <- args[2] 
colnames(dat) <- c("chr","pos","depth")

#make a histogram of the coverage for each file
pdf(paste(title, ".pdf", sep="")) 
hist(dat$depth, xlim=c(0,500), breaks=500) 
dev.off() 
``` 

![](ProjectFigures/CMO_75R_smallbins2_cov.pdf)

This showed that all samples have a proportion of regions with >200x depth. Because this is much deeper than the genome was sequenced to, this points to PCR duplicates or reads mapped to repetitive regions of the genome in the data set.

I can also look at the average read depth and look for the presence of duplicates in the data set using the following scripts. These will produce a histogram of read depth. We predict our average read depth to be about 100x so regions with higher coverage indicate the presence of PCR duplicates or highly repetitive regions where reads can not be reliably mapped.


###7. GatK realignment 
The regions around indels can produce errors in mapping. Therefore, it is important to relaign reads around these sites. To do this, I will use the gatk tool (https://software.broadinstitute.org/gatk/).

Before mapping, RG information must be added to all files. As mentioned above, in the future, I would do this step earlier in the pipeline so that I can add more accurate and information. 

``` 
#! /bin/bash

#Variable for project:
project_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/

#Path to Picard
picdir=/usr/local/picard-tools/picard.jar


files=(${project_dir}*.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`

java -jar ${picdir} AddOrReplaceReadGroups I=${project_dir}/${base}.bam \
  O=${project_dir}/${base}_RG.bam \
  RGID=L001_L002 \
  RGLB=library1 \
  RGPL=illumina \
  RGPU=None \
  RGSM=${base}

done
``` 
Then indexed files need to be created for each file (Li *et al.* 2009) . This must be run within the directory with the files to index. This will spawn many jobs (one for each file) simultaneously, but the jobs run fairly quickly. 
``` 
#!/bin/bash /

files=(*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam` 
samtools index ${base}_RG.bam &
done
``` 
Interval files also need to be created for all samples. This also needs to be run in the directory with all the files and will spawn many jobs. 
```
#!/bin/bash/
#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar

files=(*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx32g -jar ${gatk} -I ${base}_RG.bam \
		-R /home/katie/flyGenome/dmel_r6.23/gatk/dmel-all-chromosome-r6.23.fasta \
 		-T RealignerTargetCreator \
  	    -o ${base}.intervals &
done
```
Flags: "-I" = input files. "-R" = reference genome. "-T" = what to do with the files. "-o" = output file. 

Next the genome needs to be indexed with samtools (Li *et al.* 2009) .
``` 
samtools faidx dmel-all-chromosome-r6.23.fasta.gz
``` 
Now I can use gatk to realign reads around indels. 
``` 
#!/bin/bash

#Path to input directory
final_bam=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/

#Path to output directory
gatk_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/

#Variable for reference genome (non-zipped)
ref_genome=/home/katie/flyGenome/dmel_r6.23/gatk/dmel-all-chromosome-r6.23.fasta

#Path to GATK
gatk=/usr/local/gatk/GenomeAnalysisTK.jar


files=(${final_bam}/*_RG.bam)
for file in ${files[@]}
do
name=${file}
base=`basename ${name} _RG.bam`

java -Xmx32g -jar ${gatk} -I ${final_bam}/${base}_RG.bam -R ${ref_genome} \
  -T IndelRealigner -targetIntervals ${gatk_dir}/${base}.intervals \
  -o ${gatk_dir}/${base}_realigned.bam

done
``` 
Flags: "-I" = input files. "-R" = reference genome. "-T" = what to do with the files. "-o" = output file. 

A number files failed this step. The error indicates that some multi mapped reads were missing their read pair. I will need to correct this and add these groups into my analysis in the future. The missing samples are: all CMO samples, control samples for FVW12, FVW13/14, PHO and control selection linages for the emc selection experiment.

###8. Remove  duplicates
Following the realignment, we want to remove and PCR duplicates from our data. These can artificially inflate allele counts as the duplicated sequence is not representative of an unique allele. To do this, I will use picard (http://broadinstitute.github.io/picard/index.html).
``` 
#! /bin/bash

mapped_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/
outdir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/wildpop/
picdir=/usr/local/picard-tools/picard.jar

files=(${mapped_dir}/*.bam)

for file in ${files[@]}
do
name=${file}
base=`basename ${name} .bam`
#echo ${name}
java -Xmx2g -jar ${picdir} MarkDuplicates I=${mapped_dir}/${base}.bam O=${outdir}/${base}_rmd.bam M=${outdir}/{base}_dupstat.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true
done
``` 
Flags: "-I" = input files. "-R" = reference genome.  "-O" = output file. "-M" = log file. "VALIDATION_STRINGENCY" = Paul says this needs to be set to silent otherwise picard won't run. "REMOVE_DUPLICATES" = removes duplicates from files. 

I can now look at the depth plots again using the same scripts as above to see if the regions of high coverage have been removed.

I do see the trend that I expect, fewer positions with >200x coverage. Although there is still some high coverage regions. I predict that this is due to repetitive regions with poor mapping. 

###9. Creating mpileup and sync files 
I not want to create a sigle multiple pileup (mpileup) file for each of these experiments so I can use this for future analysis. 
``` 
#! /bin/bash

mapped_dir=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/wildpop
picdir=/usr/local/picard-tools/picard.jar
genome=/home/katie/flyGenome/dmel_r6.23/bwa/dmel-all-chromosome-r6.23.fasta
out_dir=/2/scratch/Katie/WingShapeBSA/pileup/wildpop/

files=(${mapped_dir}/*_rmd.bam)
for file in ${files[@]}
do
#echo ${files[@]}
name=${file}
base=`basename ${name} _rmd.bam`
#echo ${base}
samtools mpileup -B -Q 0 -f ${genome} ${base}_rmd.bam > ${out_dir}/${base}.mpileup &
done
``` 
Flags: "-B" = disable base quality alignment, we have already had a lot of checks for this. "-Q" = sets min quality to zero as we have already checked this many times.  "-f" = fasta reference file.

I also created a sync file for each mpileup, run in the same directory os the mpileup file. 
``` 
java -Xmx32g -jar /home/katie/bin/popoolation2_1201/mpileup2sync.jar --input ee_nogatk.mpileup --output ee_nogatk.sync

```

Following these checks, I can begin the analysis of my samples.

## *PoPoolation analysis*
PoPoolation is designed to calculate population parameters within populations of pool sequencing data (Kofler *et al.* 2011a). For this project, I will be using it to calculate Tiajama's pi

###Calculation of pi in samples 
First, I had to calculate pi across the genome, I used relatively large windows for this to speed up analysis since I am only interested in general trends right now. 
```
#/bin/bash

#Path to PoPoolation
pi=/home/katie/bin/popoolation_1.2.2/Variance-sliding.pl

# Path to input directory
input=/2/scratch/Katie/WingShapeBSA/pileup/ee

# Path to output Tajima Pi files
output=/2/scratch/Katie/WingShapeBSA/population/

files=(${input}/*.bam)

for file in ${files[@]}

do

name=${file}

base=`basename ${name} _merged_aligned_PE.bam`

perl ${pi} \
	--input ${input}/${base}._merged_aligned_PE.bam` \
	--output ${output}/${base}.pi \
	--measure pi \
	--window-size 1000 \
	--step-size 1000 \
	--min-count 2 \
	--min-coverage 4 \
	--max-coverage 250 \
	--min-qual 20 \
	--pool-size 150 \
	--fastq-type sanger \
	--snp-output ${output}/${base}.snps \
	--min-covered-fraction 0.5
done 
``` 
This can then be visualized with this script In my samples, pi values matched what I expected to see: a reduction of diveristy on the X chromosome, 4th chromosome and around centromeres and telomeres. 

``` 
library(tidyverse)

dat <- read.table("CMO_75L_ATTACTCG-TATAGCCT.pi")
colnames(dat) <- c('chr', 'window', 'windowCount', ' propInwindow', 'Pi')

##Taking out the NAs
dat <- dat[-which(dat$Pi=="na"),]


##cleaning to only keep what I care about 
#Because I am using v6 of the genome that has a bunch of scafolds so its easier to just select what I want. 
#The second line gives a number to each line to number along the chromosome. 
datX <- subset(dat, chr == "X")
a <- dim(datX)[1]
datX$number <- 1:a


dat2L <- subset(dat, chr == "2L")
b <- dim(dat2L)[1]
dat2L$number <- (a+1):(a+b)

dat2R <- subset(dat, chr == "2R")
c <- dim(dat2R)[1]
dat2R$number <- (a+b+1):(a+b+c)

dat3L <- subset(dat, chr == "3L")
d <- dim(dat3L)[1]
dat3L$number <- (a+b+c+1):(a+b+c+d)

dat3R <- subset(dat, chr == "3R")
e <- dim(dat3R)[1]
dat3R$number <- (a+b+c+d+1):(a+b+c+d+e)

dat4 <- subset(dat, chr == "4")
f <- dim(dat4)[1]
dat4$number <- (a+b+c+d+e+1):(a+b+c+d+e+f)

##add all together 
dat2 <- rbind(datX, dat2L, dat2R, dat3L, dat3R, dat4)

##Pi needs to be numeric 
dat2$Pi=as.numeric(levels(dat2$Pi))[dat2$Pi]


plot1 <- ggplot(dat2, aes(x = number, y= Pi, colour = chr))+
  geom_point(size=0.3, show.legend = F) +
  #scale_y_continuous(limits=c(0, 0.02), breaks=seq(0, 0.02, 0.005)) + 
  xlab("") +
  scale_x_discrete(limits=c(a, a+b, a+b+c, a+b+c+d, a+b+c+d+e, a+b+e+c+d+e+f), labels = c("X", "2L", '2R', '3L', '3R', "4")) +
  theme(text = element_text(size=20), 
        axis.text.x= element_text(size=15), axis.text.y= element_text(size=15)) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4'))
  
pdf("",width=1060,height=412,units="px")
plot1
dev.off()
  ```

### *PoPoolation2 analysis*
PoPoolation2 is a population genetics tool to compare between populations in pool-seq data (Kofler *et al.* 2011b). For this project, I chose to calculate the Fst values because this is the easiest parameter to interrupt and compare for both the wild population and artificial selection experiments.

### Calculation of Fst in samples 
First Fst was calculated with this script. This makes pairwise comparisons between all populations in the sync file 
``` 
perl /home/katie/bin/popoolation2_1201/fst-sliding.pl \
				--input /2/scratch/Katie/WingShapeBSA/pileup/ee/ee_nogatk.sync \
				--output /2/scratch/Katie/WingShapeBSA/population/ee/ee_nogatk.fst \
				--suppress-noninformative \
				--min-count 2 \
				--min-coverage 10 \
				--max-coverage 300 \
				--min-covered-fraction 1 \
				--window-size 50 \
				--step-size 50 \
				--pool-size 150
```

Then I needed to remove the comparisons that I am interested in from the Fst tile Population numbers are assigned based on the order that they were listed in when the mpileup was made (alphabetical)

```
#Each line should be run within the directory with the fst calculation file. 
#ee ds UPvsDN gatk
awk '{print $1, $2, $3, $4, $5, $47, $48, $48, $57, $58, $59, $66, $67, $68}' ee_gatk.fst > UPvDNgatk.fst

#ee ds UPvsDN nogatk
awk '{print $1, $2, $3, $4, $5, $47, $48, $48, $57, $58, $59, $66, $67, $68}' ee_nogatk.fst > UPvDNnogatk.fst

#ee ds UPvsCN gatk
awk '{print $1, $2, $3, $4, $5, $11, $12, $13, $24, $25, $26, $36, $37, $38}' ee_gatk.fst > UPvCNgatk.fst

#ee ds UPvsCN nogatk
awk '{print $1, $2, $3, $4, $5, $11, $12, $13, $24, $25, $26, $36, $37, $38}' ee_nogatk.fst > UPvCNnogatk.fst

#ee ds DNvsCN gatk
awk '{print $1, $2, $3, $4, $5, $8, $9, $10, $21, $22, $23, $33, $34, $35}' ee_gatk.fst > DNvCNgatk.fst

#ee ds DNvsCN nogatk
awk '{print $1, $2, $3, $4, $5, $8, $9, $10, $21, $22, $23, $33, $34, $35}' ee_nogatk.fst > DNvCNnogatk.fst

##ee emc UPvsDN gatk
awk '{print $1, $2, $3, $4, $5, $98, $99, $100, $102, $103, $104, $105, $106, $107}' ee_gatk.fst > emcUPvDNgatk.fst

##ee emc UPvsDN nogatk
awk '{print $1, $2, $3, $4, $5, $98, $99, $100, $102, $103, $104, $105, $106, $107}' ee_nogatk.fst > emcUPvDNnogatk.fst

#wild ds LvR gatk
awk '{print $1, $2, $3, $4, $5, $6, $14, $16, $23, $25, $66, $68, $69, $71}' wild_gatk.fst > wild_ds_gatk.fst

#wild ds LvR nogatk
awk '{print $1, $2, $3, $4, $5, $6, $14, $16, $23, $25, $66, $68, $69, $71}' wildpop_nogatk.fst > wild_ds_nogatk.fst

#wild neur LvR gatk
awk '{print $1, $2, $3, $4, $5, $27, $29, $31, $36, $38, $44, $46, $51, $57}' wild_gatk.fst > wild_neur_gatk.fst

#wild neur LvR gatk
awk '{print $1, $2, $3, $4, $5, $27, $29, $31, $36, $38, $44, $46, $51, $57}' wildpop_nogatk.fst > wild_neur_nogatk.fst 
```

Fst for the whole genome could be visualized with this script. This was performed on my local machine so that the number of collums could be easily edited.

```
#All chromosomes
library(data.table)
library(ggplot2)

ddat2 <- fread("UPvDNgatk.fst")


ccol <- ncol(ddat2)

#cleans up col to get only the fst value
for (i in 6:ccol){
  ddat2[[i]] <- gsub(".*=","", ddat2[[i]])
}


#to change the columns to numeric 
for (i in 6:ccol){
  ddat2[[i]] <- as.numeric(ddat2[[i]])
}

#cal mean Fst for all comparisons. 
ddat2$meanFst <- rowMeans(subset(ddat2, select = c(6:ccol)), na.rm = TRUE)

ddat2 <- ddat2[ddat2$meanFst!='NaN',]


#choose the last column for the fst mean
#####YOU WILL NEED TO CHANGE THE LAST NUMBER HERE
ddat2 <- ddat2[,c(1,2,3,4,5,15)]

colnames(ddat2) <- c('chr', 'window', "num", 'frac', 'meanCov','meanFst')

##need to select only the chr (get rid of all the shit)
ddat22L <- ddat2[which(ddat2$chr=='2L'),]
ddat22R <- ddat2[which(ddat2$chr=='2R'),]
ddat23L <- ddat2[which(ddat2$chr=='3L'),]
ddat23R <- ddat2[which(ddat2$chr=='3R'),]
ddat24 <- ddat2[which(ddat2$chr=='4'),]
ddat2X <- ddat2[which(ddat2$chr=='X'),]

ddat2 <- rbind(ddat2X, ddat22L, ddat22R, ddat23L, ddat23R, ddat24)

#this part below is so you change the order of the chromosomes to make a nice map going across the x-axis. 
#I need to check and see which order I had my chromosomes in by pulling out rows of data equal to letters below. 
#So I'd look at ddat2[1,] then if the first value is on X I would look at ddat2[1 + l,]. This is how you choose the order of the chromosomes as you will see below.
g <- nrow(ddat2[which(ddat2$chr=='2L'),])
h <- nrow(ddat2[which(ddat2$chr=='2R'),])
i <- nrow(ddat2[which(ddat2$chr=='3L'),])
j <- nrow(ddat2[which(ddat2$chr=='3R'),])
k <- nrow(ddat2[which(ddat2$chr=='4'),])
l <- nrow(ddat2[which(ddat2$chr=='X'),])




#NOTE: Chaging Orders:
#To change the order for X to be first:
# need to figure out the order the chromosomes are put in the data frame to give each the correct number in the sequence

#Example: I want X first but the rest the same, this will have X have numbers 1:l (last line) and then start with 2L (first line)



#X-2L-2R-3L-3R-4
ddat2$number <-  c((1:l),
                   (l+1):(l+g), 
                   (l+g+1):(l+g+h), 
                   (l+g+h+1):(l+g+h+i),
                   (l+g+h+i+1):(l+g+h+i+j),
                   (l+g+h+i+j+1):(l+g+h+i+j+k))
#Basically: each chromosome will correspond to one split, just need to move it around based initial order (assumeing the order you want is X-2L-2R-3L-3R-4)



### PLOTS:

plt <-  ggplot(data = ddat2, aes(x=number, y=meanFst, color=chr))
plt2 <- plt + 
  geom_point(size=0.5, show.legend = F) + 
  theme(panel.background = element_blank()) +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.1)) +
  xlab("Chromosome") +
  ylab(expression(F[ST])) +
  scale_x_discrete(limits=c(l/2, l+(g/2), (l+g+(h/2)), (l+g+h+(i/2)), (l+g+h+i+(j/2)), (l+g+h+i+j+(k/2))), labels = c("X","2L", "2R", '3L', '3R', "4")) +
  scale_colour_manual(values=c("#56B4E9", "#E69F00", 'grey30', 'grey46', 'wheat3', 'lemonchiffon4')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))

png("UPvDNgatk.png",width=1060,height=412,units="px")
plt2
dev.off() 
```
Finally, I wanted to look at allele frequency changes around the gene of interest that had been selected on. To do this, I subset the chromosome of interest and coloured a small region around the gene of interest blue in the image. 

```
library(data.table)
library(ggplot2)

ddat2 <- fread("wild_ds_gatk.fst")


ccol <- ncol(ddat2)

#cleans up col to get only the fst value
for (i in 6:ccol){
  ddat2[[i]] <- gsub(".*=","", ddat2[[i]])
}


#to change the columns to numeric 
for (i in 6:ccol){
  ddat2[[i]] <- as.numeric(ddat2[[i]])
}


#cal mean Fst for all comparisons. 
ddat2$meanFst <- rowMeans(subset(ddat2, select = c(6:ccol)), na.rm = TRUE)

ddat2 <- ddat2[ddat2$meanFst!='NaN',]


#choose the last column for the fst mean
#####YOU WILL NEED TO CHANGE THE LAST NUMBER HERE
ddat2 <- ddat2[,c(1,2,3,4,5,15)]

colnames(ddat2) <- c('chr', 'window', "num", 'frac', 'meanCov','meanFst')

chr2L <- ddat2[which(ddat2$chr=='2L'),]
a <- nrow(chr2L)
chr2L$number <- c(1:a)

##selcect window around ds 
#looked up ds and found flanking windows 
b <- chr2L$number[which(chr2L$window==640025)]
c <- chr2L$number[which(chr2L$window==715075)]

pre <- chr2L[1:(b-1),]
dach <- chr2L[b:c,]
post <- chr2L[(c+1):a]

pre$position <- "notds"
dach$position <- "ds"
post$position <- "notds"

chr2L <- rbind(pre, dach, post)

### PLOTS:

plt <-  ggplot(data = chr2L, aes(x=number, y=meanFst, color=position))
plt2 <- plt + 
  geom_point(size=0.5, show.legend = F) + 
  theme(panel.background = element_blank()) +
  #scale_y_continuous(limits=c(0, 0.5), breaks=seq(0, 0.5, 0.1)) +
  xlab("Chromosome 2L") +
  ylab(expression(F[ST])) +
  scale_colour_manual(values=c('#56B4E9', 'grey45')) +
  theme(text = element_text(size=20),
        axis.text.x= element_text(size=15), 
        axis.text.y= element_text(size=15))

png("UPvDNgatk_2L.png",width=1060,height=412,units="px")
plt2
dev.off()
``` 
Only representative images are shown here, but because I want to discuss these results in detail below, all images are provided within my Biol722 Github directory.

As expected, Fst is much higher across the genome in the artificial selection than the wild population analysis. Visually, there is little to no change between the graphs using data that had or had not undergone indel realignment. However, as discussed below this step did have an effect on the number of variants detected in the genome. 

What is also unexpected is that there is no clear change in allele frequency associated with the gene of interest in any case. Although we know that this is a highly polygenic trait, we expected a change in allele frequency at this site in addition to a number of others. This is interesting and I will consider my Fst calculation parameters more closely to determine if something has been left out of the analysis.  

### *CRISP variant caller*
CRISP is a variant caller designed to work with pool-seq data It requires multiple pools of individuals and compares polymorphisms found between pools, with more significance assigned to variants found in more than one sample (Bansal et al 2011). This has been demonstrated to be the variant caller that is most accurate for pool-seq data with large pools of individuals.

CRISP is run with the script below. This script had to be modified for each population, this is one example. 
``` 
#/bin/bash

#Path to CRISP
crisp=/home/katie/bin/CRISP-122713/CRISP

#Variable for reference genome (non-zipped)
ref_genome=/home/katie/flyGenome/dmel_r6.23/samtools/dmel-all-chromosome-r6.23.fasta

#Path to .bam files from GATK
input=/2/scratch/Katie/WingShapeBSA/mapped/bwa/BAM_files/merge/gatkindel/ee

#Output
output=/2/scratch/Katie/WingShapeBSA/crisp/ee/


${crisp} --bams ${input}/emcUP.txt \
			--ref ${ref_genome} \
			--poolsize 150 \
 			--qvoffset 33 \
			--mbq 10 \
			--mmq 10 \
 			--minc 2 \
 			--VCF ${output}/emcUP_gatk.vcf > ${output}/emcUP_gatk.log
``` 
I counted the number of lines in these files (each line corresponds to a SNP, with the first 22 lines as a header)
Some files are also missing here because CRISP would not read some of my files, I am not sure why this happened because the same code worked with all the other files. I will need to look more into this. 
```
   1973247 dsCN_notgatk.vcf
   1931251 dsDN_gatk.vcf
   1955317 dsDN_notgatk.vcf
   1949294 dsUP_notgatk.vcf
   1925600 emcDN_gatk.vcf
   1952376 emcDN_notgatk.vcf
   1775582 emcUP_gatk.vcf
   3064330 dsL_nogatk.vcf
   2967749 dsR_nogatk.vcf
   3151222 neurL_gatk.vcf
   3162525 neurL_nogatk.vcf
```
There are >10000 fewer variants detected after indel realignment than before indel realignment. If we are interested in identifying 'true' variants for downstream analysis, the reallignment step is crucial to help remove false positives from the data set. 

***

###Additions to the pipeline in the future 
###
Although CRISP has been demonstrated to be the most accurate SNP caller for pool sequencing data with large pools of individuals, there is still a degree of identified polymorphisms that are false positives (Huang *et al* 2015). Because I want to do further analysis on the identified SNPs that would require a significant investment in time and resources, I want to be confident that the polymorphisms I am working with are real. One way that I can approach this problem is to cross reference my list of identified SNPs with a list generated from a different program. A number of other programs have been identified for their ability to call SNPs in poolSeq
data such as LoFreq and VarScan (Wilm *et al* 2012, Koboldt *et al* 2012). Although these programs all have a higher rate of false positives for poolSeq data, we can put more trust in a SNPs identified by both programs. 

In addition to a second variant caller, we can also use a second mapper. Because different mappers use unique algorithms to map reads, each program will produce unique errors in mapping. By comparing the polymorphisms identified at the end of the pipeline from the two different mapping programs, we should be able to increase our ability to identify ‘real’ polymorphisms by only considering those found in both analysis pipelines. A number of other mapping programs exist but I plan to use novallign because it is already in use in the lab (http://www.novocraft.com). 

The pipeline is also missing a repeat masker step. This is an essential consideration to include when we are doing analyses in the future. Not masking these regions will introduce additional polymorphisms nearby. This will help to further reduce the number of false positives in the data set. In the future, I would implement RepeatMasker as a tool to remove these regions from analysis (http://www.repeatmasker.org/). 

##Next steps 
With these data, I will identify the polymorphisms contributing to phenotypic change in populations. The Fst analysis can be used to identify the regions of the genome with significant allele frequency changes between populations. With that information, I can look for variants present within the selected population but absent in the opposing population and control population. In addition, I can perform a Cochran-Mantel-Haenszel (CMH) test to test for consistent allele changes among biological replicates (Kofler et al.2011b). Although this is designed for single pairs of populations, as is the case for the wild-caught populations in this analysis, it has also been applied to identify changes in multiple replicate lineages undergoing artificial selection. Once variants of interest have been identified, a number of follow up experiments can be done to study the effect size and epistatic interactions with these SNPs. 

#References 

***

Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data." Babraham Bioinformatics - FastQC A Quality Control Tool for High Throughput Sequence Data. 

Bolger, A.M., Lohse, M., & Usadel, B. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, btu170.

Flint J. & Mackay T.F. 2009. Genetic architecture of quantitative traits in mice, flies, and humans. Genome Res. 19(5):723-33

Houle D., Bolstad G.H., van der Lind K., Handsen T.F. 2017. Mutation predicts 40 million years of fly wing evolution. Nature. 548:447–450

Huang, H.W., NISC Comparative Sequencing Program, Mullikin J.C., Hansen N.F. 2015. Evaluation of variant detection software for pooled next-generation sequence data. BMC Bioinformatics. 16: 235

Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, et al. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res. 2012;22(3):568–76.

Kofler R, Orozco-terWengel P, De Maio N, Pandey R.V, Nolte V, Futschik A, Kosiol C, Schlotterer C. 2011a. PoPoolation: a toolbox for popula- tion genetic analysis of next generation sequencing data from pooled individuals. PLoS One6:e15925.

Kofler R, Pandey R.V, Schlotterer C. 2011b. PoPoolation2: identifying dif- ferentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27:3435–3436.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 25:2078–2079.

Li H. and Durbin R. 2010. Fast and accurate long-read alignment with Burrows-Wheeler Transform. Bioinformatics, Epub. [PMID: 20080505]

Mackay, T.F. *et al*. 2012.The Drosophila melanogaster Genetic Reference Panel. Nature.482: 173–178

Pitchers W., Nye J., Marquez E.J., Kowalski A., Dworkin I., Houle D. (2019) A Multivariate Genome-Wide Association Study of Wing Shape in Drosophila melanogaster. 211(4):1429-1447

Schotterer C., Trobler R.,  Kofler R, Nolte, V. 2014. Sequencing pools of individuals — mining genome-wide polymorphism data without big funding. Nature Review Genetics. 15:749–763

Wilm A, Aw PP, Bertrand D, Yeo GH, Ong SH, Wong CH, et al. LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Res. 2012;40(22):11189–201.
