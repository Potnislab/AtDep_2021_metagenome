1.	The first step after receiving the raw reads from sequencing is trimming for adapter and contaminant removal. We used BBDuk (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) to align our metagenomic reads to the adapter file and remove those adapters to obtain a clean output file without adapters. We loaded BBMap from the module as it was already installed in the cluster.

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bbmap/37.36
for fq1 in ./*_R1_001.fastq.gz
    do
    echo "working with file $fq1"
    base=$(basename $fq1 _R1_001.fastq.gz)
    echo "base name is $base"
    fq1=./${base}_R1_001.fastq.gz
    fq2=./${base}_R2_001.fastq.gz
bbduk.sh -Xmx1g in1=$fq1 in2=$fq2 out1=${base}.clean1.fq out2=${base}.clean2.fq ref=adapter.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo

bbduk.sh -Xmx1g in1=${base}.clean1.fq in2=${base}.clean2.fq out1=${base}.clean_trim1.fq out2=${base}.clean_trim2.fq qtrim=rl trimq=10
done

2.	The adapter trimmed reads were then trimmed for quality using SolexaQA++ (https://solexaqa.sourceforge.net). Dynamictrim and Lengthsort algorithm of the program was used to remove the poor-quality reads with Phred quality score <20 and length less than 50bp

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load solexaqa/3.1.7.1
for F in *.fq; do
N=$(basename $F .fq) ;
SolexaQA++ dynamictrim -h 20 $F ; done

And for length

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load solexaqa/3.1.7.1
for F in *.trimmed; do N=$(basename $F .trimmed) ;
SolexaQA++ lengthsort -l 50 $F ; done

3.	The quality trimmed reads were then used to remove host contamination using KneadData (https://github.com/biobakery/kneaddata).  KneadData performs in sillico separation of microbial reads from those of host. We used NCBI to download the complete pepper (GCA_021292125.1) genomes, then created a bowtie2 database (https://github.com/BenLangmead/bowtie2). Here, we renamed the fna file of pepper from NCBI to pepper.fna for our ease.

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bowtie2/2.2.9
bowtie2-build /scratch/aubrrb/pepper.fna /scratch/aubrrb/pepper_ref --threads 4

The above command will create pepper_ref database folder, which will be our reference for the pepper reads. We then use KneadData to remove the pepper reads.

#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load anaconda/3-2018.12
for fq1 in ./*_trim1.fq.trimmed.single
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _trim1.fq.trimmed.single )
    echo "base name is $base"

    fq1=./${base}_trim1.fq.trimmed.single
    fq2=./${base}_trim2.fq.trimmed.single

kneaddata --input $fq1 --input $fq2 --reference-db ./pepper_ref --output ./${base}.kneaddata --bypass-trim -t 6
done

KneadData will produce the clean reads, removed reads, and log files as an output of the above step. Then we copied all the kneaddata output into a final folder using the following command for both the forward and reverse reads (use the same command for reverse(*_paired_2.fastq) reads).

find . -name '*_paired_1.fastq' -exec mv --target-directory=' /scratch/aubrrb/raw/quality_reads/AtDep/AtDep_2020' '{}' +


4.	KneadData is not pair end aware, so it might have an issue with the different number of forward and reverse reads. To remove that error, we used repair.sh script from BBMap (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/). 

#!/bin/bash

source /opt/asn/etc/asn-bash-profiles-special/modules.sh

module load bbmap

for fq1 in ./*.clean_trim1.fq.trimmed_kneaddata_paired_1.fastq

    do
    echo "working with file $fq1"
    base=$(basename $fq1 .clean_trim1.fq.trimmed_kneaddata_paired_1.fastq)
    echo "base name is $base"
    fq1=./${base}.clean_trim1.fq.trimmed_kneaddata_paired_1.fastq
    fq2=./${base}.clean_trim1.fq.trimmed_kneaddata_paired_2.fastq

repair.sh overwrite=t in1=$fq1 in2=$fq2 out1=${base}.paired.1.fq out2=${base}.paired.2.fq outs=${base}.singletons.fq repair
done

5.	Our reads now are free from adapter, low quality reads and host contamination. These quality reads will now be used for our downstream analysis. We really don’t have to convert the fastq file to fasta for any of our analysis but if you need to convert, we can use the enveomics package (https://github.com/enveomics/enveomics). 

!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load enveomics/20dec2017
for F in *.fastq; do N=$(basename $F .fastq) ;
awk -f ./FastQ.toFastA.awk < $F > $N.fasta ; done

