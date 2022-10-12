
1.	Functional profiling of the metagenome reads was done using HUMAnN3.0 (https://github.com/biobakery/humann). HUMAnN profiles the presence/absence and abundance of microbial pathways from metagenome sequence data. 
#We installed HUMAnN3.0 using the conda installation 
conda create -n humann3
source activate humann3
#Install human using the following code
conda install -c biobakery human

2.	HUMAnN3.0  doesn’t work with paired-end reads, so the reads need to be either merged or concatenated. We used BBMerge (https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/) to merge these reads.
#!/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load bbmap/37.36
for fq1 in ./* .paired.1.fq
    do
    echo "working with file $fq1"
    base=$(basename $fq1 .paired.1.fq)
    echo "base name is $base"
    fq1=./${base}.paired.1.fq
    fq2=./${base}.paired.2.fq
bbmerge.sh in1=$fq1 in2=$fq2 out=${base}_merged.fastq
done
3.	Then we download the ChocoPhlAn database for HUMAnN3.0. We created a folder chocophlan database and downloaded the database.
humann_databases --download chocophlan full /home/aubrrb/humann_database/chocophlan

#We also download a translated search database for analysis of genes if required. 
human_databases –download uniref uniref90_diamond /home/aubrrb/humann_database/uniref

4.	Then we updated our HUMAnN 3.0 configuration file to include the locations of our downloaded databases using the following steps.
#Update the location of the ChocoPhlAn database 
humann_config --update database_folders nucleotide /home/aubrrb/humann_database/chocophlan/chocophlan
#Update the location of the UniRef database
 humann_config --update database_folders protein home/aubrrb/humann_database/uniref/uniref

5.	HUMAnN requires MetaPhlAn indexed files so we also installed c3.0 (https://huttenhower.sph.harvard.edu/metaphlan/). 
#It is recommended to create an isolated conda environment and install MetaPhlAn.
conda create --name mpa -c bioconda python=3.7 metaphlan
source activate mpa
#If we have installed MetaPhlAn using Anaconda, it is advised to install the database in a folder outside the conda environment. 
#!/bin/bash
module load bowtie2/2.2.9
metaphlan --install --bowtie2db /scratch/aubrrb/metaphlan


# If we need to run the index separately, we can use MetaPhlAn. But HUMAnN can do this work for us.
module load anaconda/2-4.2.0_cent
source activate mpa
for fq1 in ./*_paired_1.fastq
    do
    echo "working with file $fq1"
    base=$(basename $fq1 _paired_1.fastq)
    echo "base name is $base"
    fq1=./${base}_paired_1.fastq
    fq2=./${base}_paired_2.fastq
metaphlan $fq1,$fq2 --nproc 5 --input_type fastq -o ./${base}metagenome.txt --index mpa_v30_CHOCOPhlAn_201901 --bowtie2out ${base}.bz2
done
6.	Then we ran HUMAnN with the following command (make sure we have a metaphlan index also)
#!/bin/bash
source activate humann3
for i in *.fastq
do
  humann --input $i \
    --output hmn2_output \
--metaphlan-options "-x done


#The results will be three main output files for each input file named $SAMPLE_genefamilies.tsv, $SAMPLE_pathabundance.tsv, and $SAMPLE_pathcoverage.tsv.
7.	Then the HUMAnN output tables were normalized using cpm for both the pathabundance and gene families

#!/bin/bash
for SAMPLE in ./*_pathabundance.tsv; do
    humann_renorm_table --input $SAMPLE --output ${SAMPLE%.tsv}_cpm.tsv --units cpm --update-snames

done


#!/bin/bash
for SAMPLE in ./*_genefamilies.tsv; do
    humann_renorm_table --input $SAMPLE --output ${SAMPLE%.tsv}_cpm.tsv --units cpm --update-snames
done


8.	The we copied all gene and pathways files in a newfolder and the tables were joined together for a final single table.

humann_join_tables --input . --output pepper_genefamilies.tsv --file_name genefamilies

humann_join_tables --input . --output pepper_pathabundance.tsv --file_name pathabundance

9.	Then the tables was stratified to split the table with stratified (pathways+microbes) and unstratified (pathways only. 
humann_split_stratified_table --input pepper_genefamilies.tsv --output pepper_genefamilies

humann_split_stratified_table --input pepper_pathabundance.tsv --output pepper_pathabundance


10.	Then the table were normalized for their abundance using renorm script. The script should be submitted in queue as it requires more memory.

humann_renorm_table --input pepper_genefamilies.tsv --output 
pepper_genefamilies_cpm.tsv --units cpm

humann_renorm_table --input pepper_pathabundance.tsv --output pepper_pathabundance_relab.tsv --units relab


humann_split_stratified_table --input pepper_genefamilies_cpm.tsv --output pepper_genefamilies_cpm

humann_split_stratified_table --input pepper_pathabundance_relab.tsv --output pepper_pathabundance_relab

11.	(Optional) If we need to regroup the gene families nased on KEGG Orthogroups or EggNOG, we can use the following script.
#run regroup in queue
humann_regroup_table --input pepper_genefamilies.tsv --custom /home/aubrrb/humann_database/utility_mapping/map_ko_uniref90.txt.gz  \
--output pepper_genefamilies_ko.tsv

humann_renorm_table --input pepper_genefamilies_ko.tsv --output pepper_genefamilies_ko_cpm.tsv --units cpm

humann_split_stratified_table --input pepper_genefamilies_ko_cpm.tsv --output .

humann_split_stratified_table --input pepper_genefamilies_ko.tsv --output .


#run in quueue

humann_regroup_table --input pepper_genefamilies.tsv --groups uniref90_eggnog  \
--output pepper_genefamilies_cog.tsv


humann_renorm_table --input pepper_genefamilies_cog.tsv --output pepper_genefamilies_cog_cpm.tsv --units cpm

humann_split_stratified_table --input pepper_genefamilies_cog_cpm.tsv --output .

humann_split_stratified_table --input pepper_genefamilies_cog.tsv --output .

12.	We then imported the data into QIIME2 (https://qiime2.org/) to calculate various alpha and beta diversity indices.  We first convert the tsv files into biom files from which we will convert them to QIIME format.

#!/bin/bash
source activate qiime
biom convert -i /scratch/aubrrb/atdep_humann3/newfolder/pepper_pathabundance/pepper_pathabundance_unstratified.tsv \
-o ./pepper_pathabundance_unstratified.biom \
--table-type "Pathway table" --to-hdf5

biom convert -i /scratch/aubrrb/atdep_humann3/newfolder/pepper_genefamilies_ko_unstratified.tsv \
-o ./pepper_genefamilies_ko_unstratified.biom \
--table-type="OTU table" --to-hdf5

biom convert -i /scratch/aubrrb/atdep_humann3/newfolder/pepper_genefamilies_cog_unstratified.tsv \
-o ./pepper_genefamilies_cog_unstratified.biom \
--table-type="OTU table" --to-hdf5

13.	We then use the summarization script to have the overall sampling depth and other summary statistics for our sample. The information will be used later in the analysis.
biom summarize-table -i ./pepper_pathabundance_unstratified.biom \
-o ./pepper_pathabundance_unstratified_summary.txt

biom summarize-table -i ./pepper_genefamilies_ko_unstratified.biom \
-o ./pepper_genefamilies_ko_unstratified_summary.txt

biom summarize-table -i ./pepper_genefamilies_cog_unstratified.biom \
-o ./pepper_genefamilies_cog_unstratified_summary.txt

14.	We then converted the files into .qza format which is QIIME format.
#For pathways

#!/bin/bash
source activate qiime2
qiime tools import \
--input-path pepper_pathabundance_unstratified.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path pepper_pathabundance_unstratified.qza

#For gene families KO

#!/bin/bash
source activate qiime2
qiime tools import \
--input-path /scratch/aubrrb/atdep_humann3/newfolder/pepper_genefamilies_ko_unstratified.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path pepper_genefamilies_ko_unstratified.qza

#Gene families cog

#!/bin/bash
source activate qiime2 
qiime tools import \
--input-path /scratch/aubrrb/atdep_humann3/newfolder/pepper_genefamilies_cog_unstratified.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV210Format \
--output-path /scratch/aubrrb/atdep_humann3/newfolder/pepper_genefamilies_cog_unstratified.qza

15.	Then we rarefied the total samples based on the lowest possible read depth from our summarize table step.

#!/bin/bash
source activate qiime2
qiime diversity core-metrics \
--i-table pepper_pathabundance_unstratified.qza \
--p-sampling-depth  2454811 \
--m-metadata-file pepper metadata_pa.txt \
--output-dir pepper_pathabundance_unstratified

#Just in case we want to see what we have in the qza file we have to first convert .qza to biom and then to tsv.

qiime tools export --input-path pepper_genefamilies_ko_unstratified.qza --output-path exported-feature-table

biom convert -i feature-table.biom -o table.from_biom.txt --to-tsv


#The metadata file for each of them has to be made separately as they have different extension in their sample id name. Please check the .qza or .tsv file for each of the output and prepare the metadata accordingly.

#!/bin/bash
source activate qiime2 
qiime diversity core-metrics \
--i-table pepper_genefamilies_ko_unstratified.qza \
--p-sampling-depth 5774881 \
--m-metadata-file pepper metadata_gf.txt \
--output-dir pepper_genefamilies_unstratified



#!/bin/bash
source activate qiime2
qiime diversity core-metrics \
--i-table pepper_genefamilies_cog_unstratified.qza \
--p-sampling-depth 6103171 \
--m-metadata-file pepper metadata_gf.txt \
--output-dir pepper_genefamilies_cog_unstratified

16.	Then we export the distance matrix from qza based on Bray-Curtis and Jaccard distance so we can use them with R for further analysis.

#!/bin/bash
source activate qiime2
qiime tools export \
--input-path ./pepper_pathabundance_unstratified/bray_curtis_distance_matrix.qza \
--output-path bray-curtis-distance-matrix-pa


#!/bin/bash
source activate qiime2
qiime tools export \
--input-path ./pepper_pathabundance_unstratified/jaccard_distance_matrix.qza \
--output-path jaccard-distance-matrix-pa

#!/bin/bash
source activate qiime2
qiime tools export \
--input-path ./pepper_genefamilies_cog_unstratified/bray_curtis_distance_matrix.qza \
--output-path bray-curtis-distance-matrix-cog

#!/bin/bash
source activate qiime2
qiime tools export \
--input-path ./pepper_genefamilies_cog_unstratified/jaccard_distance_matrix.qza \
--output-path jaccard-distance-matrix-cog

#!/bin/bash
source activate qiime2
qiime tools export \
--input-path ./pepper_genefamilies_unstratified/bray_curtis_distance_matrix.qza \
--output-path bray-curtis-distance-matrix-gf

#!/bin/bash
source activate qiime2
qiime tools export \
--input-path ./pepper_genefamilies_unstratified/jaccard_distance_matrix.qza \
--output-path jaccard-distance-matrix-gf
