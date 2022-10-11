1. To identify the eukaryotic diversity across the samples, we use EukDetect (https://github.com/allind/EukDetect). EukDetect was 
installed using the conda installation and the database was downloaded.

conda env update --name eukdetect -f environment.yml

source activate eukdetect

# install eukdetect

python setup.py install

#We then created a directory for eukaryotic database and downloaded the database.

mkdir eukdb
cd eukdb
wget https://ndownloader.figshare.com/files/34885596
tar -zxvf 34885596
rm 34885596

#We then copied the default_configfile.yml file and changed all parameters in the config file as described, including the location of the database, 
raw files, output, etc. 

#Then we ran EukDetect with the following command

#!/bin/bash
source activate eukdetect
eukdetect --mode runall –configfile my_configfile.yml –cores 10

#The output from the run will have two outputs with the name *hits_table.txt and *filtered_hits_eukfrac. The RPKS (Read Per Kilobase of Sequence)  value 
from eufrac is related to absolute abundance. The library size of all the samples was used to create a scaling factor and the normalization was done. Then the RPKS was multiplied with the scaling factor to compare the abundance across samples.
#Finally, the tables are joined using comine_rpks_eukfrac.py command to give a final single table.
