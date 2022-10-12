1.	For differentially abundant pathways analysis, LEfSe (https://huttenhower.sph.harvard.edu/lefse/) was used. LEfSe was installed with the conda installation. 
We combined the pathway output files from HUMAnN by selecting only the control and inoculated in both the cultivars, ambient and elevated ozone and without 
inoculation, and ambient and inoculated plus elevated ozone to see the pathways enriched under these different stress. 

conda create -n lefse -c bioconda python=2.7
conda activate lefse
conda install -c bioconda lefse


2. Then we formatted the input the way LEfSe supports. 

lefse_format_input.py ambient-vs-ozone.txt  ambient-vs-ozone.in -c 1 -u 1 -o 1000000

lefse_format_input.py inoculated-vs-control.txt  inoculated-vs-control.in -c 1 -u 1 -o 1000000

lefse_format_input.py ambient-vs-combined-stress.txt ambient-vs-combined-stress.in -c 1 -u 1 -o 1000000

3.	Then run LEfSe with a bootstrap of 200 and LDA score of 3

lefse_run.py ambient-vs-ozone.in ambient-vs-ozone.res -b 200 -l 3

lefse_run.py inoculated-vs-control.in inoculated-vs-control.res -b 200 -l 3

lefse_run.py ambient-vs-combined-stress.in ambient-vs-combined-stress.res -b 200 -l 3

3.	Finally, plot the differential abundance 
lefse_plot_res.py --dpi 600 ambient-vs-ozone.res --format pdf --feature_font_size 14 –title ambient-vs-ozone --title_font_size 14 --class_legend_font_size 20 --max_feature_len 200 --right_space 0.1 --left_space 0.4 --width 14

lefse_plot_res.py --dpi 600 inoculated-vs-control.res --format pdf --feature_font_size 14 –title inoculated-vs-control.res --title_font_size 14 --class_legend_font_size 20 --max_feature_len 200 --right_space 0.1 --left_space 0.4 --width 14

lefse_plot_res.py --dpi 600 ambient-vs-combined-stress.res --format pdf --feature_font_size 14 –title ambient-vs-combined-stress.res --title_font_size 14 --class_legend_font_size 20 --max_feature_len 200 --right_space 0.1 --left_space 0.4 --width 14

