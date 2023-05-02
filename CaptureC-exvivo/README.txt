#############################################################################
##Data alignment and filtering was performed with Hicup as described in M&M. 
#############################################################################

##Briefly, merged fastq files were aligned to in silico DpnII-digested mm10, using the standard settings with minimum fragment set to 100 and maximum to 800.

##Digestion:
hicup_digester --genome mm10 --re1 ^GATC,DpnII <pathTo/mm10.fa>

##alinment/filtering were performed from the hicup.conf file
hicup --config *.conf --nofill

#############################################################################
##Quantification and statistical analysis was performed using Chicago, 
##as described in M&M
#############################################################################

##Paired bam files were converted to Chicago-accessible objects using the script bam2chicago_adapted_for_rep2.sh
markdowns/bam2chicago_adapted_for_rep2.sh
designDir=Raw_data
${script} <my.bam> ${designDir}*.baitmap ${designDir}*.rmap <my.bam>_adapted_for_rep2_dir

##Chicago files were prepared using .rmap and .baitmap tables
python markdowns/makeNPerBinFile.py --binsize 20000 --designDir Raw_data --minFragLen=100
python markdowns/makeProxOEFile.py --binsize 20000 --designDir Raw_data --minFragLen=100
python markdowns/makeNBaitsPerBinFile.py --binsize 20000 --designDir Raw_data --minFragLen=100
