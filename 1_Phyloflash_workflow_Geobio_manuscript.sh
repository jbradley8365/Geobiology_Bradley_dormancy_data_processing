#!/bin/bash

# Phyloflash workflow for totalRNA for BONCAT paper


# First make sure appropriate tools are installed in your conda environment. This workflow will use:
# phyloFlash (https://github.com/HRGV/phyloFlash)
# illumina-utils (https://github.com/merenlab/illumina-utils)
# FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)




# Copy raw metaT files (R1 and R2) to Raw data folder
mkdir Raw # made inside working folder


# If no QC has been done on raw reads do this first -


### QC ###

# Running FastQC (v0.11.9) from /Working/BONCAT

# Make QC folder for results
mkdir 01_Pre_QC


# Move to conda environment for FastQC and MultiQC
conda activate multiqc

# Run for all files in a folder
fastqc Raw/*.fastq.gz --outdir 01_Pre_QC --threads 16

# Visualize QC results using MultiQC together
multiqc 01_Pre_QC/ -o 01_Pre_QC/multiqc_report

# Create Filtered folder
mkdir 02_Filtered


# Generate a .txt file with sample names and locations. Do this from inside the folder with the raw data.

# To generate file path text files
ls *R1* > ../R1s.txt # for R1s
ls *R2* > ../R2s.txt # for R2s

for f in *_R1_001.fastq.gz;
    do SAMPLE=`basename ${f%%_S*_*_R1_001*}`;                 echo "$SAMPLE"; 
    done; # for sample names

# Download (if you're on a shared server) the file path text files, open both, add the headings "sample", "r1", and "r2", and merge R1 and R2 file paths into one file. Save as .txt then reupload (i.e. 'sample_names_for_filtering.txt').


# Create config files before running the filtering program.
# HAVE TO RUN THIS FROM THE DATA FOLDER OR IT WON'T CAPTURE THE CORRECT DATA INPUT FOLDER.
iu-gen-configs ../sample_names_for_filtering.txt -o ../02_Filtered/


# If you want to filter one sample at a time, use this command -
# iu-filter-quality-minoche Filtered/<sample_name>.ini 

# OR

# Instead, to filter multiple samples using a for loop (run from your working folder)
for ini in 02_Filtered/*.ini; 
    do iu-filter-quality-minoche $ini; 
done

# Check sample stats after running and output to text file (besides looking at individual stat files)
grep -E 'number of pairs analyzed|total pairs passed' ../02_Filtered/*STATS.txt > ../All_passed_stats.txt



### QC Post Filtering ###

# Make QC folder for results
mkdir 03_Post_QC

# fastqc Data/10-F-DNA_S13_L001_R1_001.fastq.gz Data/10-F-DNA_S13_L001_R2_001.fastq.gz --outdir /home/ctrivedi/Working/Meta/Genomics/01_QC --threads 16

# Run for all files in a folder
fastqc 02_Filtered/*.fastq --outdir 03_Post_QC --threads 16


# Visualize QC results using MultiQC together
multiqc 03_Post_QC/ -o 03_Post_QC/multiqc_report

# Check sample stats after running and output to text file (besides looking at individual stat files)
grep -E 'number of pairs analyzed|total pairs passed' 02_Filtered/*STATS.txt > All_passed_stats.txt



### Time for phyloFlash ###

# Create a new sample list for phyloflash (from inside the filtered data folder)
cd 02_Filtered/

for f in *_R1.fastq;         
do SAMPLE=`basename ${f%%-QUALITY*}`;                 
echo "$SAMPLE"; 
done;
# copy this output and paste this inside nano and save as 'pf_names.txt' (make sure to save in working directory)


mkdir 04_phyloflash
cd 04_phyloflash
conda deactivate
conda activate pf


# Because we have different read lengths due to different sequencing machines (MiSeq and NextSeq)

# MIT6
/home/ctrivedi/miniconda3/envs/pf/bin/phyloFlash.pl -lib GR19-MIT6 -read1 ../02_Filtered/GR19-MIT6-QUALITY_PASSED_R1.fastq -read2 ../02_Filtered/GR19-MIT6-QUALITY_PASSED_R2.fastq -CPUs 16 -readlength 150 -taxlevel 20 -id 80 -tophit -zip

# IS19-13
/home/ctrivedi/miniconda3/envs/pf/bin/phyloFlash.pl -lib IS19-13 -read1 ../02_Filtered/IS19-13-QUALITY_PASSED_R1.fastq -read2 ../02_Filtered/IS19-13-QUALITY_PASSED_R2.fastq -CPUs 16 -readlength 250 -taxlevel 20 -id 80 -tophit -zip

# IS19-14
/home/ctrivedi/miniconda3/envs/pf/bin/phyloFlash.pl -lib IS19-14 -read1 ../02_Filtered/IS19-14-QUALITY_PASSED_R1.fastq -read2 ../02_Filtered/IS19-14-QUALITY_PASSED_R2.fastq -CPUs 16 -readlength 250 -taxlevel 20 -id 80 -tophit -zip

# Have a look at the output
/home/ctrivedi/miniconda3/envs/pf/bin/phyloFlash_compare.pl --allzip --task barplot,heatmap


# Unzip all the tarballs after!
for a in *.tar.gz
do
    a_dir=`expr $a : '\(.*\).tar.gz'`
    mkdir $a_dir 2>/dev/null
    tar -xvzf $a -C $a_dir
done


# Find NTU files and copy to new directory
mkdir phyloflash_NTUabundance_files

cp */*phyloFlash.NTUabundance.csv phyloflash_NTUabundance_files/ # From within 04_phyloflash/


# Move analysis into R

