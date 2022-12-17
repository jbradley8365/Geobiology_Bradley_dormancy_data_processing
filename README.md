# Bradley, et al. (2022) Active and dormant microbes - data processing code repository.

This is the repository for the manuscript [Active and dormant microorganisms on glacier surfaces](https://onlinelibrary.wiley.com/doi/full/10.1111/gbi.12535). The raw sequencing data can be obtained from the Sequence Read Archive at [NCBI](https://www.ncbi.nlm.nih.gov/) under:

BioProject PRJNA657180<br/>
IS19-13 DNA: SAMN15830840<br/>
IS19-14 DNA: SAMN15830841<br/>
IS19-13 RNA: SAMN26570441<br/>
IS19-14 RNA: SAMN26570444

BioProject PRJNA893676<br/>
MIT5 DNA: SAMN31525315, SAMN31525316<br/>
MIT5 RNA: SAMN31525317

Within the Seq_analysis folder are files used to process the raw Illumina paired-end metagenomic and metatranscriptomic data for the above publication. 

These include:
* __1_Phyloflash_workflow_Geobio_manuscript.sh__ - using the tool [phyloFlash](http://hrgv.github.io/phyloFlash/] to construct the SSU rRNAs (for 16S and 18S) for phylogenetic inference from the raw paired-end sequencing data.
*  __2_NTU_abundance_to_Phyloseq_Geobio_manuscript.R__ to convert phyloFlash output data in [phyloseq](https://joey711.github.io/phyloseq/) objects in R for exoploration and visualization.
    * output: __NTU_combined.csv__ - the combined output of the multiple NTU_abundance files from the phyloFlash workflow.
    * output: __Fixed_ps.rds__ - the output Phyloseq object with taxonomy names fixed using [microViz](https://david-barnett.github.io/microViz/index.html)
    * output: __phyloseq_biom.csv__ - A CSV created from the phyloseq object that we used in the last script.

* __3_phyloseq_to_csvs_for_plots_Geobio_manuscript.R__  - Code using the Tidyverse set of packages to modify our dataset for calculation of relative abundances and production of plots.

We offer these scripts and in order to be as transparent as possible with our sequence data processing. Please note that these are not intended to be run as stand alone Shell or R scripts and instead are meant to be used as a guided walkthrough. Thank you!