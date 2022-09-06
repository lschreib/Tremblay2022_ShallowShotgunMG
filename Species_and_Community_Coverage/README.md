# Analysis of species and community coverage based on rpoB gene count tables

The provided R script was used to calculate the species and community coverage of the agricultural soil data set as described in the manuscript: Tremblay,J., Schreiber,L., Greer, C.W. 2022. "High resolution shotgun metagenonimcs: the more data the better?"

The analyses were carried out using the attached R packages:
* ggplot2 3.3.5
* entropart 1.6-11
* vegan 2.5-7
* lattice_0.20-44
* permute 0.9-5
* readr 1.4.0

Instructions on how to carry out the analysis are provided within the R script and it is strongly recommended to run the script interactively in R Studio.

The rpoB count tables ("ct_rpoB_..._.RData") of the agricultural soil dataset of Tremblay et al. (2022) are provided as an example dataset. The data of sequencing depth of the complete dataset ("Agricultural_soil_sequencing_depth.all_clusters.tsv") was retrieved from NCBI's Short Read Archive (SRA) for the Bioproject PRJNA448773.
