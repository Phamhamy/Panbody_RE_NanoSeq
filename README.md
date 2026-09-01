# A comprehensive survey of somatic mutation rates and mutational signatures in normal human cells

For any queries or support requests, please raise an issue on the GitHub page or email My H. Pham (Mimy) via mp29@sanger.ac.uk

## RE-NanoSeq sequencing data 

Sequencing data are available for access via the European Genome-Phenome Archive (EGA) under accession number EGAD0000101607.

## Data Processing

Variant calling and filtering were processed using the NasoSeq pipeline, available at https://github.com/cancerit/NanoSeq. 

## Mutation burden analysis

SBS, INDEL, and DBS mutation burden for each sample was reported in Supplementary Table 1. 

The Figure was generated using `SBS_ID_DBS_mutation_burden.R`

## Mutation rate analysis 

We investigated the correlation between mutation burden and age for each cell type, using data from Supplementary Table 1. 

Code for regression models can be found in `Regression models test.R`

SBS mutation rates were calculated, and Figure 2 was generated using `SBS_mutation_rates.R`

INDEL mutation rates were calculated, and Extended Figure 2 was generated using `INDEL_mutation_rates.R`

## Mutational signatures analysis

96-context SBS and 83-context INDEL mutational profiles can be found in Supplementary Table 1

Mutational signature extraction tools can be installed by following: 

https://github.com/nicolaroberts/hdp for HDP 

https://github.com/SigProfilerSuite/SigProfilerExtractor for SigProfilerExtractor

https://github.com/parklab/MuSiCal for MuSiCal 

Mutational signature attribution tools can be installed via https://github.com/SigProfilerSuite/SigProfilerAssignment

Attribution of mutational signatures with bootstrapping was implemented using `manual_sigprofiler_bootstrapping.py`

Figure 3 and Extended Figure 3 were generated using `SBS_sigs_attribution_figures.R` and `ID_sigs_attribution_figures.R`






