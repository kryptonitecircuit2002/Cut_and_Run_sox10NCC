# CUT&RUN analysis of histone variant occupancy and chromatin-state dynamics during neural crest development

This repository contains the computational workflows and downstream analysis scripts used to process and analyze CUT&RUN datasets generated from zebrafish neural crest cells across development. The analyses were performed to investigate the genomic occupancy and developmental dynamics of the histone variants H2A.Z and H3.3 and their relationship with chromatin states defined by H3K4me3, H3K27ac, and H3K27me3.

Raw sequencing data were processed using the nf-core/CUT&RUN Nextflow pipeline. Downstream analyses include peak annotation, chromatin-state classification, histone variant occupancy dynamics, variant–chromatin state integration, motif enrichment, gene ontology analysis, and data visualization.

## Software
Nextflow
nf-core/cutandrun
Bowtie2
MACS2
BEDTools
deepTools
HOMER
R

## Major Analyses
CUT&RUN preprocessing and spike-in normalization
Histone variant occupancy dynamics
Double-variant nucleosome analysis
Chromatin-state classification
Variant-associated poised chromatin analysis
Motif enrichment
Gene ontology enrichment

## Reproducibility

All analyses were performed using scripts contained in this repository. Software versions and key parameters are documented within individual workflows and scripts.
