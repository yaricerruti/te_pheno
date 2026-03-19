**Code used for the paper "Sequence composition and conservation predict the phenotypic relevance of transposable elements" 
by Yari Cerruti, Daniela Fusco, Mahdis Jahanbin, Davide Marnetto, Paolo Provero**

The entire pipeline is organized into the following snakefiles. Scripts called by these snakefiles are available in the src directory. For convenience the data used for the linear models and the complete linear model results are available in the data directory

# Snakefile_prep: prepare predictors and outcomes for linear models

The pipeline can be run with the command

snakemake -s Snakefile_prep

after specifying the path to the data files at line 17 and that to the scripts (i.e. the src directory of the repository) at line 18

The pipelines needs the following external data files:

- encodeCcreCombined.bb can be obtained from https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb
- rmsk_hg38.txt can be obtained from https://hgdownload.gi.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
- human_genomic_age.bed can be obtained from https://zenodo.org/records/19089267
- cactus241way.phyloP.bw can be obtained from https://hgdownload.soe.ucsc.edu/goldenpath/hg38/cactus241way/cactus241way.phyloP.bw
- DR_score_Stefansson_Nature_2022.bed can be obtained from https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-04965-x/MediaObjects/41586_2022_4965_MOESM3_ESM.gz
- Dfam-RepeatMasker.lib can be obtained from https://www.dfam.org/releases/current/families/Dfam-RepeatMasker.lib.gz

## Snakefile_egene: compute eGenes associated to TEs

The pipeline can be run with the command

snakemake -s Snakefile_prep

after specifying the path to the data files at line 17 and that to the scripts (i.e. the src directory of the repository) at line 18

The pipelines needs the following external data files:



## Snakefile_models: add GC content and eGene data and fit linear models

## Snakefile_tf: ChIP-Atlas enrichment

