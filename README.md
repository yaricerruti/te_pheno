**Code used for the paper "Sequence composition and conservation predict the phenotypic relevance of transposable elements" 
by Yari Cerruti, Daniela Fusco, Mahdis Jahanbin, Davide Marnetto, Paolo Provero**

The entire pipeline is organized into the following snakefiles. Scripts called by these snakefiles are available in the src directory. For convenience the data used for the linear models and the complete linear model results are available in the data 
directory

- Snakefile_prep: prepare predictors and outcomes for linear models

- Snakefile_egene: compute eGenes associated to TEs

- Snakefile_models: add GC content and eGene data and fit linear models

- Snakefile_tf: ChIP-Atlas enrichment
