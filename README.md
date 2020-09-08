# Variant-interpreter
Variant Interpreter R script for Hipathia

This script obtain the expression files needed to perform variant interpreter analysis in Hipathia tool, and to compute pathway activation values in order to evaluate the impact that a loss of function mutation in one (or several) genes would have over signaling pathways.

The user has to indicate the path to where gtex files are located, this folder (path_to_GTEx) must contain the following RData files:

**gtexInHipathia.rda** Contains the gene expression data from tissues available in GTEx, but only from the genes used by Hipathia to compute pathway activation values in order to reduce the time to compute.

**tissue_gtex.rda** Contains the header of gtexInHipathia.rda, that corresponds to a list of the tissues from GTEX and samples.

**allgenes.rda** A list of Hipathia genes and corresponding entrez.


The user has to indicate a gene or a list of genes to evaluate, and finally a tissue or a list of tissues in which the impact has to be calculated.

This script will save a RData named Matrix.rda, that will serve as input for regular Hipathia analysis, but simulating a loss of function of indicated genes over desired tissues.

There is a version of the script adapted to be run in parallel and to be inlcuded in Hipathia-web tool.
