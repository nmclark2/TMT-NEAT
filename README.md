# TMT-NEAT
Tandem Mass Tag Normalization, Expression Analysis, and statistical Testing

[![DOI](https://zenodo.org/badge/232925706.svg)](https://zenodo.org/badge/latestdoi/232925706)

# Install instructions
Download all R code to your computer, then run the file "RUN_TMT.R" to start the RShiny app. Make sure to change your working directory to the folder that contains the code before running this file.

# Known issues
1) This method will automatically write over output files with the same name in your working directory. If you would like to compare your results, make sure to move them from your working directory or change your directory to prevent them being written over.
2) Certain parameters were deprecated in version > 1.4 of the EnhancedVolcano package. If you are a new user of TMT-NEAT, or if you have updated the EnhancedVolcano package to version >1.4, please make sure to use Version 1.5.1 and beyond.
3) Please be aware that your samples are read into TMT-NEAT based on their order in the MaxQuant output table. This means that, if you do not use leading zeroes in your sample names, your samples may be in the incorrect order. Please ensure that your samples are in the correct order in the MaxQuant output table before using TMT-NEAT.

# Test data
Test data are included in the TEST.zip folder. These data are published in Zander et al, 2020, Nature Communications: 
https://doi.org/10.1038/s41477-020-0605-7

- Two separate analyses can be performed: one on the protein abundance data (proteinGroups.csv), and one for the phosphosite data (Phospho (STY)Sites.csv). The "sampledata.txt" file is used as the Metadata file for both analyses.
- For the protein abundance data, Experiment Name is "ProtAbun" and PTM is set to "None."
- For the phosphosite data, Experiment Name is "Phospho" and PTM is set to "P."
- Differential expression analysis may be performed using the "comps.xlsx" file with q-value < 0.1. 
- We include screenshots for each analysis to help facilitate parameter selection.

# Version History

# Version 1.6 - May 23, 2022
- New, published test data from Zander et al, 2020, Nature Communications are now included. The previous test data have been removed. By incorporating published, citable test data, we hope to improve reproducibility of results as more features are added in the future.
- Pacman is now used to install and load packages, with the exception of PoissonSeq, which is still directly installed from the Github repository as detailed in the Version 1.5.3 update.
- There is now a screenshot button which can be used to save an image of the app window.
- Differential expression is now optional. To run TMT-NEAT processing without differential expression, set the "Differential Expression?" button to "No." If not performing differential expression, the "comps.xlsx" file is not required.

# Version 1.5.3 - March 8, 2022
- PoissonSeq is now loaded and installed from its Github repository. It was archived on 03-07-2022 by CRAN. The PoissonSeq analysis and code still run as intended, but has not been updated since October 8, 2012, likely leading to this archival. The Github repository should be permanent and resolve these issues.

# Version 1.5.2 - January 14, 2022
- Fixes an issue with the Sample Loading Normalization parameter where SLN would always be skipped even if the user selected "Yes". The parameter works as intended now.

# Version 1.5.1 - October 19, 2021
- Fixes volcano plotting parameters for compatibility with EnhancedVolcano version > 1.4.

# Version 1.5 - October 6, 2021
- Adds option to skip sample loading normalization for samples where enrichment is expected, such as kinase assays, co-IPs, TurboID, etc.

- Clarifies what the experiment name parameter is used for, and allows users to proceed without an experiment name.

- Changes how contaminant sequences are removed: now based on the majority sequence rather than all detected sequences.

# Version 1.4 - August 23, 2021

- Same as Version 1.3, but archived for publication in Nature Communications: https://doi.org/10.1038/s41467-021-26165-3

# Version 1.3 - August 3, 2020

- Allows user to choose whether using p- or q-values for statistics.

- Pairwise comparisons are now set using a user-supplied Excel file. Please see the TEST data for an example.

# Version 1.2 - July 20, 2020

- Altered IRS code to allow for missing values between runs. These values are automatically removed during the differential expression analysis but are useful to retain in case one wishes to do comparisons with less than the maximum number of biological replicates.

# Version 1.1 - April 10, 2020 

- Fixes issues installing ggbiplot and EnhancedVolcano packages

- Automatically removes missing channels based on metadata file

- Fixes issue with assigning FC to wrong genes in summary table

- Adds ubiquitination as a possible PTM

- Fixes issues with setting working directory.

- Adds options for multiple PTMs to app.


# Version 1.0 - January 9, 2020 (first stable build)
