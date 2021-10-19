# TMT-NEAT
Tandem Mass Tag Normalization, Expression Analysis, and statistical Testing

[![DOI](https://zenodo.org/badge/232925706.svg)](https://zenodo.org/badge/latestdoi/232925706)

# Install instructions
Download all R code to your computer, then run the file "RUN_TMT.R" to start the RShiny app. Make sure to change your working directory to the folder that contains the code before running this file.

# Known issues
1) This method will automatically write over output files with the same name in your working directory. If you would like to compare your results, make sure to move them from your working directory or change your directory to prevent them being written over.
2) Certain parameters were deprecated in version > 1.4 of the EnhancedVolcano package. If you are a new user of TMT-NEAT, or if you have updated the EnhancedVolcano package to version >1.4, please make sure to use Version 1.5.1 and beyond.

# Test data
Test data are included in the TEST.zip folder. The metadata file is sampledata.txt, the MaxQuant output file is proteinGroups.csv, and the comparison file is comps.xlsx. PTM should be set to "None". The experiment name is "rep". We recommend comparing results using a p-value of 0.05 or a q-value of 0.1.

# Version History

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
