# TMT-Analysis-Pipeline
Tandem Mass Tag (TMT) Analysis Pipeline

# Install instructions
Download all R code to your computer, then run the file "RUN_TMT.R" to start the RShiny app. Make sure to change your working directory to the folder that contains the code before running this file.

# Test data
Test data are included in the TEST.zip folder. You can run these with the standard settings in the app. The metadata file is sample.txt and the MaxQuant output file is proteinGroups.csv. The experiment name is "rep"

# Version History

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
