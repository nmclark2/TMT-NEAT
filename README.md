# TMT-NEAT
Tandem Mass Tag Normalization, Expression Analysis, and statistical Testing

# Install instructions
Download all R code to your computer, then run the file "RUN_TMT.R" to start the RShiny app. Make sure to change your working directory to the folder that contains the code before running this file.

# Known issues
This method will automatically write over files with the same name in your working directory. If you would like to compare your results, make sure to move them from your working directory or change your directory to prevent them being written over.

# Test data
Test data are included in the TEST.zip folder. You can run these with the standard settings in the app. The metadata file is sampledata.txt, the MaxQuant output file is proteinGroups.csv, and the comparison file is comps.xlsx. PTM should be set to "None". The experiment name is "rep". We recommend comparing results using a p-value of 0.05 or a q-value of 0.1.

# Version History

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
