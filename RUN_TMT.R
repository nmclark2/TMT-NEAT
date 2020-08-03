#run all code in this file to start the TMT analysis pipeline

#check that all packages are installed, and load them
for (package in c('openxlsx', 'plyr','dplyr','PoissonSeq','shiny')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }else{
    library(package,character.only=T)
  }
}

if (!require('ggbiplot',quietly=T)) {
  install.packages('devtools')
  library(devtools)
  install_github("vqv/ggbiplot")
  library('ggbiplot')
}else{
  library('ggbiplot')
}

if (!require('EnhancedVolcano',quietly=T)) {
  install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
  library('EnhancedVolcano')
}else{
  library('EnhancedVolcano')
}

source("TMT_pseq_pipeline.R")
runApp('app.R')