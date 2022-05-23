#run all code in this file to start the TMT analysis pipeline

#check that all packages are installed, and load them using pacman
if (!require('pacman', character.only=T, quietly=T)) {
  install.packages('pacman')
  library('pacman', character.only=T)
}else{
  library('pacman',character.only=T)
}

p_load(openxlsx)
p_load(plyr)
p_load(dplyr)
p_load(shiny)
p_load(ggbiplot)
p_load(EnhancedVolcano)
p_load(shinyscreenshot)

#install PoissonSeq from Github if needed
if (!require('PoissonSeq',quietly=T)) {
  p_load(devtools)
  install_github("cran/PoissonSeq")
  library('PoissonSeq')
}else{
  library('PoissonSeq')
}

source("TMT_pseq_pipeline.R")
runApp('app.R')