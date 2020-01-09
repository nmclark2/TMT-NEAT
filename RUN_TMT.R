#run all code in this file to start the TMT analysis pipeline

library(shiny)
setwd("D:/OneDrive/Documents/Walley Lab Postdoc/TMT analysis") #change to your own directory that contains the pipeline!
source("TMT_pseq_pipeline.R")
runApp('app.R')