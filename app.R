

# Define UI ----
ui <- fluidPage(
  titlePanel('Pipeline for TMT data analysis'),
  sidebarLayout(
    sidebarPanel(h3("Files"),
      # Copy the line below to make a file upload manager
      textInput("text", label = "Working directory", value = "c:/path/to/main/dir"),
      fileInput("file", label = "Metadata file"),
      fileInput("file2", label = "MaxQuant output file"),
      fileInput("file3", label = "Comparisons file"),
      textInput("text2", label = "Experiment name", value = ""),
      radioButtons("radio", label = "PTM?", choices = c("None","P","U")),
      radioButtons('radio2',label="p- or q-value",choices = c('q','p')),
      numericInput("num", label = "cutoff", value = 0.1),
      actionButton("action", label = "Run")
    ),
    mainPanel(h3("Instructions"),
              p("This pipeline processes TMT proteomics data using the following steps:"),
              p("1. Proteins that are not present in at least 2 of the runs are eliminated"),
              p("1B. Imputes reference samples (only if a reference is not included)"),
              p("1C. Separates PTM sites by multiplicity (only if PTM analysis is included)"),
              p("2. Performs sample loading normalization"),
              p("3. Performs internal reference normalization"),
              p("After data are processed, this pipeline performs differential expression analysis on all
                pairwise comparisons using PoissonSeq and saves results."),
              h3("Input files"),
              strong("Working directory"),
              p("Name of your working directory which MUST contain your metadata and MaxQuant output files."),
              strong("Metadata file:"),
              p("text (.txt) file that contains the following information in tab-delimited format:"),
              p("Column 1 should be named \"sample\" and contain the lane # for each sample (1,2,3, etc)"),
              p("Column 2 should be named \"run\" and contain the run number (1,2,3,etc)"),
              p("Column 3 should be named \"rep\" and contain the replicate (1,2,3,etc)"),
              p("Column 4 should be called \"name\" and contain the names of your samples (text). NOTE that your
                reference names must contain the text \"ref\" (capitalization does not matter)."),
              p("An example metadata file can be downloaded in the sidebar. You then can fill out this file with your
                information, save it, and upload it to the pipeline."),
              strong("MaxQuant output file:"),
              p("Output file from MaxQuant that has been converted to comma separated (.csv) format using
                Notepad, Microsoft Excel, etc."),
              strong('Comparisons file'),
              p('List of pairwise comparisons you would like to make. The names of your samples MUST match what
                is provided in the metadata file, and must be named A_vs_B. Please see the TEST files for an example.'),
              strong("Experiment name"),
              p("Input what you called your experiment when loading into MaxQuant. If you cannot remember,
                you can check the summary file."),
              strong("PTM"),
              p("Select P for phosphorylation, U for ubiquitination experiments."),
              strong("cutoff for differential expression"),
              p("Indicate the p- or q-value cutoff you would like to use for differential expression analysis. We recommend
                starting with a q-value of 0.1 (or p-value of 0.05) and then adjusting as you see fit. p- or q-value histograms are included in the
                output files to help you choose a more suitable cutoff if desired."),
              h3("Output files"),
              strong("All output files are saved to your working directory."),br(),br(),
              p("Various comma separated (.csv) files that contain your expression data before and after each step of
                normalization. All differential expression analysis is performed on the values in the file \"IRS_Normalized_values.csv\""),
              p("Various Quality Control (QC) plots including a boxplot colored by run, hierarchical clustering, and
                Principal Component Analysis (PCA) plot that will show the differences between your samples. You should
                see that your samples group by biological replicate after our normalization process."),
              p("For each pairwise comparison we include a volcano plot and a q-value histogram. 
                The volcano plot shows the q-value and fold change distribution for each pairwise comparison. The q-value histogram shows the
                q-value distribution in more detail and is good for choosing a more suitable cutoff if the default cutoff
                does not work for your samples."),
              p("Two Excel files: the first contains results for all pariwise comparisons, the second (which includes
                your q-value cutoff in the file name) lists all the differentially expressed protein groups (or PTMs) at the q-value
                cutoff provided.")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  options(shiny.maxRequestSize=100*1024^2) 
  observe({
    if (input$action > 0){
      TMT_pseq_pipeline(workdir=input$text,datafile=input$file2$datapath,metadatafile=input$file$datapath,
                        exp=input$text2,PTM=input$radio,stat=input$radio2,qval=input$num,compsfile=input$file3$datapath)
    }
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)