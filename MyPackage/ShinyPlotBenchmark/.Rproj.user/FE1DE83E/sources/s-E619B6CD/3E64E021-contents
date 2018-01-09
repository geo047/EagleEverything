
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(rbokeh)

shinyUI(fluidPage(
  
  
  fileInput("fi_file", "Select a file"),
  
  # Application title
  titlePanel("Read timing file"),
  sidebarPanel(

  
  shinyFilesButton('files', label='Load Dataset', title='Please select a dataset', multiple=FALSE)
  # shinyFilesButton('file', 'File select', 'Please select a file', FALSE)
  # Sidebar with a slider input for number of bins
  ),
 
  mainPanel(
      rbokehOutput("rbokeh", width = 500, height = 540),
      textOutput('filepaths'),
      verbatimTextOutput('fi_text')
    )
  
))
