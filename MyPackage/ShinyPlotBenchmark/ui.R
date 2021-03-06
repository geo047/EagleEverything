
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)
library(rbokeh)
# library(htmlwidgets)

shinyUI(fluidPage(
  

  # Application title
  titlePanel("Read timing file"),
  sidebarPanel(
    shinyFilesButton('files', label='Load Dataset', title='Please select a dataset', multiple=FALSE),
    numericInput("itnum", label = h4("Subset by iteration number"), value = 2, min = 1),
    numericInput("gpunum", label = h4("Subset by number of GPUs"), value = 0, min = 0),
    #checkboxInput("funct_or_box_cb", label = "Plot total time box", value = FALSE)
    radioButtons("plottype_rb", label = h4("Plot type"), choices = list("function time" = 1, "total time" = 2),  selected = 1)
  ),
  
 
  mainPanel(
    textInput(inputId="dataset_path", label="dataset:",  value = "", width="1000px"),
    rbokehOutput("rbokeh", width = 800, height = 600),
    verbatimTextOutput('filepaths')
  )
  
))
