rootdir <-  c('Home' = Sys.getenv("HOME"))
if(.Platform$OS.type == "windows") {
     rootdir <-  c('Home' = paste0(rootdir, "\\..\\"))
}




home_intro <- function(){
  txt <- "
strong(Eagle)  is a software package for genome-wide association mapping.
It differs from most other association mapping packages in that it fits all marker-trait associations simultaneously, an 
returns the 'best' set of snp loci in strongest association with a trait as its findings. 
Eagle can handle data collected from populations of arbitrary structure. The populations can contain inbred or outbred individuals. 
br()
An analysis is performed by reading in the marker data (Read Genotypes), reading in the phenotypic data (Read Phenotypes), reading in the 
marker map if known (Read Map), reading in the Z matrix if needed, and performing the genome-wide analysis (Analyse).  
br()
Help is available by hovering over the widgets or by clicking on the help tab at the top of the screen.  "
    return(txt)
}

read_geno_intro <- function(){
  txt <- "
  Eagle can handle two different types of marker data; genotype data in a plain text space separated file 
  "
  return(txt)
}
read_pheno_intro <- function(){
  txt <- "adfadf"
  
  return(txt)
}  
  
    
  




## ui.R ##
library(shiny)
library(shinythemes)
library(shinyBS)
library(shinyjs)
library(shinyFiles)

FullPage <- navbarPage(title="Eagle: Genome-wide association mapping",  theme = shinytheme("flatly"),
                #       theme = shinytheme("slate"),
                #       theme = shinytheme("united"),

                      
                       ##----------------------------##
                       ##   Home Page                ##
                       ##----------------------------##
                            tabPanel("Home", icon=icon("fa-home", class="fa fa-home  fa-lg"),
                            tags$head(includeCSS("css.css")),
                            fluidPage(
                              fluidRow(
                                column(12,
                                tags$div(img(src = "images/HomeScreen.jpg", 
                                             style="width: 100% ; height: 100%"))
                                )
                              ) ## end fluidRow
                            
                                
                              
                            ) ## end fluidPage
                                ), ## end tabPanel("Home") 

                       ##-----------------------##
                       ##     Read Genotypes    ##
                       ##-----------------------##
                                 
                      tabPanel("Read Genotypes",  icon=icon("file-o"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),




                            fluidPage(
                              fluidRow(
                                column(12,  {
                                       tags$div(img(src = "images/marker_banner.jpg", 
                                                    style="width: 100% ; height: 100%;"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                            
                                              bsButton(inputId="dummy1", label="Hover here for details", 
                                                    style="warning", size="large", type="action", block=TRUE, 
                                                    icon=icon("question-circle-o")
                                                    )
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                                         fluidRow(
                                           column(12,
                                                  wellPanel(
                                                  radioButtons(inputId="filetype", label=h4("Step 1: Choose file type"), 
                                                               choices=c("PLINK"="plink","Text/ASCII"="text" )),
                                                  style="padding: 1px",
                                                  bsTooltip("filetype",
title='<font size="3" > click on file type </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  ),  ## wellPanel
                                           
                                           conditionalPanel(
                                             condition = "input.filetype == 'text'",
                                             wellPanel(
                                               fluidPage(
                                                 wellPanel(
                                                   h5("Assign marker genotypes to snp genotypes AA, AB, BB, and missing"),
                                                 fluidRow(
                                                  column(4, textInput(inputId="AA",label="AA", value="") ),
                                                  column(4, textInput(inputId="AB",label="AB", value="") ),
                                                  column(4, textInput(inputId="BB",label="BB", value="") ),
                                                  column(4, textInput(inputId="missing",label="missing", value="") ) ,
bsTooltip("AB", 
title='<font size="4" > Only a single value can be entered. If inbreds, leave blank  </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")),
bsTooltip("missing",
title='<font size="3" > Enter genotype code used in file that is to be missing. Leave blank if data contains no missing marker genotypes </font>' , 
placement="right", trigger="hover",
                                                            options=list(container="body"))





                                                ) ## end fluidRow
                                               ) ## end inner wellPanel
                                               )  ## end fluidPage

                                                 ) ## end Wellpanel
                                            
                                           ) ## end conditionalPanel
                                           
                                           
                                           
                                           
                                                  ) ## end column
                                           
                                         ), ## end fluidRow choose file type
                        
                                         fluidRow(
                                           column(12, wellPanel(
                                                  numericInput(inputId="memsize", label=h4("Step 2: Specify available memory in Gbytes"), 
                                                               value=8, min = 2, max = NA, step = NA),
                                                  style="padding: 1px",
                                                  bsTooltip("memsize",
title='<font size="3" > set to maximum available memory in Gbytes  </font>',
placement="right", trigger="hover",
                                                          options=list(container="body"))
                                                  )) ## end column
                                           
                                           
                                         ), ## end fluidRow specify amout of memory
                                        
                                         
                                      fluidRow(column(12, 
                                        wellPanel(
                                        h4("Step 3: Select marker file"),
                                        shinyFilesButton('choose_marker_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_marker_file_text", label = h5("or enter file name (including full path)"))
                                         
                                           
                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow


                                      






 
                                         
                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 4: Upload file"),




                                                          actionButton(inputId="marker_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          bsTooltip("marker_go", 
title='<font size="3" > Upload file. <br> This may take some time if the file is large.  </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))




 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                         
                                         
                                         
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets
                                column(7, 
                                        verbatimTextOutput("ReadMarker", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page
                                   ## for outputs from ReadMarker function
                                
                              ) ## end fluidRow
                              
                              
                              
                              
                              
                            ) ## end fluidPage    
                                
                                
                                
                                
                                
                      ),  ## end tabPanel("Read Genotypes")
                       
                       
                       ##----------------------##
                       ## Read Phenotypes      ##
                       ##----------------------##
                      
                                    
                      tabPanel("Read Phenotypes",  icon=icon("file-o"), 
                               tags$head(tags$style(HTML('
                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/pheno_banner.jpg", 
                                                    style="width: 100% ; height: 100%; "))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                             bsButton(inputId="dummy2", label="Hover here for details",
                                                    style="warning", size="large", type="action", block=TRUE, 
                                                    icon=icon("question-circle-o")
                                                    )
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                                         fluidRow(
                                           column(12,
                                                  wellPanel(
                                                  radioButtons(inputId="pheno_header", label=h4("Step 1: Select if file contains column names"), 
                                                               choices=c("yes"="yes","no"="no" ), selected="yes"),
                                                  style="padding: 1px",
                                                  bsTooltip("pheno_header",
title='<font size="3" > click on yes if the first row of the file contains the column names. Generic names will be assigned if no is clicked.  </font>',
placement="right", 
trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                           
                                           
                                           
                                             ) ## end column
                                           
                                         ), ## end fluidRow choose file type
                        
                                         fluidRow(
                                           column(12, wellPanel(
                                                  radioButtons(inputId="pheno_csv", label=h4("Step 2: Is the file comma separated"),
                                                               choices=c("yes"="yes","no"="no" ), selected="no" ),
                                                  style="padding: 1px",
                                                  bsTooltip("pheno_csv",
title='<font size="3" > click on yes if the file is a csv file. </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column
                                           
                                         ), ## end fluidRow specify amout of memory
                                        

                                         fluidRow(
                                           column(12, wellPanel(
                                                   textInput(inputId="pheno_missing", label=h4("Step 3: Code for missing value", value="") ),
                                                  style="padding: 1px",
                                                  bsTooltip("pheno_missing",
title='<font size="3" > Assign value that denotes a missing value. Leave blank if file does not contain missing data. </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column

                                         ), ## end fluidRow specify amout of memory

                                      fluidRow(column(12,
                                        wellPanel(
                                        h4("Step 4: Select phenotypic file"),
                                        shinyFilesButton('choose_pheno_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_pheno_file_text", label = h5("or enter file name (including full path)"))
                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow




                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 5: Upload file"),




                                                          actionButton(inputId="pheno_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          bsTooltip("pheno_go",
title='<font size="3" > Upload file. <br> This may take some time if the file is large. </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))




 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                         
                                         
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets

                               column(7,
                                        verbatimTextOutput("ReadPheno", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page
                                   ## for outputs from ReadMarker function
                              ) ## end fluidRow
                              
                              
                              
                              
                              
                            ) ## end fluidPage    
                                
                                
                                
                                
                                
                      ),  ## end tabPanel("Read Phenotypes")





#              #-----------------------TESTING ======================================= 
#                                    
#                     tabPanel("Read Phenotypes",  icon=icon("file-o"), 
#                               tags$head(tags$style(HTML('
#                                                         .popover {
#                                                         max-width: 80%;
#                                                         
#                                                         }
#                                                         '))
#                               ),
#
#
#                            fluidPage(
#                              fluidRow(
#                                column(12, {
#                                       tags$div(img(src = "images/pheno_banner.jpg", 
#                                                    style="width: 100% ; height: 100%; "))
#                               
#                                }
#                                       ) ## end column(12, )
#                              ), ## end fluidRow
#                              br(),
#                              fluidRow(column(12, 
#                                             bsButton(inputId="dummy2", label="Hover here for details",
#                                                    style="warning", size="large", type="action", block=TRUE, 
#                                                    icon=icon("question-circle-o")
#                                                    )
#                                           
#                                       ) ## end column
#                              ), ## end fluidRow
#                              
#                              
#                              br(),
#                              fluidRow(
#                                column(5, 
#                                       fluidPage(
#                                         fluidRow(
#                                           column(12,
#                                                  wellPanel(
#                                                  radioButtons(inputId="pheno_header", label=h4("Step 1: Select if file contains column names"), 
#                                                               choices=c("yes"="yes","no"="no" )),
#                                                  style="padding: 1px",
#                                                  bsTooltip("pheno_header",
#title='<font size="3" > click on yes if the first row of the file contains the column names. Generic names will be assigned if no is clicked.  </font>',
#placement="right", 
#trigger="hover",
#                                                            options=list(container="body")
#                                                      )
#                                                  )  ## wellPanel
#                                           
#                                           
#                                           
#                                             ) ## end column
#                                           
#                                         ), ## end fluidRow choose file type
#
#
# fluidRow(column(12,
#                                          wellPanel(
#                                            h4("Step 3: Select marker file"),
#
#                                           actionButton(inputId="choose_marker_file", h6("Choose File")), br(),
#                                           textOutput("choose_marker_file"),
#                                           style='padding: 1px',
#                                           bsTooltip("choose_marker_file",
#title='<font size="3" >WARNING! File browser window may open behind web browser  </font>',
#placement="right",
#trigger="hover",
#                                                     options=list(container="body"))
#
#
#                                          )
#                                         )
#                                         ) ## end fluidRow
#
#
#
#
#
#                        
#)
#)
#))),
#




                      ##-----------------------------------##
                      ## Read Z matrix (if needed)         ##
                      ##-----------------------------------##
                      
                       tabPanel("Read Z matrix (if needed)", icon=icon("file-o"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/Zmat_banner.jpg", 
                                                    style="width: 100% ; height: 100%"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                          bsButton(inputId="Zmat1", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )

  
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                        



                                      fluidRow(column(12,
                                        wellPanel(
                                        h4("Step 1: Select Z matrix file"),
                                        shinyFilesButton('choose_Zmat_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_Zmat_file_text", label = h5("or enter file name (including full path)"))


                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow









 
                                         
                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 2: Upload file"),




                                                          actionButton(inputId="Zmat_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          bsTooltip("Zmat_go",
title='<font size="3" > Upload file.    </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))
 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets

                               column(7,
                                        verbatimTextOutput("ReadZmat", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page
                                   ## for outputs from ReadMarker function

                                
                              ) ## end fluidRow
                              
                            ) ## end fluidPage    
                                
                                
                      ),  ## end tabPanel("Read Z matrix")











                      ##-------------------------##
                      ## Read Marker map         ##
                      ##-------------------------##
                      
                       tabPanel("Read Map (optional)", icon=icon("file-o"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/map_banner.jpg", 
                                                    style="width: 100% ; height: 100%"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                          bsButton(inputId="dummy3", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )

  
                                           
                                       ) ## end column
                              ), ## end fluidRow
                              
                              
                              br(),
                              fluidRow(
                                column(5, 
                                       fluidPage(
                        

                                         fluidRow(
                                           column(12,
                                                  wellPanel(
                                                  radioButtons(inputId="map_header", label=h4("Step 1: Select if file contains column names"),
                                                               choices=c("yes"="yes","no"="no" )),
                                                  style="padding: 1px",
                                                  bsTooltip("map_header",
title='<font size="3" > click on yes if the first row of the file contains the column names. Generic names will be assigned if no is clicked. </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel



                                             ) ## end column

                                         ), ## end fluidRow choose file type


                                         fluidRow(
                                           column(12, wellPanel(
                                                  radioButtons(inputId="map_csv", label=h4("Step 2: Is the file comma separated"),
                                                               choices=c("yes"="yes","no"="no" ), selected="no"),
                                                  style="padding: 1px",
                                                  bsTooltip("map_csv",
title='<font size="3" > click on yes/no  </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column
                                           
                                         ), ## end fluidRow specify amout of memory
                                        


                                      fluidRow(column(12,
                                        wellPanel(
                                        h4("Step 3: Select map file"),
                                        shinyFilesButton('choose_map_file', 'Select File', 'Please select file', FALSE),
                                        textInput("choose_map_file_text", label = h5("or enter file name (including full path)"))
                                          )  ## end wellPannel
                                         )
                                         ), ## end fluidRow



                                       
                                         
                                         fluidRow(column(12, 
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 4: Upload file"),




                                                          actionButton(inputId="map_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          bsTooltip("map_go",
title='<font size="3" > Upload file. <br> This may take some time if the file is large.   </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))
 
                                                        )
                                                  )
                                         ) ## end fluidRow
                                       ) ## end fluidPage -- widgets on left hand side
                                       
                                       
                                       
                                       
                                       ), ## end column(6,  )  -- left half of page
                                          ## for input widgets

                               column(7,
                                        verbatimTextOutput("ReadMap", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page
                                   ## for outputs from ReadMarker function


                                
                              ) ## end fluidRow
                              
                            ) ## end fluidPage    
                                
                                
                      ),  ## end tabPanel("Read Map")








                  ##-------------------------------------##
                  ##   Analysis - Genome-wide Analysis   ##
                  ##-------------------------------------##
                      
                       tabPanel("Analyse", icon=icon("fa-area-chart", class = "fa fa-area-chart fa-lg", lib = "font-awesome"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;
                                                         
                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/analyse_banner.jpg", 
                                                    style="width: 100% ; height: 100%"))
                               
                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12, 
                                          bsButton(inputId="dummy4", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )
                                           
                                       ) ## end column
                              ), ## end fluidRow
                           br(), 
 
                            fluidRow(
                                 column(6, 
                                 fluidPage(
                                
                                    fluidRow(column(12,  
                                      wellPanel(
                                            uiOutput("analyse_names"),

                                           bsTooltip("analyse_names",
title='<font size="5" > Select a single variable to be treated as the trait for the analysis  </font>',
                                         placement="right", trigger="hover", options=list(container="body"))

                                      ) ## end wellPanel
                                          ) ## end column
                                    ), ## end fluidRow                             

                                    fluidRow(column(12, 
                                        wellPanel(
                                            uiOutput("analyse_fnames"),

                                           bsTooltip("analyse_fnames",
title='<font size="5" > Select the variables, if any, to be used as fixed effects in the analysis. If no variables are selected, then only an overall mean will be fitted. </font>',
                                         placement="right", trigger="hover", options=list(container="body")),

                                            textOutput("fmodel")
                                        ) ## end wellPanel
                                     ) ## end column
                                  ), ## end fluidRow


                                   fluidRow(column(12,  wellPanel(
                                         numericInput(inputId="analyse_cpu", label=h4("Step 3: Specify number of cpu"), value=1), 

                                         style="padding: 1px",
                                         bsTooltip("analyse_cpu",
title='<font size="5" > set to the number of cpu available for distributed computing. </font>',
placement="right", trigger="hover",
                                                          options=list(container="body"))
                                    ) ## end wellpanel
                                   )),  ## end column and fluidRow



                                   fluidRow(
                                      column(12,  wellPanel(
                                          h4("Step 4: Additional Options"),
                                          actionButton(inputId="options_go", h6("Click Here")),
                                          conditionalPanel(
                                            condition="input.options_go > 0 ",
                                                 wellPanel( 
                                                    fluidPage(
                                                        fluidRow(column(12, 

                                                        sliderInput(inputId="analyse_maxits", label=h4("Specify maximum number of iterations"),
                                                               value=20, min = 1, max = 40, step = NA),
                                                        style="padding: 1px",
                                                        bsTooltip("analyse_maxits",
title='<font size="5" > set to the maximum number of detectable marker-trait associations. <br> Very rarely will this need to be adjusted. Its a safety feature to prevent analyses taking too long. </font>',
placement="right", trigger="hover",
                                                          options=list(container="body"))


                                                        )), ## end column and fluidRow

                                                        fluidRow(column(12, 
                                                         
                                                 radioButtons(inputId="analyse_quiet", label=h4("Verbose mode"),
                                                               choices=c("No"="no", "Yes"="yes"  )),
                                                  style="padding: 1px",
                                                  bsTooltip("analyse_quiet",
title='<font size="5" > Click yes if detailed output is wanted </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )

                                

                                                         ))  ## end column and fluidRow

                                                    ) ## end fluidPage
                                                 ) ## end wellPanel
                                           ) ## end conditionalPanel

                                                           ) ## end wellPanel
                                      ) ## column
                                   ),  ## fluidRow










                                   fluidRow(column(12,
                                                       wellPanel(
                                                          shinyjs::useShinyjs(),
                                                          h4("Step 5: Perform genome-wide analysis"),




                                                          actionButton(inputId="analyse_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                                       icon=icon("upload", lib="glyphicon")),



                                                          style='padding: 1px',
                                                          bsTooltip("analyse_go",
title='<font size="5" > Click here to find the set of snp in strongest association with the trait  </font>',
placement="right", trigger="hover",
                                                                     options=list(container="body"))





                                                        )
                                                  )
                                         ) ## end fluidRow

                                  ) ## end fluidPage
                             
                                ),  ## end column

                               column(6,
                                        verbatimTextOutput("AM", placeholder=TRUE)
                                )  ## end column(6, ) -- right half of page
                                   ## for outputs from ReadMarker function


                   ) ## end fluidRow 
                              
                            ) ## end fluidPage    
                                
                                
                                
                                
                                
                      ),  ## end tabPanel("Analysis ")
                      
                       tabPanel("Findings", icon=icon("fa-puzzle-piece", class="fa fa-puzzle-piece fa-lg"), 
                               tags$head(tags$style(HTML('

                                                         .popover {
                                                         max-width: 80%;

                                                         }
                                                         '))
                               ),


                            fluidPage(
                              fluidRow(
                                column(12, {
                                       tags$div(img(src = "images/findings_banner.jpg",
                                                    style="width: 100% ; height: 100%"))

                                }
                                       ) ## end column(12, )
                              ), ## end fluidRow
                              br(),
                              fluidRow(column(12,
                                          bsButton(inputId="dummy5", label="Hover here for details",
                                          style="warning", size="large", type="action", block=TRUE,
                                          icon=icon("question-circle-o")
                                          )



                                       ) ## end column
                              ), ## end fluidRow
                           br(), 
                              fluidRow(
                                column(4,
                                       fluidPage(

                                         fluidRow(
                                           column(12, wellPanel(


                                                  h4("Click to calculate additional summary information"),
                                                  actionButton(inputId="pvalue_go",label="", width='35%', style='padding:5px 5px 5px 5px; font-size:180%',
                                                  icon=icon("glyphicon glyphicon-stats", lib="glyphicon")),

                                                  style="padding: 1px",
                                                  bsTooltip("summary",
title='<font size="5" > click on yes/no </font>',
placement="right", trigger="hover",
                                                            options=list(container="body")
                                                      )
                                                  )  ## wellPanel
                                             ) ## end column

                                         ) ## end fluidRow specify amout of memory
                                ) ## ene fluidPage
                         ),  ## end column
                          column(8, 
                              fluidPage(
                               fluidRow(
                                   column(12, 
tags$div(
         HTML(paste( tags$span(style="color: #ad1d28; font-size: 32px", "Marker-trait Associations"), sep = ""))),
                                       tableOutput("findings")
                                  ) ## end column
                               ), ## end fluidRow
                             fluidRow(
                                column(12, 
                                    conditionalPanel(condition="input.pvalue_go > 0", 
tags$div(
         HTML(paste( tags$span(style="color: #ad1d28; font-size: 32px", "Size and Significance of Effects"), sep = ""))),
                                    tableOutput("size")

                                    ) ## end conditionalPanel
                               ) ## end column
                           ), ## end fluidRow



                            fluidRow(
                                column(12, 
                                    conditionalPanel(condition="input.pvalue_go > 0", 
tags$div(
         HTML(paste( tags$span(style="color: #ad1d28; font-size: 32px", "Proportion of Variance Explained as Markers Added to Model"), sep = ""))),
                                    tableOutput("R")
                                    )
                               ) ## end column
                           ) ## end fluidRow


                         ) ## ene fluidPage 
                       ) # end column
                       ) # end fluidRow
                       ) # end fluidPage

),  ## end tabPAnel pvalue


   tabPanel("Help",  icon=icon("question-circle-o", class="fa fa-question-circle-o fa-lg  "),
                fluidPage(
                      fluidRow(
                        column(12, 
                            includeHTML("help.html")
                        )  ## end column
                      ) ## end fluidRow
                 ) ## end fluidPage
            ) ## end tabPanel("Help") 









                       ) ## end navbarPage


## displays eagle log on navbar
FullPage[[3]][[1]]$children[[1]]$children[[1]]$children[[1]] <- 
  tags$img(src = 'images/logo.jpg', width = 80, height = 60)
ui <- FullPage

get_path <- function (defaultpath="/R/library/Eagle/shiny_app/shinydata/genoDemo.dat") {
            path_to_file_res <- tryCatch({
                if(.Platform$OS.type=="unix"){
                    path_to_file_res <- tk_choose.files()
                    #print(path_to_file_res)
                } else {
                    path_to_file_res <- file.choose()                   
                }                
                }, warning = function(war) {
                    print(paste("Eagle::get_path() Warning: ",war))
                    path_to_file_res<-defaultpath
                    return (path_to_file_res)
                }, error = function(err) {
                    print(paste("Eagle::get_path() Error: ",err))
                    path_to_file_res<-defaultpath
                    return (path_to_file_res)
                }, finally = {
                   # path_to_file_res<-"/R/library/Eagle/shiny_app/shinydata/genoDemo.dat"
                  #  return (path_to_file_res)
                }) # END tryCatch
    
            return (path_to_file_res)
      }


server <- function(input, output, session){
  library("Eagle")

#rootdir <-  c('Home' = Sys.getenv("HOME"))
#rootdir <-  c('rootdir'="C:\\", 'Home' = Sys.getenv("HOME"))
#rootdir <- c('roodir'=c(wd="."))
#if(.Platform$OS.type == "windows") {
#     print(" in the window part of root")
#     rootdir <<-  c('rootdir'="C:/", 'Home' = Sys.getenv("HOME"))
#     print(rootdir)
#  } else {
#     rootdir <<- c(rootdir="/", home=Sys.getenv("HOME"))
#  }
  ##------------------------------------------
  ## Intros to pages
  ##-------------------------------------------
 
#  output$home_intro <- renderText(home_intro())
  output$read_geno_intro <- renderText(read_geno_intro())
  output$read_pheno_intro <- renderText(read_pheno_intro())
  
  ##----------------------------------------
  ##  Read marker path and file name
  ##---------------------------------------- 
  ## upload path and file name
    shinyFileChoose(input=input, id='choose_marker_file', session=session, roots=rootdir)
    # path_to_marker_file <- NULL    
    observeEvent(input$choose_marker_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_marker_file)
           updateTextInput(session, "choose_marker_file_text", value =  as.character(inFile$datapath))
           path_to_marker_file  <<- as.character(inFile$datapath)
    })

    observeEvent(input$choose_marker_file_text, {
           path_to_marker_file  <<- as.character(input$choose_marker_file_text)
    })

 



   ## Read marker information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   geno <- NULL
   observeEvent(input$marker_go, {
   withProgress(message = 'Loading marker data', value = 1, {
    
     if(input$filetype == "plink"){
       withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 if (file.exists(path_to_marker_file) == TRUE) {
                   geno <<- ReadMarker(filename = path_to_marker_file, type = "PLINK", availmemGb = input$memsize, quiet = TRUE)
                 } else {
                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", path_to_marker_file))
              }
          }, ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
             })

     }

     if(input$filetype == "text"){
             withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 aa <- input$AA
                 ab <- input$AB
                 bb <- input$BB
                 missing <- input$missing
                 if(input$AA=="")
                     aa <- NULL
                 if(input$AB=="")
                     ab <- NULL
                 if(input$BB=="")
                     bb <- NULL
                 if(input$missing=="")
                     missing <- NULL
                 if (file.exists(path_to_marker_file) == TRUE) {
                 geno <<- ReadMarker(filename = path_to_marker_file, type = "text", AA = aa, 
                            AB = ab  , BB = bb, availmemGb = input$memsize,  quiet = TRUE , missing=missing) 
                } else {
                    shinyjs::html(id = "ReadMarker", html = paste0("ReadMarker", "  File does not exist:", path_to_marker_file))
                 }

              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
             })


     }  ## end if(input$filetype == "text")


  })  ## withProgress

  })  ## end observeEvent






  ##----------------------------------------
  ##  Read phenotypic path and file name
  ##---------------------------------------- 
  ## upload path and file name

#  This doesn't work - don't know why - the new rootdir is not used by shinyFileChoose
#  observeEvent(input$marker_go, {
#      
#        rootdir <<-  c(rootdir="/", rootdir2="/flush1/geo047/")
#        print(" in here ")
#        print(rootdir)
#})


        shinyFileChoose(input=input, id='choose_pheno_file', session=session, roots=rootdir )

        observeEvent(input$choose_pheno_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_pheno_file)
           updateTextInput(session, "choose_pheno_file_text", value =  as.character(inFile$datapath))
           path_to_pheno_file  <- as.character(inFile$datapath)
       })

        observeEvent(input$choose_pheno_file_text, {
           path_to_pheno_file  <<- as.character(input$choose_pheno_file_text)
         })



   ## Read phenotypic  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   pheno <- NULL
   observeEvent(input$pheno_go, {
   withProgress(message = 'Loading phenotypic data', value = 1, {

   header_flag <- FALSE
   if(input$pheno_header == "yes")
      header_flag <- TRUE
   csv_flag <- FALSE
   if(input$pheno_csv == "yes")
      csv_flag <- TRUE

   pheno_missing <- input$pheno_missing
   if(input$pheno_missing=="")
      pheno_missing <- NULL




   withCallingHandlers({
                shinyjs::html("ReadPheno", "")
                 if (file.exists(path_to_pheno_file) == TRUE) {
                 pheno  <<- ReadPheno(filename = path_to_pheno_file, header=header_flag, csv=csv_flag, missing= pheno_missing)
                 } else {
                    shinyjs::html(id = "ReadPheno", html = paste0("ReadPheno", "File does not exist:", path_to_pheno_file))
                 }
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadPheno", html = m$message, add = TRUE)
       })



  })  ## end withProgress

  })  ## end observeEvent





  ##------------------------------##
  ## Read Z matrix                ## 
  ##------------------------------ss


  ##----------------------------------------
  ##  Read Z matrix path and file name
  ##---------------------------------------- 
  ## upload path and file name
        shinyFileChoose(input=input, id='choose_Zmat_file', session=session, roots=rootdir )

        observeEvent(input$choose_Zmat_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_Zmat_file)
           updateTextInput(session, "choose_Zmat_file_text", value =  as.character(inFile$datapath))
           path_to_Zmat_file  <- as.character(inFile$datapath)
       })
        observeEvent(input$choose_Zmat_file_text, {
           path_to_Zmat_file  <<- as.character(input$choose_Zmat_file_text)
         })



   ## Read Zmat  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   Zmat <- NULL
   observeEvent(input$Zmat_go, {
     withProgress(message = 'Loading Z matrix file', value = 1, {




       withCallingHandlers({
                 shinyjs::html("ReadZmat", "")
         
                 if (file.exists(path_to_Zmat_file) == TRUE) {
                 Zmat  <<- ReadZmat(filename = path_to_Zmat_file)
                 } else {
                    shinyjs::html(id = "ReadZmat", html = paste0("ReadZmat", "File does not exist:", path_to_Zmat_file))
                 }
         
                 # Zmat  <<- ReadZmat(filename = path_to_Zmat_file
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadZmat", html = m$message, add = TRUE)
       })


  })


  })  ## end observeEvent











  ##----------------------------------------
  ##  Read map path and file name
  ##---------------------------------------- 
  ## upload path and file name
        shinyFileChoose(input=input, id='choose_map_file', session=session, roots=rootdir )

        observeEvent(input$choose_map_file, {
           inFile <- parseFilePaths(roots=rootdir, input$choose_map_file)
           updateTextInput(session, "choose_map_file_text", value =  as.character(inFile$datapath))
           path_to_map_file  <- as.character(inFile$datapath)
       })

        observeEvent(input$choose_map_file_text, {
           path_to_map_file  <<- as.character(input$choose_map_file_text)
         })
   





   ## Read map  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   map <- NULL
   observeEvent(input$map_go, {
    withProgress(message = 'Loading marker map file', value = 1, {


     csv_flag <- FALSE
     if(input$map_csv == "yes")
        csv_flag <- TRUE


   map_header_flag <- FALSE
   if(input$map_header == "yes")
      map_header_flag <- TRUE

       withCallingHandlers({
                 shinyjs::html("ReadMap", "")
         
                 if (file.exists(path_to_map_file) == TRUE) {
                    map  <<- ReadMap(filename = path_to_map_file, csv=csv_flag, header= map_header_flag)
                 }  else {
                    shinyjs::html(id = "ReadMap", html = paste0("ReadMap", "File does not exist:", path_to_map_file))
                 }
         
               #   map  <<- ReadMap(filename = path_to_map_file, csv=csv_flag, header= map_header_flag)
         
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMap", html = m$message, add = TRUE)
       })




  })
  })  ## end observeEvent







  ##-------------------
  ## Analyse Data
  ##-------------------
  traitn <- NULL
  ## gets column names of pheno file
  nms <- reactive({
     if(input$pheno_go && input$pheno_header == "yes")
        return(names(pheno))
     if(input$pheno_go && input$pheno_header == "no")
     {  ## pheno file is not named
        nms <- paste("V", 1:ncol(pheno), sep="")
        return(nms)
     } 


     })



  ## gets column names length of  pheno file
  sz <- reactive({
     if(input$pheno_go && input$pheno_header == "yes")
     {
        return(length(names(pheno)))
     }
     if(input$pheno_go && input$pheno_header == "no")
     {  ## pheno file is not named
        nms <- paste("V", 1:ncol(pheno), sep="")
        return(length(nms))
     }


     })



  sz <- 0
  output$analyse_names <- renderUI({
   if (length(nms) < 5){
       sz <- length(nms)
   } else {
       sz <- 5
   }

  selectInput(inputId="nmst", label=h4("Step 1: Choose trait"), choices=nms(), size = sz   , selectize=FALSE )     
  })  ## end renderUI


  output$analyse_fnames <- renderUI({
      #nms <- names(pheno)
      checkboxGroupInput("nmsf", h4("Step 2: Choose fixed effects"), nms(), inline=TRUE)
    })  ## end renderUI

  fform <- NULL
   output$fmodel <- renderText({
        fform <<- paste(input$nmsf, collapse="+")
   })


   ## how to get traitn and effectsn from UI for later use ??????
   traitn <- reactive({input$nmst})

 res <- NULL
   observeEvent(input$analyse_go, {
   withProgress(message = 'Analysing data', value = 1, {



       
       withCallingHandlers({
                 shinyjs::html("AM", "")
                 quietvalue <-  TRUE
                 if(input$analyse_quiet == "yes")
                    quietvalue <- FALSE
                 res <<- AM(trait=input$nmst , fformula=fform , availmemGb = input$memsize , 
                            quiet = quietvalue,
                            ncpu = input$analyse_cpu, maxit = input$analyse_maxits , pheno = pheno, geno=geno, map=map) 

              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "AM", html = m$message, add = TRUE)
       })
  })
  })  ## end observeEvent


 ##--------------------------
 ## Print findings .... 
 ##--------------------------

 ## form data frame of results 
 observeEvent(input$analyse_go, {
 dfres <- NULL
 if(length(res$Mrk)>1){
    dfres <- data.frame(snps=res$Mrk[-1], chrm=res$Chr[-1], position=res$Pos[-1])

 }

 if(is.null(dfres)){
    ## no associations found
    output$findings <- renderText("No significant associations between snp and trait were found.")
  }
  if(!is.null(dfres)){
    output$findings <- renderTable(dfres)
  }

     observeEvent(input$pvalue_go, {
   withProgress(message = ' Calculating additional summary measures', value = 1, {



      withCallingHandlers({
                 shinyjs::html("summary", "")
                  sumres <- SummaryAM(AMobj=res, pheno=pheno, geno=geno, map=map)
                  output$pvalue <- renderTable(sumres[["pvalue"]], digits=-1, hover=TRUE, bordered=TRUE)
                  output$size <- renderTable(sumres[["size"]], digits=-1, hover=TRUE, bordered=TRUE)
                  output$R <- renderTable(sumres[["R"]],  hover=TRUE, bordered=TRUE)


              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "summary", html = m$message, add = TRUE)
       })



  })

  })  ## end observeEvent




}) ## end observeEvent


##----------------------
## Help - Docs
##---------------------



## FAQ html 
#   output$faq <- renderUI({
#        HTML(markdown::markdownToHTML(knit('faq.rmd', 
#              quiet = TRUE),options=c('toc'), fragment.only=TRUE))
#    })

# observeEvent(input$quickstart, {
#      RShowDoc("QuickStart", package="Eagle")
#})

  ##---------------------------------
  ## Help files   - addPopover
  ##--------------------------------

 
  addPopover(session, "dummy1", "Details", content = HTML("
Eagle can handle two types of marker genotype file; a space separated plain text file and PLINK
ped file. We assume the marker loci are snps. 
Missing marker genotypes are allowed but the 
proportion of missing genotypes is assumed to be low. 
<br><br>
The marker genotype file should not contain column names. 
We also assume that each row of data in the file corresponds to data on a different 
individual. The ordering of the rows, by individual, must be the same for the marker genotype file and phenotypic file.<br><br> 


If the file is a plain text file, then character or numeric genotypes can be used to 
denote the snp genotypes. However, Eagle needs 
to map these genotypes to its internal snp genotypes. We do this by asking the user to assign their snp genotype codes to our AA, 
AB and BB codes. If data are collected on inbred individuals, only AA and BB need be specified. 
<br> <br>

To load the marker genotype data into Eagle, follow the four steps.  Upon cliking Upload, Eagle checks the genotype file for errors, 
and recodes the genoytpes for later analysis. If the marker genotype is large (many Gbytes), this step can take several minutes. <br><br>
Output from reading in the marker genotype file will appear in the right hand-side panel. 
                                                           "), trigger = 'hover')

  
  addPopover(session, "dummy2", "Details", content = HTML("  
  Eagle assumes the phenotypic file is either a space separated or comma separated file. The rows correspond to data 
  collected on the individuals. The first row of the file can contain column headings or not.  
  The number of rows of data in the phenotypic file must equal the number of rows in the marker genotype file otherwise an error occurs.  
  Also, Eagle assumes the phenotypic data is row ordered by individual in the same way as the marker genotype data. 
 <br> <br>
Data on multiple traits and fixed effects that may or may not be used in the analysis can be included in this file.  <br> <br>
Missing values are allowed.  <br> <br>
Output from reading in the phenotypic file will appear in the right hand-side panel. 
  "), trigger = 'hover')


addPopover(session, "Zmat1", "Details", content = HTML("
The Z matrix contains only zeros and ones. The number of rows must be greater than the number of columns. 
It is used for those situations where multiple observations of the same trait have been recorded for an individual. 
"), trigger = "hover") 





  
addPopover(session, "dummy3", "Details", content = HTML("
    Eagle does not
     require a known marker map in order to analyse the data.  
     If a map file is read into Eagle, then the
     marker names are used when results are reported in  'Findings'. If a
     map file is not supplied, generic names M1, M2, ..., are
     given to the marker loci. 
      <br> <br>
     The map file can have three or four columns. If the
     map file has three columns, then it is assmed that the three
     columns are the marker locus names, the chromosome number, and the
     map position (in any units). If the map file has four columns as
     with a PLINK map file, then the columns are assumed to be the
     marker locus names, the chromosome number, the map position in
     centimorgans, and the map position in base pairs.
      <br> <br>
     Missing values are not allowed.
      <br> <br>
    The order of the marker loci in this file is assumed to be in the
     same order as the loci (or columns) in the marker data file.
  "), trigger = "hover") 



addPopover(session, "dummy4", "Details", content = HTML(paste("
    Here,  multiple_locus association mapping is performed. The analysis
     simultaneously accounts for  familial
     relatedness and nuisance fixed effects while detecting 
     multiple marker-trait associations. Unlike other association mapping methods, 
     there are no regularization parameters to be tuned, nor significance thresholds to be set. 
     <br><br>
     Output from performing the analysis is printed to the right hand panel. A table of results is printed in 'Findings'. 
     <br><br>", tags$span(style="color:red",
     "Once an analysis has been completed, a new analysis can be performed  
     by selecting a new trait in 'Step1' or different fixed effects in 'Step2' and clicking the 'Perform genome-wide analysis' button.", sep=""))


), trigger = "hover")



addPopover(session, "dummy5", "Details", content = HTML(paste("
By default, the 'best' set of snp in strongest association with the trait are reported. 
These results are given in table form. 
<br><br>
By clicking the 'Additional Summary' button on the left, two additional tables of results are shown; a table on the size and significance of the snp 
in the model and a table for the amount of phenotypic variance explained as they are added one at a time to the model.
<br><br>
There is additional computation needed to produce these extra tables. It may take a few minutes before these tables appear. 
     ", sep="")


), trigger = "hover")





 
session$onSessionEnded(stopApp)
  
}



shinyApp(ui=ui, server=server)


