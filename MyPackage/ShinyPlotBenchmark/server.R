
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# cat am.R | OMP_NUM_THREADS=1 singularity exec mro_cuda8_eagle_acc2.img /usr/bin/runR.sh 1 > new_1_thread_0_gpu_v2.out 2>&1
# cat new*.out | grep -e "profile,:" -e " is .." > times.txt
# cat times.txt | sort  -h -r  | tail -n +6 > times2.out


library(shiny)
library(shinyFiles)
library(rbokeh)

shinyServer(function(input, output, session) {
  flush1    <- Sys.getenv("FLUSH1DIR")
  flush1 <- ".."
  output$fi_text <- renderText({
    filePath <- input$fi_file$datapath
    # fileText <- paste(readLines(filePath), collapse = "\n")
    # fileText
    filePath
  })
  
  # volumes <- c('Home directory'='/flush1/bow355')
  volumes <- c('Home directory'=flush1)
  shinyFileChoose(input, 'files', root=volumes, session=session, filetypes=c('', 'txt')) 
  
  output$filepaths  <- renderPrint({
    testfile <-  parseFilePaths(volumes,input$files)[["datapath"]]
    fileText <- ""

    if (file.exists(as.character(testfile[1]))) {
       fileText <- paste(readLines(as.character(testfile[1])), collapse = " ")
    }
    fileText
  }) 
  
  output$rbokeh <- renderRbokeh({
    figure() %>% ly_points(1:10) %>%
      x_range(callback = shiny_callback("x_range"))
  })
  
})

