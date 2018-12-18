#' @title Browser-based Graphical User Interface
#' @description Opens a web browser to act as a user-friendly interface to 'Eagle'
#' @details
#' \code{OpenGUI} is an easy to use web-based interface for 'Eagle'. By clicking on the navigation 
#' tabs at the top of a page, data can be read and analysed. By using this GUI, a user can avoid having to write R code. 
#'
#'
#' Note, that even though a web browser is being used as the user interface, everything remains local to the computer. 
#' @examples
#'\dontrun{
#'# opens a web browser 
#' OpenGUI()
#'}
#'
OpenGUI <- function() {
  appDir <- system.file('shiny_app', package = 'Eagle')
  if (appDir == "") {
    message("Could not find shiny-app directory. Try re-installing `Eagle` package.")
    return(NULL)
  }

  shinyAppDir(appDir, options=list(port = 3838))
}

