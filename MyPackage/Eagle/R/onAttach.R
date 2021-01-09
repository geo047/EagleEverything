.onAttach <- function(libname, pkgname) {
  # Runs when attached to search() path such as by library() or require()
  if (interactive()) {


##-------------------------------
## Opening message
##--------------------------------
packageStartupMessage(' ')
packageStartupMessage('   +--------------------------------------------------------------------------+ ')
packageStartupMessage('   |                  Eagle Version 2.4.3 has been loaded                     | ')  
packageStartupMessage('   |                                                                          | ')             
packageStartupMessage('   |    Visit eagle.r-forge.r-project.org for Quick Start guide and FAQ       | ')               
packageStartupMessage('   |                                                                          | ')
packageStartupMessage('   |                Type OpenGUI() to open web-based GUI                      | ')
packageStartupMessage('   +--------------------------------------------------------------------------+ ')
packageStartupMessage(' ')
}

}

