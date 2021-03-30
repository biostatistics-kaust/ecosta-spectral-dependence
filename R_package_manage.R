# remotes::install_github("wviechtb/mathjaxr")
# devtools::install_github("thomasp85/patchwork")
library(devtools)
library(roxygen2)

startDirectory <- getwd()
packageName <- "ECOSTASpecDepRef"
packageDir <- file.path(paste0(startDirectory, "/", packageName))

print(startDirectory)
print(packageDir)
#print(getwd())
#
if(!file.exists(packageDir)){
  #package.skeleton(name=packageName)
  defaultAuthors <- 'c(
    person("Pinto", "Marco", email="pinto.marco@live.com", role=c("aut", "cre")),
    person("Ombao", "Hernano", email="hernando.ombao@kaust.edu.sa", role=c("aut"))
  )'
  defaultProperties <- list(
    "Title"=paste(packageName, "Toolbox"),
    #"Authors@R"=list(other, myself),
    "Authors@R"=defaultAuthors,
    #"Maintainer"=defaultAuthors,
    "License"="MIT/GPL"
  )
  create(packageDir, defaultProperties)
  sampleFunctionPath <- paste0(packageDir, "/R/cat_function.R")
  sampleFunction <- "
#' A Cat Function 
#' 
#' This function allows you to express your love of cats. 
#' @param love Do you love cats? Defaults to TRUE. 
#' @keywords cats 
#' @export 
#' @examples 
#' cat_function()  
cat_function <- function(love=TRUE) {
  if(love==TRUE) {
    print('I love cats!')
  } else {
    print('I am not a cool person.')
  }
}
  "
  fp <- file(sampleFunctionPath, "w+")
  write(sampleFunction, fp)
  unlink(fp)
}

setwd(packageDir)
#devtools::document()
devtools::document(roclets = c('rd', 'collate', 'namespace'))
devtools::build(quiet=TRUE)
devtools::install(quiet=TRUE, upgrade="never")
setwd(startDirectory)


# Additional resources for handling R packages
# hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# https://www.rdocumentation.org/packages/devtools/versions/1.13.6/topics/create
# https://r-pkgs.org/man.html