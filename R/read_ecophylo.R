################################################################################
####                               ABOUT                                    ####
################################################################################
#                                                                              #
# This script is intended to be used along with the ecophylo library in Python #
# language, to import simulated datasets.You can redistribute it and/or modify #
# it under the terms of the CeCill licence.                                    #
#                                                                              #
# Authors: Elizabeth Barthelemy & Maxime Jaunatre                              #
#                                                                              #
################################################################################
require(readr)
require(stringr)
require(ape)

read_ecophylo <- function(file = NULL, type = c("params","sumstats","trees")) {
  cat(paste("\n Operating on file",file))
  # RTFM
  if(!is.character(file)) stop("file must be a string with the location of a file.")
  
  # import the file 
  file <- readr::read_file(file = file)
  
  # split into different output types
  parts <- strsplit(file, "\r*\n\\#+\r*\n")[[1]]
  parts <- parts[nchar(parts)> 0]
  
  if(length(parts)!= length(type)){
    stop("The file must contain as many output types as supplied in 'type'. 
         Default is parameters, abundance table and tree list in Newick format")
  }
  
  parameters <- NULL
  # Parameter table
  if("params" %in% type){
    params_raw <- parts[1]
    #process params_raw
    names <- unlist(strsplit(sub('\r*\n.+','', params_raw), " "))
    names <- names[nchar(names)> 0]
    parameters <- matrix(unlist(read.table(text= gsub('^.+[[:alpha:]]+[[:digit:]]*[[:space:]]{2,}','',gsub('\r*\n[[:digit:]]+','', parts[1])), 
                                           stringsAsFactors = F)), ncol = length(names), byrow = T)
    colnames(parameters) <- names
  }
  
  sumstats <- NULL
  # Summary statistics table
  if("sumstats" %in% type){
    if(length(type) == 1){
      sumstats_raw <- parts[1]
    }else{
      sumstats_raw <- parts[2]
    }
    #process sumstats_raw
    sumstats <- as.matrix(read.table(text=  sub("^([[:space:]]+[[:digit:]]+){1,}\r*\n","",sumstats_raw)))[,-1]
    colnames(sumstats) <- stringr::str_replace(colnames(sumstats), "X", "sp")
    
  }
  
  trees <- NULL
  # Phylogenies
  if("trees" %in% type){
    if(length(type)==1){
      trees_raw <- parts[1]
    }
    if(length(type) == 2){
      trees_raw <- parts[2]
    }
    if(length(type) == 3){
      trees_raw <- parts[3]
    }
    tmp <- unlist(strsplit(trees_raw,"\r*\n"))
    trees <- list()
    for(i in seq(length(tmp))){
      trees[[i]] <- ape::read.tree(text = tmp[i])
    }
    
  }
  list(parameters,sumstats, trees)
}
