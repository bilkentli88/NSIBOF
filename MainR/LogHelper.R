writeToLogFile<- function(textToWrite){
  log_filename<-"output-TSBF_main_grid_search.txt"
  for(line in textToWrite)
  {
    line_to_write=toString(line)
    write(line_to_write, file=log_filename,append=TRUE)
    print(line_to_write)
  }
}

writeToCsvLogFile<- function(list_to_write){
  log_filename<-"output.csv"

  line_to_write = ""  
  for (name in names(list_to_write)) {
    word = list_to_write[[name]]
    line_to_write= paste(line_to_write,toString(word),sep=",")
  }
  line_to_write <- sub(',', '', line_to_write)
  write(line_to_write, file=log_filename,append=TRUE)
  print(line_to_write)
  

}


writeToCsvLogFileHeader<- function(list_to_write){
  log_filename<-"output.csv"
  
  line_to_write = ""  
  for (name in names(list_to_write)) {
    word = list_to_write[[name]]
    line_to_write= paste(line_to_write,toString(word),sep=",")
  }
  line_to_write <- sub(',', '', line_to_write)
  write(line_to_write, file=log_filename,append=FALSE)
  print(line_to_write)
  
  
}
