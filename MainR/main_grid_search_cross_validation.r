
source("TSBF_main_grid_search_inner_loop.r")
source("DatasetImportCode/list_of_dataset_name.R")
source("Helpers/LogHelper.R")





list_of_width_of_intervals_list <- list(c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10))

number_of_subsequences = 1
cross_validation_count = 4
#list_of_dataset_name = seq(20)

list_of_dataset_name = c(
#  "Car"
  "Beef"
  
  )

list_to_write = list()
list_to_write$main_dataset_name <- "main_dataset_name"
list_to_write$width_of_intervals_list <- "width_of_intervals_list"
list_to_write$cv_index <- "cv_index"
list_to_write$current_dataset_name <- "dataset_name"
list_to_write$test_error_rate_mean <- "test_error_rate_mean"

writeToCsvLogFileHeader(list_to_write)


for(main_dataset_name in list_of_dataset_name)
{
  
  for(width_of_intervals_list in list_of_width_of_intervals_list)
  {
      for(cv_index in seq(cross_validation_count))
      {
         dataset_name = paste(main_dataset_name,"train-cross-validation", cross_validation_count,cv_index,sep="-")
          dataset_name_inner = paste(dataset_name,"_",sep="")
          train_filename = paste ("DataSets/",dataset_name, "_TRAIN", sep = "", collapse = NULL) 
          text_width_of_intervals_list = paste("width_of_intervals_list=c(",toString(width_of_intervals_list),")",sep="")
          text_number_of_subsequences = paste("number_of_subsequences=",toString(number_of_subsequences),sep="")
          
          
          run_result = TSBF_main_grid_search_inner_loop(dataset_name_inner,width_of_intervals_list,number_of_subsequences)
          
          list_to_write = list()
          list_to_write$main_dataset_name <- main_dataset_name
          list_to_write$width_of_intervals_list <- width_of_intervals_list
          list_to_write$cv_index <- cv_index
          list_to_write$current_dataset_name <- dataset_name
          list_to_write$test_error_rate_mean <- run_result$test_error_rate_mean
          writeToCsvLogFile(list_to_write)
          
      }
    }
  
  
}





