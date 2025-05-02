
source("TSBF_main_grid_search_inner_loop.r")
source("DatasetImportCode/list_of_dataset_name.R")
source("Helpers/LogHelper.R")

#set.seed(1234567899)


randomize =FALSE

list_of_dataset_name = c("UWaveGestureLibraryZ")

list_of_width_of_intervals_list <- list(c(7),c(8),c(9))
#list_of_width_of_intervals_list <- list(c(39),c(40),c(41),c(42),c(43),c(44),c(45),c(46),c(47),c(48),c(49),c(50),c(51),c(52))

#list_of_width_of_intervals_list <- list(c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11),c(12),c(13),c(14),c(15),c(16), c(17),c(18),c(19),c(20),c(21),c(22),c(23),c(24),c(25),c(26),c(28),c(29),c(30),c(31),c(32),c(39),c(40),c(41),c(42),c(43),c(44),c(45),c(46),c(47),c(48),c(49),c(50),c(51),c(52),c(53),c(54),c(55),c(56),c(60),c(61),c(99))

#list_of_width_of_intervals_list <- list(c(2),c(3),c(4),c(5),c(6),c(7),c(8),c(9),c(10),c(11),c(12),c(13),c(14),c(15),c(16), c(17),c(18),c(19),c(20),c(21),c(30),c(40),c(50))

#list_of_width_of_intervals_list <- list(c(2),c(3),c(4),c(5),c(6),c(7),c(8))

number_of_subsequences = 1

for(dataset_name in list_of_dataset_name)
{
    for(width_of_intervals_list in list_of_width_of_intervals_list)
    {
        print(dataset_name)
        dataset_name2 = paste(dataset_name,"_",sep="")
        train_filename = paste ("DataSets/",dataset_name, "_TRAIN", sep = "", collapse = NULL) 
        text_width_of_intervals_list = paste("width_of_intervals_list=c(",toString(width_of_intervals_list),")",sep="")
        text_number_of_subsequences = paste("number_of_subsequences=",toString(number_of_subsequences),sep="")
        
        
        result = tryCatch({
            run_result = TSBF_main_grid_search_inner_loop(dataset_name2,width_of_intervals_list,number_of_subsequences,randomize)
            writeToLogFile(c("SUCCESS",train_filename,text_width_of_intervals_list,text_number_of_subsequences,run_result$string_to_return))
        }, warning = function(w) {
            writeToLogFile(c("WARNING",train_filename,text_width_of_intervals_list,text_number_of_subsequences,w))
        }, error = function(e) {
            writeToLogFile(c("ERROR",train_filename,text_width_of_intervals_list,text_number_of_subsequences,e))
        }, finally = {
            # NONE
        })
    }
    
}





