source("Helpers/OperatingSystemHelpers.R")

require(randomForest)
library(foreach)
nofthreads=2 #number of cores to be used for training random forests



if(!isWindows())
{
    require("doMC")
    registerDoMC(nofthreads)
} 





source("BaydoganCode/TSBF_functions.r") #include the functions for subsequence and codebook generation
source("KeypointDescriptorsClusteringFeatures/find_interval_list_using_SiZer.R")
source('KeypointDescriptorsClusteringFeatures/find_descriptors_from_dataset.R')


TSBF_main_grid_search_inner_loop <- function(dataset_prefix,width_of_intervals_list,number_of_subsequences,randomize)
{
    gentype=4 #subsequence generation scheme -> 
    # 1:random (default) 2:uniform 3:totally random 4: Sizer List Deneme
    #Data characteristics
    
    


    train_filename=paste("DataSets/",dataset_prefix,"TRAIN",sep = "")
    test_filename=paste("DataSets/",dataset_prefix,"TEST",sep = "")
    print(train_filename)
    print(test_filename)
    
    #TSBF parameters
    minimum_interval_length=6;   #minimum interval length
    
    
    
    
    zlevels=c(1.0,0.1,0.25,0.5,0.75) #minimum subsequence length factors (z) to be evaluated
    zlevels=c(1.0) #minimum subsequence length factors (z) to be evaluated
    binsize=10      #bin size for codebook generation   
    
    #Experiment parameters
    number_of_replications=10       #number of replications
    noftree_step=50 #step size for tree building process
    tolerance=0.05  #marginal improvement in OOB error rates required for growing more trees in a forest
    verbose=1       #verbose=1 for detailed info about subsequences, 
    
    #define arrays to store information about computation times, OOB and test error rates
    test_error_rate<-array(1,number_of_replications) 
    OOB_error_rate<-array(1,number_of_replications)	  
    computation_time<-array(0,number_of_replications) 
    
    
    #Start TSBF
    traindata=as.matrix(read.table(train_filename))	
    testdata=as.matrix(read.table(test_filename))
    
    #standardize data if necessary
    traindata[,2:ncol(traindata)]=t(apply(traindata[,2:ncol(traindata)], 1, function(x) (x-mean(x))/sd(x)))
    testdata[,2:ncol(testdata)]=t(apply(testdata[,2:ncol(testdata)], 1, function(x) (x-mean(x))/sd(x)))
    
    trainclass=as.factor(traindata[,1])
    testclass=as.factor(testdata[,1])
    number_of_train=nrow(traindata)
    number_of_test=nrow(testdata)
    
    for(nrep in 1:number_of_replications){ 
        print(paste("replication:",nrep,"/",number_of_replications,sep=""))
    	ptm <- proc.time()	#record current time for replication start
    	test_error=array(1,length(zlevels))	#create test error rate array for each z level
    	OOB_error=array(1,length(zlevels))	#create OOB error rate array for each z level
    	for(z in 1:length(zlevels)){
    		ptm1 <- proc.time() #record current time for each zlevel
            if(gentype == 4)
            {
              interval_pair_list = find_interval_list_using_SiZer(train_filename,width_of_intervals_list,randomize)
              #x=generate_subsequences_new(traindata, testdata, interval_pair_list, verbose,number_of_subsequences)
              x=find_descriptors_from_dataset(dataset_prefix,width_of_intervals_list)
            }
            else
            {
                x=generate_subsequences(traindata, testdata, number_of_subsequences, zlevels[z], minimum_interval_length, gentype, verbose)
                
            }
    		
    		if(x$nsub<1) next #if minimum subsequence length is smaller than minimum_interval_length skip this z setting
    
    		if(isWindows())
    		{
    		    RFsub <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %do% randomForest(x$trainsub,x$classtr,ntree=ntree,na.action=na.roughfix)
    		}
    		else
    		{
    		    RFsub <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %dopar% randomForest(x$trainsub,x$classtr,ntree=ntree,na.action=na.roughfix)
    		}
    		prev_OOBerror=1; cur_OOBerror=1-sum(predict(RFsub,type='response')==x$classtr)/nrow(x$trainsub)				 
    		iter=1
    		while(iter<20&&cur_OOBerror<(1-tolerance)*prev_OOBerror){    
    			prev_OOBerror=cur_OOBerror
    			if(isWindows())
    			{
    			    RFsubmid <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %do% randomForest(x$trainsub,x$classtr,ntree=ntree,na.action=na.roughfix)
    			}
    			else
    			{
    			    RFsubmid <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %dopar% randomForest(x$trainsub,x$classtr,ntree=ntree,na.action=na.roughfix)
    			}
    			RFsub <- combine(RFsub, RFsubmid)
    			cur_OOBerror=1-sum(predict(RFsub,type='response')==x$classtr)/nrow(x$trainsub)
    			iter=iter+1
    		}
    
    		codetrain=generate_codebook(predict(RFsub,type='prob'),binsize,x$membertr)
    		codetest=generate_codebook(predict(RFsub,x$testsub,type='prob'),binsize,x$membertst)
    		
    		if(exists("RFsub"))
    		{
    		  rm(RFsub)
    		}
    		if(exists("RFsubmid"))
    		{
    		  rm(RFsubmid)
    		}
    		if(exists("x"))
    		{
    		  rm(x)
    		}
    		gc(FALSE) #clear memory
    
    		if(isWindows())
    		{
    		    RFts <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %do% randomForest(codetrain,trainclass,ntree=ntree,,na.action=na.roughfix)
    		}
    		else
    		{
    		    RFts <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %dopar% randomForest(codetrain,trainclass,ntree=ntree,,na.action=na.roughfix)
    		}
    		prev_OOBerror=1; cur_OOBerror=1-sum(predict(RFts,type='response')==trainclass)/number_of_train
    		iter=1
    		while(iter<20&&cur_OOBerror<(1-tolerance)*prev_OOBerror){    
    			prev_OOBerror=cur_OOBerror
    			if(isWindows())
    			{
    			    RFtsmid <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %do% randomForest(codetrain,trainclass,ntree=ntree,,na.action=na.roughfix)
    			}
    			else
    			{
    			    RFtsmid <- foreach(ntree=rep(noftree_step/nofthreads, nofthreads), .combine=combine, .packages='randomForest') %dopar% randomForest(codetrain,trainclass,ntree=ntree,,na.action=na.roughfix)
    			}
    			RFts <- combine(RFts, RFtsmid)
    			cur_OOBerror=1-sum(predict(RFts,type='response')==trainclass)/number_of_train
    			iter=iter+1
    		}	
    	
    		OOB_error[z]=1-sum(predict(RFts,type='response')==trainclass)/number_of_train
    		test_error[z]=1-sum(predict(RFts,codetest,type='response')==testclass)/number_of_test
    
    		rm(RFts,RFtsmid); gc(FALSE); #clear memory
    		passedZ=proc.time() - ptm1	
    
    		print(sprintf("Level z=%.3f is over for rep. %d in %.2f sec, OOB error=%.3f Test error=%.3f",zlevels[z],nrep,passedZ[3],OOB_error[z],test_error[z]))
    	}
    
    	passed=proc.time() - ptm
    	selected=which(OOB_error==min(OOB_error))
    	selectedtest=which.max(test_error[selected])
    	print(sprintf("Selected z=%.3f, test error=%.3f",zlevels[selected[selectedtest]],test_error[selected[selectedtest]]))
    	print(sprintf("Total run time (train+test) %.2f secs for rep. %d",passed[3],nrep))
    
    	test_error_rate[nrep]=test_error[selected[selectedtest]]
    	OOB_error_rate[nrep]=OOB_error[selected[selectedtest]]
    	computation_time[nrep]=passed[3]
    }

    result = list()
    string_to_return = sprintf ("RESULTS SUMMARY OVER %d REPLICATIONS\n Average test error rate: %.3f \n Min and Max test error rates: %.3f and %.3f \nAverage computation time (train+test):%.2f secs\n",
            nrep ,mean(test_error_rate),min(test_error_rate),max(test_error_rate),mean(computation_time))
    
    result$string_to_return <- string_to_return
    result$test_error_rate_mean <- mean(test_error_rate)
    result$test_error_rate_min <- min(test_error_rate)
    result$test_error_rate_max <- max(test_error_rate)
    result$test_error_rate <- test_error_rate
    result$replication_count <- nrep
    
    
    return(result)



}
