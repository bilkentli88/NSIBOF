// in command line compile with R CMD SHLIB extract_sub_new2.c 
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>

double find_interval_mean(double *input, int id, int len, int start, int end);
double find_interval_variance(double *input, int id, int len, double interval_mean, int start, int end);
double find_interval_max(double *input, int id, int len, int start, int end);
int find_interval_convexity(double *input, int id, int len, int start, int end,double interval_mean);


void extract_sub_new(
          double *trainseries
        , double *testseries
        , int *length_of_series
        , int *number_of_train
        , int *number_of_test
        , int *number_of_subsequences
        , int *min_sub_length
        , int *number_of_intervals
        , int *type
        , int *verbose
        , double *train_subsequences
        , double *testart_subsequencesequences
        , double *interval_pair_list
        ){
	if(*verbose>0){
        Rprintf("Calling extract_sub_new with length_of_series %d \n",*length_of_series); 
        Rprintf("Calling extract_sub_new with number_of_train %d \n",*number_of_train); 
        Rprintf("Calling extract_sub_new with number_of_test %d \n",*number_of_test); 
        Rprintf("Calling extract_sub_new with number_of_subsequences %d \n",*number_of_subsequences); 
        Rprintf("Calling extract_sub_new with min_sub_length %d \n",*min_sub_length); 
        Rprintf("Calling extract_sub_new with number_of_intervals %d \n",*number_of_intervals); 
        Rprintf("Calling extract_sub_new with type %d \n",*type); 
        Rprintf("Calling extract_sub_new with verbose %d \n",*verbose); 
	}  
            
    int i,j,k,start_subsequence,current_interval_length,feature_count_train,feature_count_test;
    double interval_convexity,interval_mean,interval_variance;

	GetRNGstate();
	
	if(*type==4 ){ // if generation scheme is not new 
		feature_count_train=0;
		feature_count_test=0;			
		for(i=0;i<*number_of_subsequences;i++){  
      		Rprintf("type %d. index i - number_of_subsequences %d - %d : \n",*type,i,*number_of_subsequences);
            int seq_start = *(interval_pair_list + 2*i) ;
            int seq_end = *(interval_pair_list + 2*i +1);
            start_subsequence = seq_start;
            current_interval_length = seq_end-seq_start;
            Rprintf("seq_start %d seq_end %d start_subsequence %d. current_interval_length %d \n",seq_start,seq_end,start_subsequence,current_interval_length);
            
			if(*verbose>0){
  				Rprintf("%d. Partition: Generated length %d, %d partitions, step size %d, start point %d, end point %d\n",i+1,*number_of_intervals*current_interval_length,*number_of_intervals,current_interval_length,start_subsequence+1,start_subsequence+*number_of_intervals*current_interval_length);
			}       			
			for(j=0;j<*number_of_train;j++){                      
				for(k=0;k<*number_of_intervals;k++){
                    interval_mean=find_interval_mean(trainseries,j,*length_of_series,start_subsequence+k*current_interval_length, start_subsequence+(k+1)*current_interval_length);
				    interval_convexity=find_interval_convexity(trainseries,j,*length_of_series,start_subsequence+k*current_interval_length, start_subsequence+(k+1)*current_interval_length,interval_mean);                    
				    interval_variance=find_interval_variance(trainseries,j,*length_of_series,interval_mean,start_subsequence+k*current_interval_length, start_subsequence+(k+1)*current_interval_length);  
				    train_subsequences[feature_count_train++]=interval_convexity;
				    train_subsequences[feature_count_train++]=interval_mean;
				    train_subsequences[feature_count_train++]=interval_variance;                                    
				}
				interval_mean=find_interval_mean(trainseries,j,*length_of_series, start_subsequence, start_subsequence+*number_of_intervals*current_interval_length);
				interval_variance=find_interval_variance(trainseries,j,*length_of_series,interval_mean,start_subsequence,start_subsequence+*number_of_intervals*current_interval_length);
				train_subsequences[feature_count_train++]=interval_mean;
			        train_subsequences[feature_count_train++]=interval_variance;
				train_subsequences[feature_count_train++]=start_subsequence+1;   
				train_subsequences[feature_count_train++]=start_subsequence+*number_of_intervals*current_interval_length;                                                                   
			}
			for(j=0;j<*number_of_test;j++){                      
				for(k=0;k<*number_of_intervals;k++){
                    interval_mean=find_interval_mean(testseries,j,*length_of_series,start_subsequence+k*current_interval_length, start_subsequence+(k+1)*current_interval_length);
                    interval_convexity=find_interval_convexity(testseries,j,*length_of_series,start_subsequence+k*current_interval_length, start_subsequence+(k+1)*current_interval_length,interval_mean);				    
				    interval_variance=find_interval_variance(testseries,j,*length_of_series,interval_mean,start_subsequence+k*current_interval_length, start_subsequence+(k+1)*current_interval_length);  
				    testart_subsequencesequences[feature_count_test++]=interval_convexity;
				    testart_subsequencesequences[feature_count_test++]=interval_mean;
				    testart_subsequencesequences[feature_count_test++]=interval_variance;                                    
				}
				interval_mean=find_interval_mean(testseries,j,*length_of_series, start_subsequence, start_subsequence+*number_of_intervals*current_interval_length);
				interval_variance=find_interval_variance(testseries,j,*length_of_series,interval_mean,start_subsequence,start_subsequence+*number_of_intervals*current_interval_length);
				testart_subsequencesequences[feature_count_test++]=interval_mean;
			        testart_subsequencesequences[feature_count_test++]=interval_variance;
				testart_subsequencesequences[feature_count_test++]=start_subsequence+1;   
				testart_subsequencesequences[feature_count_test++]=start_subsequence+*number_of_intervals*current_interval_length;                                                                        
			}
		}
	}

	else{
		Rprintf("Invalid option selected: type=%d\n",*type);
	}

	PutRNGstate();
}


double find_interval_mean(double *input, int id, int len, int start, int end){
    int k,stin;
    double sum,average;
    
    stin=id*len;
    sum=0;
    for(k=(start+stin);k<(end+stin);k++){                        
       sum=sum+input[k];                         
    }
    average=sum/(end-start);
    //Rprintf("find_interval_mean average %f \n",average);            
    return average;   
}

double find_interval_max(double *input, int id, int len, int start, int end){
    int k,stin;
    double max = DBL_MIN;
    double current;    
    stin=id*len;
    for(k=(start+stin);k<(end+stin);k++){                        
        current = input[k];                         
        if(current > max)
        {
            max = current;
        }
    }
    return max;   
}

int find_interval_convexity(double *input, int id, int len, int start, int end,double interval_mean)
{
    double max = find_interval_max(input, id, len, start, end);
    if(max > interval_mean)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}


double find_interval_variance(double *input, int id, int len, double interval_mean, int start, int end){

    int k,stin;
    double current,variance;
    
    stin=id*len;
    variance=0;
    for(k=(start+stin);k<(end+stin);k++){ 
        current = input[k];
        variance = variance + (current - interval_mean)*(current - interval_mean);                      
    }
    variance=variance/(end-start);
    //Rprintf("find_interval_variance variance %f \n",variance);            
    return variance;   
 
}

