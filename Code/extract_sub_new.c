// in command line compile with 
// R CMD SHLIB extract_sub_new.c 
#include <R.h>
#include <Rmath.h>
#include <math.h>

double findmean(double *input, int id, int len, int start, int end);
double findvariance(double *input, int id, int len, double mean, int start, int end);
double findslope(double *input, int id, int len, int start, int end);


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
        , double *test_subsequences
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
            
    int i,j,k,st_sub,min_intlen,max_intlen,cur_intlen,feat_cnt,feat_cnt_test;
    double rnd,slp,mean,vrn;

	GetRNGstate();
	
	min_intlen=*min_sub_length/(*number_of_intervals);	
	if(*type==4 ){ // if generation scheme is not new 
		feat_cnt=0;
		feat_cnt_test=0;			
		for(i=0;i<*number_of_subsequences;i++){  
      		Rprintf("type %d. index i - number_of_subsequences %d - %d : \n",*type,i,*number_of_subsequences);
            int seq_start = *(interval_pair_list + 2*i) ;
            int seq_end = *(interval_pair_list + 2*i +1);
            st_sub = seq_start;
            cur_intlen = seq_end-seq_start;
            Rprintf("seq_start %d seq_end %d st_sub %d. cur_intlen %d \n",seq_start,seq_end,st_sub,cur_intlen);
            
			if(*verbose>0){
  				Rprintf("%d. Partition: Generated length %d, %d partitions, step size %d, start point %d, end point %d\n",i+1,*number_of_intervals*cur_intlen,*number_of_intervals,cur_intlen,st_sub+1,st_sub+*number_of_intervals*cur_intlen);
			}       			
			for(j=0;j<*number_of_train;j++){                      
				for(k=0;k<*number_of_intervals;k++){
				    slp=findslope(trainseries,j,*length_of_series,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    mean=findmean(trainseries,j,*length_of_series,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);
				    vrn=findvariance(trainseries,j,*length_of_series,mean,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    train_subsequences[feat_cnt++]=slp;
				    train_subsequences[feat_cnt++]=mean;
				    train_subsequences[feat_cnt++]=vrn;                                    
				}
				mean=findmean(trainseries,j,*length_of_series, st_sub, st_sub+*number_of_intervals*cur_intlen);
				vrn=findvariance(trainseries,j,*length_of_series,mean,st_sub,st_sub+*number_of_intervals*cur_intlen);
				train_subsequences[feat_cnt++]=mean;
			        train_subsequences[feat_cnt++]=vrn;
				train_subsequences[feat_cnt++]=st_sub+1;   
				train_subsequences[feat_cnt++]=st_sub+*number_of_intervals*cur_intlen;                                                                   
			}
			for(j=0;j<*number_of_test;j++){                      
				for(k=0;k<*number_of_intervals;k++){
				    slp=findslope(testseries,j,*length_of_series,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    mean=findmean(testseries,j,*length_of_series,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);
				    vrn=findvariance(testseries,j,*length_of_series,mean,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    test_subsequences[feat_cnt_test++]=slp;
				    test_subsequences[feat_cnt_test++]=mean;
				    test_subsequences[feat_cnt_test++]=vrn;                                    
				}
				mean=findmean(testseries,j,*length_of_series, st_sub, st_sub+*number_of_intervals*cur_intlen);
				vrn=findvariance(testseries,j,*length_of_series,mean,st_sub,st_sub+*number_of_intervals*cur_intlen);
				test_subsequences[feat_cnt_test++]=mean;
			        test_subsequences[feat_cnt_test++]=vrn;
				test_subsequences[feat_cnt_test++]=st_sub+1;   
				test_subsequences[feat_cnt_test++]=st_sub+*number_of_intervals*cur_intlen;                                                                        
			}
		}
	}

	else{
		Rprintf("Invalid option selected: type=%d\n",*type);
	}

	PutRNGstate();
}


double findmean(double *input, int id, int len, int start, int end){
    int k,stin;
    double sum,average;
    
    stin=id*len;
    sum=0;
    for(k=(start+stin);k<(end+stin);k++){                        
       sum=sum+input[k];                         
    }
    average=sum/(end-start);            
    return average;   
}

double findvariance(double *input, int id, int len, double mean, int start, int end){
    int k,stin;
    double avy;
    
    stin=id*len;
    avy=0;
    for(k=(start+stin);k<(end+stin);k++){                         
       avy=avy+pow((input[k]-mean),2);                         
    }
    avy=avy/(end-start-1);
            
    return avy;   
}

double findslope(double *input, int id, int len, int start, int end){
    int k,stin;
    double sx,sxy,avx,avy;
    
    stin=id*len;
    avx=0;
    avy=0;
    for(k=(start+stin);k<(end+stin);k++){    
       avx=avx+k;                         
       avy=avy+input[k];                         
    }
    avx=avx/(end-start);
    avy=avy/(end-start);
    
    sx=0;
    sxy=0;
    for(k=(start+stin);k<(end+stin);k++){  
       sx=sx+pow((k-avx),2);
       sxy=sxy+((k-avx)*(input[k]-avy));                        
    }

    return sxy/sx;    
}

