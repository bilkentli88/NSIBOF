// in command line compile with R CMD SHLIB extract_sub.c 
#include <R.h>
#include <Rmath.h>
#include <math.h>

double findmean(double *input, int id, int len, int start, int end);
double findvariance(double *input, int id, int len, double mean, int start, int end);
double findslope(double *input, int id, int len, int start, int end);

void extract_sub(double *trainseries, double *testseries, int *lenseries, int *noftrain, int *noftest, int *nsub, int *minsublength, int *nofint, int *type, int *verbose, double *trainsub, double *testsub){
	int i,j,k,st_sub,min_intlen,max_intlen,cur_intlen,feat_cnt,feat_cnt_test;
	double rnd,slp,mean,vrn;

	GetRNGstate();
	
	min_intlen=*minsublength/(*nofint);	
	if(*type==1||*type==2){ // if generation scheme is random or uniform over the same time points for train and test	
		feat_cnt=0;
		feat_cnt_test=0;			
		for(i=0;i<*nsub;i++){  //extract subsequences of random length and location, i is subsequence index  
			if(*type==1){			
				st_sub=floor(runif(0,*lenseries-*minsublength));
				max_intlen=((*lenseries)-st_sub)/(*nofint);
				cur_intlen=floor(runif(min_intlen,max_intlen));
			}
			else{
				st_sub=min_intlen*i;
				cur_intlen=min_intlen;
			}
			if(*verbose>0){
  				Rprintf("%d. Partition: Generated length %d, %d partitions, step size %d, start point %d, end point %d\n",i+1,*nofint*cur_intlen,*nofint,cur_intlen,st_sub+1,st_sub+*nofint*cur_intlen);
			}       			
			for(j=0;j<*noftrain;j++){                      
				for(k=0;k<*nofint;k++){
				    slp=findslope(trainseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    mean=findmean(trainseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);
				    vrn=findvariance(trainseries,j,*lenseries,mean,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    trainsub[feat_cnt++]=slp;
				    trainsub[feat_cnt++]=mean;
				    trainsub[feat_cnt++]=vrn;                                    
				}
				mean=findmean(trainseries,j,*lenseries, st_sub, st_sub+*nofint*cur_intlen);
				vrn=findvariance(trainseries,j,*lenseries,mean,st_sub,st_sub+*nofint*cur_intlen);
				trainsub[feat_cnt++]=mean;
			        trainsub[feat_cnt++]=vrn;
				trainsub[feat_cnt++]=st_sub+1;   
				trainsub[feat_cnt++]=st_sub+*nofint*cur_intlen;                                                                   
			}
			for(j=0;j<*noftest;j++){                      
				for(k=0;k<*nofint;k++){
				    slp=findslope(testseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    mean=findmean(testseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);
				    vrn=findvariance(testseries,j,*lenseries,mean,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    testsub[feat_cnt_test++]=slp;
				    testsub[feat_cnt_test++]=mean;
				    testsub[feat_cnt_test++]=vrn;                                    
				}
				mean=findmean(testseries,j,*lenseries, st_sub, st_sub+*nofint*cur_intlen);
				vrn=findvariance(testseries,j,*lenseries,mean,st_sub,st_sub+*nofint*cur_intlen);
				testsub[feat_cnt_test++]=mean;
			        testsub[feat_cnt_test++]=vrn;
				testsub[feat_cnt_test++]=st_sub+1;   
				testsub[feat_cnt_test++]=st_sub+*nofint*cur_intlen;                                                                        
			}
		}
	}
	else if(*type==3){     // if generation scheme is totally random	
		feat_cnt=0;
		feat_cnt_test=0;			
		for(j=0;j<*noftrain;j++){   
			for(i=0;i<*nsub;i++){  
				st_sub=floor(runif(0,*lenseries-*minsublength));
				max_intlen=((*lenseries)-st_sub)/(*nofint);
				cur_intlen=floor(runif(min_intlen,max_intlen));                
				for(k=0;k<*nofint;k++){
				    slp=findslope(trainseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    mean=findmean(trainseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);
				    vrn=findvariance(trainseries,j,*lenseries,mean,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    trainsub[feat_cnt++]=slp;
				    trainsub[feat_cnt++]=mean;
				    trainsub[feat_cnt++]=vrn;                                    
				}
				mean=findmean(trainseries,j,*lenseries, st_sub, st_sub+*nofint*cur_intlen);
				vrn=findvariance(trainseries,j,*lenseries,mean,st_sub,st_sub+*nofint*cur_intlen);
				trainsub[feat_cnt++]=mean;
			        trainsub[feat_cnt++]=vrn;
				trainsub[feat_cnt++]=st_sub+1;   
				trainsub[feat_cnt++]=st_sub+*nofint*cur_intlen;                                                                   
			}
		}
		for(j=0;j<*noftest;j++){   
			for(i=0;i<*nsub;i++){  
				st_sub=floor(runif(0,*lenseries-*minsublength));
				max_intlen=((*lenseries)-st_sub)/(*nofint);
				cur_intlen=floor(runif(min_intlen,max_intlen));                   
				for(k=0;k<*nofint;k++){
				    slp=findslope(testseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    mean=findmean(testseries,j,*lenseries,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);
				    vrn=findvariance(testseries,j,*lenseries,mean,st_sub+k*cur_intlen, st_sub+(k+1)*cur_intlen);  
				    testsub[feat_cnt_test++]=slp;
				    testsub[feat_cnt_test++]=mean;
				    testsub[feat_cnt_test++]=vrn;                                    
				}
				mean=findmean(testseries,j,*lenseries, st_sub, st_sub+*nofint*cur_intlen);
				vrn=findvariance(testseries,j,*lenseries,mean,st_sub,st_sub+*nofint*cur_intlen);
				testsub[feat_cnt_test++]=mean;
			        testsub[feat_cnt_test++]=vrn;
				testsub[feat_cnt_test++]=st_sub+1;   
				testsub[feat_cnt_test++]=st_sub+*nofint*cur_intlen;                                                                        
			}
		}
	}
	else{
		Rprintf("Invalid option selected: type=%d\n",*type);
	}

	PutRNGstate();
}


void generate_codebook(double *votes, int *nallsub, int *member, int *nofclass, int *nofbag, int *nofbin, double *codebook){
	int i,j,k,cnt=0,ind;
	double temp;
	for(i=0;i<*nofbag;i++){ 
		for(j=0;j<*nofclass;j++){ 
			for(k=0;k<*nofbin;k++){
				codebook[cnt++]=0;
			}
			codebook[cnt++]=0;
		}
	}
	
	for(i=0;i<*nallsub;i++){ 
		k=member[i]-1;
		temp=-1;
		for(j=0;j<*nofclass;j++){
			cnt=floor(*nofbin*votes[i+*nallsub*j]); 
			codebook[*nofbag*(cnt+j*(*nofbin))+k]++;
			if(votes[i+*nallsub*j]>temp){
				ind=j;
				temp=votes[i+*nallsub*j];
			}
		}
		cnt=*nofbin*(*nofclass); 
		codebook[*nofbag*(cnt+ind)+k]++;
	}
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
