
/* based on gibbs_pottsv3_seqeunce_saq_v1.c Has numParallel MCMC chains 
   that interact. 
   outputs:
            the output_samples of first 2 chains
            energy of all nosim steps for all numParallel chains
 
          
  different in PT_final_chain_1 and PT_final_chain_4 is that the 1st one
  ony stores samples of the 1st chain in memeory while the 2econd one 
  stroes sampes of 4 chains */




#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

#define NEW(x) (x*)mxMalloc(sizeof(x))




void initProc(int total_lengthIn, double * curr_vectorIn, int * oneposIn, int * onepos_lengthIn, double * energy_currIn, double * J_MINFLOW_mat_array){

    int ivalue_onepos;
    int jvalue_onepos;
    int countIn, aba, indi, indj;    
    double currvalueIn;


    countIn=0;

    for (aba=0;aba<total_lengthIn;aba++) {
        currvalueIn = curr_vectorIn[aba];
        /*printf("currvalueIn is %f and aba is %d\n",currvalueIn,aba);*/
        if (currvalueIn==1){
            oneposIn[countIn]=aba;
            /*printf("oneposIn is %d\n",oneposIn[countIn]);*/
            countIn=countIn+1;

            }
    }
    *onepos_lengthIn=countIn;
    /*printf("oneposIn_length is %d\n",oneposIn_length);*/


    *energy_currIn=0;
    for (indi=0;indi<*onepos_lengthIn;indi++) {
        ivalue_onepos = oneposIn[indi];
        for (indj=0;indj<*onepos_lengthIn;indj++) {
            jvalue_onepos= oneposIn[indj];
            *energy_currIn = *energy_currIn + J_MINFLOW_mat_array[total_lengthIn*ivalue_onepos+jvalue_onepos];
        }
    }

}




void runMCMC(int aba, int total_length_In, int * onepos_In, int * dOnePos_In, int * dvalue_In, int * onepos_lengthPtr, double * random_site_array_In, double * rand_amino_array_In, double * phi_cumulative, double * phi_curr, double * energy_curr_main, double * J_MINFLOW_mat_array, double * random_array_In, double beta){
        

    int rand_site, rand_amino, curr_start, curr_end, curr_change_pos;
    int flag, onePosChange, indi, indj;
    double energy_new_temp, energy_curr_temp, trans_value, trans_prob;
    
    
    
    
	rand_site = random_site_array_In[aba]-1; /* random site starting from zero*/
	rand_amino=rand_amino_array_In[aba]-1; /* random amino in the site*/
	
	/*length_curr = phi_curr[rand_site];*/
	if (rand_site>0) /* find the start of the random site*/
		curr_start = phi_cumulative[rand_site-1];
	else
		curr_start=0;
	
	curr_end = curr_start + phi_curr[rand_site]-1;
	curr_change_pos = curr_start + rand_amino;
	
	/*printf("start is %d and end is %d and curr_change_pos is %d\n",curr_start,curr_end,curr_change_pos);*/
	/* *********************************** */
	/* form d value*/
	
	/* cancel d part */
	flag=0;
	dOnePos_In[0]=-99;
	onePosChange=-99;
	for (indi=0;indi<*onepos_lengthPtr;indi++) {
		if ((onepos_In[indi]<=curr_end) && (onepos_In[indi]>=curr_start)) /* if  the current vector is NOT the wildtype in the amino acid locations*/
		{
			dOnePos_In[0]= onepos_In[indi]; /* if there is a one, store the location of the one*/
			dvalue_In[0]= -1; /* have to subtract the one to get new vector*/
			flag=1;
			onePosChange=indi; /* the position in onepos_In which is a one and needs changing*/
		}
	}
	if (flag==0) /* if  the current vector is the wildtype in the amino acid locations*/
		dvalue_In[0]=0; /* don't have to substract anything*/
	
	/* add d part */
	dOnePos_In[1] = curr_change_pos; /* location where the bit is going to be flipped */
	if(dOnePos_In[1]!=dOnePos_In[0])
		dvalue_In[1] = 1; /* if flip from zero to one */
	else
		dvalue_In[1] = 0; /*        if flip from one to zero  */
	
	/*printf("energy curr is %f\n",energy_curr_temp);
	printf("aba is %d\n,",aba);
	printf("donepos is %d %d and dvalue_In is %d %d and oneposchange is %d and oneposlength is %d\n",dOnePos_In[0],dOnePos_In[1],dvalue_In[0],dvalue_In[1],onePosChange,*onepos_lengthPtr); */
	
	/*for (indj=0;indj<*onepos_lengthPtr;indj++) {
	printf("onepos %d is %d \n",indj,onepos_In[indj]);	
	}*/
	energy_curr_temp = *energy_curr_main;
	energy_new_temp=energy_curr_temp;
	if(  dOnePos_In[0]==-99)
	{
		
		indi=1;
		indj=1;
		
		energy_new_temp = energy_new_temp + J_MINFLOW_mat_array[total_length_In*dOnePos_In[indi]+dOnePos_In[indj]]*dvalue_In[indi]*dvalue_In[indj];
		
		for (indj=0;indj<*onepos_lengthPtr;indj++) {
			energy_new_temp = energy_new_temp + 2*J_MINFLOW_mat_array[total_length_In*dOnePos_In[indi]+onepos_In[indj]]*dvalue_In[indi];
		}
		
	}
	else
	{
		
		/*printf("here=1"); */
		for (indi=0;indi<2;indi++) {
		/*	printf("energy_new_temp is %f ",energy_new_temp);*/
			for (indj=0;indj<2;indj++) {
				energy_new_temp = energy_new_temp + J_MINFLOW_mat_array[(total_length_In*dOnePos_In[indi])+dOnePos_In[indj]]*dvalue_In[indi]*dvalue_In[indj];
			}
			for (indj=0;indj<*onepos_lengthPtr;indj++) {
				energy_new_temp = energy_new_temp + 2*(J_MINFLOW_mat_array[(total_length_In*dOnePos_In[indi])+onepos_In[indj]]*dvalue_In[indi]);
			}
			
		}
	}
	
	trans_value = 1/(1+exp(beta*(-energy_curr_temp+energy_new_temp)));
	
	if (1<trans_value)
		trans_prob=1;
	else
		trans_prob = trans_value;
	
		/*printf("energy_curr_temp is %f, energy_new_temp is %f and transprob is %f\n\n",energy_curr_temp,energy_new_temp,trans_prob);
		printf("energy_new_temp is %f\n",energy_new_temp);*/

	/*printf("aba is %d and thin is %d and burnin is %d\n",aba,thin,burnin); */

    
    
    

    if (random_array_In[aba]<trans_prob) /* change  */
    {
        /*curr_vector_bin[dOnePos_In(0)] = curr_vector_bin[dOnePos_In[0]]+dvalue_In[0]; */
        /*curr_vector_bin[dOnePos_In[1]] = curr_vector_bin[dOnePos_In[1]]+dvalue_In[1];*/

        /* change OnePos */
        if (dOnePos_In[1]!=dOnePos_In[0]) /* if not change to the wildtype */
        {

            if (dvalue_In[0]!=0) /* if the current vector is not the wildtype */

                onepos_In[onePosChange]=curr_change_pos; /* replace entry in onePos with new entry */
            else { /* if the current vector is the wildtype */

                onepos_In[*onepos_lengthPtr]=curr_change_pos; /* add entry to end */
                *onepos_lengthPtr=*onepos_lengthPtr+1;
            }

        }
        else { /* if change to the wildtype */
            /* shift entry back one place */
            for (indi=onePosChange;indi<*onepos_lengthPtr-1;indi++) {
                onepos_In[indi] = onepos_In[indi+1];
            }
            onepos_In[*onepos_lengthPtr]=-99;
            *onepos_lengthPtr=*onepos_lengthPtr-1;

        }
        energy_curr_temp=energy_new_temp;
        *energy_curr_main = energy_curr_temp;
    }

/*    printf("\n onepos_length = %d", *onepos_lengthPtr); */
    /*for(indj=0;indj<*onepos_lengthPtr;indj++) */
        

}






void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]) {

int count, aba, indi, indj,  rand_site, rand_amino, curr_start, curr_end, curr_change_pos, flag, onePosChange;
int burnin, thin, total_length, nosim, numParallel, mcSweepLength;
double energy_curr, energy_new, trans_value, trans_prob, beta;
double swap_value, swap_prob;
int * dOnePos;
int * dvalue;
double * random_array;
double * random_site_array;
double * rand_amino_array;
double * random_array_swap;
double * double_sum;
double * curr_vector_In;
double * J_MINFLOW_mat_array;
double * randomArray_sweep;
double currvalue; /*delete*/
int * onepos;
int onepos_length;
int ivalue_onepos;
int jvalue_onepos;
double * phi_cumulative;
double * phi_curr;
double * num_samples;
double * samples_output;
double * samples_output_2;
double * samples_output_3;
double * samples_output_4;
double * betaArray;
double * Vall;
double * numMutVecAll;
int sample_length;

int * countVec;
int * onepos_lengthVec;
double * energy_curr_main;
double * energyVecAll;


int  oneposTemp;
int  dOnePosTemp;
int  dvalueTemp;
double energy_curr_mainTemp;
int onepos_lengthVecTemp;



/* inputs*/
random_array = mxGetPr(prhs[0]);
random_site_array = mxGetPr(prhs[1]);
rand_amino_array = mxGetPr(prhs[2]);
curr_vector_In = mxGetPr(prhs[3]);
J_MINFLOW_mat_array = mxGetPr(prhs[4]);
total_length = mxGetScalar(prhs[5]);
nosim = mxGetScalar(prhs[6]); 
phi_cumulative= mxGetPr(prhs[7]);
phi_curr= mxGetPr(prhs[8]);
burnin= mxGetScalar(prhs[9]);
thin= mxGetScalar(prhs[10]);
curr_vector_In = mxGetPr(prhs[11]);
sample_length= mxGetScalar(prhs[12]);
numParallel = mxGetScalar(prhs[13]);
betaArray = mxGetPr(prhs[14]);
mcSweepLength = mxGetScalar(prhs[15]);
random_array_swap = mxGetPr(prhs[16]);

/* ************************************ */
plhs[0]= mxCreateDoubleMatrix(1,(mwSize) 1,mxREAL);
double_sum =  mxGetPr(plhs[0]);

plhs[1]= mxCreateDoubleMatrix(1,1,mxREAL);
num_samples =  mxGetPr(plhs[1]);

plhs[2]= mxCreateDoubleMatrix(1,(mwSize) numParallel*nosim,mxREAL);
energyVecAll =  mxGetPr(plhs[2]);

plhs[3]= mxCreateDoubleMatrix(1,(mwSize) numParallel*nosim,mxREAL);
numMutVecAll =  mxGetPr(plhs[3]);

plhs[4]= mxCreateDoubleMatrix(1,(mwSize) nosim/mcSweepLength*(numParallel - 1) ,mxREAL);
Vall =  mxGetPr(plhs[4]);

plhs[5]= mxCreateDoubleMatrix(1,(mwSize) total_length*sample_length,mxREAL);
samples_output =  mxGetPr(plhs[5]);
/*
plhs[6]= mxCreateDoubleMatrix(1,(mwSize) total_length*sample_length,mxREAL);
samples_output_2 =  mxGetPr(plhs[6]);

plhs[7]= mxCreateDoubleMatrix(1,(mwSize) total_length*sample_length,mxREAL);
samples_output_3 =  mxGetPr(plhs[7]);

plhs[8]= mxCreateDoubleMatrix(1,(mwSize) total_length*sample_length,mxREAL);
samples_output_4 =  mxGetPr(plhs[8]);
*/


/* *********************************** */
onepos=mxMalloc(sizeof(int)*total_length*numParallel+total_length);
dOnePos=mxMalloc(sizeof(int)*3*numParallel);
dvalue=mxMalloc(sizeof(int)*3*numParallel);
energy_curr_main = mxMalloc(sizeof(double)*numParallel);
countVec = mxMalloc(sizeof(int)*numParallel);
onepos_lengthVec = mxMalloc(sizeof(int)*numParallel);

/*
oneposTemp=mxMalloc(sizeof(int)*1);
dOnePosTemp=mxMalloc(sizeof(int));
dvalueTemp=mxMalloc(sizeof(int));
energy_curr_mainTemp = mxMalloc(sizeof(double));
 no need to swap countVec as it same for all js 

onepos_lengthVecTemp = mxMalloc(sizeof(int));

*/

/* *********************************** */


double curr_vector[total_length];
int i, j;




double * curr_vectorPtrD;
int * oneposPtrI;
int * onepos_lengthPtrI;
int * dOnePosPtrI;
int * dvaluePtrI;
double * energy_currPtrD;
double * random_site_arrayPtrD;
double * rand_amino_arrayPtrD;
double * random_arrayPtrD;
double * energyVecPtrD;
double * numMutVecPtrD;

int MCswapCount, MCswapCount_2;
MCswapCount = 0;
MCswapCount_2 = 0;

for(i=0;i<numParallel;i++)
    energy_curr_main[i] = 0;

for (j=0;j<numParallel;j++){
    

    for(i=0;i<total_length;i++)
        curr_vector[i] = curr_vector_In[i + j*total_length];
    
    
    curr_vectorPtrD = &curr_vector_In[j*total_length]; /* curr_vectorPtrD = &curr_vector_In[j*total_length]; */
    oneposPtrI = &onepos[j*total_length];
    onepos_lengthPtrI = &onepos_lengthVec[j];  /*onepos_lengthPtrI = &onepos_length;*/
    energy_currPtrD = &energy_curr_main[j]; /*&energy_curr;*/

/* printf("this is the %dth run...", j); */
    /*initProc(int total_lengthIn, int * curr_vectorIn, int * oneposIn, int * onepos_lengthIn, double * energy_curr, double * J_MINFLOW_mat_array) */

    initProc(total_length, curr_vectorPtrD, oneposPtrI,onepos_lengthPtrI, energy_currPtrD, J_MINFLOW_mat_array);



}

/*
dOnePos[0]=1;
dOnePos[1]=1;
dvalue[0]=0;
dvalue[0]=1;
count = 0;*/

for (j=0;j<numParallel;j++){
    countVec[j]=0;
    
    dOnePosPtrI = &dOnePos[j*3];
    dvaluePtrI = &dvalue[j*3];
    
    dOnePosPtrI[0]=1;
    dOnePosPtrI[1]=1;
    dvaluePtrI[0]=0;
    dvaluePtrI[0]=1;
    
}




for (aba=0;aba<nosim;aba++) {

    for (j=0;j<numParallel;j++){
        
        oneposPtrI = &onepos[j*total_length];
        dOnePosPtrI = &dOnePos[j*3];
        dvaluePtrI = &dvalue[j*3];
        onepos_lengthPtrI = &onepos_lengthVec[j];  /*onepos_lengthPtrI = &onepos_length;*/
        random_site_arrayPtrD = &random_site_array[j*nosim]; /* random_site_arrayPtrD = &random_site_array[0]; */
        rand_amino_arrayPtrD = &rand_amino_array[j*nosim]; /* rand_amino_arrayPtrD = &rand_amino_array[0]; */
        random_arrayPtrD = &random_array[j*nosim]; /* random_arrayPtrD = &random_array[0]; */
        energy_currPtrD = &energy_curr_main[j]; /*energy_currPtrD = &energy_curr; */
        beta = betaArray[j];
        energyVecPtrD = &energyVecAll[j*nosim + aba];
        numMutVecPtrD = &numMutVecAll[j*nosim + aba];
                
                
    /*runMCMC(int aba, int total_length_In, int * onepos_In, int * dOnePos_In, int * dvalue_In, int * onepos_lengthPtr, double * random_site_array_In, double * rand_amino_array_In, double * phi_cumulative, double * phi_curr, double * energy_curr_main, double * J_MINFLOW_mat_array, double * random_array, double beta));

        runMCMC(int aba, int total_length, 
        int * onepos_In, int * dOnePos_In, int * dvalue_In, int * onepos_lengthPtr, 
        double * random_site_array_In, double * rand_amino_array_In, 
        double * phi_cumulative, double * phi_curr, double * energy_curr_main, 
        double * J_MINFLOW_mat_array, double * random_array)    */


        /*runMCMC(aba, total_length, oneposPtrI, dOnePosPtrI, dvaluePtrI, onepos_lengthPtrI, random_site_arrayPtrD, rand_amino_arrayPtrD, phi_cumulative, phi_curr, energy_currPtrD, J_MINFLOW_mat_array, random_arrayPtrD, beta);*/
        runMCMC(aba, total_length, oneposPtrI, dOnePosPtrI, dvaluePtrI, onepos_lengthPtrI, random_site_arrayPtrD, rand_amino_arrayPtrD, phi_cumulative, phi_curr, energy_currPtrD, J_MINFLOW_mat_array, random_arrayPtrD, beta);
        *energyVecPtrD = *energy_currPtrD;
        *numMutVecPtrD = (double) *onepos_lengthPtrI;

    /*    printf("onepos_lengthPtrI= %d", *onepos_lengthPtrI);*/



        if (aba>=burnin){

            if ( (aba+1) % thin == 0){
           
                /*  energyVecPtrD = &energyVecAll[j*sample_length + countVec[j]];
                    *energyVecPtrD = *energy_currPtrD;
                    change size of energyVecAll
                 */
                
                 
                /* printf("aba is %d and thin is %d\n",aba,thin); */

                /*             samples2(count,:) =curr_vector_bin;  */
                for (indi=0;indi<*onepos_lengthPtrI;indi++) {
                    /* double_sum[total_length*onepos[indi]+onepos[indi]] = double_sum[total_length*onepos[indi]+onepos[indi]]+1;*/

                    /*for (indj=indi+1;indj<onepos_length;indj++) {
                        double_sum[total_length*onepos[indi]+onepos[indj]] = double_sum[total_length*onepos[indi]+onepos[indj]]+1;
                    }*/
                    if(j==0)
                        samples_output[countVec[j]*total_length+oneposPtrI[indi]]=1;
/*                    if(j == 1)
                        samples_output_2[countVec[j]*total_length+oneposPtrI[indi]]=1;
                    if(j == 2)
                        samples_output_3[countVec[j]*total_length+oneposPtrI[indi]]=1;
                    if(j == 3)
                        samples_output_4[countVec[j]*total_length+oneposPtrI[indi]]=1;
                    
  */                  
                     
                }

                countVec[j]=countVec[j]+1;
            }   
        }
        

    }
    
    
    



    if ( (aba+1) % mcSweepLength == 0){
        for (j=numParallel-2; j>-1 ;j--){
     
            /* swap replicas every MCseep 
               for loop works backward, swaping M-1, and M-2 chain first all
               the way down to chain no. 0  */
            
            
            /* check if want to swap or not */
            swap_value = 1/(1+exp((betaArray[j+1] - betaArray[j])*(-energy_curr_main[j+1]+energy_curr_main[j])));
	
            if (1<swap_value)
                swap_prob=1;
            else
                swap_prob = swap_value;
            
     
            
            Vall[MCswapCount_2] = swap_prob;
            MCswapCount_2 = MCswapCount_2 + 1;
            if (random_array_swap[MCswapCount]<swap_prob){ /* change  */
                
                for(i=0;i<total_length;i++){
                    oneposTemp = onepos[(j+1)*total_length + i];
                    onepos[(j+1)*total_length + i] = onepos[j*total_length + i];
                    onepos[j*total_length + i] = oneposTemp;
                }
                for(i=0;i<3;i++){
                    dOnePosTemp = dOnePos[(j+1)*3 + i];
                    dOnePos[(j+1)*3 + i] = dOnePos[j*3 + i];
                    dOnePos[j*3 + i] = dOnePosTemp;
                    
                    dvalueTemp = dvalue[(j+1)*3 + i];
                    dvalue[(j+1)*3 + i] = dvalue[j*3 + i];
                    dvalue[j*3 + i] = dvalueTemp;
                }
                
                energy_curr_mainTemp = energy_curr_main[j+1];
                energy_curr_main[j+1] = energy_curr_main[j];
                energy_curr_main[j] = energy_curr_mainTemp;
                
                onepos_lengthVecTemp = onepos_lengthVec[j+1];
                onepos_lengthVec[j+1] = onepos_lengthVec[j];
                onepos_lengthVec[j] = onepos_lengthVecTemp;
                
                
            }
            
        }
        
        MCswapCount = MCswapCount + 1;
        
    }
    
}
    
/*printf("Onepos_length is %d and size is %f",onepos_length,sizeof(onepos)); */

/*mxFree(onepos);
free(dOnePos);
free(dvalue); */
      num_samples[0] = countVec[0];
      
/*      for (j=0;j<numParallel;j++)
          printf("countVec of j %d", countVec[j]);
 */     
      /*  double_sum[0]=1; */


}
	
	
        
        

