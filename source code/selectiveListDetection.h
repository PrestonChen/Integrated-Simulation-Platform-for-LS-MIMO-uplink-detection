/*
 * selectiveListDetection.h
 *
 *  Created on: May 10, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef SELECTIVELISTDETECTION_H_
#define SELECTIVELISTDETECTION_H_
#include "SENIA.h"
#include "IU.h"
#include "linkedList.h"
#include "fullfact.h"
#include "diversity_max_selection.h"
#include "channel_independent_selection.h"
#include "MMSE_SIC_ordering.h"
#include "MMSE_SIC_ordering_fast.h"
void selectiveListDetection(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
		double snr, double pav, int N, int k, int L, gsl_vector_complex *pSymConstell,
		gsl_vector_complex *psymOut);
void selectiveListDetection(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
		double snr, double pav, int N, int k, int L, gsl_vector_complex *pSymConstell,
		gsl_vector_complex *psymOut){
	const int Nr=pH->size1;
	const int Nt=pH->size2;
	const int M=pSymConstell->size;
	const int listSize=pow(M, N);   //the number of vector candidates in the list
	gsl_rng *pr;    //random number generator
	int seed = time (NULL);
	const gsl_rng_type *pT;
	pT = gsl_rng_default;
	pr = gsl_rng_alloc (pT);
	gsl_rng_set (pr,seed);
	gsl_matrix_complex *pH1, *pH2, *WInv, *G;
    gsl_vector_complex *preceiveSubtmp, *psymOut_Sub1tmp, *psymOut_Sub2,
    *psymOut_Sub2tmp;
    gsl_matrix_complex *psymOut_Sub2_M;
    gsl_matrix *List;
    gsl_vector *subList;
    gsl_vector *EuclideanV;
    struct Node *index1_head, *index2_head;
    index1_head=create(1);
    index2_head=create(1);
    WInv=gsl_matrix_complex_calloc(Nt-N, Nt-N);
    G=gsl_matrix_complex_calloc(Nt-N, Nr);
	pH1=gsl_matrix_complex_calloc(Nr, N);
	pH2=gsl_matrix_complex_calloc(Nr, Nt-N);
    preceiveSubtmp=gsl_vector_complex_calloc(Nr);
    psymOut_Sub1tmp=gsl_vector_complex_calloc(N);
    psymOut_Sub2=gsl_vector_complex_calloc(Nt-N);
    psymOut_Sub2tmp=gsl_vector_complex_calloc(Nt-N);
    psymOut_Sub2_M=gsl_matrix_complex_calloc(Nt-N, listSize);
    subList=gsl_vector_calloc(N);
    List=gsl_matrix_calloc(N,listSize);
    EuclideanV=gsl_vector_calloc(listSize);
    double Euclidean;
    int count, count1;
    gsl_complex alpha, beta1, beta2, betasnr;
    GSL_SET_COMPLEX(&alpha, 1, 0);
    GSL_SET_COMPLEX(&beta1, 0, 0);
    GSL_SET_COMPLEX(&beta2, -1, 0);
    GSL_SET_COMPLEX(&betasnr, pow(snr, -1), 0);
    int m;
    int index;
    fullfact(N, M, List);   //generate all the possible sub symbol vector hypotheses
#if defined (DMS)   //diversity maximization selection
    diversity_max_selection(pH, N, k, L, snr, index1_head, index2_head, pH1, pH2, WInv);
#elif defined(CIS) //channel independent selection
    channel_independent_selection(pH, N, k, L, snr , pr, index1_head,
    		index2_head, pH1, pH2, WInv);
#endif

//    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, WInv, pH2, beta1, G);
    for (count=0;count<listSize; count++){
    	//construct sub symbol vector s1
    	gsl_matrix_get_col(subList, List, count);
    	for (count1=0;count1<N;count1++){
    		gsl_vector_complex_set(psymOut_Sub1tmp, count1,
    				gsl_vector_complex_get(pSymConstell,
    						gsl_vector_get(subList, count1)));
    	}
    	//interference cancellation
    	gsl_vector_complex_memcpy(preceiveSubtmp, preceived);
    	gsl_blas_zgemv(CblasNoTrans, beta2, pH1, psymOut_Sub1tmp, alpha, preceiveSubtmp);
        //MMSE-SIC sub-detection
    	MMSE_SIC_ordering_fast(preceiveSubtmp, pH2, WInv, snr, pav, M, psymOut_Sub2tmp);
    	gsl_matrix_complex_set_col(psymOut_Sub2_M, count, psymOut_Sub2tmp);
    	gsl_blas_zgemv(CblasNoTrans, beta2, pH2, psymOut_Sub2tmp, alpha, preceiveSubtmp);
    	Euclidean=gsl_blas_dznrm2(preceiveSubtmp);
    	gsl_vector_set(EuclideanV, count, Euclidean);

    }
    //reconstruct psymOut
    m=gsl_vector_min_index(EuclideanV);
   	gsl_matrix_get_col(subList, List, m);
    	for (count1=0;count1<N;count1++){  //construct sub symbol vector
    		gsl_vector_complex_set(psymOut_Sub1tmp, count1,
    				gsl_vector_complex_get(pSymConstell,
    						gsl_vector_get(subList, count1)));
    	}

    	gsl_matrix_complex_get_col(psymOut_Sub2, psymOut_Sub2_M, m);

    	for (count1=0;count1<N;count1++){
         index=get(0, index1_head);
         gsl_vector_complex_set(psymOut, index,
        		 gsl_vector_complex_get(psymOut_Sub1tmp,count1));
    	}
    	for (count1=0;count1<(Nt-N);count1++){
    		index=get(0, index2_head);
            gsl_vector_complex_set(psymOut, index,
           		 gsl_vector_complex_get(psymOut_Sub2,count1));
    	}
    gsl_rng_free(pr);
    gsl_matrix_complex_free(WInv);
    gsl_matrix_complex_free(G);
	gsl_matrix_complex_free(pH1);
	gsl_matrix_complex_free(pH2);
	gsl_vector_complex_free(preceiveSubtmp);
	gsl_vector_complex_free(psymOut_Sub1tmp);
	gsl_vector_complex_free(psymOut_Sub2);
	gsl_vector_complex_free(psymOut_Sub2tmp);
	gsl_matrix_complex_free(psymOut_Sub2_M);
	gsl_matrix_free(List);
	gsl_vector_free(subList);
	gsl_vector_free(EuclideanV);
	freeLinkedList(index1_head);
	freeLinkedList(index2_head);

	return;
}


#endif /* SELECTIVELISTDETECTION_H_ */
