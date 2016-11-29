/*
 * channel_independent_selection.h
 *
 *  Created on: May 10, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef CHANNEL_INDEPENDENT_SELECTION_H_
#define CHANNEL_INDEPENDENT_SELECTION_H_
#include "linkedList.h"
#include "SENIA.h"
#include "IU.h"
void channel_independent_selection(gsl_matrix_complex *pH, int N, int k, int L, double snr, gsl_rng *pr,
		struct Node *index1_head, struct Node *index2_head,
		gsl_matrix_complex *pH1, gsl_matrix_complex *pH2, gsl_matrix_complex *G_preInv);
void channel_independent_selection(gsl_matrix_complex *pH, int N, int k, int L, double snr, gsl_rng *pr,
		struct Node *index1_head, struct Node *index2_head,
		gsl_matrix_complex *pH1, gsl_matrix_complex *pH2, gsl_matrix_complex *G_preInv
		){
    int Nr = pH->size1;
    int Nt = pH->size2;
    int Nu = gsl_sf_fact(Nt)/(gsl_sf_fact(Nt-N)*gsl_sf_fact(N));   //the number of subsets
    gsl_combination *subset;
    gsl_vector_complex *colReserve;
    colReserve = gsl_vector_complex_calloc(Nr);
    int m = gsl_rng_uniform_int(pr, Nu);   //randomly choose a combination subset
    int index;
    int count = 0;
    int count1,count2;
    int pilot;
    subset = gsl_combination_calloc(Nt, N);
    gsl_complex alpha, beta1, beta2;
    GSL_SET_COMPLEX(&alpha, 1,0);
    GSL_SET_COMPLEX(&beta1, pow(snr, -1),0);
    GSL_SET_COMPLEX(&beta2, 0,0);
    while (count < m){
    	gsl_combination_next(subset);
    	count++;
    }
    index2_head->next = create(Nt);
    index1_head->next = split(index2_head, subset);
    for(count = 0; count < N; count++){
    	index = gsl_combination_get(subset, count);
    	gsl_matrix_complex_get_col(colReserve, pH, index);
    	gsl_matrix_complex_set_col(pH1, count, colReserve);
    }
    count2 = 0;
    for (count = 0; count < Nt; count++){
    	pilot = 0;
    	for(count1 = 0; count1 < N; count1++){
    		if(count =  = gsl_combination_get(subset, count1)){
    			pilot = 1;
    			break;
    		}
    	}
    	if(pilot == 0){
    	 gsl_matrix_complex_get_col(colReserve, pH, count);
    	 gsl_matrix_complex_set_col(pH2, count2, colReserve);
    	 count2++;
    	}
    }
    //compute the initial matrix inverse
#if defined (SLDSENIAIU)
    gsl_matrix_complex *pH_sub1, *pH_sub2, *WIni, *WIniInv, *WIninonofuse;
    	   gsl_vector_complex *pH_col_tmp;
    	   pH_sub1 = gsl_matrix_complex_calloc(Nr, L);
    	   pH_sub2 = gsl_matrix_complex_calloc(Nr, Nt - N - L);
    	   WIni = gsl_matrix_complex_calloc(L, L);
    	   WIniInv = gsl_matrix_complex_calloc(L, L);
    	   WIninonofuse = gsl_matrix_complex_calloc(L, L);
    	   pH_col_tmp = gsl_vector_complex_calloc(Nr);
    		for(count1 = 0;count1 < Nt - N;count1++){
    			if(count1<L){
    			gsl_matrix_complex_get_col(pH_col_tmp, pH2, count1);
    			gsl_matrix_complex_set_col(pH_sub1, count1, pH_col_tmp);
    			}else{
    			gsl_matrix_complex_get_col(pH_col_tmp, pH2, count1);
    			gsl_matrix_complex_set_col(pH_sub2, count1-L, pH_col_tmp);
    			}
    		}
    		gsl_matrix_complex_set_identity(WIni);
    		gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH_sub1, pH_sub1, beta1, WIni);
    		SENIA(WIni, k, WIniInv, WIninonofuse); //SENIA
    		IU(WIniInv, pH_sub1, pH_sub2, snr, G_preInv);  //inflate update
    		gsl_matrix_complex_free(pH_sub1);
    		gsl_matrix_complex_free(pH_sub2);
    		gsl_matrix_complex_free(WIni);
    		gsl_matrix_complex_free(WIniInv);
    		gsl_matrix_complex_free(WIninonofuse);
    		gsl_vector_complex_free(pH_col_tmp);
#elif defined (SLD)
    	gsl_permutation *p = gsl_permutation_calloc(Nt-N);
    	int *signum = (int*)calloc(1, sizeof(int));
        *signum = 1;
       gsl_matrix_complex *G_pre = gsl_matrix_complex_calloc(Nt-N, Nt-N);
	   gsl_matrix_complex_set_identity(G_pre);
	   gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH2, pH2, beta1, G_pre);
	   gsl_linalg_complex_LU_decomp(G_pre, p, signum);
	   gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);
	   gsl_permutation_free(p);
	   free(signum);
	   gsl_matrix_complex_free(G_pre);
#endif
    gsl_combination_free(subset);
    gsl_vector_complex_free(colReserve);
	return;
}


#endif /* CHANNEL_INDEPENDENT_SELECTION_H_ */
