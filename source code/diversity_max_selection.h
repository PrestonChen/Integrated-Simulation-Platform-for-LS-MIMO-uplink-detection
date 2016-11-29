/*
 * diversity_max_selection.h
 *
 *  Created on: Dec 11, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "linkedList.h"
#include "SENIA.h"
#include "IU.h"
#ifndef DIVERSITY_MAX_SELECTION_H_
#define DIVERSITY_MAX_SELECTION_H_
#define LEN sizeof(struct Node)
void diversity_max_selection(gsl_matrix_complex *pH, int N, int k, int L, double snr,
		struct Node *index1_head, struct Node *index2_head,
		gsl_matrix_complex *pH1, gsl_matrix_complex *pH2,
		gsl_matrix_complex *WInv);
void diversity_max_selection(gsl_matrix_complex *pH, int N, int k, int L, double snr,
		struct Node *index1_head, struct Node *index2_head,
		gsl_matrix_complex *pH1, gsl_matrix_complex *pH2,
		gsl_matrix_complex *WInv){
	int Nr = pH->size1;
	int Nt = pH->size2;
	gsl_combination *subset;
	gsl_matrix_complex *G_pre, *G_preInv;
	gsl_permutation *p;
	int *signum = (int*)calloc(1, sizeof(int));
	*signum = 1;
	gsl_matrix_complex *H1Tmp, *H2Tmp;
	gsl_vector_complex_view diag_viewComplex;
	gsl_vector_view diag_viewReal;
	gsl_vector *diagMax, *diag;
	gsl_vector_complex *colReserve;
	subset = gsl_combination_calloc(Nt, N);
	int Nu = gsl_sf_fact(Nt)/(gsl_sf_fact(Nt-N)*gsl_sf_fact(N));   //the number of subsets
	int count, count1, count2, count3, count4;
	double minValue;
	p = gsl_permutation_calloc(Nt-N);
	G_pre = gsl_matrix_complex_calloc(Nt-N, Nt-N);
	G_preInv = gsl_matrix_complex_calloc(Nt-N, Nt-N);
	H1Tmp = gsl_matrix_complex_calloc(Nr, N);
	H2Tmp = gsl_matrix_complex_calloc(Nr, Nt-N);
	diagMax = gsl_vector_calloc(Nu);   //store the minimum diagonal value of all the subset
	diag = gsl_vector_calloc(Nt-N);   //the diagonal vector of one subset
	colReserve = gsl_vector_complex_calloc(Nr);
	gsl_complex alpha, beta1, beta2;
	GSL_SET_COMPLEX(&alpha, 1, 0);
	GSL_SET_COMPLEX(&beta1, pow(snr, -1), 0);
	GSL_SET_COMPLEX(&beta2, 0,0);
//	int m;
	int pilot = 0;
	double maxValueTmp = 0;
   for (count = 0;count<Nu;count++){
      //construct H1Tmp, H2Tmp
	   count3 = 0;
	   count4 = 0;
	   for (count1 = 0; count1 < Nt; count1++){
		   pilot = 0;
		   for (count2 = 0; count2 < N; count2++){
			   if (count1 =  = gsl_combination_get(subset,count2)){
				   pilot = 1;
				   gsl_matrix_complex_get_col(colReserve, pH, count1);
				   gsl_matrix_complex_set_col(H1Tmp, count3, colReserve);
				   count3++;
				   break;
			   }
		   }
		   if(pilot =  = 1){
			   continue;
		   }else{
			   gsl_matrix_complex_get_col(colReserve, pH, count1);
			   gsl_matrix_complex_set_col(H2Tmp, count4, colReserve);
			   count4++;
		   }
	   }
#if defined (SLDSENIAIU)

#elif defined (SLD)
	   gsl_matrix_complex_set_identity(G_pre);
	   gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, H2Tmp, H2Tmp, beta1, G_pre);
#endif
/*
* compute the inverse of G_pre G_preInv (consider EMI and AMI cases)
*/
#if defined (SLD)
	   gsl_linalg_complex_LU_decomp(G_pre, p, signum);
	   gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);
#elif defined (SLDSENIAIU)
	   gsl_matrix_complex *pH_sub1, *pH_sub2, *WIni, *WIniInv, *WIninonofuse;
	   gsl_vector_complex *pH_col_tmp;
	   pH_sub1 = gsl_matrix_complex_calloc(Nr, L);
	   pH_sub2 = gsl_matrix_complex_calloc(Nr, Nt-N-L);
	   WIni = gsl_matrix_complex_calloc(L, L);
	   WIniInv = gsl_matrix_complex_calloc(L, L);
	   WIninonofuse = gsl_matrix_complex_calloc(L, L);
	   pH_col_tmp = gsl_vector_complex_calloc(Nr);
		for(count1 = 0;count1<Nt-N;count1++){
			if(count1<L){
			gsl_matrix_complex_get_col(pH_col_tmp, H2Tmp, count1);
			gsl_matrix_complex_set_col(pH_sub1, count1, pH_col_tmp);
			}else{
			gsl_matrix_complex_get_col(pH_col_tmp, H2Tmp, count1);
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
#endif
	   diag_viewComplex = gsl_matrix_complex_diagonal(G_preInv);
	   diag_viewReal = gsl_vector_complex_real(&diag_viewComplex.vector);
	   gsl_vector_memcpy(diag, &diag_viewReal.vector);
	   maxValueTmp = gsl_vector_max(diag);
	   //update the information of the subset with strongest weakest sub-data stream
	   if(count =  = 0){
		   minValue = maxValueTmp;
		   gsl_matrix_complex_memcpy(pH1, H1Tmp);
		   gsl_matrix_complex_memcpy(pH2, H2Tmp);
		   gsl_matrix_complex_memcpy(WInv, G_preInv);
		   index2_head->next = create(Nt);
		   index1_head->next = split(index2_head, subset);
	   }else if(maxValueTmp<minValue){
		   minValue = maxValueTmp;
		   gsl_matrix_complex_memcpy(pH1, H1Tmp);
		   gsl_matrix_complex_memcpy(pH2, H2Tmp);
		   gsl_matrix_complex_memcpy(WInv, G_preInv);
		   freeLinkedList(index2_head->next);
		   freeLinkedList(index1_head->next);
		   index2_head->next = create(Nt);
		   index1_head->next = split(index2_head, subset);
	   }
	   gsl_combination_next(subset);
   }

  	gsl_combination_free(subset);
  	gsl_matrix_complex_free(G_pre);
  	gsl_matrix_complex_free(G_preInv);
  	gsl_permutation_free(p);
  	free(signum);
  	gsl_matrix_complex_free(H1Tmp);
  	gsl_matrix_complex_free(H2Tmp);
  	gsl_vector_free(diagMax);
  	gsl_vector_free(diag);
  	gsl_vector_complex_free(colReserve);
	return;
}



#endif /* DIVERSITY_MAX_SELECTION_H_ */
