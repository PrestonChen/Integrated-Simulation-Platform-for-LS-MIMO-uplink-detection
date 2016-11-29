/*
 * MMSE_SIC_ordering.h
 *  This routine considers the MMSE-SIC algorithm with various ordering strategies: V-BLAST,
 *   LLR reliability ordering, simplified LLR reliability ordering
 *  Created on: May 9, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#include "linkedList.h"
#ifndef MMSE_SIC_ORDERING_H_
#define MMSE_SIC_ORDERING_H_
void MMSE_SIC_ordering(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
       double snr, double pav, int M, gsl_vector_complex *psymOut);
void MMSE_SIC_ordering(gsl_vector_complex *preceived, gsl_matrix_complex *pH, double snr,
		double pav, int M, gsl_vector_complex *psymOut){
   int Nr=pH->size1;
   int Nt=pH->size2;
   struct Node *head=create(1);
   head->next=create(Nt);
   gsl_matrix_complex *pH_inter=gsl_matrix_complex_calloc(Nr, Nt);
   gsl_matrix_complex_memcpy(pH_inter, pH);
   gsl_vector_complex *preceive_tmp=gsl_vector_complex_calloc(Nr);
   gsl_vector_complex_memcpy(preceive_tmp, preceived);
   gsl_vector_complex_view diag_viewComplex;
   gsl_vector_view diag_viewReal;
   gsl_vector *diag;
   gsl_matrix_complex *pHtemp, *G_pre, *G_preInv, *row_M;
   gsl_permutation *p;
   int *signum=(int*)calloc(1, sizeof(int));
   *signum=1;
   gsl_matrix_complex  *Gmmse=gsl_matrix_complex_calloc(1, Nr);
   gsl_vector_complex  *GmmseR, *colNulling, *colReserve;
   gsl_vector_complex *row;
   gsl_vector_complex *psymOut_tmpTmp;
   GmmseR=gsl_vector_complex_calloc(Nr);
   colNulling=gsl_vector_complex_calloc(Nr);
   colReserve=gsl_vector_complex_calloc(Nr);
   gsl_complex symCurrent;
   gsl_vector_complex *symCurrent_V=gsl_vector_complex_calloc(1);
   gsl_complex alpha, beta1, beta2, beta3, psymOut_tmpEle;
   GSL_SET_COMPLEX(&alpha, 1, 0);
   GSL_SET_COMPLEX(&beta1, pow(snr, -1), 0);
   GSL_SET_COMPLEX(&beta2, 0, 0);
   GSL_SET_COMPLEX(&beta3, -1, 0);
   int k;
   gsl_matrix_complex *Gequal;
   gsl_vector_complex *psymOut_tmp;
#if defined (MMSEIOSIC)
   double diag_tmp;
   Gequal=gsl_matrix_complex_calloc(Nt, Nr);
   psymOut_tmp=gsl_vector_complex_calloc(Nt);
#endif
   int count, count1, count2;
   int index;
//   double temp;
   for (count=0;count<Nt;count++){
	   G_pre=gsl_matrix_complex_calloc(Nt-count, Nt-count);
	   G_preInv=gsl_matrix_complex_calloc(Nt-count, Nt-count);
	   gsl_matrix_complex_set_identity(G_pre);
	   p=gsl_permutation_calloc(Nt-count);
	   diag=gsl_vector_calloc(Nt-count);
	   row=gsl_vector_complex_calloc(Nt-count);
	   row_M=gsl_matrix_complex_calloc(1, Nt-count);
	   gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH_inter,
			   pH_inter, beta1, G_pre);
	   gsl_linalg_complex_LU_decomp(G_pre, p, signum);
	   gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);//calculate the inverse matrix
	   diag_viewComplex=gsl_matrix_complex_diagonal(G_preInv);
	   diag_viewReal=gsl_vector_complex_real(&diag_viewComplex.vector);
	   gsl_vector_memcpy(diag, &diag_viewReal.vector);
#if defined(MMSESIC)
	   k=0;   //without any ordering
	   gsl_matrix_complex_get_row(row, G_preInv, k);
	   gsl_matrix_complex_set_row(row_M, 0, row);
	   gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, row_M, pH_inter, beta2, Gmmse);
	   gsl_matrix_complex_get_row(GmmseR, Gmmse,0);
	   gsl_blas_zdotu(GmmseR, preceive_tmp, &symCurrent);
	   gsl_vector_complex_set(symCurrent_V, 0, symCurrent);  //estimation
	   RectangularQAMSlicer(symCurrent_V, pav, M);
//#elif defined(MMSEOSIC)||defined(MGDOSIC)||defined(MGDSENIAIUOSIC)
#elif defined(MMSEVBLASTSIC)
	   k=gsl_vector_min_index(diag);   //the current index of the strongest post processing SINR
	   gsl_matrix_complex_get_row(row, G_preInv, k);
	   gsl_matrix_complex_set_row(row_M, 0, row);
	   gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, row_M, pH_inter, beta2, Gmmse);
       gsl_matrix_complex_get_row(GmmseR, Gmmse,0);
       gsl_blas_zdotu(GmmseR, preceive_tmp, &symCurrent);
	   gsl_vector_complex_set(symCurrent_V, 0, symCurrent);  //estimation
	   RectangularQAMSlicer(symCurrent_V, pav, M);
#elif defined (MMSEIOSIC)
// 	   LLR=gsl_vector_calloc(Nt-count);
 	   gsl_matrix_complex_free(Gequal);
 	   gsl_vector_complex_free(psymOut_tmp);
 	   Gequal=gsl_matrix_complex_calloc(Nt-count, Nr);
 	   psymOut_tmp=gsl_vector_complex_calloc(Nt-count);
       gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, G_preInv, pH_inter, beta2, Gequal);
       gsl_blas_zgemv(CblasNoTrans, alpha, Gequal, preceive_tmp, beta2, psymOut_tmp);
	   for (count1=0;count1<(Nt-count);count1++){
		   psymOut_tmpEle=gsl_vector_complex_get(psymOut_tmp, count1);
		   diag_tmp=(fabs(GSL_REAL(psymOut_tmpEle))+fabs(GSL_IMAG(psymOut_tmpEle)))
				   /(gsl_vector_get(diag, count1));
		   gsl_vector_set(diag, count1, diag_tmp);
	   }
       //the current index of the data stream with the largest reliability
	   k=gsl_vector_max_index(diag);      //the index of the symbol to be detected
//	   gsl_vector_free(LLR);
	   gsl_vector_complex_set(symCurrent_V, 0, gsl_vector_complex_get(psymOut_tmp,k));
	   RectangularQAMSlicer(symCurrent_V, pav, M);
#endif
	   index=get((int)k, head);
	   //reconstruct the desicion symbol vector
	   gsl_vector_complex_set(psymOut, index, gsl_vector_complex_get(symCurrent_V, 0));
	   //update the observatted signal vector
	   gsl_matrix_complex_get_col(colNulling, pH_inter, k);
	   gsl_vector_complex_scale(colNulling, gsl_vector_complex_get(symCurrent_V, 0));
	   gsl_vector_complex_sub(preceive_tmp, colNulling);
	   //update channel matrix and the soft estimate psymOut_tmpEle
	   if((Nt-count)==1){
		   gsl_matrix_complex_free(pH_inter);
		   gsl_matrix_complex_free(G_pre);
		   gsl_matrix_complex_free(G_preInv);
		   gsl_permutation_free(p);
		   gsl_vector_free(diag);
		   gsl_vector_complex_free(row);
		   gsl_matrix_complex_free(row_M);
		   break;
	   }
	   pHtemp=gsl_matrix_complex_calloc(Nr, Nt-count-1);
		count2=0;
		for (count1=0;count1<(Nt-count);count1++){
		  if (count1==k){
			  continue;
		  }
		  gsl_matrix_complex_get_col(colReserve, pH_inter, count1);
		  gsl_matrix_complex_set_col(pHtemp, count2, colReserve);
		  count2++;
		  }
		gsl_matrix_complex_free(pH_inter);
		pH_inter=gsl_matrix_complex_calloc(Nr, Nt-count-1);
		gsl_matrix_complex_memcpy(pH_inter, pHtemp);
	    gsl_matrix_complex_free(pHtemp);
	   gsl_matrix_complex_free(G_pre);
	   gsl_matrix_complex_free(G_preInv);
	   gsl_permutation_free(p);
	   gsl_vector_free(diag);
	   gsl_vector_complex_free(row);
	   gsl_matrix_complex_free(row_M);



   }
   free(head);
   gsl_vector_complex_free(preceive_tmp);
   free(signum);
   gsl_vector_complex_free(GmmseR);
   gsl_vector_complex_free(colNulling);
   gsl_vector_complex_free(colReserve);
   gsl_matrix_complex_free(Gmmse);
   gsl_vector_complex_free(symCurrent_V);
#if defined (MMSEIOSIC)
   gsl_matrix_complex_free(Gequal);
   gsl_vector_complex_free(psymOut_tmp);
#endif
	return;
}


#endif /* MMSE_SIC_ORDERING_H_ */
