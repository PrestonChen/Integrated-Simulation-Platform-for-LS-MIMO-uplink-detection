/*
 * MMSE_IU.h
 * This routine performs the MMSE detection with matrix inversion inflate update
 *  Created on: Mar 16, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */


#ifndef MMSE_IU_H_
#define MMSE_IU_H_
#include "IU.h"
void MMSE_IU(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
	 double snr, double pav, int M, int L,gsl_vector_complex *psymOut);
void MMSE_IU(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
	 double snr, double pav, int M, int L,gsl_vector_complex *psymOut){
	gsl_complex alpha, beta1,beta2;
				GSL_SET_COMPLEX(&alpha, 1,0);
				GSL_SET_COMPLEX(&beta1, 1/snr, 0);
				GSL_SET_COMPLEX(&beta2, 0, 0);
				int Nr = pH->size1;
				int Nt = pH->size2;
				int count;
				gsl_matrix_complex *G_pre, *G_preInv, *G_nonofuse, *G_preInvIU, *G;
				gsl_matrix_complex *pH_sub1, *pH_sub2;
				gsl_vector_complex *pH_col_tmp;
				gsl_permutation *p = gsl_permutation_calloc(L);
				int *signum = (int*)calloc(1, sizeof(int));
				*signum = 1;
				G_pre = gsl_matrix_complex_calloc(L, L); //H1'H1
				G_preInv = gsl_matrix_complex_calloc(L, L);// (H1'H1+rho^(-1)I)^(-1)
				G_nonofuse = gsl_matrix_complex_calloc(L, L);
				G_preInvIU = gsl_matrix_complex_calloc(Nt, Nt);
				G = gsl_matrix_complex_calloc(Nt,Nr);
				pH_sub1 = gsl_matrix_complex_calloc(Nr, L);
				pH_sub2 = gsl_matrix_complex_calloc(Nr, Nt-L);
				pH_col_tmp = gsl_vector_complex_calloc(Nr);
				for(count = 0; count < L; count++){
					gsl_matrix_complex_get_col(pH_col_tmp, pH, count);
					gsl_matrix_complex_set_col(pH_sub1, count, pH_col_tmp);
				}
				for (count = 0;count < (Nt-L); count++){
					gsl_matrix_complex_get_col(pH_col_tmp, pH, count+L);
					gsl_matrix_complex_set_col(pH_sub2, count, pH_col_tmp);
				}
				gsl_matrix_complex_set_identity(G_pre);
				gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH_sub1, pH_sub1, beta1, G_pre);
				gsl_linalg_complex_LU_decomp(G_pre, p, signum);
				gsl_linalg_complex_LU_invert(G_pre, p , G_preInv);
				IU(G_preInv, pH_sub1, pH_sub2, snr, G_preInvIU);
				gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, G_preInvIU, pH, beta2, G);
				gsl_blas_zgemv(CblasNoTrans, alpha, G, preceived, beta2, psymOut);
				RectangularQAMSlicer(psymOut, pav, M);
				gsl_matrix_complex_free(G_pre);
				gsl_matrix_complex_free(G_preInv);
				gsl_matrix_complex_free(G_nonofuse);
				gsl_matrix_complex_free(G_preInvIU);
				gsl_matrix_complex_free(G);
				gsl_matrix_complex_free(pH_sub1);
				gsl_matrix_complex_free(pH_sub2);
				gsl_vector_complex_free(pH_col_tmp);
				gsl_permutation_free(p);
				free(signum);
}


#endif /* MMSE_IU_H_ */
