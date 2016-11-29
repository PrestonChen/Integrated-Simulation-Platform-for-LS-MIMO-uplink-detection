/*
 * MMSE_SENIA.h
 * This routine peforms MMSE detection with SENIA approximate matrix inversion
 *  Created on: Mar 16, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef MMSE_SENIA_H_
#define MMSE_SENIA_H_
#include "SENIA.h"
void MMSE_SENIA(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
		double snr, double pav, int M, int k, gsl_vector_complex *psymOut);
void MMSE_SENIA(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
		double snr, double pav, int M, int k, gsl_vector_complex *psymOut){
	gsl_complex alpha, beta1,beta2;
		GSL_SET_COMPLEX(&alpha, 1,0);
		GSL_SET_COMPLEX(&beta1, 1 / snr, 0);
		GSL_SET_COMPLEX(&beta2, 0, 0);
		int Nr = pH->size1;
		int Nt = pH->size2;
		gsl_matrix_complex *G_pre, *G_preInv, *G, *G_aided;
		gsl_permutation *p = gsl_permutation_calloc(Nt);
		int *signum = (int*)calloc(1, sizeof(int));
		*signum = 1;
		G_pre = gsl_matrix_complex_calloc(Nt, Nt);
		G_preInv = gsl_matrix_complex_calloc(Nt, Nt);
		G = gsl_matrix_complex_calloc(Nt,Nr);
		G_aided = gsl_matrix_complex_calloc(Nt,Nt);
		gsl_matrix_complex_set_identity(G_pre);
		gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH, pH, beta1, G_pre);
		SENIA(G_pre, k, G_preInv, G_aided);
		gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, G_preInv, pH, beta2, G);
		gsl_blas_zgemv(CblasNoTrans, alpha, G, preceived, beta2, psymOut);
		RectangularQAMSlicer(psymOut, pav, M);
		gsl_matrix_complex_free(G_pre);
		gsl_matrix_complex_free(G_preInv);
		gsl_matrix_complex_free(G);
		gsl_matrix_complex_free(G_aided);
		gsl_permutation_free(p);
		free(signum);

}


#endif /* MMSE_SENIA_H_ */
