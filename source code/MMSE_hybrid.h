/*
 * MMSE_hybrid.h
 *
 *  Created on: Mar 23, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef MMSE_HYBRID_H_
#define MMSE_HYBRID_H_
#include "SENIAEMIhybrid.h"
int MMSE_hybrid(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
		double snr, double pav, int M, int k, double threshold, gsl_vector_complex *psymOut);
int MMSE_hybrid(gsl_vector_complex *preceived, gsl_matrix_complex *pH,
		double snr, double pav, int M, int k, double threshold, gsl_vector_complex *psymOut){
	int Nr = pH->size1;
	int Nt = pH->size2;
	int hybridLabel;
	gsl_complex alpha1, alpha2, beta, OM_tmp1Complex;
	GSL_SET_COMPLEX(&alpha1, 1, 0);
	GSL_SET_COMPLEX(&alpha2, (double)1 / snr, 0);
	GSL_SET_COMPLEX(&beta, 0, 0);
	double OM;  //orthogonality measure of this channel realization
	double OM_tmp1;
    double OM_tmp2 = 1;
    int count;
    int signum[1];
	gsl_matrix_complex *G, *G_pre, *G_pre_tmp, *G_aided, *LU, *G_preInv;
	gsl_vector_complex *G_preCol;
    gsl_permutation *p;
    G = gsl_matrix_complex_calloc(Nt,Nr);
    G_pre = gsl_matrix_complex_calloc(Nt,Nt);
    G_pre_tmp = gsl_matrix_complex_calloc(Nt,Nt);
    G_aided = gsl_matrix_complex_calloc(Nt,Nt);
    LU = gsl_matrix_complex_calloc(Nt,Nt);
    G_preInv = gsl_matrix_complex_calloc(Nt,Nt);
    G_preCol = gsl_vector_complex_calloc(Nr);
    p = gsl_permutation_calloc(Nt);
    signum[0] = 1;
	gsl_matrix_complex_set_identity(G_pre);
	gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha1, pH, pH, beta, G_pre);
	//Computation of OM

    gsl_matrix_complex_memcpy(LU, G_pre);
    gsl_linalg_complex_LU_decomp(LU, p, signum);
    OM_tmp1Complex = gsl_linalg_complex_LU_det(LU, signum[0]);
    OM_tmp1 = OM_tmp1Complex.dat[0];
    for (count = 0; count < Nt; count++){
    	gsl_matrix_complex_get_col(G_preCol, pH, count);
    	OM_tmp2 *= pow(gsl_blas_dznrm2(G_preCol), 2);
    }
    OM = OM_tmp1 / OM_tmp2;
    //MMSE detection
	gsl_matrix_complex_set_identity(G_pre_tmp);
	gsl_matrix_complex_scale(G_pre_tmp, alpha2);
	gsl_matrix_complex_add(G_pre, G_pre_tmp);
	hybridLabel = SENIAEMIhybrid(G_pre, k, OM, threshold, G_preInv, G_aided);
	gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha1, G_preInv, pH, beta, G);
	gsl_blas_zgemv(CblasNoTrans, alpha1, G, preceived, beta, psymOut);
	RectangularQAMSlicer(psymOut, pav, M);
	gsl_matrix_complex_free(G);
	gsl_matrix_complex_free(G_pre);
	gsl_matrix_complex_free(G_pre_tmp);
	gsl_matrix_complex_free(G_aided);
	gsl_matrix_complex_free(LU);
	gsl_matrix_complex_free(G_preInv);
	gsl_vector_complex_free(G_preCol);
    gsl_permutation_free(p);
    return hybridLabel;

}


#endif /* MMSE_HYBRID_H_ */
