/*
 * MMSE.h
 *
 *  Created on: Dec 4, 2015
 *      Author: Preston Chen
 */

#ifndef MMSE_H_
#define MMSE_H_
void MMSE(gsl_vector_complex *preceived, gsl_matrix_complex *pH, double snr, double pav, int M, gsl_vector_complex *psymOut);
void MMSE(gsl_vector_complex *preceived, gsl_matrix_complex *pH, double snr, double pav, int M, gsl_vector_complex *psymOut)
{// notice the snr is the symbol snr
gsl_complex alpha, beta1,beta2;
GSL_SET_COMPLEX(&alpha, 1,0);
GSL_SET_COMPLEX(&beta1, 1/snr, 0);
GSL_SET_COMPLEX(&beta2, 0, 0);
int Nr = pH->size1;
int Nt = pH->size2;
gsl_matrix_complex *G_pre, *G_preInv, *G;
gsl_permutation *p = gsl_permutation_calloc(Nt);
int *signum = (int*)calloc(1, sizeof(int));
*signum = 1;
G_pre = gsl_matrix_complex_calloc(Nt, Nt);
G_preInv = gsl_matrix_complex_calloc(Nt, Nt);
G = gsl_matrix_complex_calloc(Nt,Nr);
gsl_matrix_complex_set_identity(G_pre);
gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, alpha, pH, pH, beta1, G_pre);
//gsl_matrix_complex *LU = gsl_matrix_complex_calloc(Nt, Nt);
//gsl_matrix_complex_memcpy(LU, G_pre);
//gsl_linalg_complex_LU_decomp(LU, p, signum);
gsl_linalg_complex_LU_decomp(G_pre, p, signum);
gsl_linalg_complex_LU_invert(G_pre, p, G_preInv);
//gsl_linalg_complex_LU_invert(LU, p, G_preInv);
//gsl_matrix_complex *I = gsl_matrix_complex_calloc(Nt, Nt);
//gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, G_preInv, G_pre, beta2, I);
//int count1, count2;
//printf("The Identity matrix in MMSE is\n");
//for (count1 = 0;count1<Nt; count1++){
//	for(count2 = 0;count2<Nt;count2++){
//		printf("%g+i%g, ", gsl_matrix_complex_get(I, count1, count2).dat[0],
//				gsl_matrix_complex_get(I, count1, count2).dat[1]);
//	}
//	printf("\n");
//}
//printf("\n");
gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, alpha, G_preInv, pH, beta2, G);
gsl_blas_zgemv(CblasNoTrans, alpha, G, preceived, beta2, psymOut);
RectangularQAMSlicer(psymOut, pav, M);
gsl_matrix_complex_free(G_pre);
gsl_matrix_complex_free(G_preInv);
gsl_matrix_complex_free(G);
gsl_permutation_free(p);
free(signum);
}




#endif /* MMSE_H_ */
