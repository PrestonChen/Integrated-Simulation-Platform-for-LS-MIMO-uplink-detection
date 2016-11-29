/*
 * SENIA.h
 * This routine performs SENIA-K approximate matrix inversions
 *  Created on: Mar 16, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef SENIA_H_
#define SENIA_H_
#include "SE.h"
//Declaration
void SENIA(gsl_matrix_complex *Win, int k, gsl_matrix_complex *Wout,
		gsl_matrix_complex *Waided);

void SENIA(gsl_matrix_complex *Win, int k, gsl_matrix_complex *Wout,
		gsl_matrix_complex *Waided){
	 int Nt = Win->size1;
	 int count;
     gsl_complex alpha1, alpha2,  beta1, beta2;
     GSL_SET_COMPLEX(&alpha1, 1, 0);
     GSL_SET_COMPLEX(&alpha2, -1, 0);
     GSL_SET_COMPLEX(&beta1, 2, 0);
     GSL_SET_COMPLEX(&beta2, 0, 0);
     gsl_matrix_complex *M_current, *M_Inter1, *M_Inter2, *M_Inter3;
     M_current = gsl_matrix_complex_calloc(Nt, Nt);
     M_Inter1 = gsl_matrix_complex_calloc(Nt,Nt);
     M_Inter2 = gsl_matrix_complex_calloc(Nt,Nt);
     M_Inter3 = gsl_matrix_complex_calloc(Nt,Nt);
     //get M0 and (-D^(-1)E)D^(-1)
     SE(Win,2, M_current, Waided);
     //compute Mk
     for (count = 0; count < k; count++){
    	 gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
    			 alpha2, M_current, Win, beta2, M_Inter1);
		gsl_matrix_complex_set_identity(M_Inter2);
		gsl_matrix_complex_scale(M_Inter2, beta1);
		gsl_matrix_complex_add(M_Inter2, M_Inter1);
		gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
		         alpha1, M_Inter2, M_current, beta2, M_Inter3);
		gsl_matrix_complex_memcpy(M_current, M_Inter3);

     }
     gsl_matrix_complex_memcpy(Wout, M_current);

     gsl_matrix_complex_free(M_current);
     gsl_matrix_complex_free(M_Inter1);
     gsl_matrix_complex_free(M_Inter2);
     gsl_matrix_complex_free(M_Inter3);

}


#endif /* SENIA_H_ */
