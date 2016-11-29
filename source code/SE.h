/*
 * SE.h
 * This routine performs the Neumann series expansion matrix inversion
 *  Created on: Mar 16, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef SE_H_
#define SE_H_

//Declaration
void SE(gsl_matrix_complex *Win, int k, gsl_matrix_complex *Wout,
		gsl_matrix_complex *Waided);
void SE(gsl_matrix_complex *Win, int k, gsl_matrix_complex *Wout,
		gsl_matrix_complex *Waided){
	 int Nt = Win->size1;
	 gsl_complex alpha1, alpha2,  beta, W_diagInv_tmp;;
	 GSL_SET_COMPLEX(&alpha1, 1, 0);
	 GSL_SET_COMPLEX(&alpha2, -1, 0);
	 GSL_SET_COMPLEX(&beta, 0, 0);
	 gsl_matrix_complex *W_offdiag, *W_InvTmp,
	 *W_Inter1, *Wout_Tmp, *W_Inter2;
	 gsl_vector_complex *W_diag, *W_Inter1_row, *W_diagInv, *W_Inter2_col;
	 W_offdiag = gsl_matrix_complex_calloc(Nt, Nt);  //off diagonal matrix of Win
	 W_InvTmp = gsl_matrix_complex_calloc(Nt, Nt); //temporary matrix in iteration
	 W_Inter1 = gsl_matrix_complex_calloc(Nt, Nt); // (-D^(-1)E)
	 Wout_Tmp = gsl_matrix_complex_calloc(Nt,Nt); // D^(-1) matrix
	 W_Inter2 = gsl_matrix_complex_calloc(Nt, Nt); //(-D^(-1)E)D^(-1)
	 W_diag = gsl_vector_complex_calloc(Nt); //D vector
	 W_Inter1_row = gsl_vector_complex_calloc(Nt); //each row of (-D^(-1)E)
	 W_diagInv = gsl_vector_complex_calloc(Nt);  //D^(-1) vector
	 W_Inter2_col = gsl_vector_complex_calloc(Nt); //each col of (-D^(-1)E)D^(-1)
     int count;
     //get the diagonal vector and off diagonal matrix
     gsl_matrix_complex_memcpy(W_offdiag, Win);
     for (count = 0; count < Nt; count++){
    	 gsl_vector_complex_set(W_diag, count,
    			 gsl_matrix_complex_get(Win, count, count));
    	 gsl_matrix_complex_set(W_offdiag, count, count, beta);

     }
     //get D^(-1) vector
     for (count = 0; count < Nt; count++){
    	 W_diagInv_tmp=gsl_complex_div(alpha1,
    			 gsl_vector_complex_get(W_diag,count));
//    	 W_diagInv_tmp=gsl_complex_mul(W_diagInv_tmp, alpha2);
    	 gsl_vector_complex_set(W_diagInv, count, W_diagInv_tmp);
     }
     //compute -D^(-1)E
     for (count = 0; count < Nt; count++){
    	 gsl_matrix_complex_get_row(W_Inter1_row, W_offdiag, count);
    	 gsl_vector_complex_scale(W_Inter1_row,
    			 gsl_vector_complex_get(W_diagInv,count));
    	 gsl_vector_complex_scale(W_Inter1_row, alpha2);
    	 gsl_matrix_complex_set_row(W_Inter1, count, W_Inter1_row);

     }
     gsl_matrix_complex_memcpy(Waided, W_Inter1);
    //initialize Wout=D^(-1)
     for (count = 0; count < Nt; count++){
    	gsl_matrix_complex_set(Wout_Tmp, count, count,
    			gsl_vector_complex_get(W_diagInv,count));
     }
     gsl_matrix_complex_memcpy(Wout,Wout_Tmp);
     //The case k=1
     if (k == 1){
         gsl_matrix_complex_free (W_offdiag);
         gsl_matrix_complex_free(W_InvTmp);
         gsl_matrix_complex_free(W_Inter1);
         gsl_matrix_complex_free(Wout_Tmp);
         gsl_matrix_complex_free(W_Inter2);
         gsl_vector_complex_free(W_diag);
         gsl_vector_complex_free(W_Inter1_row);
         gsl_vector_complex_free(W_diagInv);
         gsl_vector_complex_free(W_Inter2_col);
    	 return;
     }
     //compute (-D^(-1)E)D^(-1)
     for (count = 0; count < Nt; count++){
    	 gsl_matrix_complex_get_col(W_Inter2_col, W_Inter1, count);
    	 gsl_vector_complex_scale(W_Inter2_col,
    			 gsl_vector_complex_get(W_diagInv, count));
    	 gsl_matrix_complex_set_col(W_Inter2, count, W_Inter2_col);
     }
     gsl_matrix_complex_add(Wout, W_Inter2);
     //The case k=2
     if (k == 2){
         gsl_matrix_complex_free (W_offdiag);
         gsl_matrix_complex_free(W_InvTmp);
         gsl_matrix_complex_free(W_Inter1);
         gsl_matrix_complex_free(Wout_Tmp);
         gsl_matrix_complex_free(W_Inter2);
         gsl_vector_complex_free(W_diag);
         gsl_vector_complex_free(W_Inter1_row);
         gsl_vector_complex_free(W_diagInv);
         gsl_vector_complex_free(W_Inter2_col);
    	 return;
     }
     //The case k>=3
     //compute Wout iteratively
     for(count = 2; count < k; count++){
    	gsl_matrix_complex_memcpy(W_InvTmp, Wout_Tmp);
    	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha1,
    			W_Inter1, Wout, alpha1, W_InvTmp);
    	gsl_matrix_complex_memcpy(Wout,W_InvTmp);

     }
     gsl_matrix_complex_free (W_offdiag);
     gsl_matrix_complex_free(W_InvTmp);
     gsl_matrix_complex_free(W_Inter1);
     gsl_matrix_complex_free(Wout_Tmp);
     gsl_matrix_complex_free(W_Inter2);
     gsl_vector_complex_free(W_diag);
     gsl_vector_complex_free(W_Inter1_row);
     gsl_vector_complex_free(W_diagInv);
     gsl_vector_complex_free(W_Inter2_col);
}
#endif /* SE_H_ */
