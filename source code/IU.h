/*
 * IU.h
 * This routine performs the inflate update of matrix inversion
 *  Created on: Mar 16, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef IU_H_
#define IU_H_
void IU(gsl_matrix_complex *Win, gsl_matrix_complex *pH1,
		gsl_matrix_complex *pH2, double snr, gsl_matrix_complex *Wout);
void IU(gsl_matrix_complex *Win, gsl_matrix_complex *pH1,
		gsl_matrix_complex *pH2, double snr, gsl_matrix_complex *Wout){
	int Nr=pH1->size1;
	int L=pH1->size2;
	int ReL=pH2->size2;
	int count, count1;
	gsl_complex alpha, beta, eta, c_sub2_complex, c_complex;
	GSL_SET_COMPLEX(&alpha, 1, 0);
	GSL_SET_COMPLEX(&beta, 0, 0);
	GSL_SET_COMPLEX(&eta, -1, 0);
	gsl_matrix_complex *pH_tmp, *W_current;
	gsl_matrix_complex *F11_Inv, *W_next, *W_next_sub1,*pH_tmpUp;
	gsl_vector_complex *T1, *T2, *T3, *T3_conj, *W_next_sub2;
	gsl_vector_complex *W_next_sub1_row, *W_next_col;
	gsl_vector_complex *H_col, *pH_tmp_col;
	double c, c_sub1, c_sub2;
	pH_tmp=gsl_matrix_complex_calloc(Nr, L);
	W_current=gsl_matrix_complex_calloc(L, L);
	gsl_matrix_complex_memcpy(pH_tmp, pH1);
	gsl_matrix_complex_memcpy(W_current, Win);
	H_col=gsl_vector_complex_calloc(Nr); //column h_(n) extracted from H2
    pH_tmp_col=gsl_vector_complex_calloc(Nr);

	//(Nt-L) times of iterative inflate update
	for (count=0;count<ReL;count++){
		T1=gsl_vector_complex_calloc(L+count);
		T2=gsl_vector_complex_calloc(L+count);
		T3=gsl_vector_complex_calloc(L+count);
		T3_conj=gsl_vector_complex_calloc(L+count);
		F11_Inv=gsl_matrix_complex_calloc(L+count , L+count);
		W_next=gsl_matrix_complex_calloc(L+count+1, L+count+1);
		W_next_sub1=gsl_matrix_complex_calloc(L+count+1, L+count);
		W_next_sub2=gsl_vector_complex_calloc(L+count+1);
		W_next_sub1_row=gsl_vector_complex_calloc(L+count);
		W_next_col=gsl_vector_complex_calloc(L+count+1);
		pH_tmpUp=gsl_matrix_complex_calloc(Nr, L+count);
		gsl_matrix_complex_get_col(H_col, pH2, count);
		gsl_blas_zgemv(CblasConjTrans, alpha, pH_tmp, H_col, beta, T1);
		gsl_blas_zgemv(CblasNoTrans, alpha, W_current, T1, beta, T2);
	    c_sub1=gsl_blas_dznrm2(H_col);
	    c_sub1=pow(c_sub1,2);
	    gsl_blas_zdotc(T1, T2, &c_sub2_complex);
	    c_sub2=GSL_REAL(c_sub2_complex);
	    c=(double)1/(c_sub1+((double)1/snr)-c_sub2);
	    gsl_vector_complex_memcpy(T3, T2);
	    GSL_SET_COMPLEX(&c_complex, -c, 0);
	    gsl_vector_complex_scale(T3, c_complex);
        gsl_matrix_complex_memcpy(F11_Inv, W_current);
        gsl_blas_zgerc(eta, T3, T2, F11_Inv);
        //get T3 conjugate
        for(count1=0;count1<(L+count);count1++){
        	gsl_vector_complex_set(T3_conj, count1,
        			gsl_complex_conjugate(gsl_vector_complex_get(T3, count1)));
        }
        //get W_next_sub1
        for(count1=0;count1<(L+count+1);count1++){
        	if(count1==(L+count)){
        	gsl_matrix_complex_set_row(W_next_sub1, count1, T3_conj);
        	}else{
        	gsl_matrix_complex_get_row(W_next_sub1_row, F11_Inv,count1);
        	gsl_matrix_complex_set_row(W_next_sub1, count1, W_next_sub1_row);
        	}
        }
        //get W_next_sub2
        for (count1=0;count1<(L+count+1);count1++){
        	if(count1==(L+count)){
        		gsl_vector_complex_set(W_next_sub2, count1,
        				gsl_complex_mul(eta, c_complex));
        	}else{
        		gsl_vector_complex_set(W_next_sub2, count1,
        				gsl_vector_complex_get(T3, count1));
        	}
        }

        //construct W_next
        for (count1=0;count1<(L+count+1);count1++){
        	if(count1==(L+count)){
        		gsl_matrix_complex_set_col(W_next,count1, W_next_sub2);
        	}else{
               gsl_matrix_complex_get_col(W_next_col, W_next_sub1,count1);
               gsl_matrix_complex_set_col(W_next, count1, W_next_col);
        	}
        }
        //update W_current
        gsl_matrix_complex_free(W_current);
        W_current=gsl_matrix_complex_calloc(L+count+1, L+count+1);
        gsl_matrix_complex_memcpy(W_current, W_next);
        //update pH_tmp

        gsl_matrix_complex_memcpy(pH_tmpUp, pH_tmp);
        gsl_matrix_complex_free(pH_tmp);
        pH_tmp=gsl_matrix_complex_calloc(Nr, L+count+1);
        for(count1=0;count1<(L+count+1);count1++){
        	if(count1==(L+count)){
        		gsl_matrix_complex_set_col(pH_tmp, count1, H_col);
        	}else{
        		gsl_matrix_complex_get_col(pH_tmp_col, pH_tmpUp, count1);
        		gsl_matrix_complex_set_col(pH_tmp, count1, pH_tmp_col);
        	}
        }

      gsl_vector_complex_free(T1);
      gsl_vector_complex_free(T2);
      gsl_vector_complex_free(T3);
      gsl_vector_complex_free(T3_conj);
      gsl_matrix_complex_free(F11_Inv);
      gsl_matrix_complex_free(W_next);
      gsl_matrix_complex_free(W_next_sub1);
      gsl_vector_complex_free(W_next_sub2);
      gsl_vector_complex_free(W_next_sub1_row);
      gsl_vector_complex_free(W_next_col);
      gsl_matrix_complex_free(pH_tmpUp);



	}
    gsl_matrix_complex_memcpy(Wout, W_current);
	gsl_matrix_complex_free(pH_tmp);
	gsl_matrix_complex_free(W_current);
	gsl_vector_complex_free(H_col);
	gsl_vector_complex_free(pH_tmp_col);

}


#endif /* IU_H_ */
