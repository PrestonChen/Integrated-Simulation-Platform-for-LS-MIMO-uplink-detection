void matrix_sqrt (gsl_matrix_complex *pstart, gsl_matrix_complex *presult);
void matrix_sqrt (gsl_matrix_complex *pstart, gsl_matrix_complex *presult)
{
    int n, count=0;
    n = pstart -> size2;


	gsl_complex z;
	gsl_complex alpha;
	gsl_complex beta;

	GSL_SET_COMPLEX (&alpha, 1.0, 0.0);
	GSL_SET_COMPLEX (&beta, 0.0, 0.0);
	GSL_SET_COMPLEX (&z, 0.0, 0.0);

	gsl_matrix_complex *ptemp, *ptemp2, *ptemp3;
	ptemp = gsl_matrix_complex_calloc (n,n);
	ptemp2 = gsl_matrix_complex_calloc (n,n);
	ptemp3 = gsl_matrix_complex_calloc (n,n);

	gsl_matrix_complex_memcpy (ptemp, pstart);

	gsl_vector *eval = gsl_vector_alloc (n);
       gsl_matrix_complex *evec = gsl_matrix_complex_alloc (n, n);
     
       gsl_eigen_hermv_workspace *w = 
         gsl_eigen_hermv_alloc (n);
       
       gsl_eigen_hermv (ptemp, eval, evec, w);

	for (count=0; count<n; count++)
	   {
		GSL_SET_COMPLEX (&z, sqrt(gsl_vector_get (eval, count)), 0.0);
	    	gsl_matrix_complex_set (ptemp2, count, count, z);
	   }
     
	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, evec, ptemp2, beta, ptemp3);
	gsl_blas_zgemm (CblasNoTrans, CblasConjTrans, alpha, ptemp3, evec, beta, ptemp);

	gsl_matrix_complex_memcpy (presult, ptemp);

	gsl_blas_zgemm (CblasNoTrans, CblasNoTrans, alpha, ptemp, ptemp, beta, ptemp3);

	//matrixviewer("pstart",pstart);
	//matrixviewer("presult", presult);
	//matrixviewer("ptemp3", ptemp3);

	gsl_matrix_complex_free (ptemp); 
	gsl_matrix_complex_free (ptemp2);
	gsl_matrix_complex_free (ptemp3);

	  gsl_eigen_hermv_free (w);
	     gsl_vector_free (eval);
       gsl_matrix_complex_free (evec);

	return;
    
}



