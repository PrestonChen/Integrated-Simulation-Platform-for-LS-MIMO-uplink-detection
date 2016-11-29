void error_channel_generator (gsl_matrix_complex *pH, gsl_rng *pr, double sigmaerr);
void error_channel_generator (gsl_matrix_complex *pH, gsl_rng *pr, double sigmaerr)
{
	int Nt, Nr, count1=0, count2=0;
	Nr=pH->size1;
	Nt=pH->size2;
	
	gsl_complex z;
	GSL_SET_COMPLEX (&z, 0.0, 0.0);

	for (count1=0 ; count1<Nr ; count1++)
	    {
		for (count2=0; count2<Nt; count2++)
		    { 
			GSL_SET_COMPLEX (&z, gsl_ran_gaussian (pr, sigmaerr), gsl_ran_gaussian (pr, sigmaerr));
			gsl_matrix_complex_set (pH, count1, count2, z);
			//printf("z = %g \n", gsl_ran_gaussian (pr, 1.0));
		     }
	    }
	return;
		   
}
