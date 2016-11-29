void noise_generator (gsl_vector_complex *pnoise, gsl_rng *pr, double noiseV);
void noise_generator (gsl_vector_complex *pnoise, gsl_rng *pr, double noiseV)
{
	int Nr, count=0;
	Nr=pnoise->size;
	
	gsl_complex z;
	for (count=0 ; count<Nr ; count++)
	    {
			GSL_SET_COMPLEX (&z, gsl_ran_gaussian (pr, sqrt(noiseV/2)), gsl_ran_gaussian (pr, sqrt(noiseV/2)));
			gsl_vector_complex_set (pnoise, count, z);
			//printf("z = %g \n", gsl_ran_gaussian (pr, 1.0));
	    }
	return;
		   
}
