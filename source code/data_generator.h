void data_generator(gsl_vector_ulong *pdata, gsl_rng *pr, int M);
void data_generator(gsl_vector_ulong *pdata, gsl_rng *pr, int M)
{
	int Nt;
	Nt = pdata -> size;	
	int count = 0;
	for (count = 0; count < Nt ; count++)
	  {
		gsl_vector_ulong_set(pdata, count, gsl_rng_uniform_int (pr,M));
	  }
	return;
}
