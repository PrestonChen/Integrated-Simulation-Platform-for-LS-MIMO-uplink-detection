void modulator (gsl_vector_ulong *pdata, gsl_vector_complex *psymbolconstellation, gsl_vector_ulong *pgraydata,
		gsl_vector_complex *ptransmitted, gsl_vector_ulong *pgrayInput);
void modulator (gsl_vector_ulong *pdata, gsl_vector_complex *psymbolconstellation, gsl_vector_ulong *pgraydata,
		gsl_vector_complex *ptransmitted, gsl_vector_ulong *pgrayInput)
{
	int Nt;
	int count=0;
	Nt = pdata->size;
	unsigned long temp;
	for (count=0;count<Nt; count++){
		temp=gsl_vector_ulong_get(pdata, count);
		gsl_vector_complex_set (ptransmitted, count, gsl_vector_complex_get (psymbolconstellation,temp));
	    gsl_vector_ulong_set(pgrayInput, count, gsl_vector_ulong_get(pgraydata, temp));

	}

	return;
}

