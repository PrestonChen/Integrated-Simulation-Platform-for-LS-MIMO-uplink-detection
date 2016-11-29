void grayencoder (gsl_vector_ulong *pgraydata);
void grayencoder (gsl_vector_ulong *pgraydata)
{
	int Q=pgraydata->size;
    int count, count1,count2;
    int pgraydata_tmp;
	int M=0;
	M = (int) ceil(sqrt(Q));
	gsl_vector_ulong *pPamdata;
	pPamdata = gsl_vector_ulong_calloc (M);
	for (count=0; count<M; count++){
		count1=count;
		count2 = count;
		count2 = count2 >> 1;
		count1 = count1^count2;
		gsl_vector_ulong_set (pPamdata, count, count1);
	  }
	for (count=0; count<Q; count++){
		pgraydata_tmp = M*gsl_vector_ulong_get(pPamdata, count/M) + gsl_vector_ulong_get (pPamdata, count%M);
		gsl_vector_ulong_set (pgraydata, count, pgraydata_tmp);

	  }


	gsl_vector_ulong_free (pPamdata);

	return;
}
