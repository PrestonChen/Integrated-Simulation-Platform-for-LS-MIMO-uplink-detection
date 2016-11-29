void binaryerrors (gsl_vector_ulong *pgrayOut,  int *ErrorIndex, gsl_vector_ulong *pgrayInput,  int Q, int *errors)
{
	int  count1, count2;
	int M  = (int)log2(Q);
	int l = pgrayOut->size;  // the number of the error gray code
	int temp1, temp2;
	int currentBit;
	(*errors) = 0;
	for (count1 = 0; count1 < l; count1++){
	  temp1 = gsl_vector_ulong_get(pgrayOut, count1);
	  temp2 = gsl_vector_ulong_get(pgrayInput, ErrorIndex[count1]);
	  temp1 = temp1 ^ temp2;
	  for (count2 = 0; count2 < M; count2++){
		  currentBit = (int)temp1 % 2;
		  temp1 = temp1 >> 1;
		  if (currentBit == 1){
			  (*errors)++;
		  }
	  }
	}
	return;
}
