void demodulator(gsl_vector_complex *psymOut, gsl_vector_complex *psymbolconstellation, gsl_vector_ulong *pgraydata,  int *ErrorIndex,  gsl_vector_ulong *pgrayOut);
void demodulator(gsl_vector_complex *psymOut, gsl_vector_complex *psymbolconstellation, gsl_vector_ulong *pgraydata,  int *ErrorIndex,  gsl_vector_ulong *pgrayOut)
{
//notice: put the condition judgement ErrorIndex in the main function if ErrorIndex =  = NULL is true, that jump over demodulator and binaryerrors
     int count1, count2;
     int l = pgrayOut->size; //the length of the error symbol vector
     int M = psymbolconstellation->size;
     gsl_complex temp, psymbolconstellation_tmp;
     for (count1 = 0; count1 < l; count1++){
    	 temp = gsl_vector_complex_get(psymOut, ErrorIndex[count1]);
    	 for (count2 = 0; count2 < M; count2++){
    		 psymbolconstellation_tmp = gsl_vector_complex_get(psymbolconstellation, count2);
    		 if(fabs(temp.dat[0] - psymbolconstellation_tmp.dat[0]) < 1e-5 && fabs(temp.dat[1] - psymbolconstellation_tmp.dat[1]) < 1e-5){
    			 gsl_vector_ulong_set(pgrayOut, count1, gsl_vector_ulong_get(pgraydata, count2));
    			 break;
    		 }
    	 }
     }
}

