/*
 * symErrorCheck.h
 *
 *  Created on: Dec 4, 2015
 *      Author: Preston Chen
 */
/*
 * check the number of symbol errors and the corresponding indexes
 */
#ifndef SYMERRORCHECK_H_
#define SYMERRORCHECK_H_
void symErrorCheck(gsl_vector_complex *ptransmit, gsl_vector_complex *psymOut, int *ErrorIndex_V,  int *symError, int *frameError);
void symErrorCheck(gsl_vector_complex *ptransmit, gsl_vector_complex *psymOut, int *ErrorIndex_V,  int *symError, int *frameError)
{
	int count1;
	int count;
	*frameError = 0;
	*symError = 0;
	int Nt = ptransmit->size;
	gsl_complex pt, ps;
	count = 0;
	for (count1 = 0; count1 < Nt; count1++){
		pt = gsl_vector_complex_get(ptransmit,count1);
		ps = gsl_vector_complex_get(psymOut,count1);
		if (fabs(pt.dat[0] - ps.dat[0]) > 1e-5 || fabs(pt.dat[1] - ps.dat[1]) > 1e-5){
			(*symError)++;
			ErrorIndex_V[count] = count1;
			count++;
		}
	}
	if((*symError) > 0){
		(*frameError)++;
	}
    return;
}




#endif /* SYMERRORCHECK_H_ */
