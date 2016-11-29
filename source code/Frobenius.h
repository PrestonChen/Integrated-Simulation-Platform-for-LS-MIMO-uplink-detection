/*
 * Frobenius.h
 *
 *  Created on: Mar 18, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */
#ifndef FROBENIUS_H_
#define FROBENIUS_H_
double Frobenius(gsl_matrix_complex *W);
double Frobenius(gsl_matrix_complex *W){
	double result;
	int m = W->size1;
	int n = W->size2;
	int count1, count2;
	for (count1 = 0; count1 < m; count1++){
		for (count2 = 0; count2 < n; count2++){
			result += gsl_complex_abs2(gsl_matrix_complex_get(W,count1,count2));
		}
	}
	result = sqrt(result);
	return result;
}


#endif /* FROBENIUS_H_ */
