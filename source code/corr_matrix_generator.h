#include "matrix_sqrt.h"
void corr_matrix_generator (gsl_matrix_complex *pRt, gsl_matrix_complex *pRr, double Rt, double Rr);
void corr_matrix_generator (gsl_matrix_complex *pRt, gsl_matrix_complex *pRr, double Rt, double Rr)
{
	int Nt, Nr, count1 = 0, count2 = 0;
	Nr = pRr->size1;
	Nt = pRt->size1;
	gsl_matrix_complex *ptemp, *ptemp2;
	ptemp = gsl_matrix_complex_calloc(Nr, Nr);
	ptemp2 = gsl_matrix_complex_calloc(Nt, Nt);
	double temp = 0.0;
	gsl_complex z;
	GSL_SET_COMPLEX (&z, 0.0, 0.0);
	for (count1 = 0 ; count1 < Nr ; count1++)
	    {
		GSL_SET_COMPLEX (&z, 1.0, 0.0);
		gsl_matrix_complex_set (pRr, count1, count1, z);
		for (count2 = count1 + 1; count2 < Nr; count2++)
		   {
			temp = count2 - count1;
			temp = temp * temp;
			GSL_SET_COMPLEX (&z, pow (Rr, temp), 0.0);
			gsl_matrix_complex_set (pRr, count1, count2, z);
			gsl_matrix_complex_set (pRr, count2, count1, z);
		   }

	    }
	//matrixviewer ("pRr", pRr);
	gsl_matrix_complex_memcpy (ptemp, pRr);
	matrix_sqrt(ptemp, pRr);
	//matrixviewer ("pRr", pRr);
	for (count1 = 0; count1 < Nt; count1++)
	    {
		GSL_SET_COMPLEX (&z, 1.0, 0.0);
		gsl_matrix_complex_set (pRt, count1, count1, z);
		for (count2 = count1 + 1 ; count2 < Nt; count2++)
		   {
			temp = count2 - count1;
			temp = temp * temp;
			GSL_SET_COMPLEX (&z, pow (Rt, temp), 0.0);
			gsl_matrix_complex_set (pRt, count1, count2, z);
			gsl_matrix_complex_set (pRt, count2, count1, z);
		   }

	    }
	//matrixviewer ("pRt", pRt);
	gsl_matrix_complex_memcpy (ptemp2, pRt);
	matrix_sqrt(ptemp2, pRt);
	//matrixviewer ("pRt", pRt);
	gsl_matrix_complex_free (ptemp);
	gsl_matrix_complex_free (ptemp2);
	return;
}
