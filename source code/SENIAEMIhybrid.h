/*
 * SENIAEMIhybrid.h
 *
 *  Created on: Mar 23, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef SENIAEMIHYBRID_H_
#define SENIAEMIHYBRID_H_
#include "gslCommon.h"
#include "SENIA.h"
int SENIAEMIhybrid(gsl_matrix_complex *Win, int k, double OM, double threshold,
		gsl_matrix_complex *Wout, gsl_matrix_complex *Waided){
   int Nt = Win->size1;
   gsl_matrix_complex *G;
   gsl_permutation *p;
   int signum[1];
   int hybridLabel;
   G = gsl_matrix_complex_calloc(Nt, Nt);
   p = gsl_permutation_calloc(Nt);
   signum[0] = 1;
   gsl_matrix_complex_memcpy(G, Win);
   if (OM > threshold){
	  SENIA(Win, k, Wout, Waided);
	  hybridLabel = 1;
   }else{
	   gsl_linalg_complex_LU_decomp(G, p, signum);
	   gsl_linalg_complex_LU_invert(G, p, Wout);
	   hybridLabel = 0;
   }
   gsl_matrix_complex_free(G);
   gsl_permutation_free(p);
   return hybridLabel;
}


#endif /* SENIAEMIHYBRID_H_ */
