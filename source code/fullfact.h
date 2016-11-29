/*
 * fullfact.h
 *
 *  Created on: Dec 12, 2015
 *      Author: tianpei.chen@mail.mcgill.ca
 */

#ifndef FULLFACT_H_
#define FULLFACT_H_
#include "gslCommon.h"
void fullfact(
		int N,  //the number of elements that use full expansion
		int M,    //constellation size
		gsl_matrix *s_sub_index    //the indexes of the full expansion matrix  (pow(M,N))
		){
	int ssize, ncycles, columns, nreps;
	int count1, count2, count3, count4;
	int *settings2,*settings1;
	settings2 = (int*)calloc(1, (int)(pow(M,N))*sizeof(int));
	ssize = (int)(pow(M,N));
	ncycles = ssize;
	for(count1 = 0; count1 < N; count1++)
	{
//		memset(settings2,0,int(pow(M,N))*sizeof(int));

		nreps = (int)(ssize / ncycles);
		ncycles = ncycles / M;
		settings1 = (int*)calloc(1, M * nreps * sizeof(int));
		for (count2 = 0; count2 < M; count2++)
		{
			for (count3 = 0; count3 < nreps; count3++)
			{
				settings1[count2 * nreps + count3] = count2;
			}
		}


	for(count2 = 0; count2 < ncycles; count2++)
	{
		for(count3 = 0; count3 < M * nreps; count3++)
		{
			settings2[count2 * M * nreps + count3] = settings1[count3];
		}
	}
	for(count4 = 0; count4 < pow(M,N);count4++)
	{
	gsl_matrix_set(s_sub_index, count1, count4, settings2[count4]);
	}

    free(settings1);
	}
//	memset(settings2,0,int(pow(M,N))*sizeof(int));
	free(settings2);
//	free(settings1);
	return;
}




#endif /* FULLFACT_H_ */
