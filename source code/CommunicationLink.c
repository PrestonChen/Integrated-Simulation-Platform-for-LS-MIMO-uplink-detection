/*
 ============================================================================
 Name        : CommunicationLink.c
 Author      : Tianpei Chen
 Version     : final 3.0
 Copyright   : copyright (c) 2016 Tianpei Chen-All Rights Reserved 
 Description : This is the integrated simulation platform for the various detection
 techniques discussed in Tianpei Chen's Master thesis.

This platform is the general software testbed for the performance test of
LS-MIMO detectors proposed. The simulation results are presented in the form
of frame error rate (FER), symbol error rate (SER) and bit error rate (BER),
symbol error rate (SER) and frame error rate (FER) versus the average
received signal to noise ratio (SNR) in dB.

At each SNR level, Monte-Carlo (MC) simulation is performed in the step
of channel realizations. The MC simulation stops until a minimum number of channel
realization is achieved as well as a minimum number of frame errors are accumulated.

 ============================================================================
 */

#include "gslCommon.h"
#include "commonSettings.h"
#include "systemConfigurations.h"
#include "RectangularQAMSlicer.h"
#include "grayencoder.h"
#include "symbolconstellation.h"
#include "data_generator.h"
#include "modulator.h"
#include "channel_generator.h"
#include "corr_matrix_generator.h"
#include "error_channel_generator.h"
#include "noise_generator.h"
#include "MMSE.h"
#include "MMSE_IU.h"
#include "MMSE_SENIA.h"
#include "MMSE_SENIA_IU.h"
#include "MMSE_SE.h"
#include "MMSE_SE_IU.h"
#include "MMSE_hybrid.h"
#include "MMSE_SIC_ordering.h"
#include "selectiveListDetection.h"
#include "symErrorCheck.h"
#include "demodulator.h"
#include "binaryerrors.h"




int main(void) {
	/**
	 * Basic platform initialization
	 */
	char CWD[1024];   //the current work directory
	getcwd(CWD, sizeof(CWD));  //get the current directory
	const int Nr=receiveAntennas;   // number of receive antennas
	const int Nt=transmitAntennas; // number of transmit antennas
	const int M=symConstellationSize; //symbol constellation size
	const double pav=(double)1/(double)Nt;  //average symbol power
	clock_t start, end; //time clock
	int SNR_tmp;
	FILE *pfile; //the output file
	pfile=fopen(fileName, "a");
    fprintf(pfile, "==============================================================================\n");
    fprintf(pfile, "Output file of the integrated simulation platform of LS-MIMO detectors\n");
    fprintf(pfile, "Current work directory is %s\n", CWD);
    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "SYSTEM CONFIGURATIONS\n");
    fprintf(pfile, "%d times %d MIMO with %d QAM modulation\n", Nr, Nt, M);
    fprintf(pfile, "the average transmit symbol power is %g\n", pav);
    fprintf(pfile, "SNR (dB) are\n");
     SNR_tmp=Start_SNR;
     while(SNR_tmp<=End_SNR){
     	fprintf(pfile, "%d, ", SNR_tmp );
     	SNR_tmp+=Step_SNR;
     }
     fprintf(pfile, "\n");
     fprintf(pfile, "the minimum channel realization is %g\n", (double) minChannelRealizations);
     fprintf(pfile, "the minimum frame errors accumulated is %g\n", (double) minFrameErrors);

/**
 * Command Window Interface output
 */
     printf("===================================================\n");
     printf("The program begin.\n");
     printf("Integrated simulation platform of LS-MIMO detections\n");
     printf("%d X %d %d QAM system\n", Nr, Nt, M);
     printf("Current work directory is %s\n", CWD);

printf("The SNR (dB) considered are\n");
     SNR_tmp=Start_SNR;
     while(SNR_tmp<=End_SNR){
     	printf("%d, ", SNR_tmp);
     	SNR_tmp+=Step_SNR;
     }
printf("\n");

/**
 * print the name of the algorithm tested
 */
#if defined(MMSEEMI)
    printf("The detector tested is MMSE-EMI\n");
    fprintf(pfile, "The detector tested is MMSE-EMI\n");
#elif defined(MMSESE3)
    printf("The detector tested is MMSE-SE3\n");
    fprintf(pfile, "The detector tested is MMSE-SE3\n");
#elif defined(MMSESENIA)
    printf("The detector tested is MMSE-SENIA\n");
    fprintf(pfile, "The detector tested is MMSE-SENIA\n");
#elif defined(MMSEEMIIU)
    printf("The detector tested is MMSE-EMI-IU\n");
    fprintf(pfile, "The detector tested is MMSE-EMI-IU\n");
#elif defined(MMSESE3IU)
    printf("The detector tested is MMSE-SE3-IU\n");
    fprintf(pfile, "The detector tested is MMSE-SE3-IU\n");
#elif defined(MMSESENIAIU)
    printf("The detector tested is MMSE-SENIA-IU\n");
    fprintf(pfile, "The detector tested is MMSE-SENIA-IU\n");
#elif defined(MMSESIC)
    printf("The detector tested is MMSE-SIC\n");
    fprintf(pfile, "The detector tested is MMSE-SIC\n");
#elif defined(MMSEVBLASTSIC)
    printf("The detector tested is MMSE-VBLASTSIC\n");
    fprintf(pfile,"The detector tested is MMSE-VBLASTSIC\n");
#elif defined(MMSEIOSIC)
    printf("The detector tested is MMSE-IOSIC\n");
    fprintf(pfile,"The detector tested is MMSE-IOSIC\n");
#elif defined(MMSEHYB)
    printf("The detector tested is MMSE-HYB\n");
    fprintf(pfile, "The detector tested is MMSE-HYB\n");
#elif defined(DMS)&defined(SLD)&defined(SIC)
    printf("The detector tested is DMS-Ei-SIC\n");
    fprintf(pfile, "The detector tested is DMS-Ei-SIC\n");
#elif defined(DMS)&defined(SLD)&defined(VBLASTSIC)
    printf("The detector tested is DMS-Ei-VBLASTSIC\n");
    fprintf(pfile, "The detector tested is DMS-Ei-VBLASTSIC\n");
#elif defined(DMS)&defined(SLD)&defined(IOSIC)
    printf("The detector tested is DMS-Ei-IOSIC\n");
    fprintf(pfile, "The detector tested is DMS-Ei-IOSIC\n");
#elif defined(DMS)&defined(SLDSENIAIU)&defined(SIC)
    printf("The detector tested is DMS-Ai-SIC\n");
    fprintf(pfile, "The detector tested is DMS-Ai-SIC\n");
#elif defined(DMS)&defined(SLDSENIAIU)&defined(VBLASTSIC)
    printf("The detector tested is DMS-Ai-VBLASTSIC\n");
    fprintf(pfile, "The detector tested is DMS-Ai-VBLASTSIC\n");
#elif defined(DMS)&defined(SLDSENIAIU)&defined(IOSIC)
    printf("The detector tested is DMS-Ai-IOSIC\n");
    fprintf(pfile, "The detector tested is DMS-Ai-IOSIC\n");
#elif defined(CIS)&defined(SLD)&defined(SIC)
    printf("The detector tested is CIS-Ei-SIC\n");
    fprintf(pfile, "The detector tested is CIS-Ei-SIC\n");
#elif defined(CIS)&defined(SLD)&defined(VBLASTSIC)
    printf("The detector tested is CIS-Ei-VBLASTSIC\n");
    fprintf(pfile, "The detector tested is CIS-Ei-VBLASTSIC\n");
#elif defined(CIS)&defined(SLD)&defined(IOSIC)
    printf("The detector tested is CIS-Ei-IOSIC\n");
    fprintf(pfile, "The detector tested is CIS-Ei-IOSIC\n");
#elif defined(CIS)&defined(SLDSENIAIU)&defined(SIC)
    printf("The detector tested is CIS-Ai-SIC\n");
    fprintf(pfile, "The detector tested is CIS-Ai-SIC\n");
#elif defined(CIS)&defined(SLDSENIAIU)&defined(VBLASTSIC)
    printf("The detector tested is CIS-Ai-VBLASTSIC\n");
    fprintf(pfile, "The detector tested is CIS-Ai-VBLASTSIC\n");
#elif defined(CIS)&defined(SLDSENIAIU)&defined(IOSIC)
    printf("The detector tested is CIS-Ai-IO-SIC\n");
    fprintf(pfile, "The detector tested is CIS-Ai-IOSIC\n");
#endif

/**
 * Print the specifications of the tested algorithms
 */

//Applications of approximate matrix inversion in linear MMSE detection
#if defined(MMSEEMI)||defined (MMSESENIA)||defined(MMSESE3)||defined(MMSEEMIIU)||defined(MMSESENIAIU)||defined(MMSESE3IU)
#if defined (MMSEEMI)||defined(MMSEEMIIU)
    printf("MMSE with exact matrix inversion\n");
    fprintf(pfile, "MMSE with exact matrix inversion\n");
#if defined (MMSEEMIIU)
    printf("with IU\n");
    printf("The size of the initial matrix inversion (approximate) is %d\n", L);
    fprintf(pfile, "with IU\n");
    fprintf(pfile, "The size of the initial matrix inversion (approximate) is %d\n", L);
#endif
#else
	printf("approximate matrix inversion aided MMSE mode\n");
	printf("The iteration time of SENIA or SE3 is %d\n", k);
	fprintf(pfile, "approximate matrix inversion aided MMSE mode\n");
	fprintf(pfile,"The iteration time of SENIA or SE3 is %d\n", k);
#if defined(MMSESENIAIU)||defined(MMSESE3IU)
	printf("with IU\n");
	printf("The size of the initial matrix inversion (approximate) is %d\n", L);
	fprintf(pfile, "with IU\n");
	fprintf(pfile, "The size of the initial matrix inversion (approximate) is %d\n", L);
#endif
#endif
//hybrid MMSE detection
#elif defined(MMSEHYB)
    printf("Hybrid linear detector mode\n");
    printf("The estimated percentage of SENIA is %g\n", perSENIA);
    printf("The threshold of orthogonality measure is %g\n", OM);
    printf("The iteration time of SENIA is %d\n", k);
    fprintf(pfile, "Hybrid linear detector mode\n");
    fprintf(pfile, "The estimated percentage of SENIA is %g\n", perSENIA);
    fprintf(pfile, "The threshold of orthogonality measure is %g\n", OM);
    fprintf(pfile, "The iteration time of SENIA is %d\n", k);
//selection based list detection
#elif defined (SLD)||defined(SLDSENIAIU)
    printf("Selection based list detection mode\n");
    fprintf(pfile, "Selection based list detection mode\n");
	const int N=Nsel; //number of antennas selected at channel partition stage
    printf("The number of antennas selected at the channel partition stage is %d\n", N);
    fprintf(pfile, "The number of antennas selected at the channel partition stage is %d\n", N);
#if defined (DMS)
    printf("The selection scheme is diversity maximization selection\n");
    fprintf(pfile, "The selection scheme is diversity maximization selection\n");
#endif
#if defined (CIS)
    printf("The channel selection scheme is channel independent selection\n");
    fprintf(pfile, "The channel selection scheme is channel independent selection\n");
#endif
#if defined (SLD)
    printf("Ei mode\n");
    fprintf(pfile, "Ei mode\n");
#endif
#if defined (SLDSENIAIU)
    printf("Ai mode\n");
    fprintf(pfile, "Ai mode\n");
    printf("The number of iterations of SENIA is %d\n", k);
    fprintf(pfile, "The number of iterations of SENIA is %d\n", k);
    printf("The size of the initial matrix inversion for IU is %d\n", L);
    fprintf(pfile, "The size of the initial matrix inversion for IU is  %d\n", L);
#endif
#if defined (SIC)
    printf("SIC sub-detection without ordering\n");
    fprintf(pfile, "SIC sub-detection without ordering\n");
#endif
#if defined (VBLASTSIC)
    printf("SIC sub-detection with V-BLAST ordering\n");
    fprintf(pfile, "SIC sub-detection with V-BLAST ordering\n");
#endif
#if defined (IOSIC)
    printf("SIC sub-detection with improved ordering\n");
    fprintf(pfile, "SIC sub-detection with improved ordering\n");
#endif


#endif
printf("===================================================\n");
fprintf(pfile, "*************************************\n");
fprintf(pfile, "\n");
fprintf(pfile, "\n");


/**
 * Communication Link data setup
 */
    int count;
    gsl_vector_ulong *pgraydata = gsl_vector_ulong_calloc (M); //gray code book
    gsl_vector_complex *pSymConstell = gsl_vector_complex_calloc (M); //symbol constellation alphabet
    gsl_vector_ulong *pdata = gsl_vector_ulong_calloc (Nt); //the indexes of transmitted data
    gsl_vector_ulong *pgrayInput = gsl_vector_ulong_calloc (Nt);// the gray code of the transmitted data
    gsl_vector_complex *ptransmitted = gsl_vector_complex_calloc (Nt); //the transmitted symbol vector
    gsl_matrix_complex *pH = gsl_matrix_complex_calloc (Nr, Nt);   //the physical channel matrix
    gsl_matrix_complex *pH_tmp=gsl_matrix_complex_calloc(Nr, Nt); //the temporary channel matrix
    gsl_matrix_complex *pHest = gsl_matrix_complex_calloc (Nr, Nt);   //the estimated channel matrix (with error)
    gsl_matrix_complex *pRr = gsl_matrix_complex_calloc (Nr, Nr);  //receive spatial correlation matrix
    gsl_matrix_complex *pRt = gsl_matrix_complex_calloc (Nt, Nt);  //transmit spatial correlation matrix
    gsl_vector_complex *pnoise = gsl_vector_complex_calloc (Nr);  //AWGN vector
    gsl_vector_complex *preceived = gsl_vector_complex_calloc (Nr);  //the received signal vector
    gsl_vector_complex *psymOut=gsl_vector_complex_calloc(Nt); //the output of the detectors
    int *ErrorIndex_V=(int*)calloc(1, sizeof(int)*Nt);    //the indexes of the erroneous symbol
    gsl_vector_ulong *pgrayOut;   //output gray code for erroneous symbols
    gsl_rng *pr;    //random number generator
	int seed = time (NULL);
	const gsl_rng_type *pT;
	pT = gsl_rng_default;
	pr = gsl_rng_alloc (pT);
    gsl_rng_set (pr,seed);
    gsl_complex alpha, beta;   //auxiliary complex numbers
    GSL_SET_COMPLEX(&alpha, 1,0);
    GSL_SET_COMPLEX(&beta, 0, 0);
    int frameError, symError, bitError;  //accumulated frame errors, symbol errors and bit errors in one SNR point
    int *frameError_sub, *symError_sub, *bitError_sub; //frame errors, symbol errors and bit errors in one channel realization
    frameError_sub=(int*)calloc(1, sizeof(int));
    symError_sub=(int*)calloc(1, sizeof(int));
    bitError_sub=(int*)calloc(1, sizeof(int));
    double  noiseV, snr, sigmaerr;
    double FER, SER, BER;   //frame error rate, symbol error rate and bit error rate
    int Realizations=0;   //the number of channel realizations
    SNR_tmp=Start_SNR;
    //generate gray code book
    grayencoder (pgraydata);
#if defined(DEBUG)
    printf("the gray data are\n");
    for (count=0;count<M;count++){
    	printf("%d, ", gsl_vector_ulong_get(pgraydata, count));
    }
    printf("\n");
#endif
    //generate symbol constellation alphabet
    symbolconstellation(pSymConstell, pav);
#if  defined(DEBUG)
    printf("the symbol constellation are:\n");
    for (count=0;count<M;count++){
    	printf("%f+i%f, ", gsl_vector_complex_get(pSymConstell, count).dat[0],
    			gsl_vector_complex_get(pSymConstell, count).dat[1]);
    }
    printf("\n");
#endif
    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "SYSTEM OUTPUT\n");
#if defined(MMSEHYB)
    fprintf(pfile, "SNR \t channel realizations \t frame errors \t symbol errors \t bit errors \t FER \t SER \t BER \t Percentage of AMI \t Operation time\n");
#else
    fprintf(pfile, "SNR \t channel realizations \t frame errors \t symbol errors \t bit errors \t FER \t SER \t BER \t Operation time\n");
#endif
    fclose(pfile);
/**
 * Monte-Carlo Simulation process
 */
    SNR_tmp=Start_SNR;
    while (SNR_tmp<=End_SNR){
    	frameError=0;
    	symError=0;
    	bitError=0;
    	FER=0;
    	SER=0;
    	BER=0;
    	Realizations=0;
    	pfile=fopen(fileName, "a");
		snr=pow(10,((double)SNR_tmp/(double)10)); //SNR in decimal
		noiseV=(double)1/snr;     //noise variance
#if defined(MMSEHYB)
		int hybridLabel=0;  //set the counter that records the number of K-SENIA-IU aided MMSE used in the hybrid MMSE detector
#endif
    	start=clock();
    	while (frameError<minFrameErrors||Realizations<minChannelRealizations){
    		//generate the random index of the data to be transmitted
    		data_generator(pdata, pr, M);
#if defined(DEBUG)
    		printf("Test of pdata\n");

    		for ( count=0;count<Nt; count++){
    			printf("%u, ", gsl_vector_ulong_get(pdata, count));
    		}
    		printf("\n");
#endif
            //generate transmitted symbol vector and the corresponding gray code
    		modulator(pdata,  pSymConstell, pgraydata, ptransmitted, pgrayInput);
    		//generate channel
    		channel_generator (pH, pr);
    		//generate AWGN vector
    		noise_generator (pnoise, pr, noiseV);   //generate noise vector
    		gsl_vector_complex_memcpy(preceived, pnoise);
    		gsl_blas_zgemv(CblasNoTrans, alpha, pH, ptransmitted, alpha, preceived);  //generate  receive signal vector
/**
 * Detectors
 */
//Applications of approximate matrix inversion in linear MMSE detection
#if defined(MMSEEMI)||defined(MMSESE3)||defined(MMSESENIA)||defined(MMSEEMIIU)||defined(MMSESE3IU)||defined(MMSESENIAIU)
#if defined (MMSEEMI)
    		MMSE(preceived, pH, snr/(double)(Nt), pav,  M,  psymOut);
#endif
#if defined (MMSESE3)
    		 MMSE_SE(preceived, pH, snr/(double)Nt,  pav, M, k, psymOut);
#endif
#if defined (MMSESENIA)
    		 MMSE_SENIA(preceived, pH, snr/(double)Nt, pav, M,  k, psymOut);
#endif
#if defined (MMSEEMIIU)
    		 MMSE_IU(preceived, pH, snr/(double)Nt, pav,  M, L, psymOut);
#endif
#if defined (MMSESE3IU)
    		 MMSE_SE_IU (preceived, pH, snr/(double)Nt,  pav, M,  k,  L, psymOut);
#endif
#if defined (MMSESENIAIU)
    		 MMSE_SENIA_IU (preceived, pH, snr/(double)Nt,  pav, M,  k,  L, psymOut);
#endif
//hybrid linear detector
#elif defined(MMSEHYB)
  hybridLabel+=MMSE_hybrid(preceived, pH, snr/(double)Nt, pav, M, k, OM, psymOut);
//MMSE-SIC with different ordering strategies
#elif defined(MMSESIC)||defined(MMSEVBLASTSIC)||defined(MMSEIOSIC)
  MMSE_SIC_ordering(preceived, pH, snr/(double)Nt, pav, M, psymOut);
//selection based list detections
#elif (defined (CIS)||defined (DMS))&(defined (SLD)||defined(SLDSENIAIU))&(defined(SIC)||defined(VBLASTSIC)||defined(IOSIC))
   selectiveListDetection(preceived, pH, snr/(double)Nt, pav, N, k, L, pSymConstell,psymOut);
#endif
    		symErrorCheck(ptransmitted, psymOut, ErrorIndex_V, symError_sub, frameError_sub);   //check symbol and frame errors
    		Realizations++;
    		if (*symError_sub==0){
    			continue;
    		}
    		int *ErrorIndex=(int*)calloc(1, sizeof(int)*(*symError_sub));
    		for (count=0;count<(*symError_sub); count++){
    			ErrorIndex[count]=ErrorIndex_V[count];
    		}
    		pgrayOut=gsl_vector_ulong_calloc(*symError_sub);
#if defined(DEBUG)
    		printf("the indexes of the erroneous symbols are\n");
    		for(count=0;count<(pgrayOut->size);count++){
    			printf("%d, ", ErrorIndex[count]);
    		}
    		printf("\n");

#endif
    		demodulator(psymOut, pSymConstell, pgraydata, ErrorIndex,  pgrayOut);   //decode the estimated symbol vector into gray code
#if defined(DEBUG)
    		printf("the output gray codes are\n");
    		for(count=0;count<(pgrayOut->size);count++){
    			printf("%d, ", gsl_vector_ulong_get(pgrayOut, count));
    		}
    		printf("\n");

#endif
    		binaryerrors (pgrayOut, ErrorIndex, pgrayInput,  M,  bitError_sub); //count the number of bit errors in this channel realization
    		symError+=*symError_sub;
    		frameError+=*frameError_sub;
    		bitError+=*bitError_sub;
    		free(ErrorIndex);
    		gsl_vector_ulong_free(pgrayOut);

    	}
        end=clock();
        FER=(double)frameError/((double)(Realizations));  //calculate frame error rate
        SER=((double)symError/((double)Nt))/((double)(Realizations)); //calculate symbol error rate
        BER=((double)bitError/((double)Nt*(double)ceil(log2(M))))/((double)(Realizations)); //calculate bit error rate
#if defined(MMSEHYB)
    	printf("SNR=%d, Realization=%d, FER=%g, SER=%g, BER=%g, Percentage of AMI=%g, OperationTime=%g s\n", SNR_tmp, Realizations, FER, SER, BER,
    			(double)hybridLabel/(double)Realizations, (end-start)/(double)CLOCKS_PER_SEC);
#else
    	printf("SNR=%d, Realization=%d, FER=%g, SER=%g, BER=%g, OperationTime=%g s\n", SNR_tmp, Realizations, FER, SER, BER,
    			(end-start)/(double)CLOCKS_PER_SEC);
#endif
#if defined(MMSEHYB)
    	fprintf(pfile, "%d \t %d \t %d \t %d \t %d \t %g \t %g \t %g \t %g \t %g\n", \
    			SNR_tmp, Realizations, frameError, symError, bitError, FER, SER, BER,\
    			(double)hybridLabel/(double)Realizations, (end-start)/(double)CLOCKS_PER_SEC);

#else
    	fprintf(pfile, "%d \t %d \t %d \t %d \t %d \t %g \t %g \t %g \t %g\n", \
           SNR_tmp, Realizations, frameError, symError, bitError, FER, SER, BER, \
           (end-start)/(double)CLOCKS_PER_SEC);
#endif
    	fclose(pfile);
    	SNR_tmp+=Step_SNR;

    }
    pfile=fopen(fileName, "a");
    fprintf(pfile, "*************************************\n");
    fprintf(pfile, "The whole program ends successfully!!\n");
    fprintf(pfile, "==============================================================================\n");
    fclose(pfile);
    printf("the whole program ends successfully!!\n");
    gsl_vector_ulong_free(pgraydata);
    gsl_vector_complex_free(pSymConstell);
    gsl_vector_ulong_free(pdata);
    gsl_vector_ulong_free(pgrayInput);
    gsl_vector_complex_free(ptransmitted);
    gsl_matrix_complex_free(pH);
    gsl_matrix_complex_free(pH_tmp);
    gsl_matrix_complex_free(pHest);
    gsl_matrix_complex_free(pRr);
    gsl_matrix_complex_free(pRt);
    gsl_vector_complex_free(pnoise);
    gsl_vector_complex_free(preceived);
    gsl_vector_complex_free(psymOut);
    gsl_rng_free(pr);
    free(frameError_sub);
    free(symError_sub);
    free(bitError_sub);
    free(ErrorIndex_V);
	return 0;
}
