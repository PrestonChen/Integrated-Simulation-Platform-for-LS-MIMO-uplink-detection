/*
 * RectangularQAMSlicer.h
 *
 *  Created on: Dec 6, 2015
 *      Author: Preston Chen
 */

#ifndef RECTANGULARQAMSLICER_H_
#define RECTANGULARQAMSLICER_H_
void RectangularQAMSlicer(gsl_vector_complex *psymOut, double pav, int M);
void RectangularQAMSlicer(gsl_vector_complex *psymOut, double pav, int M){
	int Nt=psymOut->size;
	double d;
	if(M==2){
		d=sqrt(pav);
	}else{
		d=sqrt(3*pav/(2*(M-1)));
	}
	int count;
	gsl_complex temp;
	for (count=0;count<Nt; count++){
		temp=gsl_vector_complex_get(psymOut, count);
	if (M==2){
         if (temp.dat[0]<0){
        	 GSL_SET_COMPLEX(&temp, -d,0);
         }else{
        	 GSL_SET_COMPLEX(&temp, d, 0);
         }
	}else if (M==4){
		if (temp.dat[0]<0){
			GSL_SET_REAL(&temp, -d);
		}else{
			GSL_SET_REAL(&temp, d);
		}

		if (temp.dat[1]<0){
			GSL_SET_IMAG(&temp, -d);
		}else{
			GSL_SET_IMAG(&temp, d);
		}


	}else if (M==16){
		if (temp.dat[0]<-2*d){
			GSL_SET_REAL(&temp, -3*d);
		}else if(temp.dat[0]>=-2*d&&temp.dat[0]<0){
			GSL_SET_REAL(&temp, -d);
		}else if(temp.dat[0]>=0&&temp.dat[0]<2*d){
			GSL_SET_REAL(&temp, d);
		}else if(temp.dat[0]>=2*d){
			GSL_SET_REAL(&temp, 3*d);
		}

		if (temp.dat[1]<-2*d){
			GSL_SET_IMAG(&temp, -3*d);
		}else if(temp.dat[1]>=-2*d&&temp.dat[1]<0){
			GSL_SET_IMAG(&temp, -d);
		}else if(temp.dat[1]>=0&&temp.dat[1]<2*d){
			GSL_SET_IMAG(&temp, d);
		}else if(temp.dat[1]>=2*d){
			GSL_SET_IMAG(&temp, 3*d);
		}

	}else if (M==64){
		if (temp.dat[0]<-6*d){
			GSL_SET_REAL(&temp, -7*d);
		}else if(temp.dat[0]>=-6*d&&temp.dat[0]<-4*d){
			GSL_SET_REAL(&temp, -5*d);
		}else if (temp.dat[0]>-4*d&&temp.dat[0]<=-2*d){
			GSL_SET_REAL(&temp, -3*d);
		}else if (temp.dat[0]>-2*d&&temp.dat[0]<=0){
			GSL_SET_REAL(&temp, -d);
		}else if (temp.dat[0]>0&&temp.dat[0]<=2*d){
			GSL_SET_REAL(&temp, d);
		}else if (temp.dat[0]>2*d&&temp.dat[0]<=4*d){
			GSL_SET_REAL(&temp, 3*d);
		}else if (temp.dat[0]>4*d&&temp.dat[0]<=6*d){
			GSL_SET_REAL(&temp, 5*d);
		}else {
			GSL_SET_REAL(&temp, 7*d);
		}

		if (temp.dat[1]<-6*d){
			GSL_SET_IMAG(&temp, -7*d);
		}else if(temp.dat[1]>=-6*d&&temp.dat[1]<-4*d){
			GSL_SET_IMAG(&temp, -5*d);
		}else if (temp.dat[1]>-4*d&&temp.dat[1]<=-2*d){
			GSL_SET_IMAG(&temp, -3*d);
		}else if (temp.dat[1]>-2*d&&temp.dat[1]<=0){
			GSL_SET_IMAG(&temp, -d);
		}else if (temp.dat[1]>0&&temp.dat[1]<=2*d){
			GSL_SET_IMAG(&temp, d);
		}else if (temp.dat[1]>2*d&&temp.dat[1]<=4*d){
			GSL_SET_IMAG(&temp, 3*d);
		}else if (temp.dat[1]>4*d&&temp.dat[1]<=6*d){
			GSL_SET_IMAG(&temp, 5*d);
		}else {
			GSL_SET_IMAG(&temp, 7*d);
		}


	}else {
		printf("The modulation scheme is not supported!\n");
		return;
	}
	gsl_vector_complex_set(psymOut, count, temp);
	}


   return;
}




#endif /* RECTANGULARQAMSLICER_H_ */
