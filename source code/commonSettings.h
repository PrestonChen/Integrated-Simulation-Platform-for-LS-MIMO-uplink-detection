/*
 * commonSettings.h
 * This head file sets the detetion algorithm to be tested
 */

#ifndef COMMONSETTINGS_H_
#define COMMONSETTINGS_H_
//debugging mode
//#define DEBUG
//MMSE with exact matrix inversion (EMI)
//#define MMSEEMI
//MMSE with Neumann series expansion Newton iteration approximate matrix inversion (K-SENIA)
//#define MMSESENIA
//MMSE with 3-term Neumann series expansion (SE3) approximate matrix inversion
//#define MMSESE3
//MMSE with exact matrix inversion and inflate update (IU)
//#define MMSEEMIIU
////MMSE with K-SENIA and IU
//#define MMSESENIAIU
////MMSE with SE3 and IU
//#define MMSESE3IU
//hybrid linear MMSE detection
//#define MMSEHYB
//MMSE-SIC without ordering
//#define MMSESIC
//MMSE-SIC with conventional V-BLAST ordering
//#define MMSEVBLASTSIC
//MMSE-SIC with improved ordering
//#define MMSEIOSIC
//selection based list detection with EMI
#define SLD
//selection based list detection with AMI (SENIAIU)
//#define SLDSENIAIU
//diversity maximization selection
//#define DMS
//channel independent selection
#define CIS
//SIC sub-detection without ordering
//#define SIC
//SIC with V-BLAST ordering
//#define VBLASTSIC
//SIC with improved ordering
#define IOSIC




#endif /* COMMONSETTINGS_H_ */
