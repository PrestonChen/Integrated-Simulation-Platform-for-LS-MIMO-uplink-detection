/*
 * systemConfigurations.h
 *  Created on: Apr 3, 2016
 *      Author: tianpei.chen@mail.mcgill.ca
 *      This file sets all the system configurations for this integrated simulation platform.
 */
#ifndef SYSTEMCONFIGURATIONS_H_
#define SYSTEMCONFIGURATIONS_H_
//number of receive antenna
const int receiveAntennas = 128;
//number of transmit antenna
const int transmitAntennas = 64;
//symbol alphabet size (M-QAM (M = 4,16,64))
const int symConstellationSize = 4;
//the minimum symbol error accumulated
const int minFrameErrors = 100;
//the minimum channel realizations
const int minChannelRealizations = 1e4;
//the start SNR (dB)
const int Start_SNR = 0;
//the end SNR (dB)
const int End_SNR = 8;
//the step of SNR (dB)
const int Step_SNR = 2;
//the iteration time of K-SENIA
const int k = 3;
//the size of the initial matrix inversion for IU (EMI-IU, K-SENIA-IU and SE3-IU)
const int L = 16;
//the threshold of orthogonality measure (MMSE-HYB)
const double OM = 0.192;
//Estimated percentage of usage of K-SENIA in MMSE-HYB
const double perSENIA = 0.75;
//number of antennas selected at the channel partition stage of selection based list detections
const int Nsel = 1;
//Output file name
#define fileName "SimulationResults.txt"  //the output file
#endif /* SYSTEMCONFIGURATIONS_H_ */
