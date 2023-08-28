#ifndef __CHANNEL_H__
#define __CHANNEL_H__

#define ERASE_MARK	2

/* CHANNEL.H  */

double	getStd_dev(double EbNo, double rate);
void	channel_AWGN(char *codeword, int N, double * recv, double *received_LR, double * received_LLR, double std_dev);
void   channel_BSC(char* codeword, int N, int* recv, double* received_LR, double* received_LLR, double std_dev,double p);
void	channel_BEC(char *codeword, int N, int * recv, double eps);

#endif