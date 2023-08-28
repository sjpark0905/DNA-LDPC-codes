#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "rand.h"
#include "channel.h"

double getStd_dev(double EbNo, double rate)
{
	double ENL;
	double ENLxRATE;
	ENL = pow(10.0, (EbNo * 0.1));
	ENLxRATE = 2 * rate * ENL;
	return 1/sqrt(ENLxRATE);
}

// codeword		: transmitted codeword
// N			: length of codeword
// recv			: received codeword
// received_LR	: received likelihood ratio
// std_dev		: noise std_dev
void channel_AWGN(char *codeword, int N, double * recv, double *received_LR, double * received_LLR, double std_dev)
{
	int i;
	rand_gaussian(recv,N,std_dev);
	for(i = 0 ;i < N ; i++)
	{
		if(codeword[i] == 0)	{	recv[i] += 1;	}
		else					{	recv[i] += -1;	}

		received_LLR[i] = 2.0*recv[i]/(std_dev*std_dev);
		received_LR[i]  = exp(received_LLR[i]);
	}
}

void channel_BSC(char* codeword, int N, int* recv, double* received_LR, double* received_LLR, double std_dev, double p) //bpsk에서 std_dev를 통해 crossover probability p를 얻을 수 있을 것. 
{
	//static double* noise = NULL;
	double temp_error;
	static int      noise_len = 0;

	/*
	if (noise == NULL)
	{
		noise = (double*)calloc(N, sizeof(double));
		noise_len = N;
	}
	else
	{
		if (noise_len != N)
		{
			free(noise);
			noise = (double*)calloc(N, sizeof(double));
			noise_len = N;
		}
	}
	*/
	//rand_gaussian(noise, N, std_dev);
	//*p = 0.5 * erfc(1 / std_dev / sqrt(2)); //crossover probability
	//printf("%5f\n", p);
	for (int i = 0; i < N; i++)
	{

		rand_uniform(&temp_error, 1);

		//if (codeword[i] == 0)   noise[i] += 1;
		//else               noise[i] += -1;

		// hard decision
		if (temp_error >= p)      recv[i] = 1;
		else               recv[i] = -1;

		//printf("%5f", p);
		if (recv[i] == 1)
		{
			received_LR[i] = (1 - p) / p;
			received_LLR[i] = log(received_LR[i]);
		}
		else
		{
			received_LR[i] = p / (1 - p);
			received_LLR[i] = log(received_LR[i]);
		}

	}


}



// array_erased = TRUE : erasure event 
// array_erased = FALSE : 
void channel_BEC(char *codeword, int N, int * recv, double eps)
{
	//
	double temp;

	

	// 0 <= temp <= 1 
	// uniform dist.

	for(int i = 0 ; i < N ; i++)
	{
		rand_uniform(&temp, 1);

		if(temp  < eps) // 0<= x < eps : erasure event 
		{

			recv[i] = ERASE_MARK;
		}
		else
		{
			recv[i] = codeword[i];
		}

	}
}