/* DEC.C - Decoding procedures. */

/* Copyright (c) 2000, 2001 by Radford M. Neal
*
* Permission is granted for anyone to copy, use, modify, or distribute this
* program and accompanying programs and documents for any purpose, provided
* this copyright notice is retained and prominently displayed, along with
* a note saying that the original programs are available from Radford Neal's
* web page, and note is made of any changes made to the programs.  The
* programs and documents are distributed without any warranty, express or
* implied.  As the programs were written for research purposes only, they have
* not been tested to the degree that would be advisable in any important
* application.  All use of these programs is entirely at the user's own risk.
*/
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))  
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y)) 

/* NOTE:  See decoding.html for general documentation on the decoding methods */
#define FORWARD 0
#define BACKWARD 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <cmath>

#include "alloc.h"
#include "mod2sparse.h"
#include "rcode.h"
#include "check.h"
#include "dec.h"
#include "rand.h"
#include "common.h"



int		max_iter;
int		max_iter2;
int		D_v;
int		D_c;
int		bRegular_dv;
int		bRegular_dc;

int cst_row, cst_col;
int b_SC, c_SC, w_SC, L_conv, WIN_SC;
int M_win, N_win;




int activated_chk;
int activated_var;
int activated_dec;

int PL_activated_chk[39];
int PL_activated_var[70];
int PL_activated_dec[70];



char* sdr_check;


int		precision;
int		min_value;
int		max_value;
double	step_size;
int		offset_beta;

int		type_FAID;
int		type_FAID_weight;


int		 bSave_word_state;
int** g_word_state;
int** g_check_state;
int** g_variable_state;
int* g_unsat_check_num;
int* g_err_num;
int		 save_state_num;
int		 tot_save_state_num;

//===============================================================================
// code의 길이 n
// I_max까지 iteration이 진행되었을 때 마지막 iter개의 상태를 저장하려고 한다
// [save_state_num] 번째에는 unsatisfied check수가 가장 적은 상태를 저장한다
//===============================================================================
void Init_Decoder(int n, int m, int iter)
{
	save_state_num = iter + 1; //	0 ~ iter_max
	tot_save_state_num = save_state_num + 1; //가장 unsatisfied check 수가 적은 상태

	g_word_state = (int**)calloc(tot_save_state_num, sizeof(int*));
	for (int i = 0; i < tot_save_state_num; i++)
	{
		g_word_state[i] = (int*)calloc(n, sizeof(int));
	}
	g_check_state = (int**)calloc(tot_save_state_num, sizeof(int*));
	for (int i = 0; i < tot_save_state_num; i++)
	{
		g_check_state[i] = (int*)calloc(m, sizeof(int));
	}

	for (int i = 0; i < tot_save_state_num; i++)
	{
		for (int j = 0; j < n; j++)
		{
			g_word_state[i][j] = 0;
		}

		for (int j = 0; j < m; j++)
		{
			g_check_state[i][j] = 0;
		}
	}

	g_unsat_check_num = (int*)calloc(tot_save_state_num, sizeof(int));
	g_err_num = (int*)calloc(tot_save_state_num, sizeof(int));

	for (int i = 0; i < tot_save_state_num; i++)
	{
		g_unsat_check_num[i] = 0;
		g_err_num[i] = 0;
	}

	g_variable_state = (int**)calloc(2, sizeof(int*));
	for (int i = 0; i < 2; i++)
	{
		g_variable_state[i] = (int*)calloc(n, sizeof(int));
	}
}

//===============================================================================
// code의 regular 여부를 체크한다
//===============================================================================
void CheckRegular(mod2sparse* H)
{
	int i, j;
	int N, M;
	int temp;
	mod2entry* e;

	N = mod2sparse_cols(H);
	M = mod2sparse_rows(H);

	D_v = -1;
	D_c = -1;
	bRegular_dv = true;
	bRegular_dc = true;
	for (j = 0; j < N; j++)
	{
		temp = 0;
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			temp++;
		}

		if (D_v == -1)
		{
			D_v = temp;
		}
		else
		{
			if (temp != D_v)		bRegular_dv = false;
			if (temp > D_v)		D_v = temp;
		}
	}

	for (i = 0; i < M; i++)
	{
		temp = 0;
		for (e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			temp++;
		}

		if (D_c == -1)
		{
			D_c = temp;
		}
		else
		{
			if (temp != D_c)		bRegular_dc = false;
			if (temp > D_c)		D_c = temp;
		}
	}
}


int Run_Belief_Propagation_Decoder_SAVE
(
	mod2sparse* H,			/* Parity check matrix */
	double* lratio,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_L,
	int SC_N
)
{
	int n, c;
	Init_Belief_Propagation(H, lratio, dblk);
	//	test_BER(0,dblk,position_BER,SC_L,SC_N);
	for (n = 0; ; n++)
	{

		//c = check(H, dblk, pchk);
		//if ((n == max_iter) || (c == 0))	{ break;}
		if ((n == max_iter)) { break; }

		Iter_Belief_Propagation(H, lratio, dblk);
		//		test_BER(n+1,dblk,position_BER,SC_L,SC_N);
	}

	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	return n;
}

void test_BER(int iter, char* dblk, double** position_BER, int pos, int V_start, int V_end)
{
	int count = 0;


	for (int j = V_start; j < V_end; j++)
	{
		if (dblk[j] != 0) count++;
	}

	position_BER[iter][pos] = (double)((double)count / (double)(V_end - V_start + 1));

}

//===============================================================================
// standard BEC CHANNEL Decoder
//===============================================================================

int Run_BEC_Decoder
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	int* Target_VN
)
{
	int n, c;
	int N = mod2sparse_cols(H);
	char* dblk_prev;
	dblk_prev = (char*)calloc(N, sizeof(int));

	Init_BEC_Decoder(H, recv, dblk);

	for (n = 0; ; n++)
	{
		for (int i = 0; i < N; i++)
		{
			dblk_prev[i] = dblk[i];
		}
		c = check(H, dblk, pchk);

		Iter_BEC_Decoder(H, dblk);

		int exit = FALSE;
		for (int i = 0; i < N; i++)
		{
			if (dblk_prev[i] != dblk[i])
			{
				break;
			}
			if (i == N - 1)
			{
				exit = TRUE;
			}
		}

		for (int i = Target_VN[0] - 1; i < Target_VN[1]; i++)
		{
			if (dblk[i] != 0)
				break;
			if (i == Target_VN[1] - 1)
				exit = TRUE;
		}

		if ((n == max_iter) || (c == 0) || exit) { break; }
	}



	if (c == 0)
		* bIsCodeword = true;

	free(dblk_prev);
	return n;

}
int Run_BEC_Decoder_TARGET
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,	/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_L,
	int* Mv,
	int* Target_VN
)
{
	int n, c;
	int N = mod2sparse_cols(H);
	char* dblk_prev;
	dblk_prev = (char*)calloc(N, sizeof(int));

	Init_BEC_Decoder(H, recv, dblk);

	for (n = 0; ; n++)
	{
		for (int i = 0; i < N; i++)
		{
			dblk_prev[i] = dblk[i];
		}
		c = check(H, dblk, pchk);

		Iter_BEC_Decoder(H, dblk);

		int exit = FALSE;
		for (int i = 0; i < N; i++)
		{
			if (dblk_prev[i] != dblk[i])
			{
				break;
			}
			if (i == N - 1)
			{
				exit = TRUE;
			}
		}

		for (int i = Target_VN[0] - 1; i < Target_VN[1]; i++)
		{
			if (dblk[i] != 0)
				break;
			if (i == Target_VN[1] - 1)
				exit = TRUE;
		}

		if ((n == max_iter) || (c == 0) || exit) { break; }
	}

	//int temp=0;
	//for(int i=Target_VN[0]-1;i<Target_VN[1];i++)
	//{
	//	if(dblk[i]!=0)
	//	{
	//		printf("%d ",i+1);
	//		temp = 1;
	//	}
	//}
	//if(temp==1)
	//	printf("\n");


	if (c == 0)
		* bIsCodeword = true;
	free(dblk_prev);
	return n;

}


int Run_BEC_Decoder_SAVE
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,	/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_L,
	int* Mv
)
{
	int n, c;
	int N = mod2sparse_cols(H);
	char* dblk_prev;
	dblk_prev = (char*)calloc(N, sizeof(int));

	int D_Start = 0;
	int D_End = Mv[0];

	Init_BEC_Decoder(H, recv, dblk);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	for (int t = 1; t < SC_L; t++)
	{
		D_Start += Mv[t - 1];
		D_End += Mv[t];

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
	}
	for (n = 0; ; n++)
	{

		for (int i = 0; i < N; i++)
		{
			dblk_prev[i] = dblk[i];
		}
		c = check(H, dblk, pchk);

		Iter_BEC_Decoder(H, dblk);

		int exit = FALSE;
		for (int i = 0; i < N; i++)
		{
			if (dblk_prev[i] != dblk[i])
			{
				break;
			}
			if (i == N - 1)
			{
				exit = TRUE;
			}
		}
		if (n < 199)
		{
			D_Start = 0;
			D_End = Mv[0];
			test_BER(n + 1, dblk, position_BER, 0, D_Start, D_End);
			for (int t = 1; t < SC_L; t++)
			{
				D_Start += Mv[t - 1];
				D_End += Mv[t];

				test_BER(n + 1, dblk, position_BER, t, D_Start, D_End);
			}
		}


		if ((n == max_iter) || (c == 0) || exit) { break; }
	}

	//D_Start = 0;
	//D_End = Mv[0];
	//test_BER(1,dblk,position_BER,0,D_Start,D_End);
	//for(int t = 1; t<SC_L; t++)
	//{
	//	D_Start += Mv[t-1];
	//	D_End += Mv[t];

	//	test_BER(1,dblk,position_BER,t,D_Start,D_End);
	//}




	if (c == 0)
		* bIsCodeword = true;
	free(dblk_prev);
	return n;

}
void Init_BEC_Decoder(mod2sparse* H, int* recv, char* dblk)
{
	mod2entry* e;
	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->msg_v_to_c = recv[j];
		}

		dblk[j] = recv[j];
	}// end for(j)

}
void Iter_BEC_Decoder(mod2sparse* H, char* dblk)
{
	Check_Update_BEC(H);
	Variable_Update_BEC(H, dblk);
	//	Decision_BEC(H, dblk);
}

void Check_Update_BEC(mod2sparse* H)
{
	int temp;
	int M = mod2sparse_rows(H);
	for (int i = 0; i < M; i++)
	{
		for (mod2entry* e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			temp = 0;
			// e를 제외한 edge로부터 전달된 incoming msg를 보고 outgoing msg를 결정
			for (mod2entry* edge = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_row(edge))
			{
				if (e->col != edge->col)
				{
					if (edge->msg_v_to_c == ERASE_MARK)
					{
						temp = ERASE_MARK;
						break;
					}
					else
					{
						temp = (bool)temp + edge->msg_v_to_c;
					}
				}
			}//end for(edge)

			e->msg_c_to_v = temp;
		}//end for(e)
	}//end for(i)
}

void Variable_Update_BEC(mod2sparse* H, char* dblk)
{

	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		if (dblk[j] == ERASE_MARK)
		{
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				if (e->msg_c_to_v != ERASE_MARK)
				{
					dblk[j] = e->msg_c_to_v;
					for (mod2entry* edge = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
					{
						edge->msg_v_to_c = e->msg_c_to_v;
					}
				}
				// e를 제외한 edge로부터 전달된 incoming msg를 보고 outgoing msg를 결정
				//for(mod2entry * edge = mod2sparse_first_in_col(H,j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
				//{
				//	if(e->row != edge->row && edge->msg_c_to_v !=ERASE_MARK)
				//	{
				//		e->msg_v_to_c = edge->msg_c_to_v;
				//		break;
				//	}
				//}//end for(edge)

			}//end for(e)
		}
	}//end for (j)
}

void Decision_BEC(mod2sparse* H, char* dblk)
{

	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		if (dblk[j] == ERASE_MARK)
		{
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				if (e->msg_c_to_v != ERASE_MARK)
				{
					dblk[j] = e->msg_c_to_v;
					break;
				}
			}//end for(e)
		}

	}//end for(j)
}






//===============================================================================
// standard Belief Propagation Decoder
//===============================================================================
int Run_Belief_Propagation_Decoder
(
	mod2sparse* H,			/* Parity check matrix */
	double* lratio,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword		/* decode 결과 codeword라고 판단되는지 */
)
{
	int n, c;
	Init_Belief_Propagation(H, lratio, dblk);
	for (n = 0; ; n++)
	{
		c = check(H, dblk, pchk);
		if ((n == max_iter) || (c == 0)) { break; }
		Iter_Belief_Propagation(H, lratio, dblk);
	}

	if (c == 0)
		* bIsCodeword = TRUE;

	return n;
}

/* INITIALIZE PROBABILITY PROPAGATION.  Stores initial ratios, probabilities, and guess at decoding. */
void Init_Belief_Propagation(mod2sparse* H, double* lratio, char* dblk)
{
	mod2entry* e;
	int N;
	int j;

	N = mod2sparse_cols(H);
	for (j = 0; j < N; j++)
	{
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->pr = lratio[j];				//	p0/p1
			e->lr = 1;

			//p0 = 1 - 1/(1+lratio[j])
			//p1 = 1/(1+lration[j])
		}//end for(e)

		dblk[j] = (lratio[j] < 1);	// p0/p1 > 1  : 0
		// p0/p1 <= 1 : 1
	}//end for(j)
}

/* DO ONE ITERATION OF PROBABILITY PROPAGATION. */
void Iter_Belief_Propagation(mod2sparse* H, double* lratio, char* dblk)
{
	double pr, dl, t;
	mod2entry* e;
	int N, M;
	int i, j;



	M = mod2sparse_rows(H);
	N = mod2sparse_cols(H);

	/* Recompute likelihood ratios. */
	/* Compute check node message.  */
	for (i = 0; i < M; i++)
	{
		dl = 1;
		for (e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			e->lr = dl;
			dl *= 1 - 2 / (1 + e->pr);	// (수정)
		}//end for(e)

		dl = 1;
		for (e = mod2sparse_last_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_prev_in_row(e))
		{
			t = e->lr * dl;
			e->lr = (1 + t) / (1 - t);		// (수정)			
			dl *= 1 - 2 / (1 + e->pr);	// (수정)
		}//end for(e)
	}//end for(i)

	/* Recompute probability ratios.  Also find the next guess based on the
	individually most likely values. */
	/* Compute variable node message */
	for (j = 0; j < N; j++)
	{
		pr = lratio[j];
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->pr = pr;
			pr *= e->lr;
		}//end for(e)

		if (_isnan(pr))
			pr = 1;

		//p0 = 1 - 1/(1+pr)

		dblk[j] = (pr <= 1);	// (수정)
		pr = 1;

		for (e = mod2sparse_last_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_prev_in_col(e))
		{
			e->pr *= pr;
			if (_isnan(e->pr))
			{
				e->pr = 1;
			}
			pr *= e->lr;
		}//end for(e)
	}//end for(j)
}

//===============================================================================
// Message Alphabet = {-1,1}
//===============================================================================
int Run_Gallager_Decoder
(
	mod2sparse* H,		// parity check matrix
	int* recv,	// hard decision received value
	char* decoded,// decoded codeword
	char* pchk,	// syndrome
	int				type,	// 0 : A algorithm, 1 : B algorithm
	int* bIsCodeword
)
{
	int n, c;
	int N = mod2sparse_cols(H);
	Init_Gallager(H, recv, decoded);

	for (n = 0; ; n++)
	{
		c = check(H, decoded, pchk);
		if ((n == max_iter) || (c == 0)) { break; }
		Iter_Gallager(H, recv, decoded, type);
	}

	if (c == 0)
		* bIsCodeword = true;
	return n;
}

void Init_Gallager(mod2sparse* H, int* recv, char* decoded)
{
	mod2entry* e;
	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->msg_v_to_c = recv[j];		// m0를 그대로 전달
		}

		if (recv[j] >= 0)	decoded[j] = 0;
		else				decoded[j] = 1;
	}// end for(j)
}

void Iter_Gallager(mod2sparse* H, int* recv, char* decoded, int type)
{
	Check_Update_Gallager(H, recv);
	Variable_Update_Gallager(H, recv, type);
	Decision_Gallager(H, recv, decoded, type);
}

void Check_Update_Gallager(mod2sparse* H, int* recv)
{
	int temp;
	int M = mod2sparse_rows(H);
	for (int i = 0; i < M; i++)
	{
		for (mod2entry* e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			temp = 1;
			// e를 제외한 edge로부터 전달된 incoming msg를 보고 outgoing msg를 결정
			for (mod2entry* edge = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_row(edge))
			{
				if (e->col != edge->col)
				{
					temp = temp * edge->msg_v_to_c;
				}
			}//end for(edge)

			e->msg_c_to_v = temp;
		}//end for(e)
	}//end for(i)
}

void Variable_Update_Gallager(mod2sparse* H, int* recv, int type)
{
	int num;
	int message;
	int b;

	if (type == 0)		b = D_v - 1;
	else if (type == 1)	b = D_v - 2;
	else				b = (D_v / 2) + (D_v % 2);

	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		message = -recv[j];	// -m0
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			num = 0;
			// e를 제외한 edge로부터 전달된 incoming msg를 보고 outgoing msg를 결정
			for (mod2entry* edge = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
			{
				if (e->row != edge->row)
				{
					if (edge->msg_c_to_v == message)
						num++;
				}
			}//end for(edge)

			if (num >= b)	e->msg_v_to_c = message;
			else			e->msg_v_to_c = recv[j];
		}//end for(e)
	}//end for (j)
}

void Decision_Gallager(mod2sparse* H, int* recv, char* decoded, int type)
{
	int num;
	int message;
	int b;
	int temp;

	if (type == 0)		b = D_v;
	else if (type == 1)	b = D_v - 1;
	else				b = (D_v / 2) + 1;

	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		num = 0;
		message = -recv[j];	// -m0
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			if (e->msg_c_to_v == message)
				num++;
		}//end for(e)

		if (num >= b)	temp = message;
		else			temp = recv[j];

		if (temp >= 0)	decoded[j] = 0;
		else			decoded[j] = 1;
	}//end for(j)
}

//===============================================================================
//  Shiva Kumar Planjery, Bane Vasic
//===============================================================================
int Run_Finite_Alphabet_Iterative_Decoder
(
	mod2sparse* H,			// parity check matrix
	int* recv,		// hard decision received value
	char* decoded,	// decoded codeword
	char* pchk,		// syndrome
	int* bIsCodeword	// check codeword
)
{
	int n, c;
	int N = mod2sparse_cols(H);
	Init_FAID(H, recv, decoded);

	for (n = 0; ; n++)
	{
		c = check(H, decoded, pchk);

		if ((n == max_iter2) || (c == 0)) { break; }
		Iter_FAID(H, recv, decoded);
	}

	if (c == 0)
		* bIsCodeword = true;
	return n;
}


//===============================================================================
//
//===============================================================================
void Init_FAID(mod2sparse* H, int* recv, char* decoded)
{
	int N = mod2sparse_cols(H);

	for (int j = 0; j < N; j++)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			if (recv[j] >= 0)	e->msg_v_to_c = 1;
			else				e->msg_v_to_c = -1;
		}//end for(e)

		if (recv[j] >= 0)	decoded[j] = 0;
		else				decoded[j] = 1;
	}//end for(j)
}

//===============================================================================
//
//===============================================================================
void Iter_FAID(mod2sparse* H, int* recv, char* decoded)
{
	Check_Update_FAID(H);
	Variable_Update_FAID(H, recv);
	Decision_FAID(H, recv, decoded);
}

//===============================================================================
//
//===============================================================================
void Check_Update_FAID(mod2sparse* H)
{
	int sign;
	int mag_min;
	int M = mod2sparse_rows(H);

	for (int i = 0; i < M; i++)
	{
		for (mod2entry* e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			sign = 1;
			mag_min = 999;

			for (mod2entry* edge = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_row(edge))
			{
				if (e->col != edge->col)
				{
					if (edge->msg_v_to_c >= 0)	sign *= 1;
					else						sign *= -1;

					if (mag_min > abs(edge->msg_v_to_c))
					{
						mag_min = abs(edge->msg_v_to_c);
					}
				}
			} //end for(edge)

			e->msg_c_to_v = sign * mag_min;
		} // end for(e)
	} //end for(i)
}

//===============================================================================
//
//===============================================================================
void Variable_Update_FAID(mod2sparse* H, int* recv)
{
	int m1, m2;
	int N = mod2sparse_cols(H);
	for (int j = 0; j < N; j++)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			m1 = 999;
			m2 = 999;
			for (mod2entry* edge = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
			{
				if (e->row != edge->row)
				{
					if (m1 == 999)			m1 = edge->msg_c_to_v;
					else					m2 = edge->msg_c_to_v;
				}
			} //end for(edge)

			e->msg_v_to_c = Variable_FAID_LUT(m1, m2, recv[j]);
		}// end for (e)
	}//end for (j)
}

//===============================================================================
//
//===============================================================================
void Decision_FAID(mod2sparse* H, int* recv, char* decoded)
{
	double sum;
	int		N = mod2sparse_cols(H);
	double	C;
	int		mag;
	double	sign;
	double	weight[4];

	if (type_FAID_weight == 0)
	{
		C = 0.5;
	}
	else
	{
		C = 1.5;
	}

	weight[0] = 1;
	weight[1] = 1;
	weight[2] = 1;
	weight[3] = 1;

	for (int j = 0; j < N; j++)
	{
		sum = ((double)recv[j]) * C;
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			mag = abs(e->msg_c_to_v);
			if (e->msg_c_to_v >= 0)	sign = 1;
			else					sign = -1;

			sum += (sign * weight[mag]);
		} // end for(e)

		if (sum > 0)					decoded[j] = 0;
		else if (sum < 0)			decoded[j] = 1;
		else						decoded[j] = recv[j];
	}//end for (j)
}

//===============================================================================
//
//===============================================================================
void Set_FAID(int type)
{
	type_FAID = type;
}

void Set_FAID_Weight(int type)
{
	type_FAID_weight = type;
}


int FAID_LUT_ZIZE[N_DIVERSITY] =
{
	//2,
	2,
	2,
	//2,
	//2,
	3,
	3,
	3
};

int FAID_LUT_2[N_LUT_2_SIZE][5][5] =
{
	//Finite alphabet iterative decoders for LDPC codes surpassing floating-point iterative decoders
	//TABLE 1
	{
		//        -2  -1   0   1   2
		/* -2 */ {-2, -2, -2, -2,  0},
		/* -1 */ {-2, -2, -2, -1,  0},
		/*  0 */ {-2, -2, -1,  0,  1},
		/*  1 */ {-2, -1,  0,  0,  1},
		/*  2 */ { 0,  0,  1,  1,  2}
},

////Iterative Decoding Beyond Belief Propagation
//{
//	//        -2  -1   0   1   2
//	/* -2 */ {-2, -2, -2, -2, -1},
//	/* -1 */ {-2, -2, -2, -1,  1},
//	/*  0 */ {-2, -2, -2, -1,  1},
//	/*  1 */ {-2, -1, -1,  1,  2},
//	/*  2 */ {-1,  1,  1,  2,  2}
//},


//Finite Alphabet Iterative Decoding of the (155,64,20) Tanner Code
//TABLE V
{
	//        -2  -1   0   1   2
	/* -2 */ {-2, -2, -2, -2,  0},
	/* -1 */ {-2, -2, -1, -1,  1},
	/*  0 */ {-2, -1, -1,  0,  1},
	/*  1 */ {-2, -1,  0,  1,  2},
	/*  2 */ { 0,  1,  1,  2,  2}
},

////Finite Alphabet Iterative Decoding of the (155,64,20) Tanner Code
////TABLE VI
//{
//	//        -2  -1   0   1   2
//	/* -2 */ {-2, -2, -2, -2,  0},
//	/* -1 */ {-2, -2, -1, -1,  1},
//	/*  0 */ {-2, -1, -1,  0,  2},
//	/*  1 */ {-2, -1,  0,  1,  2},
//	/*  2 */ { 0,  1,  2,  2,  2}
//},

////Finite Alphabet Iterative Decoding of the (155,64,20) Tanner Code
////TABLE VII
//{
//	//        -2  -1   0   1   2
//	/* -2 */ {-2, -2, -2, -2,  0},
//	/* -1 */ {-2, -1, -1, -1,  2},
//	/*  0 */ {-2, -1, -1,  0,  2},
//	/*  1 */ {-2, -1,  0,  2,  2},
//	/*  2 */ { 0,  2,  2,  2,  2}
//},


};

int FAID_LUT_3[N_LUT_3_SIZE][7][7] =
{
	//Finite alphabet iterative decoders for LDPC codes surpassing floating-point iterative decoders
	//TABLE 2
	{
		//        -3  -2  -1   0   1   2   3
		/* -3 */ {-3, -3, -3, -3, -3, -3, -1},
		/* -2 */ {-3, -3, -3, -3, -2, -1,  1},
		/* -1 */ {-3, -3, -2, -2, -1, -1,  1},
		/*  0 */ {-3, -3, -2, -1,  0,  0,  1},
		/*  1 */ {-3, -2, -1,  0,  0,  1,  2},
		/*  2 */ {-3, -1, -1,  0,  1,  1,  3},
		/*  3 */ {-1,  1,  1,  1,  2,  3,  3}
},

//Finite Alphabet Iterative Decoding of the (155,64,20) Tanner Code
//TABLE VIII
{
	//        -3  -2  -1   0   1   2   3
	/* -3 */ {-3, -3, -3, -3, -3, -3, -1},
	/* -2 */ {-3, -3, -3, -3, -2, -1,  1},
	/* -1 */ {-3, -3, -2, -2, -1,  0,  1},
	/*  0 */ {-3, -3, -2, -1, -1,  1,  2},
	/*  1 */ {-3, -2, -1, -1,  0,  1,  2},
	/*  2 */ {-3, -1,  0,  1,  1,  1,  2},
	/*  3 */ {-1,  1,  1,  2,  2,  2,  3}
},

//
{
	//        -3  -2  -1   0   1   2   3
	/* -3 */ {-3, -3, -3, -3, -3, -3, -1},
	/* -2 */ {-3, -3, -2, -2, -1, -1,  1},
	/* -1 */ {-3, -2, -2, -1, -1,  1,  1},
	/*  0 */ {-3, -2, -1, -1, -1,  1,  2},
	/*  1 */ {-3, -1, -1, -1,  0,  1,  2},
	/*  2 */ {-3, -1,  1,  1,  1,  2,  2},
	/*  3 */ {-1,  1,  1,  2,  2,  2,  3}
}
};


//===============================================================================
//
//===============================================================================
int	Variable_FAID_LUT(int m1, int m2, int y)
{
	int temp;
	int size;
	int v = 0;
	int offset = 0;
	int t;

	if (y >= 0)
	{
		m1 = -m1;
		m2 = -m2;
	}

	size = FAID_LUT_ZIZE[type_FAID];
	if (size == 2)
	{
		t = type_FAID;
		offset = 2;
		v = FAID_LUT_2[t][m1 + offset][m2 + offset];
	}
	else if (size == 3)
	{
		t = type_FAID - N_LUT_2_SIZE;
		offset = 3;
		v = FAID_LUT_3[t][m1 + offset][m2 + offset];
	}

	if (y >= 0)
	{
		v = -v;
	}

	return v;
}




//===========================================================================
// 유한한 message size를 가지는 MSA
// g_precision > 0
//===========================================================================
int Run_MSA_Decoder(mod2sparse* H, int* qLLR, int* L_Q, char* decoded, char* pchk, int* bIsCodeword, FILE* fp1, FILE** fp_w_v, int* w_v_idx)
{
	int N = mod2sparse_cols(H);
	int n = 0;
	int c = 0;
	int erasure = 0;
	int min_c = 0;

	erasure = Init_MSA(H, qLLR, decoded);
	for (n = 0; ; n++)
	{
		c = check(H, decoded, pchk);

		if (bSave_word_state == TRUE)
		{
			Save_State(n, c, decoded, pchk, min_c);
			Print_Variable_State(H, n, c, decoded, qLLR, L_Q, pchk, fp_w_v, w_v_idx, fp1);
		}

		if (n == max_iter)
		{
			break;
		}

		if (c == 0)
		{
			break;
		}

		erasure = Iter_MSA(H, qLLR, L_Q, decoded, n, pchk);
	} // end for(n)

	if (c == 0)
		* bIsCodeword = TRUE;

	return n;
}

//===========================================================================
// 무한한 message size를 가지는 MSA
// g_precision = 0
//===========================================================================
int	Run_MSA_Decoder_INF(mod2sparse* H, double* LLR, double* L, char* decoded, char* pchk, int* bIsCodeword)
{
	int N = mod2sparse_cols(H);
	int n = 0;
	int c = 0;
	int min_c = 0;

	Init_MSA_INF(H, LLR, decoded);
	for (n = 0; ; n++)
	{
		c = check(H, decoded, pchk);

		if (bSave_word_state == TRUE)
		{
			Save_State(n, c, decoded, pchk, min_c);
		}

		if (n == max_iter)
		{
			break;
		}

		if (c == 0)
		{
			break;
		}

		Iter_MSA_INF(H, LLR, L, decoded);
	} // end for(n)

	if (c == 0)
		* bIsCodeword = TRUE;

	return n;
}


//===========================================================================
// 유한한 message size를 가지는 MSA : 0 번째 iteration
//===========================================================================
int Init_MSA(mod2sparse* H, int* qLLR, char* decoded)
{
	int erasure = 0;
	int N = mod2sparse_cols(H);
	int r;
	for (int j = 0; j < N; j++)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->msg_v_to_c = qLLR[j];
		}

		if (qLLR[j] > 0)			decoded[j] = 0;
		else if (qLLR[j] < 0)	decoded[j] = 1;
		else
		{
			erasure++;
			r = rand_int(2);
			if (r == 0)		decoded[j] = 0;
			else			decoded[j] = 1;
		}
	}// end for(j)


	for (int i = 0; i < tot_save_state_num; i++)
	{
		for (int j = 0; j < N; j++)
		{
			g_word_state[i][j] = 0;
		}
	}

	for (int i = 0; i < tot_save_state_num; i++)
	{
		g_unsat_check_num[i] = 0;
		g_err_num[i] = 0;
	}

	return erasure;
}

//===========================================================================
// 무한한 message size를 가지는 MSA : 0 번째 iteration
//===========================================================================
void Init_MSA_INF(mod2sparse* H, double* LLR, char* decoded)
{
	int N = mod2sparse_cols(H);
	int r;
	for (int j = 0; j < N; j++)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->msg_v_to_c_fp = LLR[j];
		}

		if (LLR[j] > 0)			decoded[j] = 0;
		else					decoded[j] = 1;
	}// end for(j)

	for (int i = 0; i < tot_save_state_num; i++)
	{
		for (int j = 0; j < N; j++)
		{
			g_word_state[i][j] = 0;
		}
	}

	for (int i = 0; i < tot_save_state_num; i++)
	{
		g_unsat_check_num[i] = 0;
		g_err_num[i] = 0;
	}

}


//===========================================================================
// 유한한 message size를 가지는 MSA : l 번째 iteration (l > 0)
//===========================================================================
int Iter_MSA(mod2sparse* H, int* qLLR, int* L_Q, char* decoded, int iter, char* pchk)
{
	int v = 0;
	Check_Update_MSA(H);
	Variable_Update_MSA(H, qLLR, iter, pchk);
	v = Decision_MSA(H, qLLR, L_Q, decoded);
	return v;
}

//===========================================================================
// 무한한 message size를 가지는 MSA : l 번째 iteration (l > 0)
//===========================================================================
void Iter_MSA_INF(mod2sparse* H, double* LLR, double* L, char* decoded)
{
	Check_Update_MSA_INF(H);
	Variable_Update_MSA_INF(H, LLR);
	Decision_MSA_INF(H, LLR, L, decoded);
}

//===========================================================================
// 유한한 message size를 가지는 MSA : check node update
//===========================================================================
void Check_Update_MSA(mod2sparse* H)
{
	int M = mod2sparse_rows(H);
	int x1, x2;
	int mag_min;
	int sign;
	int beta = offset_beta;

	for (int i = 0; i < M; i++)
	{
		for (mod2entry* e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			mag_min = -1;
			sign = 1;
			for (mod2entry* edge = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_row(edge))
			{
				if (e->col != edge->col)
				{
					if ((mag_min == -1) || (mag_min > abs(edge->msg_v_to_c)))
					{
						//update
						mag_min = abs(edge->msg_v_to_c);
					}

					if (edge->msg_v_to_c >= 0)	sign *= 1;
					else						sign *= -1;
				}
			} //end for(edge)

			//offset
			mag_min = mag_min - beta;
			if (mag_min < 0)		mag_min = 0;

			e->msg_c_to_v = sign * mag_min;
		} //end for(e)
	} //end for(i)
}

//===========================================================================
// 무한한 message size를 가지는 MSA : check node update
//===========================================================================
void Check_Update_MSA_INF(mod2sparse* H)
{
	int M = mod2sparse_rows(H);
	double	x1, x2;
	double	mag_min;
	int		sign;

	for (int i = 0; i < M; i++)
	{
		for (mod2entry* e = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
		{
			mag_min = -1;
			sign = 1;
			for (mod2entry* edge = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_row(edge))
			{
				if (e->col != edge->col)
				{
					if ((mag_min == -1) || (mag_min > abs(edge->msg_v_to_c_fp)))
					{
						//update
						mag_min = abs(edge->msg_v_to_c_fp);
					}

					if (edge->msg_v_to_c_fp >= 0)	sign *= 1;
					else							sign *= -1;
				}
			} //end for(edge)

			//offset
			mag_min = mag_min;
			if (mag_min < 0)		mag_min = 0;

			e->msg_c_to_v_fp = sign * mag_min;
		} //end for(e)
	} //end for(i)	
}

//===========================================================================
// 유한한 message size를 가지는 MSA : variable node update
//===========================================================================
void Variable_Update_MSA(mod2sparse* H, int* qLLR, int iter, char* pchk)
{
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);
	int sum;
	int b;
	int bPostProc = FALSE;
	int sat_count;
	int selected;
	int row_idx;
	int col_idx;
	int selected_vn_num = 0;
	int condition;


	if (bSave_word_state == TRUE)
	{
		if (iter == 30 || iter == 60 || iter == 90)
			bPostProc = TRUE;
	}


	if (bPostProc == FALSE)
	{
		//=======================================================================
		// 일반적인 VN update 실행
		//=======================================================================
		for (int j = 0; j < N; j++)
		{
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				sum = qLLR[j];	//prior L_pr
				for (mod2entry* edge = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
				{
					if (e->row != edge->row)
					{
						// extrinsic messages L_ext
						sum += edge->msg_c_to_v;
					}
				} //end for(edge)

				e->msg_v_to_c = Cal_MSA_Clip(sum, b);
				e->msg_v_to_c_sat = b;
			}// end for (e)
		}//end for (j)	



	}
	else
	{
		//=======================================================================
		// 제안하는 VN update 실행
		//=======================================================================


		//=======================================================================
		// 일단 기존 방식으로 모든 VN를 update한다
		//=======================================================================
		for (int j = 0; j < N; j++)
		{
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				sum = qLLR[j];	//prior L_pr
				for (mod2entry* edge = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
				{
					if (e->row != edge->row)
					{
						sum += edge->msg_c_to_v;
					}
				} //end for(edge)

				e->msg_v_to_c = Cal_MSA_Clip(sum, b);
				e->msg_v_to_c_sat = b;
			}// end for (e)
		}//end for (j)

		//=======================================================================
		// Trapping Set에 포함된다고 생각되는 해당하는 VN를 찾는다 
		//=======================================================================
		for (int i = 0; i < M; i++)
		{
			if (pchk[i] != 0)
			{
				selected_vn_num = 0;	// 하나의 CN에 몇 개의 VN의 값을 바꾸는지 
				for (mod2entry* e1 = mod2sparse_first_in_row(H, i); !mod2sparse_at_end(e1); e1 = mod2sparse_next_in_row(e1))
				{
					col_idx = e1->col;

					condition = FALSE;
					sat_count = 0;

					// message 재조정		
					for (mod2entry* e = mod2sparse_first_in_col(H, col_idx); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
					{
						sum = qLLR[col_idx];	// 수신값
						for (mod2entry* edge = mod2sparse_first_in_col(H, col_idx); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
						{
							if (e->row != edge->row)
							{
								sum += edge->msg_c_to_v;
							}
						} //end for(edge)

						e->msg_v_to_c = Cal_MSA_Clip(sum, b);
						e->msg_v_to_c_sat = b;

						if (b == TRUE)
							sat_count++;
					}// end for (e)

					if (sat_count == 1)
					{
						condition = TRUE;
					}


					//=======================================================================
					// 값을 변경해 주어야 할 VN에 대한 처리를 한다
					//=======================================================================
					if (condition == TRUE)
					{
						selected_vn_num++;
						for (mod2entry* e = mod2sparse_first_in_col(H, col_idx); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
						{
							if (e->msg_v_to_c_sat == TRUE)		e->msg_c_to_v = Set_Message_Maxvalue(e->msg_c_to_v);
							else								e->msg_c_to_v = Dec_Message_value(e->msg_c_to_v);
						}// end for (e)	

						//qLLR[col_idx] = -qLLR[col_idx];
						qLLR[col_idx] = 0;

						for (mod2entry* e = mod2sparse_first_in_col(H, col_idx); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
						{
							sum = qLLR[col_idx];	//prior L_pr
							for (mod2entry* edge = mod2sparse_first_in_col(H, col_idx); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
							{
								if (e->row != edge->row)
								{
									sum += edge->msg_c_to_v;
								}
							} //end for(edge)

							e->msg_v_to_c = Cal_MSA_Clip(sum, b);
							e->msg_v_to_c_sat = b;
						}// end for (e)
					}
				}
				//printf("c - %d : selected_vn_num :%d\n", i, selected_vn_num);
			}
		}
		//printf("==========================================================================\n");
	}

}

//===========================================================================
// 무한한 message size를 가지는 MSA : variable node update
//===========================================================================
void Variable_Update_MSA_INF(mod2sparse* H, double* LLR)
{
	int N = mod2sparse_cols(H);
	double sum;	// belief

	for (int j = 0; j < N; j++)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			sum = LLR[j];	//prior L_pr
			for (mod2entry* edge = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
			{
				if (e->row != edge->row)
				{
					// extrinsic messages L_ext
					sum += edge->msg_c_to_v_fp;
				}
			} //end for(edge)

			e->msg_v_to_c_fp = sum;
		}// end for (e)
	}//end for (j)
}

//===========================================================================
// 유한한 message size를 가지는 MSA : decision
//===========================================================================
int Decision_MSA(mod2sparse* H, int* qLLR, int* L_Q, char* decoded)
{
	int N = mod2sparse_cols(H);
	int sum;
	int r;
	int erasure = 0;

	for (int j = 0; j < N; j++)
	{
		sum = qLLR[j];//prior L_pr
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			// extrinsic messages L_ext
			sum += e->msg_c_to_v;
		} // end for(e)

		L_Q[j] = sum;	// posterior probability L_ps

		if (sum > 0)		decoded[j] = 0;
		else if (sum < 0)	decoded[j] = 1;
		else
		{
			erasure++;
			r = rand_int(2);
			if (r == 0)		decoded[j] = 0;
			else			decoded[j] = 1;
		}
	}//end for (j)

	return erasure;
}

//===========================================================================
// 무한한 message size를 가지는 MSA : decision
//===========================================================================
void Decision_MSA_INF(mod2sparse* H, double* LLR, double* L, char* decoded)
{
	int N = mod2sparse_cols(H);
	double sum;

	for (int j = 0; j < N; j++)
	{
		sum = LLR[j];//prior L_pr
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			// extrinsic messages L_ext
			sum += e->msg_c_to_v_fp;
		} // end for(e)

		L[j] = sum;	// posterior probability L_ps

		if (sum > 0)		decoded[j] = 0;
		else				decoded[j] = 1;
	}//end for (j)
}

//===========================================================================
// 유한한 message size를 가지는 MSA : 설정
//===========================================================================
void Set_MSA(int q, double step, int beta)
{
	precision = q;
	step_size = step;

	if (precision > 0)
	{
		min_value = (int)-(pow(2.0, (precision - 1)) - 1);
		max_value = (int)(pow(2.0, (precision - 1)) - 1);
	}
	else
	{
		//floating point min-sum decoder
		min_value = 0;
		max_value = 0;
	}

	offset_beta = beta;
}

//===========================================================================
// 유한한 message size를 가지는 MSA : message quantization
// type = 0 : uniform quantization
//        1 : quasi-uniform quantization "Quantized Min-Sum Decoders with Low Error Floor for LDPC codes" - ISIT 2012
//===========================================================================
int Cal_MSA_Q(double x, int type)
{
	int k = 0;
	int sign = 1;
	double mag = abs(x);

	if (x >= 0)		sign = 1;
	else			sign = -1;

	if (type == 0)
	{
		// uniform quantization
		k = (int)(mag / step_size + 0.5);
		if (sign == 1)
		{
			if (k > max_value)
			{
				k = max_value;
			}
		}
		else
		{
			if (k > (-min_value))
			{
				k = (-min_value);
			}
			k *= sign;
		}
	}
	else if (type == 1)
	{
		// quasi-uniform quantization
	}

	return k;
}

//===========================================================================
// 유한한 message size를 가지는 MSA : message clip
//===========================================================================
int Cal_MSA_Clip(int x, int& b)
{
	int k = x;
	b = FALSE;
	if (x > max_value)
	{
		k = max_value;
		b = TRUE;
	}
	else if (x < min_value)
	{
		k = min_value;
		b = TRUE;
	}

	return k;
}

int Set_Message_Maxvalue(int x)
{
	int mag = abs(x);
	int sign = 1;
	if (x < 0)
		sign = -1;

	//mag = max_value;
	mag = mag / 2;

	return (sign * mag);
}

int Dec_Message_value(int x)
{
	int mag = abs(x);
	int sign = 1;
	if (x < 0)
		sign = -1;

	mag = mag / 2;

	return (sign * mag * (-1));
}



//===============================================================================
//
//===============================================================================
void Save_State(int n, int c, char* decoded, char* check, int& min_c)
{
	int idx_min;
	int err_cnt;
	int update;

	idx_min = save_state_num;

	err_cnt = 0;
	update = FALSE;

	//error가 생긴 variable node 수
	for (int j = 0; j < N; j++)
	{
		if (decoded[j] != 0)
			err_cnt++;
	}

	if (n == 0)
	{
		min_c = c;
		update = TRUE;
	}
	else
	{
		if ((min_c > c) && (c > 0))
		{
			min_c = c;
			update = TRUE;
		}
	}

	g_unsat_check_num[n] = c;
	g_err_num[n] = err_cnt;

	for (int j = 0; j < N; j++)
	{
		g_word_state[n][j] = decoded[j];
	}

	for (int j = 0; j < M; j++)
	{
		g_check_state[n][j] = check[j];
	}
}

void Print_Variable_State(mod2sparse* H, int n, int c, char* decoded, int* qLLR, int* L_Q, char* pchk, FILE** fp, int* w_v_idx, FILE* fp1)
{
	int N = mod2sparse_cols(H);
	int msg;
	int b;
	int check_idx;
	for (int j = 0; j < N; j++)
	{
		if (w_v_idx[j] != 0)
		{
			if (n == 0)	fprintf(fp[j], "<iter>\tM1\tM2\tM3\tL\tMM1\tMM2\tMM3\tB\n");

			fprintf(fp[j], "<%d>\t", n);

			//=============================
			// 나가는 msg 출력
			//=============================
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				msg = e->msg_v_to_c;
				fprintf(fp[j], "%d\t", msg);
			}


			//=============================
			// 수신값 출력
			//=============================
			fprintf(fp[j], "%d\t", qLLR[j]);

			//=============================
			// 들어오는 msg 출력
			//=============================
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				msg = e->msg_c_to_v;
				fprintf(fp[j], "%d\t", msg);
			}
			fprintf(fp[j], "%d\t\n", L_Q[j]);
		}
	}



	//===================================================================================
	// unsatisfied check에 연결된 실제 오류가 난 variable node를 출력한다
	//===================================================================================
	fprintf(fp1, "%d\t/\t", n);
	for (int j = 0; j < N; j++)
	{
		//오류가 생긴 variable node 중에서 unsatisfied check과 연결된 것을 찾는다.
		if (decoded[j] != 0)
		{
			for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
			{
				check_idx = e->row;
				if (pchk[check_idx] != 0)
				{
					fprintf(fp1, " (C : %d, V : %d) ", check_idx, j);
					break;
				}
			}
		}
	}
	fprintf(fp1, "\n");
}



int Run_Pipeline_Decoder(int b, int c, int m, int L, mod2sparse* H, double* lratio, char* dblk, char* pchk, int* bIsCodeword)
{
	b_SC = b;
	c_SC = c;
	w_SC = m;
	L_conv = L;
	int chk;
	Init_Pipeline(H, lratio, dblk);
	Pipeline(H, lratio, dblk);

	chk = check(H, dblk, pchk);
	if (chk == 0)
		* bIsCodeword = TRUE;


	return max_iter;
}

/* INITIALIZE PROBABILITY PROPAGATION.  Stores initial ratios, probabilities, and guess at decoding. */
void Init_Pipeline(mod2sparse* H, double* lratio, char* dblk)
{
	mod2entry* e;
	int N;
	int j;

	N = mod2sparse_cols(H);
	for (j = 0; j < N; j++)
	{
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->pr = lratio[j];				//	p0/p1
			e->lr = 1;

			//p0 = 1 - 1/(1+lratio[j])
			//p1 = 1/(1+lration[j])
		}//end for(e)

		dblk[j] = (lratio[j] <= 1);	// p0/p1 > 1  : 0
		// p0/p1 <= 1 : 1
	}//end for(j)

}

void Pipeline(mod2sparse* H, double* lratio, char* dblk)
{

	int N, M;
	int lasttime;
	cst_row = (w_SC + 1) * c_SC;
	cst_col = (w_SC + 1) * (c_SC - b_SC);
	M = mod2sparse_rows(H);
	N = mod2sparse_cols(H);
	lasttime = L_conv + (w_SC + 1) * max_iter;



	//FILE * temp1 = fopen("PL_activated_chk.txt","wt");
	//FILE * temp2 = fopen("PL_activated_var.txt","wt");
	//FILE * temp3 = fopen("PL_activated_dec.txt","wt");

	//===============================================================================
	// pipeline decoding을 수행한다
	//===============================================================================


	// t = time index, decoding window에서 가장 아래 time index에 해당
	for (int t = 0; t < lasttime; t++)
	{
		Decoding_window(t, H, lratio, dblk);
	}

	//int i;
	//for(i=0;i<39;i++) fprintf(temp1,"%d ",PL_activated_chk[i]);
	//for(i=0;i<70;i++) fprintf(temp2,"%d ",PL_activated_var[i]);
	//for(i=0;i<70;i++) fprintf(temp3,"%d ",PL_activated_dec[i]);

	//fclose(temp1);
	//fclose(temp2);
	//fclose(temp3);


}

void Decoding_window(int time, mod2sparse* H, double* lratio, char* dblk)
{
	int i, j;
	int index_chk;
	int index_var;



	for (i = 0; i < max_iter; i++)
	{
		for (j = 0; j < (c_SC - b_SC); j++)
		{
			index_chk = (c_SC - b_SC) * (time - i * (w_SC + 1)) + j;
			if ((index_chk >= 0) && (index_chk < M))
			{
				chk_process(index_chk, H, lratio, dblk);
				//	PL_activated_chk[index_chk] ++;
			}
		}

		for (j = 0; j < c_SC; j++)
		{
			index_var = c_SC * (time - i * (w_SC + 1) - w_SC) + j;
			if ((index_var >= 0) && (index_var < N))
			{
				var_process(i + 1, index_var, H, lratio, dblk);
				//	PL_activated_var[index_var] ++;
			}
		}

	}

}

void chk_process(int ind_chk, mod2sparse* H, double* lratio, char* dblk)
{
	double pr, dl, t;
	mod2entry* e;
	int i, j;

	dl = 1;
	for (e = mod2sparse_first_in_row(H, ind_chk); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
	{
		e->lr = dl;
		dl *= 1 - 2 / (1 + e->pr);
	}

	dl = 1;
	for (e = mod2sparse_last_in_row(H, ind_chk); !mod2sparse_at_end(e); e = mod2sparse_prev_in_row(e))
	{
		t = e->lr * dl;
		e->lr = (1 + t) / (1 - t);		// (수정)			
		dl *= 1 - 2 / (1 + e->pr);	// (수정)
	}


}

void var_process(int iter, int ind_var, mod2sparse* H, double* lratio, char* dblk)
{
	double pr, dl, t;
	mod2entry* e;
	int i, j;

	pr = lratio[ind_var];
	for (e = mod2sparse_first_in_col(H, ind_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
	{
		e->pr = pr;
		pr *= e->lr;
	}//end for(e)

	if (_isnan(pr))
		pr = 1;

	//p0 = 1 - 1/(1+pr)

	if (iter == max_iter)
	{
		dblk[ind_var] = (pr <= 1);
		//	PL_activated_dec[ind_var] ++;
	}// (수정)
	pr = 1;

	for (e = mod2sparse_last_in_col(H, ind_var); !mod2sparse_at_end(e); e = mod2sparse_prev_in_col(e))
	{
		e->pr *= pr;
		if (_isnan(e->pr))
		{
			e->pr = 1;
		}
		pr *= e->lr;
	}//end for(e)

}

////=====================================
//// Sliding window decoder
////=====================================
//
int Run_SW_Decoder
(
	mod2sparse* H,			/* Parity check matrix */
	double* lratio,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{

	int n, c;
	int N = mod2sparse_cols(H);
	double avg_iter = 0;


	int SC_D;
	if (code_type == 0)
		SC_D = SC_L + SC_w - 1;
	else
		SC_D = SC_L + (SC_w - 1) / 2;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;

	int V_Check_End = 0;
	int C_Check_End = 0;
	for (int i = 0; i < SC_w; i++)
	{
		V_Check_End += Mv[i];
		C_Check_End += Mc[i];
	}

	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}

	int D_Start = 0;
	int D_End = Mv[0];


	Init_SW_Decoder(H, lratio, V_Start, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	avg_iter += Iter_SW_Decoder(H, dblk, pchk, lratio, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Check_End, C_Check_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);


	for (int t = 1; t < SC_L; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		if (t + SC_WIN >= SC_L)
			V_End = N;
		else
			V_End += Mv[t + SC_WIN - 1];

		if (t + SC_WIN >= SC_D)
			C_End = M;
		else
			C_End += Mc[t + SC_WIN - 1];


		D_Start += Mv[t - 1];
		D_End += Mv[t];

		if (t + SC_w >= SC_L)
		{
			V_Check_End = N;
			C_Check_End = M;
		}
		else
		{
			V_Check_End += Mv[t + SC_w - 1];
			C_Check_End += Mc[t + SC_w - 1];
		}

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
		if (t + SC_WIN <= SC_L)
		{

			Init_SW_Decoder(H, lratio, V_End - Mv[t + SC_WIN - 1], V_End);
		}

		avg_iter += Iter_SW_Decoder(H, dblk, pchk, lratio, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Check_End, C_Check_End);
		test_BER(1, dblk, position_BER, t, D_Start, D_End);

	}



	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;


	avg_iter = floor((double)avg_iter / SC_L);
	return avg_iter;

}

int Run_SW_Decoder_SKU //proposed by SKU
(
	mod2sparse* H,			/* Parity check matrix */
	double* lratio,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{

	int n, c;
	int N = mod2sparse_cols(H);
	double avg_iter = 0;


	int SC_D;
	if (code_type == 0)
		SC_D = SC_L + SC_w - 1;
	else
		SC_D = SC_L + (SC_w - 1) / 2;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;

	int V_Check_End = 0;
	int C_Check_End = 0;
	for (int i = 0; i < SC_w; i++)
	{
		V_Check_End += Mv[i];
		C_Check_End += Mc[i];
	}

	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}

	int D_Start = 0;
	int D_End = Mv[0];


	Init_SW_Decoder(H, lratio, V_Start, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	avg_iter += Iter_SW_Decoder(H, dblk, pchk, lratio, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Check_End, C_Check_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);


	for (int t = 1; t < SC_L - SC_WIN; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		V_End += Mv[t + SC_WIN - 1];
		C_End += Mc[t + SC_WIN - 1];


		D_Start += Mv[t - 1];
		D_End += Mv[t];
		V_Check_End += Mv[t + SC_w - 1];
		C_Check_End += Mc[t + SC_w - 1];

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
		Init_SW_Decoder(H, lratio, V_End - Mv[t + SC_WIN - 1], V_End);
		avg_iter += Iter_SW_Decoder(H, dblk, pchk, lratio, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Check_End, C_Check_End);
		test_BER(1, dblk, position_BER, t, D_Start, D_End);
	}

	V_Start += Mv[SC_L - SC_WIN - 1];
	C_Start += Mc[SC_L - SC_WIN - 1];
	V_End = N;
	C_End = M;
	D_Start += Mv[SC_L - SC_WIN - 1];
	D_End = N;
	V_Check_End = N;
	C_Check_End = M;

	test_BER(0, dblk, position_BER, SC_L - SC_WIN, D_Start, D_End);
	Init_SW_Decoder(H, lratio, V_End - Mv[SC_L - 1], V_End);
	avg_iter += Iter_SW_Decoder(H, dblk, pchk, lratio, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Check_End, C_Check_End);
	test_BER(1, dblk, position_BER, SC_L - SC_WIN, D_Start, D_End);


	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;


	avg_iter = floor((double)avg_iter / SC_L);
	return avg_iter;

}

int Run_SWW_Decoder
(
	mod2sparse* H,			/* Parity check matrix */
	double* lratio,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_b,
	int SC_c,
	int SC_L,
	int SC_w,
	int SC_z,
	int SC_WIN,
	int SC_cr
)
{

	int n, c;
	int N = mod2sparse_cols(H);


	char* dblk2 = (char*)calloc(N, sizeof(int));
	double* lratio2 = (double*)calloc(N, sizeof(double));


	int V_Start;
	int V_End;
	int C_Start;
	int C_End;
	int D_Start;
	int D_End;


	for (int t = 0; t < SC_L; t++)
	{
		V_Start = t * SC_c * SC_z;
		V_End = MIN((t + SC_WIN) * SC_c * SC_z, N);
		C_Start = t * SC_b * SC_z;
		C_End = MIN((t + SC_WIN) * SC_b * SC_z, M);
		D_Start = t * SC_c * SC_z;
		D_End = MIN((t + 1) * SC_c * SC_z, N);
		Init_SW_Decoder(H, lratio, V_Start, V_End);
		Iter_SWW_Decoder(H, dblk, lratio, V_Start, V_End, C_Start, C_End, D_Start, D_End);
		//			Fixing_Message_SW(H,lratio,dblk,t,SC_c,SC_z);
	}



	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(lratio2);
	return 0;

}


//
///* INITIALIZE PROBABILITY PROPAGATION.  Stores initial ratios, probabilities, and guess at decoding. */
void Init_SW_Decoder(mod2sparse* H, double* lratio, int V_Start, int V_End)
{
	mod2entry* e;


	for (int j = V_Start; j < V_End; j++)
	{
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->pr = lratio[j];				//	p0/p1
			e->lr = 1;

			//p0 = 1 - 1/(1+lratio[j])
			//p1 = 1/(1+lration[j])
		}//end for(e)

		//		dblk[j] = (lratio[j] <= 1);	// p0/p1 > 1  : 0
		// p0/p1 <= 1 : 1
	}//end for(j)

}

int Iter_SW_Decoder(mod2sparse* H, char* dblk, char* pchk, double* lratio, int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End, int V_Check_End, int C_Check_End)
{
	int n;
	int exit = FALSE;

	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);


	int c;


	for (n = 0; ; n++)
	{
		for (int i = C_Start; i < C_End; i++)
		{

			Check_Update_SW(i, H, lratio);
		}

		for (int i = V_Start; i < V_End; i++)
		{
			Variable_Update_SW(i, C_Start, C_End, H, lratio);
		}
		for (int i = V_Start; i < V_End; i++)
		{
			Decision_SW(i, C_Start, C_End, H, lratio, dblk);
		}

		c = check_bound(H, dblk, pchk, V_Start, V_Check_End, C_Start, C_Check_End);
		//for(int i=D_Start;i<D_End;i++)
		//{
		//	if(dblk[i]!=0)
		//		break;
		//	if(i==D_End-1)
		//		exit = TRUE;
		//}

		if (n == max_iter || c == 0)
			break;

	}

	return n;
}

void Iter_SWW_Decoder(mod2sparse* H, char* dblk, double* lratio, int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End)
{
	int n;
	int exit = FALSE;

	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);





	for (n = 0; ; n++)
	{
		for (int i = C_Start; i < C_End; i++)
		{

			Check_Update_SW(i, H, lratio);
		}

		for (int i = V_Start; i < V_End; i++)
		{
			Variable_Update_SWW(i, C_Start, C_End, H, lratio, dblk, 1);
		}
		for (int i = V_Start; i < V_End; i++)
		{
			Decision_SWW(i, C_Start, C_End, H, lratio, dblk, 1.05);
		}
		for (int i = D_Start; i < D_End; i++)
		{
			if (dblk[i] != 0)
				break;
			if (i == D_End - 1)
				exit = TRUE;
		}

		if (n == max_iter || exit)
			break;

	}


}

void Check_Update_SW(int index_chk, mod2sparse* H, double* lratio)
{
	double pr, dl, t;
	mod2entry* e;
	int i, j;

	dl = 1;
	for (e = mod2sparse_first_in_row(H, index_chk); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
	{
		e->lr = dl;
		dl *= 1 - 2 / (1 + e->pr);
	}

	dl = 1;
	for (e = mod2sparse_last_in_row(H, index_chk); !mod2sparse_at_end(e); e = mod2sparse_prev_in_row(e))
	{
		t = e->lr * dl;
		e->lr = (1 + t) / (1 - t);		// (수정)			
		dl *= 1 - 2 / (1 + e->pr);	// (수정)
	}


}

void Variable_Update_SW(int index_var, int index_chk_start, int index_chk_end, mod2sparse* H, double* lratio)
{
	double pr, dl, t;
	mod2entry* e;
	mod2entry* temp;
	int i, j;

	pr = lratio[index_var];
	for (e = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
	{
		if (mod2sparse_row(e) < index_chk_end && mod2sparse_row(e) >= index_chk_start)
		{
			e->pr = pr;
			pr *= e->lr;
		}
		//if(mod2sparse_row(e) >= index_chk_end || mod2sparse_row(e) < index_chk_start ) break;
		//e->pr = pr;
		//pr *= e->lr;
	}

	if (_isnan(pr))
		pr = 1;



	pr = 1;
	temp = mod2sparse_last_in_col(H, index_var);
	while (mod2sparse_row(temp) >= index_chk_end)
	{
		temp = mod2sparse_prev_in_col(temp);
	}


	for (e = temp; !mod2sparse_at_end(e); e = mod2sparse_prev_in_col(e))
	{
		if (mod2sparse_row(e) < index_chk_end && mod2sparse_row(e) >= index_chk_start)
		{
			e->pr *= pr;
			if (_isnan(e->pr))
			{
				e->pr = 1;
			}
			pr *= e->lr;
		}
	}//end for(e)

}
//

void Variable_Update_SWW(int index_var, int index_chk_start, int index_chk_end, mod2sparse* H, double* lratio, char* dblk, double alpha)
{
	double pr, dl, t;
	mod2entry* e;
	mod2entry* temp;
	int i, j;

	pr = lratio[index_var];
	i = 1;
	for (e = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
	{
		if (mod2sparse_row(e) >= index_chk_end || mod2sparse_row(e) < index_chk_start) break;

		e->pr = pr;
		if (i == 1)
		{
			pr *= pow(e->lr, alpha);
		}
		else
		{
			pr *= e->lr;
		}
		i++;
	}

	if (_isnan(pr))
		pr = 1;



	pr = 1;
	temp = mod2sparse_last_in_col(H, index_var);
	while (mod2sparse_row(temp) >= index_chk_end)
	{
		temp = mod2sparse_prev_in_col(temp);
	}

	i = 1;
	for (e = temp; !mod2sparse_at_end(e); e = mod2sparse_prev_in_col(e))
	{
		if (mod2sparse_row(e) >= index_chk_end || mod2sparse_row(e) < index_chk_start) break;
		e->pr *= pr;
		if (_isnan(e->pr))
		{
			e->pr = 1;
		}
		if (i == 1)
		{
			pr *= pow(e->lr, 1 / alpha);
		}
		else
		{
			pr *= e->lr;
		}

		i++;
	}//end for(e)

}



void Fixing_Message_SW(mod2sparse* H, double* lratio, char* dblk, int VP, int SC_c, int SC_z)
{
	mod2entry* e;

	for (int j = VP * SC_c * SC_z; j < (VP + 1) * SC_c * SC_z; j++)
	{
		for (e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			if (dblk[j] == 0)
				e->pr = exp(40.0);
			else
				e->pr = 0;
		}//end for(e)
	}//end for(j)

}
//
void Decision_SW(int index_var, int index_chk_start, int index_chk_end, mod2sparse* H, double* lratio, char* dblk)
{
	double pr;
	mod2entry* e;
	pr = lratio[index_var];
	for (e = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
	{
		if (mod2sparse_row(e) < index_chk_end && mod2sparse_row(e) >= index_chk_start)
		{
			pr *= e->lr;
		}
	}

	dblk[index_var] = (pr <= 1);
}


void Decision_SWW(int index_var, int index_chk_start, int index_chk_end, mod2sparse* H, double* lratio, char* dblk, double alpha)
{
	double pr;
	mod2entry* e;
	pr = lratio[index_var];
	int i = 1;
	for (e = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
	{
		if (mod2sparse_row(e) >= index_chk_end || mod2sparse_row(e) < index_chk_start) break;
		if (i == 1)
		{
			pr *= pow(e->lr, alpha);
		}
		else if (i == 3)
		{
			pr *= pow(e->lr, 1 / alpha);
		}

		pr *= e->lr;
		i++;
	}

	dblk[index_var] = (pr <= 1);
}



//===============================================================================
// BEC Sliding window Decoder
//===============================================================================
int Run_BEC_SW_Decoder
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{
	int n, c;
	int sum = 0;
	int N = mod2sparse_cols(H);
	int avg_iter = 0;

	char* dblk_prev;
	dblk_prev = (char*)calloc(N, sizeof(int));

	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;
	int SC_D;
	if (code_type == 0)
		SC_D = SC_L + SC_w - 1;
	else if (code_type == 1)
		SC_D = SC_L + (SC_w - 1) / 2;


	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;

	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}


	int D_Start = 0;
	int D_End = Mv[0];

	Init_BEC_SW_Decoder(H, recv, dblk, V_Start, V_End);
	avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);

	for (int t = 1; t < SC_D; t++)
	{
		if (Mv[t] != 0)
		{
			V_Start += Mv[t - 1];
			C_Start += Mc[t - 1];
			if (t + SC_WIN - 1 >= SC_D)
			{
				V_End = N;
				C_End = M;
			}
			else
			{
				V_End += Mv[t + SC_WIN - 1];
				C_End += Mc[t + SC_WIN - 1];
			}

			D_Start += Mv[t - 1];
			D_End += Mv[t];

			Init_BEC_SW_Decoder(H, recv, dblk, V_Start, V_End);
			avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
		}
		else
			break;
	}





	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);

	avg_iter = floor((double)avg_iter / SC_L);
	return avg_iter;

}

int Run_BEC_SW_Decoder_OC
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_Ls,
	int SC_w,
	int SC_WIN,
	int Mv,
	int Mc,
	int etha
)
{
	int n, c;
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);

	SC_Ls = (SC_Ls - 1) / etha;
	char* dblk_prev = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int Ns = SC_Ls * Mv;
	int Ms = (SC_Ls + (SC_w - 1) / 2) * Mc;
	int SC_Ds = SC_Ls + (SC_w - 1) / 2;


	int* V_Start = (int*)calloc(etha, sizeof(int));
	int* V_End = (int*)calloc(etha, sizeof(int));
	int* C_Start = (int*)calloc(etha, sizeof(int));
	int* C_End = (int*)calloc(etha, sizeof(int));
	int* D_Start = (int*)calloc(etha, sizeof(int));
	int* D_End = (int*)calloc(etha, sizeof(int));
	int* Target_pos = (int*)calloc(etha, sizeof(int));

	int RV_Start = etha * Ns;
	int RV_End = N;


	for (int p = 0; p < etha; p++)
	{
		V_Start[p] = p * Ns;
		V_End[p] = p * Ns + Mv * (SC_WIN - 2);
		C_Start[p] = p * Ms;
		C_End[p] = p * Ms + Mc * (SC_WIN - 2);
	}

	Init_BEC_SW_Decoder(H, recv, dblk, RV_Start, RV_End);
	for (int t = 0; t < etha; t++)
		Init_BEC_SW_Decoder(H, recv, dblk, V_Start[t], V_End[t]);


	Iter_BEC_OC_Init_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, RV_Start, RV_End, etha);


	test_BER(1, dblk, position_BER, etha * SC_Ls, RV_Start, RV_End);

	for (int p = 0; p < etha; p++)
	{
		V_Start[p] = p * Ns;
		V_End[p] = p * Ns + Mv * SC_WIN;
		C_Start[p] = p * Ms;
		C_End[p] = p * Ms + Mc * SC_WIN;
		D_Start[p] = p * Ns;
		D_End[p] = p * Ns + Mv;
		Target_pos[p] = p * SC_Ls;
	}

	for (int p = 0; p < etha; p++)
		Init_BEC_SW_Decoder(H, recv, dblk, V_Start[p], V_End[p]);

	for (int p = 0; p < etha; p++)
	{
		Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start[p], V_End[p], C_Start[p], C_End[p], D_Start[p], D_End[p]);
		test_BER(1, dblk, position_BER, Target_pos[p], D_Start[p], D_End[p]);
	}



	for (int t = 1; t < SC_Ls; t++)
	{
		for (int p = 0; p < etha; p++)
		{
			V_Start[p] += Mv;

			C_Start[p] += Mc;

			D_Start[p] += Mv;
			D_End[p] += Mv;
			Target_pos[p] ++;
		}

		if (t + SC_WIN <= SC_Ls)
		{
			for (int p = 0; p < etha; p++)
			{
				V_End[p] += Mv;
				Init_BEC_SW_Decoder(H, recv, dblk, V_End[p] - Mv, V_End[p]);
			}
		}
		if (t + SC_WIN <= SC_Ds)
		{
			for (int p = 0; p < etha; p++)
				C_End[p] += Mc;
		}

		for (int p = 0; p < etha; p++)
		{
			Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start[p], V_End[p], C_Start[p], C_End[p], D_Start[p], D_End[p]);
			test_BER(1, dblk, position_BER, Target_pos[p], D_Start[p], D_End[p]);
		}
	}

	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(recv2);
	free(dblk_prev);
	return 0;

}

int Run_BEC_SW_Decoder_Two
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{
	int n, c;
	int sum = 0;
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);

	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int SC_Ls = SC_L / 2;
	int SC_D = SC_L + SC_w - 1;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;
	int R_Start = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}
	for (int i = 0; i < SC_WIN - 1; i++)
		R_Start += Mv[i];
	int D_Start = 0;
	int D_End = Mv[0];
	int R_End = V_End;

	Init_BEC_SW_Decoder(H, recv, dblk, 0, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);


	int	V_Start2 = N - V_End;
	int	V_End2 = N - V_Start;
	int	C_Start2 = M - C_End;
	int	C_End2 = M - C_Start;
	int	D_Start2 = N - D_End;
	int	D_End2 = N - D_Start;
	int R_Start2 = N - R_End;
	int R_End2 = N - R_Start;

	Init_BEC_SW_Decoder(H, recv, dblk, V_Start2, N);
	test_BER(0, dblk, position_BER, SC_L - 1, D_Start2, D_End2);
	Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start2, V_End2, C_Start2, C_End2, D_Start2, D_End2);
	test_BER(1, dblk, position_BER, SC_L - 1, D_Start2, D_End2);

	for (int t = 1; t < SC_Ls; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		V_End += Mv[t + SC_WIN - 1];
		C_End += Mc[t + SC_WIN - 1];
		D_Start += Mv[t - 1];
		D_End += Mv[t];

		V_Start2 = N - V_End;
		V_End2 = N - V_Start;
		C_Start2 = M - C_End;
		C_End2 = M - C_Start;
		D_Start2 = N - D_End;
		D_End2 = N - D_Start;


		if (t + SC_WIN <= SC_Ls)
		{
			Init_BEC_SW_Decoder(H, recv, dblk, V_End - Mv[t + SC_WIN - 1], V_End);
			Init_BEC_SW_Decoder(H, recv, dblk, V_Start2, V_Start2 + Mv[SC_L - t - SC_WIN - 1]);
		}

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
		Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
		test_BER(1, dblk, position_BER, t, D_Start, D_End);


		test_BER(0, dblk, position_BER, SC_L - t - 1, D_Start2, D_End2);
		Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start2, V_End2, C_Start2, C_End2, D_Start2, D_End2);
		test_BER(1, dblk, position_BER, SC_L - t - 1, D_Start2, D_End2);
	}

	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);
	return 0;

}

int Run_BEC_SW_Decoder_Two_Cross
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{
	int n, c;
	int sum = 0;
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);

	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int SC_Ls = SC_L / 2;
	int SC_D = SC_L + SC_w - 1;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;
	int R_Start = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}
	for (int i = 0; i < SC_WIN - 1; i++)
		R_Start += Mv[i];
	int D_Start = 0;
	int D_End = Mv[0];
	int R_End = V_End;

	Init_BEC_SW_Decoder(H, recv, dblk, 0, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);


	int	V_Start2 = N - V_End;
	int	V_End2 = N - V_Start;
	int	C_Start2 = M - C_End;
	int	C_End2 = M - C_Start;
	int	D_Start2 = N - D_End;
	int	D_End2 = N - D_Start;
	int R_Start2 = N - R_End;
	int R_End2 = N - R_Start;

	Init_BEC_SW_Decoder(H, recv, dblk, V_Start2, N);
	test_BER(0, dblk, position_BER, SC_L - 1, D_Start2, D_End2);
	Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start2, V_End2, C_Start2, C_End2, D_Start2, D_End2);
	test_BER(1, dblk, position_BER, SC_L - 1, D_Start2, D_End2);

	for (int t = 1; t < SC_L; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		V_End += Mv[t + SC_WIN - 1];
		C_End += Mc[t + SC_WIN - 1];
		D_Start += Mv[t - 1];
		D_End += Mv[t];

		if (V_End > N)
			V_End = N;
		if (C_End > M)
			C_End = M;

		V_Start2 = N - V_End;
		V_End2 = N - V_Start;
		C_Start2 = M - C_End;
		C_End2 = M - C_Start;
		D_Start2 = N - D_End;
		D_End2 = N - D_Start;


		if (t + SC_WIN <= SC_Ls)
		{
			Init_BEC_SW_Decoder(H, recv, dblk, V_End - Mv[t + SC_WIN - 1], V_End);
			Init_BEC_SW_Decoder(H, recv, dblk, V_Start2, V_Start2 + Mv[SC_L - t - SC_WIN - 1]);
		}

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
		Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
		test_BER(1, dblk, position_BER, t, D_Start, D_End);


		test_BER(0, dblk, position_BER, SC_L - t - 1, D_Start2, D_End2);
		Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start2, V_End2, C_Start2, C_End2, D_Start2, D_End2);
		test_BER(1, dblk, position_BER, SC_L - t - 1, D_Start2, D_End2);
	}

	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);
	return 0;

}

int Run_BEC_SW_Decoder_Two_Indi
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{

	int n, c;
	int sum = 0;
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);

	mod2sparse* H2 = mod2sparse_allocate(M, N);
	mod2sparse_copy(H, H2);

	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int SC_Ls = SC_L / 2;
	int SC_D = SC_L + SC_w - 1;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;
	int R_Start = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}
	for (int i = 0; i < SC_WIN - 1; i++)
		R_Start += Mv[i];
	int D_Start = 0;
	int D_End = Mv[0];
	int R_End = V_End;

	Init_BEC_SW_Decoder(H, recv, dblk, 0, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);


	int	V_Start2 = N - V_End;
	int	V_End2 = N - V_Start;
	int	C_Start2 = M - C_End;
	int	C_End2 = M - C_Start;
	int	D_Start2 = N - D_End;
	int	D_End2 = N - D_Start;
	int R_Start2 = N - R_End;
	int R_End2 = N - R_Start;

	Init_BEC_SW_Decoder(H2, recv, dblk2, V_Start2, N);
	test_BER(0, dblk2, position_BER, SC_L - 1, D_Start2, D_End2);
	Iter_BEC_SW_Decoder(H2, dblk2, dblk_prev, V_Start2, V_End2, C_Start2, C_End2, D_Start2, D_End2);
	test_BER(1, dblk2, position_BER, SC_L - 1, D_Start2, D_End2);

	for (int t = 1; t < SC_Ls; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		V_End += Mv[t + SC_WIN - 1];
		C_End += Mc[t + SC_WIN - 1];
		D_Start += Mv[t - 1];
		D_End += Mv[t];

		if (V_End > N)
			V_End = N;
		if (C_End > M)
			C_End = M;

		V_Start2 = N - V_End;
		V_End2 = N - V_Start;
		C_Start2 = M - C_End;
		C_End2 = M - C_Start;
		D_Start2 = N - D_End;
		D_End2 = N - D_Start;

		Init_BEC_SW_Decoder(H, recv, dblk, V_End - Mv[t + SC_WIN - 1], V_End);
		Init_BEC_SW_Decoder(H2, recv, dblk2, V_Start2, V_Start2 + Mv[SC_L - t - SC_WIN - 1]);

		if (t + SC_WIN <= SC_Ls)
		{

		}

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
		Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
		test_BER(1, dblk, position_BER, t, D_Start, D_End);


		test_BER(0, dblk2, position_BER, SC_L - t - 1, D_Start2, D_End2);
		Iter_BEC_SW_Decoder(H2, dblk2, dblk_prev, V_Start2, V_End2, C_Start2, C_End2, D_Start2, D_End2);
		test_BER(1, dblk2, position_BER, SC_L - t - 1, D_Start2, D_End2);
	}

	for (int i = N / 2; i < N; i++)
		dblk[i] = dblk2[i];

	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);
	mod2sparse_free(H2);
	return 0;

}


int Run_BEC_SW_Decoder_SAVE
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{
	///////////////////////////////////////
/*
	int count = 0;
	int idx_s;
	int idx_e;
	int num_e[20];
*/

////////////////////////////


	int n, c;
	int sum = 0;
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);
	double avg_iter = 0;

	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int SC_D;
	if (code_type == 0)
		SC_D = SC_L + SC_w - 1;
	else
		SC_D = SC_L + (SC_w - 1) / 2;

	int V_Start = 0;
	int V_End = 0;
	int V_End_prev;
	int C_Start = 0;
	int C_End = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}
	V_End_prev = V_End;
	int D_Start = 0;
	int D_End = Mv[0];

	//	Init_BEC_SW_Decoder(H,recv,dblk,0,N);



		///////////////////////

		//V_End = N;
		//C_End = M;
		//D_End = N;
		///////////////////////

	Init_BEC_SW_Decoder(H, recv, dblk, V_Start, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);
	//

	//
	for (int t = 1; t < SC_L; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		if (t + SC_WIN - 1 >= SC_L - 1)
			V_End = N;
		else
			V_End += Mv[t - 1 + SC_WIN];
		if (t + SC_WIN - 1 >= SC_D - 1)
			C_End = M;
		else
			C_End += Mc[t - 1 + SC_WIN];

		D_Start += Mv[t - 1];
		D_End += Mv[t];


		if (t + SC_WIN <= SC_L)
		{
			Init_BEC_SW_Decoder(H, recv, dblk, V_End_prev, V_End);
		}
		V_End_prev = V_End;
		//		Init_BEC_SW_Decoder(H,recv,dblk,V_Start,V_End);

		/////////////////////////////
		//idx_s = 0;
		//idx_e = Mv[0];
		//printf("Final result");
		//for(int i=0;i<SC_L;i++)
		//{
		//	count = 0;
		//	if(i!=0)
		//	{
		//		idx_s += Mv[i-1];
		//		idx_e += Mv[i];
		//	}
		//	for(int j=idx_s;j<idx_e;j++)
		//	{
		//		if(dblk[j]==ERASE_MARK)
		//			count ++;
		//	}
		//	printf("(%d,%d) ",i+1,count);
		//}
		//printf("\n\n");
		/////////////////////////////////////

		test_BER(0, dblk, position_BER, t, D_Start, D_End);
		avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
		test_BER(1, dblk, position_BER, t, D_Start, D_End);

	}





	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);

	avg_iter = floor(avg_iter / SC_L);
	return avg_iter;

}


int Run_BEC_SW_Decoder_TARGET
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{
	int N = mod2sparse_cols(H);
	double avg_iter = 0;

	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}

	int D_Start = 0;
	int D_End = Mv[0];

	Init_BEC_SW_Decoder(H, recv, dblk, V_Start, V_End);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);



	free(dblk2);
	free(recv2);
	free(dblk_prev);

	avg_iter = floor(avg_iter);
	return avg_iter;

}


int Run_BEC_SW_Decoder_RA
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc
)
{
	int n, c;
	int sum = 0;
	int N1 = 0;
	int N = mod2sparse_cols(H);
	double avg_iter = 0;

	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int SC_D;
	if (code_type == 0)
		SC_D = SC_L + SC_w - 1;
	else
		SC_D = SC_L + (SC_w - 1) / 2;

	int V_Start = 0;
	int V_End = 0;
	int C_Start = 0;
	int C_End = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}

	int D_Start = 0;
	int D_End = Mv[0];

	int V_Start2 = 0;
	int V_End2 = 0;
	int D_Start2 = 0;
	int D_End2 = 0;

	for (int i = 0; i < SC_D; i++)
	{
		V_Start2 += Mv[i];
	}
	N1 = V_Start2;
	V_End2 = V_Start2;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End2 += Mc[i];
	}
	D_Start2 = V_Start2;
	D_End2 = V_Start2 + Mc[0];



	Init_BEC_SW_Decoder(H, recv, dblk, V_Start, V_End);
	Init_BEC_SW_Decoder(H, recv, dblk, V_Start2, V_End2);
	test_BER(0, dblk, position_BER, 0, D_Start, D_End);
	avg_iter += Iter_BEC_RA_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Start2, V_End2, D_Start2, D_End2);
	test_BER(1, dblk, position_BER, 0, D_Start, D_End);

	for (int t = 1; t < SC_D; t++)
	{
		V_Start += Mv[t - 1];
		C_Start += Mc[t - 1];
		if (t <= SC_L - SC_WIN)
		{
			V_End += Mv[t + SC_WIN - 1];
		}
		if (t <= SC_D - SC_WIN)
			C_End += Mc[t + SC_WIN - 1];

		D_Start += Mv[t];
		D_End += Mv[t];


		D_Start2 += Mc[t - 1];
		D_End2 += Mc[t];
		V_Start2 += Mc[t - 1];

		if (t <= SC_D - SC_WIN)
		{
			V_End2 += Mc[t + SC_WIN - 1];
		}

		if (t <= SC_L - SC_WIN)
		{
			Init_BEC_SW_Decoder(H, recv, dblk, V_End - Mv[t + SC_WIN - 1], V_End);
		}
		if (t <= SC_D - SC_WIN)
		{
			Init_BEC_SW_Decoder(H, recv, dblk, V_End2 - Mc[t + SC_WIN - 1], V_End2);
		}

		if (t < SC_L)
			test_BER(0, dblk, position_BER, t, D_Start2, D_End2);
		avg_iter += Iter_BEC_RA_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End, V_Start2, V_End2, D_Start2, D_End2);
		if (t < SC_L)
			test_BER(1, dblk, position_BER, t, D_Start2, D_End2);
	}




	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);

	avg_iter = floor(avg_iter / SC_L);
	return avg_iter;

}


int Run_BEC_SW_Decoder_Step
(
	mod2sparse* H,			/* Parity check matrix */
	int* recv,		/* Likelihood ratios for bits */
	char* dblk,		/* Place to store decoding */
	char* pchk,		/* Place to store parity checks */
	int* bIsCodeword,		/* decode 결과 codeword라고 판단되는지 */
	double** position_BER,
	int code_type,
	int SC_L,
	int SC_w,
	int SC_WIN,
	int* Mv,
	int* Mc,
	int etha
)
{
	int n, c;
	int sum = 0;
	int N = mod2sparse_cols(H);
	double avg_iter = 0;


	char* dblk_prev = (char*)calloc(N, sizeof(int));
	char* dblk2 = (char*)calloc(N, sizeof(int));
	int* recv2 = (int*)calloc(N, sizeof(int));
	for (int i = 0; i < N; i++)
		recv2[i] = ERASE_MARK;

	int SC_D;
	if (code_type == 0)
		SC_D = SC_L + SC_w - 1;
	else
		SC_D = SC_L + (SC_w - 1) / 2;

	int V_Start = 0;
	int V_End = 0;
	int V_End_prev;
	int C_Start = 0;
	int C_End = 0;
	for (int i = 0; i < SC_WIN; i++)
	{
		V_End += Mv[i];
		C_End += Mc[i];
	}

	int D_Start = 0;
	int D_End = 0;
	for (int p = 0; p < etha; p++)
		D_End += Mv[p];

	Init_BEC_SW_Decoder(H, recv, dblk, V_Start, V_End);
	for (int p = 0; p < etha; p++)
	{
		test_BER(0, dblk, position_BER, etha + p, D_Start, D_End);
	}
	avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
	for (int p = 0; p < etha; p++)
	{
		test_BER(1, dblk, position_BER, etha + p, D_Start, D_End);
	}

	V_End_prev = V_End;
	for (int t = 1; t < SC_L / etha; t++)
	{
		if (Mv[t] != 0)
		{
			for (int p = 0; p < etha; p++)
			{
				V_Start += Mv[p + t - 1];
				C_Start += Mc[p + t - 1];
			}
			for (int p = 0; p < etha; p++)
			{
				if (p + t * etha + SC_WIN - 1 >= SC_D)
				{
					V_End = N;
					C_End = M;
					break;
				}
				else
				{
					V_End += Mv[p + t * etha + SC_WIN - 1];
					C_End += Mc[p + t * etha + SC_WIN - 1];
				}
			}
			for (int p = 0; p < etha; p++)
			{
				D_Start += Mv[p + t * etha - 1];
				D_End += Mv[p + t * etha];
			}


			if (t * etha + SC_WIN <= SC_L + etha)
			{
				Init_BEC_SW_Decoder(H, recv, dblk, V_End_prev, V_End);
			}

			for (int p = 0; p < etha; p++)
			{
				test_BER(0, dblk, position_BER, t * etha + p, D_Start, D_End);
			}

			avg_iter += Iter_BEC_SW_Decoder(H, dblk, dblk_prev, V_Start, V_End, C_Start, C_End, D_Start, D_End);
			for (int p = 0; p < etha; p++)
			{
				test_BER(1, dblk, position_BER, t * etha + p, D_Start, D_End);
			}

			V_End_prev = V_End;
		}
		else
			break;
	}




	c = check(H, dblk, pchk);
	if (c == 0)
		* bIsCodeword = TRUE;

	free(dblk2);
	free(recv2);
	free(dblk_prev);

	avg_iter = floor((avg_iter / SC_L) * etha);
	return avg_iter;

}

void Init_BEC_SW_Decoder(mod2sparse* H, int* recv, char* dblk, int V_Start, int V_End)
{

	for (int j = V_Start; j < V_End; j++)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, j); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			e->msg_v_to_c = recv[j];
		}
		dblk[j] = recv[j];
	}// end for(j)

}


int Iter_BEC_SW_Decoder(mod2sparse* H, char* dblk, char* dblk_prev, int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End)
{
	int n;
	int exit = FALSE;


	for (n = 0; ; n++)
	{
		for (int i = V_Start; i < V_End; i++) dblk_prev[i] = dblk[i];

		for (int i = C_Start; i < C_End; i++)
		{
			Check_Update_BEC_SW(i, H);
		}

		for (int i = V_Start; i < V_End; i++)
		{
			Variable_Update_BEC_SW(i, C_Start, C_End, H, dblk);
		}

		//for(int i=V_Start;i<V_End;i++)
		//{
		//	Decision_BEC_SW(i,H,dblk);
		//}

		for (int i = D_Start; i < D_End; i++)
		{
			if (dblk[i] == ERASE_MARK)
				break;
			if (i == D_End - 1)
			{
				exit = TRUE;
			}
		}
		for (int i = V_Start; i < V_End; i++)
		{
			if (dblk[i] != dblk_prev[i])
				break;
			if (i == V_End - 1)
			{
				exit = TRUE;
			}
		}

		if ((n == max_iter) || exit) {
			break;
		}
	}

	return n;

}

int Iter_BEC_RA_SW_Decoder(mod2sparse* H, char* dblk, char* dblk_prev, int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End, int V_Start2, int V_End2, int D_Start2, int D_End2)
{
	int n;
	int exit = FALSE;


	for (n = 0; ; n++)
	{
		for (int i = V_Start; i < V_End; i++) dblk_prev[i] = dblk[i];
		for (int i = V_Start2; i < V_End2; i++) dblk_prev[i] = dblk[i];

		for (int i = C_Start; i < C_End; i++)
		{
			Check_Update_BEC_SW(i, H);
		}

		for (int i = V_Start; i < V_End; i++)
		{
			Variable_Update_BEC_SW(i, C_Start, C_End, H, dblk);
		}

		for (int i = V_Start2; i < V_End2; i++)
		{
			Variable_Update_BEC_SW(i, C_Start, C_End, H, dblk);
		}

		for (int i = D_Start; i < D_End; i++)
		{
			if (dblk[i] == ERASE_MARK)
				break;
			if (i == D_End - 1)
			{
				for (int j = D_Start2; j < D_End2; j++)
				{
					if (dblk[j] == ERASE_MARK)
						break;

					if (j == D_End2 - 1)
						exit = TRUE;
				}
			}
		}
		for (int i = V_Start; i < V_End; i++)
		{
			if (dblk[i] != dblk_prev[i])
				break;
			if (i == V_End - 1)
			{
				for (int j = V_Start2; j < V_End2; j++)
				{
					if (dblk[j] != dblk_prev[j])
						break;
					if (j == V_End2 - 1)
						exit = TRUE;
				}
			}
		}

		if ((n == max_iter) || exit) {
			break;
		}
	}

	return n;

}

int Iter_BEC_OC_Init_Decoder(mod2sparse* H, char* dblk, char* dblk_prev, int* V_Start, int* V_End, int* C_Start, int* C_End, int RV_Start, int RV_End, int etha)
{
	int n;
	int exit = FALSE;
	int Stopping_rule = FALSE;

	for (n = 0; ; n++)
	{
		for (int t = 0; t < etha; t++)
		{
			for (int i = V_Start[t]; i < V_End[t]; i++) dblk_prev[i] = dblk[i];
		}
		for (int i = RV_Start; i < RV_End; i++) dblk_prev[i] = dblk[i];

		for (int t = 0; t < etha; t++)
		{
			for (int i = C_Start[t]; i < C_End[t]; i++)
			{
				Check_Update_BEC_SW(i, H);
			}

			for (int i = V_Start[t]; i < V_End[t]; i++)
			{
				Variable_Update_BEC_SW(i, C_Start[t], C_End[t], H, dblk);
			}

			for (int i = RV_Start; i < RV_End; i++)
			{
				Variable_Update_BEC_SW(i, C_Start[t], C_End[t], H, dblk);
			}
		}

		for (int i = RV_Start; i < RV_End; i++)
		{
			if (dblk[i] == ERASE_MARK)
				break;
			if (i == RV_End - 1) exit = TRUE;
		}

		int exit_forloop = FALSE;
		for (int t = 0; t < etha; t++)
		{
			for (int j = V_Start[t]; j < V_End[t]; j++)
			{
				if (dblk[j] != dblk_prev[j])
				{
					exit_forloop = TRUE;
					break;
				}
			}
			if (exit_forloop)
				break;
		}

		if (!exit_forloop)
		{
			for (int i = RV_Start; i < RV_End; i++)
			{
				if (dblk[i] != dblk_prev[i])
					break;
				if (i == RV_End - 1)
					exit = TRUE;
			}
		}



		if ((n == max_iter) || exit) {
			break;
		}
	}


	return n;

}


void Check_Update_BEC_SW(int index_chk, mod2sparse* H)
{
	int temp;
	for (mod2entry* e = mod2sparse_first_in_row(H, index_chk); !mod2sparse_at_end(e); e = mod2sparse_next_in_row(e))
	{
		temp = 0;


		// e를 제외한 edge로부터 전달된 incoming msg를 보고 outgoing msg를 결정
		for (mod2entry* edge = mod2sparse_first_in_row(H, index_chk); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_row(edge))
		{
			if (e->col != edge->col)
			{

				if (edge->msg_v_to_c == ERASE_MARK)
				{
					temp = ERASE_MARK;
					break;
				}
				else
				{
					temp = (bool)temp + edge->msg_v_to_c;
				}
			}
		}//end for(edge)


		e->msg_c_to_v = temp;
	}//end for(e)

}

void Variable_Update_BEC_SW(int index_var, int index_chk_start, int index_chk_end, mod2sparse* H, char* dblk)
{

	if (dblk[index_var] == ERASE_MARK)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{
			if (mod2sparse_row(e) < index_chk_end && mod2sparse_row(e) >= index_chk_start && e->msg_c_to_v != ERASE_MARK)
			{
				dblk[index_var] = e->msg_c_to_v;
				for (mod2entry* edge = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
				{
					edge->msg_v_to_c = e->msg_c_to_v;
				}
				//// e를 제외한 edge로부터 전달된 incoming msg를 보고 outgoing msg를 결정
				//for(mod2entry * edge = mod2sparse_first_in_col(H,index_var); !mod2sparse_at_end(edge); edge = mod2sparse_next_in_col(edge))
				//{
				//	if(mod2sparse_row(edge) < index_chk_end && mod2sparse_row(edge)>=index_chk_start && e->row != edge->row && edge->msg_c_to_v !=ERASE_MARK)
				//	{
				//		e->msg_v_to_c = edge->msg_c_to_v;
				//		break;
				//	}
				//}//end for(edge)
			}
		}//end for(e)
	}

}



void Decision_BEC_SW(int index_var, mod2sparse* H, char* dblk)
{
	if (dblk[index_var] == ERASE_MARK)
	{
		for (mod2entry* e = mod2sparse_first_in_col(H, index_var); !mod2sparse_at_end(e); e = mod2sparse_next_in_col(e))
		{

			if (e->msg_c_to_v != ERASE_MARK)
			{
				dblk[index_var] = e->msg_c_to_v;
				break;
			}

		}//end for(e)
	}

}