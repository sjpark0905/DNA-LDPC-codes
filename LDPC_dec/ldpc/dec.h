#ifndef __DEC_H__
#define __DEC_H__

#define ERASE_MARK	2

/* DEC.H - Interface to decoding procedures. */

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


/* DECODING METHOD, ITS PARAMETERS, AND OTHER VARIABLES.  The global variables 
   declared here are located in dec.c. */

extern int max_iter;	// Maximum number of iteratons of decoding to do 
extern int max_iter2;	// Maximum number of iteratons of decoding to do 
extern int D_v;			// max degree of variable node edge
extern int D_c;			// max degree of check node edge
extern int bRegular_dv;	// dv가 regular이면 true
extern int bRegular_dc;	// dc가 regular이면 true

extern int		bSave_word_state;
extern int **	g_word_state;
extern int **	g_check_state;
extern int **	g_variable_state;
extern int *	g_unsat_check_num;
extern int *	g_err_num;
extern int		save_state_num;
extern int		tot_save_state_num;


void		Init_Decoder(int n, int m, int iter);
void		CheckRegular(mod2sparse * H);

//===================================================
// BEC Channel Decoder
//===================================================
int	Run_BEC_Decoder(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword,int *Target_VN);
void Init_BEC_Decoder(mod2sparse * H, int *recv, char * dblk);
void Iter_BEC_Decoder(mod2sparse * H, char * dblk);
void Check_Update_BEC(mod2sparse *H);
void Variable_Update_BEC(mod2sparse *H, char * dblk);
void Decision_BEC(mod2sparse * H, char * dblk);
int	Run_BEC_Decoder_SAVE(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER,int SC_L,int *Mv);
int	Run_BEC_Decoder_TARGET(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER,int SC_L,int *Mv,int *Target_VN);
//===================================================
// BEC Channel SW Decoder
//===================================================
int	Run_BEC_SW_Decoder(mod2sparse*	H, int *recv, char *dblk,char *pchk,int * bIsCodeword, int code_type, 	int SC_L,int SC_w, int SC_WIN, int *M, int *Mc);
int	Run_BEC_SW_Decoder_SAVE(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER, int code_type, 	int SC_L,int SC_w, int SC_WIN, int *Mv, int *Mc);
int	Run_BEC_SW_Decoder_RA(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER, int code_type, 	int SC_L,int SC_w, int SC_WIN, int *Mv, int *Mc);
int Run_BEC_SW_Decoder_OC(mod2sparse*	H,int *	recv,char *	dblk,char *	pchk,int * bIsCodeword,double ** position_BER,int SC_Ls,int SC_w, int SC_WIN, int Mv, int Mc,int etha);
int Iter_BEC_OC_Init_Decoder(mod2sparse * H, char * dblk, char *dblk_prev, int* V_Start, int* V_End, int* C_Start, int* C_End, int RV_Start, int RV_End, int etha);
int Run_BEC_SW_Decoder_Two(mod2sparse*	H,int *	recv,char *	dblk,char *	pchk,int * bIsCodeword,double ** position_BER,int SC_L,int SC_w, int SC_WIN, int *M, int *Mc);
int Run_BEC_SW_Decoder_Two_Cross(mod2sparse*	H,int *	recv,char *	dblk,char *	pchk,int * bIsCodeword,double ** position_BER,int SC_L,int SC_w, int SC_WIN, int *M, int *Mc);
int Run_BEC_SW_Decoder_Two_Indi(mod2sparse*	H,int *	recv,char *	dblk,char *	pchk,int * bIsCodeword,double ** position_BER,int SC_L,int SC_w, int SC_WIN, int *M, int *Mc);
int	Run_BEC_SW_Decoder_Step(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER, int code_type, int SC_L,int SC_w, int SC_WIN, int *Mv, int *Mc, int etha);
int	Run_BEC_SW_Decoder_TARGET(mod2sparse * H, int * recv, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER, int code_type, 	int SC_L,int SC_w, int SC_WIN, int *Mv, int *Mc);
void Init_BEC_SW_Decoder(mod2sparse * H, int * recv, char * dblk, int V_Start, int V_End);
int Iter_BEC_SW_Decoder(mod2sparse * H, char * dblk, char *dblk_prev,int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End);
int Iter_BEC_RA_SW_Decoder(mod2sparse * H, char * dblk, char *dblk_prev,int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End, int V_Start2, int V_End2, int D_Start2, int D_End2);
void Check_Update_BEC_SW(int index_chk, mod2sparse *H);
void Variable_Update_BEC_SW(int index_var, int index_chk_start, int index_chk_end, mod2sparse *H, char *dblk);
void Decision_BEC_SW(int index_var,mod2sparse * H, char * dblk);


//===================================================
// Sum Product Algorithm Decoder
//===================================================
int			Run_Belief_Propagation_Decoder(mod2sparse * H, double * lratio, char * dblk, char * pchk, int * bIsCodeword);
void		Init_Belief_Propagation(mod2sparse * H, double * lratio, char * dblk);
void		Iter_Belief_Propagation(mod2sparse * H, double * lratio, char * dblk);

int			Run_Belief_Propagation_Decoder_SAVE(mod2sparse * H, double * lratio, char * dblk, char * pchk, int * bIsCodeword, double ** position_BER,int SC_L, int SC_N);
void		test_BER(int iter, char * dblk, double ** position_BER, int pos, int V_start, int V_end);
//===================================================
// Gallager A, B Algorithm Decoder
//===================================================
int			Run_Gallager_Decoder(mod2sparse * H, int * recv, char * decoded, char * pchk, int type, int * bIsCodeword);
void		Init_Gallager(mod2sparse * H, int * recv, char * decoded);
void		Iter_Gallager(mod2sparse * H, int * recv, char * decoded, int type);
void		Check_Update_Gallager(mod2sparse * H, int * recv);
void		Variable_Update_Gallager(mod2sparse * H, int * recv,  int type);
void		Decision_Gallager(mod2sparse * H, int * recv, char * decoded, int type);

//===================================================
// Finite Alphabet Iterative Decoder
//===================================================

#define		N_LUT_2_SIZE	/*5*/ /*3*/ 2
#define		N_LUT_3_SIZE	3 /*1*/


#define		N_DIVERSITY		(N_LUT_2_SIZE + N_LUT_3_SIZE)
#define		N_WEIGHT_TYPE	/*3*/ 2


#define		STATE_IDX_MIN	0
#define		STATE_IDX_CUR	1


int			Run_Finite_Alphabet_Iterative_Decoder(mod2sparse * H, int * recv, char * decoded, char * pchk, int * bIsCodeword);
void		Init_FAID(mod2sparse * H, int * recv, char * decoded);
void		Iter_FAID(mod2sparse * H, int * recv, char * decoded);
void		Check_Update_FAID(mod2sparse * H);
void		Variable_Update_FAID(mod2sparse * H, int * recv);
void		Decision_FAID(mod2sparse * H, int * recv, char * decoded);
void		Set_FAID(int type);
void		Set_FAID_Weight(int type);
int			Variable_FAID_LUT(int m1, int m2, int y);


//===================================================
// Min Sum Algorithm Decoder 
//===================================================
int			Run_MSA_Decoder(mod2sparse * H, int * qLLR, int * L_Q,  char * decoded, char * pchk, int * bIsCodeword, FILE * fp1, FILE **  fp_w_v, int * w_v_idx);
int			Init_MSA(mod2sparse * H, int * qLLR, char * decoded);
int			Iter_MSA(mod2sparse * H, int * qLLR, int * L_Q,  char * decoded, int iter, char * pchk);
void		Check_Update_MSA(mod2sparse * H);
void		Variable_Update_MSA(mod2sparse * H, int * qLLR, int iter, char * pchk);
int			Decision_MSA(mod2sparse * H, int * qLLR, int * L_Q, char * decoded);

void		Set_MSA(int q, double step, int beta);
int			Cal_MSA_Q(double x, int type);
int			Cal_MSA_Clip(int x, int & b);

int			Run_MSA_Decoder_INF(mod2sparse * H, double * LLR, double * L,  char * decoded, char * pchk, int * bIsCodeword);
void		Init_MSA_INF(mod2sparse * H, double * LLR, char * decoded);
void		Iter_MSA_INF(mod2sparse * H, double * LLR, double * L,  char * decoded);
void		Check_Update_MSA_INF(mod2sparse * H);
void		Variable_Update_MSA_INF(mod2sparse * H, double * LLR);
void		Decision_MSA_INF(mod2sparse * H, double * LLR, double * L, char * decoded);


void		Save_State(int n, int c, char * decoded, char * check, int & min_c);
void		Print_Variable_State(mod2sparse * H, int n, int c, char * decoded, int * qLLR, int * L_Q, char * pchk, FILE ** fp, int * w_v_idx, FILE * fp1);

int Set_Message_Maxvalue(int x);
int Dec_Message_value(int x);


#endif


//===================================================
// Pipeline Decoder
//===================================================
int			Run_Pipeline_Decoder(int b, int c, int m, int L, mod2sparse * H, double * lratio, char * dblk, char * pchk, int * bIsCodeword);
void		Init_Pipeline(mod2sparse * H, double * lratio, char * dblk);
void		Pipeline(mod2sparse * H, double * lratio, char * dblk);
void		Decoding_window(int time, mod2sparse * H, double * lratio, char * dblk);
void		chk_process(int ind_chk, mod2sparse * H, double * Iratio, char * dblk);
void		var_process(int iter, int ind_var, mod2sparse * H, double * Iratio, char * dblk);



//===================================================
// Sliding window Decoder
//===================================================
int Run_SW_Decoder(mod2sparse*	H, double *	lratio,char * dblk, char *	pchk, int * bIsCodeword, double ** position_BER, int code_type, 	int SC_L,int SC_w, int SC_WIN, int *Mv, int *Mc);
int Run_SW_Decoder_SKU(mod2sparse*	H, double *	lratio,char * dblk, char *	pchk, int * bIsCodeword, double ** position_BER, int code_type, 	int SC_L,int SC_w, int SC_WIN, int *Mv, int *Mc);
int Run_SWW_Decoder(mod2sparse*	H, double *	lratio,char * dblk, char *	pchk, int * bIsCodeword, double ** position_BER, int code_type, int SC_b, int SC_c, int SC_L, int SC_w, int SC_M, int SC_WIN, int SC_Cr);
void Init_SW_Decoder(mod2sparse * H, double * lratio, int V_Start, int V_End);
int Iter_SW_Decoder(mod2sparse * H, char * dblk, char * pchk, double * lratio, int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End, int V_Check_End, int C_Check_End);
void Iter_SWW_Decoder(mod2sparse * H, char * dblk,  double * lratio, int V_Start, int V_End, int C_Start, int C_End, int D_Start, int D_End);
void Check_Update_SW(int index_chk, mod2sparse * H,double * lratio);
void Variable_Update_SW(int index_var, int index_chk_start, int index_chk_end, mod2sparse * H,double * lratio);
void Variable_Update_SWW(int index_var, int index_chk_start, int index_chk_end, mod2sparse * H,double * lratio, char * dblk,double alpha);
void Fixing_Message_SW( mod2sparse * H, double * lratio, char * dblk, int VP, int SC_c, int SC_M);
void Decision_SW(int index_var, int index_chk_start, int index_chk_end, mod2sparse * H,double * lratio, char * dblk);
void Decision_SWW(int index_var, int index_chk_start, int index_chk_end, mod2sparse * H,double * lratio, char * dblk, double alpha);

//===================================================
// Product Code
//===================================================
int			Run_PD_Decoder(int init,int DSW, int repeat, int &avg_check, int prod_M, int prod_N, int gap, int verk, mod2sparse * H,mod2sparse* for_H,
	mod2sparse* rev_H, double * lratio, char * dblk, char * pchk, int * bIsCodeword);
void		Init_PD(mod2sparse * H, double * lratio, char * dblk);
void		PD_window(int rev, int num, mod2sparse * H, mod2sparse * chk_H, double * lratio, char * dblk);
void		PD_chk_process(int ind_chk, mod2sparse * H, double * Iratio, char * dblk);
void		PD_var_process(int firstchk, int lastchk, int indic_dec, int ind_var, mod2sparse * H, double * Iratio, char * dblk);
void		PD_Initsymbols(int num, int rev, int decision, mod2sparse * H, double *lratio, char *dblk);
void		PD_decision(int number,int ind_var, mod2sparse * H, double *lratio,char *dblk );
int			decision_check(int num, int reverse, mod2sparse * chk_H, char *dblk);
//===================================================
// Product Code(POST process decoding)
//===================================================

int			Run_POST_PD_Decoder(int init, int DSW,int repeat, int prod_M, int prod_N, int gap, int verk, mod2sparse * H, mod2sparse* H1, double * lratio, char * dblk, char * pchk, int * bIsCodeword);
int			PD_window1(int num, mod2sparse * H, mod2sparse * H1, double * lratio, char * dblk) ;
void		PD_H1_var_process(int firstchk, int lastchk, int ind_var, mod2sparse * H, double * Iratio, char * dblk);
void		POST_PD_window(int rev, int num, mod2sparse * H, mod2sparse *H1, double * lratio, char * dblk);
void		POST_PD_Init(int num, mod2sparse * H, char *dblk, double *lratio);
void		H1_check(int num, mod2sparse * H1, char *dblk);


