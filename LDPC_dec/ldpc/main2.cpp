
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



//====================================
// MPICH2
//====================================
//#include "mpi.h"

//====================================
// MKL
//====================================
#include "rand.h"
#include "channel.h"

//====================================
// LDPC source by Radford M. Neal
//====================================
#include "alloc.h"
#include "blockio.h"
#include "check.h"
#include "dec.h"
#include "intio.h"
#include "mod2sparse.h"
#include "open.h"
#include "rcode.h"
#include "common.h"


#define ROOT_ID				0
#define ERR_FILE_LEN		316
#define ERR_TYPE_MAX		3		//수정

#define ERASE_MARK	2				//BEC channel

enum DECODER_TYPE
{
	DECODER_BP					= 0,	// belief propagation
	DECODER_GALLAGER_A			= 1,	// Gallager algorithm A
	DECODER_GALLAGER_B_1		= 2,	// Gallager algorithm B_1
	DECODER_GALLAGER_B_2		= 3,	// Gallager algorithm B_2

	DECODER_BP_POST1			= 10,

	DECODER_MSA					= 20,
	DECODER_MSA_OFFSET			= 21,
	DECODER_MSA_QUASI_UNIFORM	= 22,
	DECODER_PIPELINE			= 50,
	DECODER_SW					= 60,
	DECODER_SWW					= 61,
	DECODER_SW_SKU				= 62,
	DECODER_BEC					= 80,
	DECODER_BEC_SAVE			= 81,
	DECODER_BEC_TARGET			= 82,
	DECODER_BEC_SW				= 90,
	DECODER_BEC_SW_SAVE			= 91,
	DECODER_BEC_SW_TWO			= 92,
	DECODER_BEC_SW_OC			= 93,
	DECODER_BEC_SW_TWO_CROSS	= 94,
	DECODER_BEC_SW_TWO_INDI		= 95,
	DECODER_BEC_SW_STEP			= 96,
	DECODER_BEC_SW_RA			= 97,
	DECODER_BEC_SW_OPTION		= 98
};

enum CHANNEL_TYPE
{
	CHANNEL_AWGN	= 0,		// AWGN	(soft decision)
	CHANNEL_BSC		= 1,		// BSC	(hard decision)
	CHANNEL_BEC		= 2
};

enum RAND_TYPE
{
	RAND_SEED			= 0,	// seed 값을 주어서 한 번 초기화
	RAND_LOAD_FILE_ONCE	= 1,	// file을 읽어서 한 번 초기화
	RAND_LOAD_FILE_ALL	= 2		// error를 모은 file을 읽어서 계속 초기화
};

enum SAVE_TYPE
{
	SAVE_NONE	= 0,		// rand state를 저장하지 않음.
	SAVE_ERROR	= 1,		// frame error가 발생할때 마다 rand state를 저장.
	SAVE_LAST	= 2,
	SAVE_POSITION_BER = 10
};

enum ERR_TYPE
{
	ERR_TYPE_RAW = 0,
	ERR_TYPE_1,
	ERR_TYPE_FINAL = ERR_TYPE_MAX - 1
};

enum CODE_TYPE
{
	Two_side = 0,
	One_side = 1,
};

void	SetUp(int argc, char *argv[]);
void	Init();
void	Init_Variable();
void	Set_Code();
void	Set_Precision();
void	Set_FrameNum();
void	Alloc_Mem();
void	Set_Random_Data();
int		Check_End();
void	Run_Simulation();
void	Close_File();
void	Print_All_Result();
void	Print_One_Result();
//void	COLLECT_MPI();



int		Load_Rand_File(char * fileNameTemp, char * fileNameTot, FILE * fp);
void	Save_Rand_File(char * fileNameTemp, char * fileNameTot, FILE * fp);
void	Print_word_state(long long index, FILE * fp);
void	Print_supp(FILE * fp, int * decoded);

void	LDPC_Set_Decoder();
void	LDPC_Set_Message();
void	LDPC_Encode();
int		LDPC_Channel();
int		LDPC_Decode(int * b);
int		LDPC_BIT_Check();
int		LDPC_Codeword_Check(int raw);
int		LDPC_Raw_Error_Check();
int		LDPC_Decode_2ND_Stage(int * bIsCodeword);


void	Read_Wrong_Variable_Node(FILE * fp);

int		LDPC_Channel_proposed(int * bEnd);

//===============================
// GLOBAL VARIABLE
//===============================
int		g_CODE_N;
int		g_CODE_M;
int		g_CODE_K;

int		g_spatial;
int		g_product;




int		g_SC_CODE_b;
int		g_SC_CODE_c;
int		g_SC_CODE_cr;
int		g_SC_CODE_L;
int		g_SC_CODE_w;
int		g_SC_CODE_z;	
int		g_SC_CODE_WIN=0;
int		g_SC_CODE_D;
int		g_SC_CODE_etha;
int		*g_SC_CODE_M;
int		*g_SC_CODE_Mc;
int		*g_SC_CODE_CN_ORDER;



int		g_punctuation;
int		g_punctuation_start;
int		g_punctuation_end;
int		g_punctuation_start1;
int		g_punctuation_end1;
int		g_punctuation_start2;
int		g_punctuation_end2;
int		*g_punctuation_set;

int		g_shortening;
int		g_shortening_start;
int		g_shortening_end;

int		g_shortening_etha;
int		g_shortening_num;

int		g_targeting;
int		g_target_VN[2];

double	g_CODE_RATE;

char	FileName_Temp2[255];
char	*g_message;				// length : CODE_K
char	*g_codeword;			// length : CODE_N
double  *g_recv_codeword_soft;	// length : COND_N
int		*g_recv_codeword_hard;	// lenght : CODE_N
double	*g_received_LR;			// length : CODE_N
double  *g_received_LLR;		// length : CODE_N
int		*g_received_qLLR;		// length : CODE_N
int		*g_L_Q;					// lenght : CODE_N
double  *g_L;					// length : CODE_N
char	*g_bit_stream_trans;	// length : CODE_N
char	*g_chks;				// length : CODE_M
double	*g_bprb;				// length : CODE_N

double	**g_position_BER;
double  **g_position_total_BER;
double	g_p;

DECODER_TYPE	g_decoder_type;
CHANNEL_TYPE	g_channel_type;	
RAND_TYPE		g_rand_type;
SAVE_TYPE		g_save_type;
CODE_TYPE		g_code_type;

int			g_maxIteration;		// decoding에서 max iteration
double		g_EbNo;				// Eb/No
double		g_std_dev;			// AWGN	sigma value
double		g_eps;
int			g_bSystematic;		// systematic code이면 1, 아닐 경우 0
int			g_initial_seed;		// 입력받은 seed value
long		g_frame_num;		// 테스트 할 총 frame 수
long		g_target_frame_err; // target frame error
int			g_precision;
double		g_step_length;
int			g_one_side_level;
int			g_total_level;

long long	g_bit_err[ERR_TYPE_MAX];			// number of bit error for one process
long long	g_frame_err[ERR_TYPE_MAX];			// number of frame error for one process
long long	g_tot_bit_err[ERR_TYPE_MAX];		// number of total bit error 
long long	g_tot_frame_err[ERR_TYPE_MAX];		// number of total frame error

long long	g_undetected_frame_err[ERR_TYPE_MAX];
long long	g_tot_undetected_frame_err[ERR_TYPE_MAX];

double		g_BER[ERR_TYPE_MAX];
double		g_FER[ERR_TYPE_MAX];


//시뮬레이션 시간 계산 
time_t		g_t_start;
time_t		g_t_end;


// 시뮬레이션할 frame 수
long long	g_my_frame_num;						// 하나의 process가 시행할 총 frame 수
long long	g_frame_i[ERR_TYPE_MAX];			// 현재 진행 중인 frame 번호
long long	g_tot_frame_i[ERR_TYPE_MAX];
long long	g_total_iter;
long long	g_total_chk_update;
int   g_avg_chk_update;

int			g_myid		= 0;
int			g_numprocs	= 1;


// FILE 입출력에 사용
FILE	*fp_load_rand;
FILE	*fp_save_total;
FILE	*fp_save_temp;
FILE	*fp_result;
FILE	*fp_word_state;
FILE	*fp_wrong_check;
FILE	*fp_wrong_variable;

char	FileName_pchk[255];				// 읽을 pchk 파일 이름 (.pchk)
char	FileName_input[255];
char	FileName_input2[255];
char	FileName_dec[255];
char	FileName_pchk1[255];
char	FileName_M[255];				// SC-LDPC에서 position별 M값 (.txt)
char	FileName_rand_load_temp[255];	// 읽을 rand 임시 파일 이름 (.err)
char	FileName_rand_load[255];		// 읽을 rand 파일 이름 (.err)
char	FileName_save_temp[255];		// rand.cpp에서 save_state를 저장할 임시 파일 이름 (.err)
char	FileName_save_total[255];		// rand.cpp에서 save_state를 저장할 전체 파일 이름 (.err)
char	FileName_save_last[255];		// rand.cpp에서 state를 저장할 파일 이름(.err)

char	FileName_result[255];			
char	FileName_word_state[255];
char	FileName_wrong_check[255];

char	FileName_wrong_variable[255];

int*	g_wrong_variable_idx;	//확인해야 할 variable node이면 1, 아닐 경우 0
int		g_wrong_variable_num;	//확인해야 할 variable node의 수

FILE	**g_fp_W_V;

// 전보환 추가 부분
double g_total_rate;


// SC-LDPC_WIN_TEST
double g_eps_p;


//===============================================================================
// 입력 받은 parameter를 읽어서 알맞게 저장한다
//===============================================================================
void SetUp(int argc, char *argv[])
{
	int pos = 1;
	g_bSystematic	= atoi(argv[pos++]);
	g_decoder_type	= (DECODER_TYPE)	atoi(argv[pos++]);
	g_channel_type	= (CHANNEL_TYPE)	atoi(argv[pos++]);
	g_initial_seed	= atoi(argv[pos++]);
	g_maxIteration	= atoi(argv[pos++]);
	g_frame_num		= atol(argv[pos++]);
	if(g_frame_num == 0) g_target_frame_err	= atoi(argv[pos++]);
	
	sprintf(FileName_Temp2, "%s", argv[pos++]);
	sprintf(FileName_input, "%s.txt", FileName_Temp2);

	char FileName_Temp3[255];

	sprintf(FileName_Temp3, "%s", argv[pos++]);
	sprintf(FileName_input2, "%s.txt", FileName_Temp3);

	char FileName_Temp[255];

	sprintf(FileName_Temp, "%s", argv[pos++]);
	sprintf(FileName_pchk, "%s.pchk", FileName_Temp);


	if(g_channel_type==CHANNEL_BEC) g_eps = atof(argv[pos++]);
	else if (g_channel_type == CHANNEL_BSC) g_p = atof(argv[pos++]);
	else g_EbNo	= atof(argv[pos++]);	// target SNR

	g_punctuation	= atoi(argv[pos++]);
	g_shortening	= atoi(argv[pos++]);
	g_targeting		= atoi(argv[pos++]);

	// Punctuation
	if(g_punctuation>0 && g_punctuation!=5)
	{
		g_punctuation_start = atoi(argv[pos++]);
		g_punctuation_end	= atoi(argv[pos++]);
	}

	if(g_shortening==1)
	{
		g_shortening_start = atoi(argv[pos++]);
		g_shortening_end  = atoi(argv[pos++]);
	}

	if(g_shortening==2)
	{
		g_shortening_etha = atoi(argv[pos++]);
		g_shortening_num  = atoi(argv[pos++]);
	}

	if(g_targeting)
	{
		g_target_VN[0] = atoi(argv[pos++]);
		g_target_VN[1] = atoi(argv[pos++]);
	}

	if(g_punctuation==2)
	{
		g_SC_CODE_w			= atoi(argv[pos++]);
		g_SC_CODE_L			= atoi(argv[pos++]);

		g_SC_CODE_D = g_SC_CODE_L + g_SC_CODE_w - 1;
		FILE* FP_temp;
		sprintf(FileName_M,"%s.txt",FileName_Temp);
		FP_temp = fopen(FileName_M, "rb");
		g_SC_CODE_M = (int*) calloc(g_SC_CODE_D, sizeof(int));
		g_SC_CODE_Mc = (int*) calloc(g_SC_CODE_D, sizeof(int));
		for(int i=0;i<g_SC_CODE_D;i++)
		{
			fscanf(FP_temp,"%d\n",&g_SC_CODE_M[i]);
		}
		for(int i=0;i<g_SC_CODE_D;i++)
		{
			fscanf(FP_temp,"%d\n",&g_SC_CODE_Mc[i]);
		}

		fclose(FP_temp);
	}

	if(g_punctuation==3)
	{
		g_SC_CODE_w			= atoi(argv[pos++]);
		g_SC_CODE_L			= atoi(argv[pos++]);

		g_SC_CODE_D = g_SC_CODE_L + g_SC_CODE_w - 1;
		FILE* FP_temp;
		sprintf(FileName_M,"%s.txt",FileName_Temp);
		FP_temp = fopen(FileName_M, "rb");
		g_SC_CODE_M = (int*) calloc(g_SC_CODE_D, sizeof(int));
		for(int i=0;i<g_SC_CODE_D;i++)
		{
			fscanf(FP_temp,"%d\n",&g_SC_CODE_M[i]);
		}

		fclose(FP_temp);
		g_punctuation_start2 = atoi(argv[pos++]);
		g_punctuation_end2	= atoi(argv[pos++]);
	}

	if(g_punctuation==4)
	{
		g_punctuation_start1 = atoi(argv[pos++]);
		g_punctuation_end1	= atoi(argv[pos++]);
		g_SC_CODE_w			= atoi(argv[pos++]);
		g_SC_CODE_L			= atoi(argv[pos++]);

		g_SC_CODE_D = g_SC_CODE_L + g_SC_CODE_w - 1;
		FILE* FP_temp;
		sprintf(FileName_M,"%s.txt",FileName_Temp);
		FP_temp = fopen(FileName_M, "rb");
		g_SC_CODE_M = (int*) calloc(g_SC_CODE_D, sizeof(int));
		for(int i=0;i<g_SC_CODE_D;i++)
		{
			fscanf(FP_temp,"%d\n",&g_SC_CODE_M[i]);
		}

		fclose(FP_temp);
		g_punctuation_start2 = atoi(argv[pos++]);
		g_punctuation_end2	= atoi(argv[pos++]);
	}

	// Sliding window
	if(g_decoder_type == DECODER_SW  || g_decoder_type == DECODER_SW_SKU ||g_decoder_type == DECODER_BEC_SW || g_decoder_type == DECODER_SWW || g_decoder_type == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SAVE || g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_STEP || g_decoder_type == DECODER_BEC_SW_RA || g_decoder_type == DECODER_BEC_SW_OPTION)
	{
		g_code_type = (CODE_TYPE) atoi(argv[pos++]);
		g_SC_CODE_w			= atoi(argv[pos++]);
		g_SC_CODE_L			= atoi(argv[pos++]);
		g_SC_CODE_WIN	= atoi(argv[pos++]);
		sprintf(FileName_M,"%s.txt",FileName_Temp);


		if(g_code_type == Two_side)
			g_SC_CODE_D = g_SC_CODE_L + g_SC_CODE_w - 1;
		else if(g_code_type == One_side)
			g_SC_CODE_D = g_SC_CODE_L + (g_SC_CODE_w - 1)/2;

		FILE* FP_temp;
		FP_temp = fopen(FileName_M, "rb");
		g_SC_CODE_M = (int*) calloc(g_SC_CODE_D, sizeof(int));
		g_SC_CODE_Mc = (int*) calloc(g_SC_CODE_D, sizeof(int));

		if(g_decoder_type == DECODER_BEC_SW_OC)
		{	
			g_SC_CODE_etha = atoi(argv[pos++]);
			g_SC_CODE_L = g_SC_CODE_L*g_SC_CODE_etha+1;
			g_SC_CODE_M[0] = atoi(argv[pos++]);
			g_SC_CODE_Mc[0] = atoi(argv[pos++]);
		}
		else
		{
			if(g_decoder_type == DECODER_BEC_SW_STEP)
			{
				g_SC_CODE_etha = atoi(argv[pos++]);
			}
			for(int i=0;i<g_SC_CODE_D;i++)
			{
				fscanf(FP_temp,"%d\n",&g_SC_CODE_M[i]);
			}
			for(int i=0;i<g_SC_CODE_D;i++)
			{
				fscanf(FP_temp,"%d\n",&g_SC_CODE_Mc[i]);
			}
			fclose(FP_temp);
		}
	}

	if(g_punctuation==5)
	{
		g_punctuation_set = (int*) calloc(g_SC_CODE_L,sizeof(int));
		FILE* FP_temp;
		sprintf(FileName_M,"%s_punct.txt",FileName_Temp);
		FP_temp = fopen(FileName_M, "rb");
		for(int i=0;i<g_SC_CODE_L;i++)
			fscanf(FP_temp,"%d\n",&g_punctuation_set[i]);

		fclose(FP_temp);
	}

	if(g_decoder_type == DECODER_BEC_SW_OPTION)
	{
		g_SC_CODE_CN_ORDER = (int*) calloc(g_SC_CODE_D,sizeof(int));
		FILE* FP_temp;
		sprintf(FileName_M,"%s_order.txt",FileName_Temp);
		FP_temp = fopen(FileName_M, "rb");
		for(int i=0;i<g_SC_CODE_L;i++)
			fscanf(FP_temp,"%d\n",&g_punctuation_set[i]);

		fclose(FP_temp);
	}



	if(argc != pos) 
	{
		if(g_myid == ROOT_ID)
		{
			fprintf(stderr,"\n\nargc error!\n\n");
		}
		//MPI_Finalize();
		exit(1);
	}


}

//===============================================================================
// 시뮬레이션에 필요한 값을 초기화 한다
//===============================================================================
void Init()
{
	Set_Code();
	Init_Variable();
	Set_Random_Data();
}

//===============================================================================
// 사용되는 변수르 초기화 한다
//===============================================================================
void Init_Variable()
{
	Set_Precision();
	Alloc_Mem();
	Set_FrameNum();

	for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
	{
		g_bit_err[i]					= 0;
		g_frame_err[i]					= 0;
		g_tot_bit_err[i]				= 0;
		g_tot_frame_err[i]				= 0;
		g_undetected_frame_err[i]		= 0;
		g_tot_undetected_frame_err[i]	= 0;
		g_frame_i[i]					= 0;
		g_tot_frame_i[i]				= 0;
	}

	g_total_iter = 0;
}

//===============================================================================
// simulation할 LDPC code에 관한 변수를 초기화한다
//===============================================================================
void Set_Code()
{
	read_pchk(FileName_pchk);

	g_CODE_M		= M;
	g_CODE_N		= N;
	g_CODE_K		= g_CODE_N - g_CODE_M;

	printf("\n");
	printf("g_CODE_N : %d\n", g_CODE_N);
	printf("g_CODE_K : %d\n", g_CODE_K);
	printf("g_CODE_M : %d\n", g_CODE_M);
	printf("\n");

	if(g_targeting==0)
	{
		g_target_VN[0] = 1;
		g_target_VN[1] = N;
	}


	if(g_rand_type == RAND_LOAD_FILE_ALL)
	{
		g_CODE_RATE		= g_total_rate; // 변경
	}
	else if(g_save_type == SAVE_POSITION_BER)
	{
		if(g_total_rate != 0 ) g_CODE_RATE		= g_total_rate; //변경(2014-05-12)
		else g_CODE_RATE		= 1.0 - (double) ((double)M/(double)N);
	}
	else if(g_punctuation==1)
	{
		g_CODE_RATE		= 1.0 - (double)(((double)M-(double)(g_punctuation_end-g_punctuation_start+1))/((double)N-(double)(g_punctuation_end-g_punctuation_start+1)));
	}
	else if(g_punctuation==2)
	{
		int num_punct = (g_punctuation_end-g_punctuation_start+1)*g_SC_CODE_L;
		g_CODE_RATE		= 1.0 - (double)((double)(M-num_punct)/(N-num_punct));
	}
	else if(g_punctuation==3)
	{
		int num_punct = (g_punctuation_end-g_punctuation_start+1)*g_SC_CODE_L+(g_punctuation_end2-g_punctuation_start2+1);
		g_CODE_RATE		= 1.0 - (double)((double)(M-num_punct)/(N-num_punct));
	}
	else if(g_punctuation==4)
	{
		int num_punct = (g_punctuation_end-g_punctuation_start+1)*g_SC_CODE_L+(g_punctuation_end2-g_punctuation_start2+1)+(g_punctuation_end1-g_punctuation_start1+1)*(g_SC_CODE_w-1);
		g_CODE_RATE		= 1.0 - (double)((double)(M-num_punct)/(N-num_punct));
	}
	else if(g_shortening==1)
	{
		int num_shortening = (g_shortening_end-g_shortening_start+1);
		g_CODE_RATE		= 1.0 - (double)((double)(M)/(N-num_shortening));
	}
	else if(g_shortening==2)
	{
		int num_shortening = g_shortening_num * floor((double) g_SC_CODE_L/(double)g_shortening_etha);
		g_CODE_RATE		= 1.0 - (double)((double)(M)/(N-num_shortening));
	}
	else 
	{
		g_CODE_RATE		= 1.0 - (double) ((double)M/(double)N);
	}

	g_std_dev		= getStd_dev(g_EbNo, g_CODE_RATE);
}

//===============================================================================
// message passing에 finite bit을 사용할 경우 알맞게 정한다
//===============================================================================
void Set_Precision()
{
	g_one_side_level	= 0;

	// 유한한 bit를 사용할 경우
	if(g_precision > 0)
	{
		g_total_level		= (int) pow(2.0, g_precision) - 1;
		g_one_side_level	= (int) pow(2.0, (g_precision - 1)) - 1;
	}
}

//===============================================================================
// simulation할 frame 수를 정한다
//===============================================================================
void Set_FrameNum()
{
	long long	temp1;
	long long	temp2;

	if(g_frame_num == 0)
	{	
		temp1 = 100 % g_numprocs;
		temp2 = 100 / g_numprocs;

		if(g_myid < g_numprocs - 1)		g_my_frame_num = temp2;
		else							g_my_frame_num = temp2 + temp1;
	}

	else
	{
		temp1 = g_frame_num % g_numprocs;
		temp2 = g_frame_num / g_numprocs;

		if(g_myid < g_numprocs - 1)		g_my_frame_num = temp2;
		else							g_my_frame_num = temp2 + temp1;
	}
}

//===============================================================================
// simulation에 사용되는 array를 만든다
//===============================================================================
void Alloc_Mem()
{
	// memory allocation for variables	  
	g_message				= (char*)	calloc(g_CODE_K, sizeof(char));
	g_codeword				= (char*)	calloc(g_CODE_N, sizeof(char));
	g_received_LR			= (double*)	calloc(g_CODE_N, sizeof(double));
	g_received_LLR			= (double*) calloc(g_CODE_N, sizeof(double));
	g_received_qLLR			= (int*)	calloc(g_CODE_N, sizeof(int));
	g_bit_stream_trans		= (char *)	calloc(g_CODE_N,sizeof(char)); //2016.12.08 곽희열 수정, Calloc만 사용하면 Iteration이 돌아가지 않는 변수 노드 또한 옳은 비트로 계산을 하기 때문에 수정함
	for(int i=0;i<g_CODE_N;i++)
		g_bit_stream_trans[i] = 2;
	//	g_bit_stream_trans		= (char*)	calloc(g_CODE_N, sizeof(char));
	g_chks					= (char*)	calloc(g_CODE_M, sizeof(char));
	g_recv_codeword_soft	= (double*) calloc(g_CODE_N, sizeof(double));
	g_recv_codeword_hard	= (int*)	calloc(g_CODE_N, sizeof(int));
	g_bprb					= (double*) calloc(g_CODE_N, sizeof(double));
	g_L_Q					= (int*)	calloc(g_CODE_N, sizeof(int));
	g_L						= (double*) calloc(g_CODE_N, sizeof(double));
	g_wrong_variable_idx	= (int *)	calloc(g_CODE_N, sizeof(int));
	//g_p						= (double*)calloc(1, sizeof(double));
	if(g_decoder_type == DECODER_SW || g_decoder_type == DECODER_SW_SKU || g_decoder_type == DECODER_SWW || g_decoder_type == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SAVE || g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_STEP || g_decoder_type == DECODER_BEC_SW_RA) {
		g_position_BER			= (double **) calloc(200, sizeof(double*));
		for(int i=0; i<200; i++)
		{
			g_position_BER[i] = (double *) calloc(g_SC_CODE_L+1, sizeof(double));
		}

		g_position_total_BER			= (double **) calloc(200, sizeof(double*));
		for(int i=0; i<200; i++)
		{
			g_position_total_BER[i] = (double *) calloc(g_SC_CODE_L+1, sizeof(double));
		}
	}
}

//===============================================================================
// random을 어떤 방식으로 초기화 할지 결정한다
//===============================================================================
void Set_Random_Data()
{
	int seed;
	bSave_word_state = FALSE;

	if (g_rand_type == RAND_SEED)
	{					
		seed = g_myid + g_initial_seed;
		rand_seed(seed);
	}
	else if(g_rand_type == RAND_LOAD_FILE_ONCE)
	{
		if(rand_load(FileName_rand_load))
		{
			fprintf(stderr,"can't open rand file - (%s)\n",FileName_rand_load);
			//MPI_Finalize();
			exit(1);
		}
	}
	else if(g_rand_type == RAND_LOAD_FILE_ALL)
	{
		// error를 한 파일에 다 모은 방식
		fp_load_rand = fopen(FileName_rand_load, "rb");
		if(fp_load_rand == NULL)
		{
			fprintf(stderr,"can't open rand file - (%s)\n", FileName_rand_load);
			//MPI_Finalize();
			exit(1);
		}

		//bSave_word_state = TRUE;
	}

	if(g_save_type == SAVE_ERROR)
	{
		fp_save_total = fopen(FileName_save_total, "wb");
		if(fp_save_total == NULL)
		{
			fprintf(stderr,"can't open rand file - (%s)\n", FileName_save_total);
			//MPI_Finalize();
			exit(1);
		}
	}

	if(bSave_word_state == TRUE)
	{
		fp_word_state		= fopen(FileName_word_state, "w");
		fp_wrong_check		= fopen(FileName_wrong_check, "w");
		fp_wrong_variable	= fopen(FileName_wrong_variable, "r");

		g_fp_W_V = (FILE**) calloc(g_CODE_N, sizeof(FILE*));
	}
	else
	{
		fp_word_state		= NULL;
		fp_wrong_check		= NULL;
	}
}

//===============================================================================
// simulation이 끝났는지를 체크한다
//===============================================================================
int Check_End()
{
	int v = FALSE;
	int state;
	char file_name[255];
	int point[15] = {2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
	int i=0;

	if(g_rand_type == RAND_LOAD_FILE_ALL)
	{
		state = Load_Rand_File(FileName_rand_load_temp, FileName_rand_load,fp_load_rand);
		if(state == FALSE)
		{
			// 다 읽었으니 종료
			v = TRUE;
		}
		else
		{
			// 읽을 부분이 있기 때문에 진행
			// 방금 만든 파일이 FileName_rand_load_temp 이름으로 있기 때문에 이를 load
			state = rand_load(FileName_rand_load_temp);
			g_frame_i[ERR_TYPE_RAW]++;
		}
	}
	else if(g_frame_num == 0) //Goal: Target frame error
	{
		if(g_frame_err[ERR_TYPE_1]>=(g_target_frame_err/g_numprocs)) v = TRUE;
		else g_frame_i[ERR_TYPE_RAW]++;


	}
	else
	{
		if(g_frame_i[ERR_TYPE_RAW] < g_my_frame_num)		g_frame_i[ERR_TYPE_RAW]++;
		else												v = TRUE;
	}


	return v;
}

//===============================================================================
// simulation을 한다
//===============================================================================
void Run_Simulation()
{
	int bEnd = FALSE;
	int state;
	int bIsCodeword;	//valid codeword인지 여부 (misscorrection 포함됨)
	int temp;
	int old;

	while(1)
	{
		if(Check_End() == TRUE)
			break;


		//=====================================
		// SET MESSAGE
		//=====================================
		//LDPC_Set_Message();

		//=====================================
		// ENCODE
		//=====================================
		LDPC_Encode();

		//=====================================
		// CHANNEL
		//=====================================
		temp = LDPC_Channel();


		g_bit_err[ERR_TYPE_RAW] += temp;
		if(temp > 0)
		{
			g_frame_err[ERR_TYPE_RAW]++;
		}

		//=====================================
		// DECODE
		//=====================================
		g_frame_i[ERR_TYPE_1]++;

		if(bSave_word_state == TRUE)
		{
			fprintf(fp_wrong_check, "\n<%d>\n", g_frame_i[ERR_TYPE_1]);
		}

		temp = LDPC_Decode(&bIsCodeword);

		//=====================================
		// POSITION_BER
		//=====================================
		if(g_decoder_type == DECODER_SW || g_decoder_type == DECODER_SW_SKU || g_decoder_type == DECODER_SWW || g_decoder_type == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SAVE || g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_STEP || g_decoder_type == DECODER_BEC_SW_RA)
		{
			for(int i=0; i<200;i++)
			{
				for(int j=0; j<g_SC_CODE_L; j++)
				{
					g_position_total_BER[i][j] += g_position_BER[i][j];
				}
			}
		}


		g_total_iter += temp;
		g_total_chk_update += g_avg_chk_update;
		temp = LDPC_BIT_Check();
		if(temp > 0)
		{
			g_bit_err[ERR_TYPE_1] += temp;
			g_frame_err[ERR_TYPE_1]++;

			if(bIsCodeword == TRUE)
			{
				g_undetected_frame_err[ERR_TYPE_1]++;
			}
		}

		if(bSave_word_state == TRUE)
		{
			if(bIsCodeword == FALSE)
			{
				Print_word_state(g_frame_i[ERR_TYPE_RAW], fp_word_state);
			}
		}

		//=====================================
		// CHECK
		//=====================================
		temp = LDPC_BIT_Check();
		if(temp > 0)
		{
			g_bit_err[ERR_TYPE_FINAL] += temp;
			g_frame_err[ERR_TYPE_FINAL]++;

			if(bIsCodeword == TRUE)
				g_undetected_frame_err[ERR_TYPE_FINAL]++;

			if(g_save_type == SAVE_ERROR)
			{
				Save_Rand_File(FileName_save_temp, FileName_save_total, fp_save_total);
			}
		}

		if(g_save_type == SAVE_ERROR)
		{
			state = rand_del(1);
		}

		if(g_rand_type == RAND_LOAD_FILE_ALL)
		{
			state = rand_del(0);
		}

		
	}

	sprintf(FileName_dec, "dec_%s.txt", FileName_Temp2);

	printf("%s", FileName_dec);
	FILE* f_i;				// 2020.07.29 박성준 수정 : decoder input 저장
	f_i = fopen(FileName_dec, "w");

	for (int i = 0;i < g_CODE_N;i++)
	{
		fprintf(f_i, "%d ", g_bit_stream_trans[i]);
	}

	fclose(f_i);


}

//===============================================================================
//
//===============================================================================
void Close_File()
{
	if(g_save_type == SAVE_LAST)
	{
		rand_save(FileName_save_last);
	}
	else if(g_save_type == SAVE_ERROR)
	{
		remove(FileName_save_temp);
		fclose(fp_save_total);
	}

	if(g_rand_type == RAND_LOAD_FILE_ALL)
	{
		fclose(fp_load_rand);
		remove(FileName_rand_load_temp);
	}

	if(bSave_word_state == TRUE)
	{
		if(fp_word_state != NULL)
			fclose(fp_word_state);	
		if(fp_wrong_check != NULL)
			fclose(fp_wrong_check);
	}
}

//===============================================================================
//
//===============================================================================
void Print_All_Result()
{
	double avg_iter;
	double avg_check_update;
	double frame_num;
	int len;

	if(g_myid == ROOT_ID)
	{
		if(g_rand_type == RAND_LOAD_FILE_ALL)
		{

			for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
			{
				g_BER[i] = (double) g_tot_bit_err[i]	/ frame_num;
				g_FER[i] = (double) g_tot_frame_err[i]	/ (double)g_frame_i[ERR_TYPE_RAW];//test
			}
			printf("%d\n", g_frame_i[ERR_TYPE_RAW]);
			avg_iter	= ((double) g_total_iter)/((double) g_frame_i[ERR_TYPE_RAW]);	
			avg_check_update = ((double) g_total_chk_update)/((double) g_frame_i[ERR_TYPE_RAW]);
		}
		else if(g_decoder_type == DECODER_SW || g_decoder_type == DECODER_SW_SKU || g_decoder_type == DECODER_SWW || g_decoder_type == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SAVE || g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_STEP || g_decoder_type == DECODER_BEC_SW_RA)
		{
			if(g_bSystematic == TRUE)	len = g_CODE_K;
			else						len = (g_target_VN[1]-g_target_VN[0]+1);

			frame_num = (double)len * g_tot_frame_i[ERR_TYPE_RAW];
			for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
			{
				g_BER[i] = (double) g_tot_bit_err[i]	/ frame_num;
				g_FER[i] = (double) g_tot_frame_err[i]	/ (double) g_tot_frame_i[ERR_TYPE_RAW];
			}

			avg_iter	= ((double) g_total_iter)/((double) g_frame_i[ERR_TYPE_RAW]);
			avg_check_update = ((double) g_total_chk_update)/((double) g_frame_i[ERR_TYPE_RAW]);

			for(int i = 0; i<200; i++) 
			{
				for(int j = 0; j<g_SC_CODE_L; j++)
				{
					g_position_total_BER[i][j] = g_position_total_BER[i][j] / g_frame_i[ERR_TYPE_RAW];
				}
			}

		}
		else
		{
			if(g_bSystematic == TRUE)	len = g_CODE_K;
			else						len = (g_target_VN[1]-g_target_VN[0]+1);

			frame_num = (double)len * g_frame_i[ERR_TYPE_RAW];
			for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
			{
				g_BER[i] = (double) g_bit_err[i]	/ frame_num;
				g_FER[i] = (double) g_frame_err[i]	/ (double)g_frame_i[ERR_TYPE_RAW];//teset
			}

			avg_iter	= ((double) g_total_iter)/((double) g_frame_i[ERR_TYPE_RAW]);
			avg_check_update = ((double) g_total_chk_update)/((double) g_frame_i[ERR_TYPE_RAW]);
		}
	}

	//=====================================
	// PRINT RESULT
	//=====================================
	if(g_myid == ROOT_ID)
	{
		int tm_hour;
		int tm_min;
		int tm_sec;
		double d;

		

		if(g_decoder_type  == DECODER_BEC_SW || g_decoder_type  == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_RA)
			sprintf(FileName_result, "result_(%s)_%s_%d_%.3f_%d_%d_%d_%d_%d.txt", FileName_input2, FileName_pchk, g_decoder_type, g_eps, g_tot_frame_i[0], g_maxIteration,g_initial_seed,g_SC_CODE_WIN,g_shortening);
		else if(g_decoder_type == DECODER_BEC_SW_STEP)
			sprintf(FileName_result, "result_(%s)_%s_%d_%.3f_%d_%d_%d_%d_%d.txt", FileName_input2, FileName_pchk, g_decoder_type, g_eps, g_tot_frame_i[0], g_maxIteration,g_initial_seed,g_SC_CODE_WIN,g_SC_CODE_etha);
		else if(g_channel_type == CHANNEL_BEC)
			sprintf(FileName_result, "result_(%s)_%s_%d_%.3f_%d_%d_%d.txt", FileName_input2, FileName_pchk, g_decoder_type, g_eps, g_tot_frame_i[0], g_maxIteration,g_initial_seed);
		else if (g_channel_type == CHANNEL_BSC)
			sprintf(FileName_result, "result_(%s)_%s_%d_%.4f_%d_%d_%d.txt", FileName_input2, FileName_pchk, g_decoder_type, g_p, g_tot_frame_i[0], g_maxIteration, g_initial_seed);
		else
			sprintf(FileName_result, "result_(%s)_%s_%d_%.3fdB_%d_%d_%d.txt", FileName_input2, FileName_pchk, g_decoder_type, g_EbNo, g_tot_frame_i[0], g_maxIteration,g_initial_seed);
		fp_result = fopen(FileName_result, "w");

		fprintf(fp_result, "code N        : %d\n", g_CODE_N);
		fprintf(fp_result, "code K        : %d\n", g_CODE_K);
		fprintf(fp_result, "code M        : %d\n", g_CODE_M);
		fprintf(fp_result, "code rate     : %.3f\n", g_CODE_RATE);
		if(g_channel_type == CHANNEL_BEC)
			fprintf(fp_result, "Eps		: %.3f \n", g_eps);
		else
		{
			fprintf(fp_result, "Eb/No         : %.2f dB\n", g_EbNo);
			fprintf(fp_result, "g_std_dev     : %.2f\n", g_std_dev);
		}
		fprintf(fp_result, "decode        : %d\n", g_decoder_type);
		fprintf(fp_result, "systematic    : %d\n", g_bSystematic);
		fprintf(fp_result, "max iteration : %d\n", g_maxIteration);
		fprintf(fp_result, "dv            : %d\n", D_v);
		fprintf(fp_result, "bRegular_dv   : %d\n", bRegular_dv);
		fprintf(fp_result, "dc            : %d\n", D_c);
		fprintf(fp_result, "bRegular_dc   : %d\n", bRegular_dc);
		if(g_decoder_type == DECODER_SW || g_decoder_type == DECODER_SW_SKU ||g_decoder_type == DECODER_SWW|| g_decoder_type  == DECODER_BEC_SW || g_decoder_type  == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_STEP || g_decoder_type == DECODER_BEC_SW_RA)
		{
			fprintf(fp_result, "w   : %d\n\n", g_SC_CODE_w);
			fprintf(fp_result, "L   : %d\n\n", g_SC_CODE_L);
			fprintf(fp_result, "W   : %d\n\n", g_SC_CODE_WIN);
		}
		if(g_targeting)
		{
			fprintf(fp_result, "target_VN:%d~%d \n\n", g_target_VN[0],g_target_VN[1]);
		}
		if(g_punctuation)
		{
			
			fprintf(fp_result, "[Type: %d] Punctuation_VN:%d~%d \n\n",g_punctuation,g_punctuation_start,g_punctuation_end);
		}
		if(g_shortening==1)
		{
			fprintf(fp_result, "[Type: %d] Shortening_VN:%d~%d \n\n",g_shortening,g_shortening_start,g_shortening_end);
		}



		fprintf(fp_result, "=============================================\n");
		fprintf(fp_result, "                 result\n");
		fprintf(fp_result, "=============================================\n");
		d	= difftime(g_t_end, g_t_start);
		tm_hour	= (int) (d / (60 * 60));
		d   = d - (tm_hour * 60 * 60);
		tm_min	= (int) (d / 60);
		d	= d - (tm_min * 60);
		tm_sec	= (int) d;

		fprintf(fp_result, "start time      : %s", ctime(&g_t_start));
		fprintf(fp_result, "end time        : %s", ctime(&g_t_end));
		fprintf(fp_result, "simulation time : %d hours %d mins %d secs\n\n", tm_hour, tm_min, tm_sec);

		fprintf(fp_result, "# of processes         : %d\n",	g_numprocs);
		fprintf(fp_result, "initial seed value     : %d\n\n",	g_initial_seed);

		for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
		{
			fprintf(fp_result, "# of Frame[%2d]          :%lld\n", i, g_frame_i[i]);
		}
		fprintf(fp_result, "\n");


		for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
		{
			fprintf(fp_result, "# of Bit Errors[%2d]     : %lld\n", i, g_bit_err[i]);
		}
		fprintf(fp_result, "\n");

		for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
		{
			fprintf(fp_result, "BER[%2d]                 : %.5e\n",i, g_BER[i]);
		}
		fprintf(fp_result, "\n");

		for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
		{
			fprintf(fp_result, "# of Frame Errors[%2d]   : %lld\n", i, g_frame_err[i]);
		}
		fprintf(fp_result, "\n");

		for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
		{
			fprintf(fp_result, "FER[%2d]                 : %.5e\n",i, g_FER[i]);
		}
		fprintf(fp_result, "\n");

		for(int i = 0 ; i < ERR_TYPE_MAX ; i++)
		{
			fprintf(fp_result, "# of undetected Frame[%2d]   : %lld\n", i, g_undetected_frame_err[i]);
		}
		fprintf(fp_result, "\n");
		fprintf(fp_result, "p				: %f\n", g_p);
		fprintf(fp_result, "maxIteration			: %d\n", g_maxIteration);
		fprintf(fp_result, "average iter			: %.5f\n", avg_iter);
		fprintf(fp_result, "average check update		: %.5f\n", avg_check_update);
		for (int i = 1; i < 30; i++)
		{
			fprintf(fp_result, "%d ", g_codeword[i]);
		}
		
		fclose(fp_result);

		//=====================================
		// SAVE POSITION BER
		//=====================================

		if(g_decoder_type == DECODER_SW || g_decoder_type == DECODER_SW_SKU || g_decoder_type == DECODER_SWW || g_decoder_type == DECODER_BEC_SW_SAVE || g_decoder_type == DECODER_BEC_SAVE|| g_decoder_type == DECODER_BEC_SW_OC || g_decoder_type == DECODER_BEC_SW_TWO || g_decoder_type == DECODER_BEC_SW_TWO_CROSS || g_decoder_type == DECODER_BEC_SW_TWO_INDI || g_decoder_type == DECODER_BEC_SW_STEP || g_decoder_type == DECODER_BEC_SW_RA) 
		{
			FILE * fp_SAVE;

			if(g_decoder_type == DECODER_SW || g_decoder_type == DECODER_SW_SKU || g_decoder_type == DECODER_SWW)
				sprintf(FileName_result, "POSITION_BER_%s_%d_%.3f_%d_%d_%d.txt", FileName_pchk,g_decoder_type, g_EbNo, g_frame_num, g_maxIteration,g_SC_CODE_WIN);
			else
				sprintf(FileName_result, "POSITION_BER_%s_%d_%.3f_%d_%d_%d.txt", FileName_pchk,g_decoder_type, g_eps, g_frame_num, g_maxIteration,g_SC_CODE_WIN);
		
			
			fp_SAVE = fopen(FileName_result, "w");
			for(int i=0; i<200 ; i++)
			{
				for(int j=0; j<g_SC_CODE_L; j++)
				{
					fprintf(fp_SAVE, "%.5e ", g_position_total_BER[i][j]);
				}
				fprintf(fp_SAVE, "\n");
			}

			fclose(fp_SAVE);
			for (int i=0;i<200;i++)
			{
				free(g_position_BER[i]);
				free(g_position_total_BER[i]);
			}
			free(g_position_BER);
			free(g_position_total_BER);
		}



	}
}

//===============================================================================
//
//===============================================================================
void Print_One_Result()
{
	printf("\n");
	if(g_channel_type == CHANNEL_BEC)
		printf("[%d]  code : (%d,%d)\trate : %.3f\tEps : %.2f \n",g_myid, g_CODE_N, g_CODE_K, g_CODE_RATE, g_eps);
	else
		printf("[%d]  code : (%d,%d)\trate : %.3f\tEb/No : %.2f dB\n",g_myid, g_CODE_N, g_CODE_K, g_CODE_RATE, g_EbNo);
	printf("[%d]  frame_num              : %d\n", g_myid, g_frame_i[ERR_TYPE_RAW]);
	printf("[%d]  bit_err                : %d\n", g_myid, g_bit_err[ERR_TYPE_FINAL]);
	printf("[%d]  frame_err              : %d\n", g_myid, g_frame_err[ERR_TYPE_FINAL]);
	printf("[%d]  without coding (bit)   : %lld\n", g_myid, g_bit_err[ERR_TYPE_RAW]);
	printf("[%d]  without coding (frame) : %lld\n\n", g_myid, g_frame_err[ERR_TYPE_RAW]);
}

//===============================================================================
// MPI에 사용된 변수를 합친다
//===============================================================================
/*void COLLECT_MPI()
{
	MPI_Reduce(g_bit_err,				g_tot_bit_err,					ERR_TYPE_MAX, MPI_LONG_LONG, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce(g_frame_err,				g_tot_frame_err,				ERR_TYPE_MAX, MPI_LONG_LONG, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce(g_frame_i,				g_tot_frame_i,					ERR_TYPE_MAX, MPI_LONG_LONG, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
	MPI_Reduce(g_undetected_frame_err,	g_tot_undetected_frame_err,		ERR_TYPE_MAX, MPI_LONG_LONG, MPI_SUM, ROOT_ID, MPI_COMM_WORLD);
}*/

//===============================================================================
// MAIN
//===============================================================================
int	main(int argc, char *argv[])  //argc : 받는 인수의 개수, argc[pos] : pos 별로 받는 순서
{


	/*MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &g_numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &g_myid);*/

	SetUp(argc, argv);
	Init();
	time(&g_t_start);
	LDPC_Set_Decoder();
	Run_Simulation();
	time(&g_t_end);
	Close_File();
	Print_One_Result();
	//COLLECT_MPI(); 
	Print_All_Result();

	//MPI_Finalize();
	//sprintf(title, "decoded.txt");
	//FILE* A = fopen(title, "wb");
	//for (i = 0; i < g_CODE_N; i++)
	//{

	//	sprintf(A, "%d ", g_bit_stream_trans[i]);

	//}
	//fclose(A);



	return 0;
}

//===============================================================================
// return value 
//			TRUE : 정상적으로 파일을 읽어서 새로운 임시 파일을 생성
//			FALSE : 더 이상 읽을 파일이 없어서 임시 파일을 생성 못함
//===============================================================================
int	Load_Rand_File(char * fileNameTemp,char * fileNameTot, FILE * fp)
{
	char buf[ERR_FILE_LEN];
	int len;
	FILE * fp_temp = NULL;
	if(fread(buf, sizeof(char), ERR_FILE_LEN, fp))
	{
		fp_temp = fopen(fileNameTemp, "wb");
		len = fwrite(buf, sizeof(char), ERR_FILE_LEN, fp_temp);
		fclose(fp_temp);
		fp_temp = NULL;
		return TRUE;
	}
	return FALSE;
}

//===============================================================================
//
//===============================================================================
void Save_Rand_File(char * fileNameTemp, char * fileNameTot, FILE * fp)
{
	int result;
	int len;
	char buf[ERR_FILE_LEN];
	FILE * src_fp = NULL;

	result = rand_save_err(fileNameTemp);
	if(result == 0)
	{
		src_fp = fopen(fileNameTemp, "rb");

		len = fread(buf, sizeof(char), ERR_FILE_LEN, src_fp);
		len = fwrite(buf, sizeof(char), ERR_FILE_LEN, fp);

		// 정리
		fclose(src_fp);
		src_fp = NULL;
	}
}

//===============================================================================
//
//===============================================================================
void LDPC_Set_Decoder()
{
	max_iter = g_maxIteration;
	max_iter2 = 20;

	CheckRegular(H);

	if(bSave_word_state == TRUE)
	{
		Init_Decoder(g_CODE_N, g_CODE_M,  max_iter);
	}

	if(g_decoder_type == DECODER_MSA || g_decoder_type == DECODER_MSA_QUASI_UNIFORM)
	{
		Set_MSA(g_precision, g_step_length, 0);
	}
	else if(g_decoder_type == DECODER_MSA_OFFSET)
	{
		Set_MSA(g_precision, g_step_length, 1);
	}
}

//===============================================================================
// 전송할 message 선택
// Assume all zero codeword
//===============================================================================
void LDPC_Set_Message()
{
	for(int i = 0 ; i < g_CODE_K ; i++)
	{
		g_message[i] = 0;
	}
}

//===============================================================================
// encode
// Assume all zero codeword
//=============================================================================== 
void LDPC_Encode()
{

	FILE* f;
	f = fopen(FileName_input, "rb");

	for (int i = 0; i < g_CODE_N; i++)
	{
		fscanf(f, "%d", &g_codeword[i]);
	}


	fclose(f);


	FILE* f2;
	f2 = fopen(FileName_input2, "rt");

	int count0 = 0;
	int count1 = 0;

	for (int i = 0; i < g_CODE_N; i++)
	{
		fscanf(f2, "%lf", &g_received_LLR[i]);
		g_recv_codeword_soft[i] = g_received_LLR[i];
		g_received_LR[i] = exp(g_received_LLR[i]);
	}


}

//===============================================================================
// channel을 통과시킨다
//===============================================================================
int	LDPC_Channel()
{
	int cnt = 0;
	int idx;
	int i;
	/*
	if(g_channel_type == CHANNEL_AWGN)
	{
		channel_AWGN(g_codeword, g_CODE_N, g_recv_codeword_soft, g_received_LR, g_received_LLR, g_std_dev);		//received_LLR, LR 값이 출력됨
	}
	else if (g_channel_type == CHANNEL_BSC)
	{
		channel_BSC(g_codeword, g_CODE_N, g_recv_codeword_hard, g_received_LR, g_received_LLR, g_std_dev, g_p);
	}
	else if(g_channel_type == CHANNEL_BEC)
	{
		channel_BEC(g_codeword, g_CODE_N, g_recv_codeword_hard, g_eps);
	}
	*/


	cnt = LDPC_Raw_Error_Check();
	if(g_shortening==1)
	{
		if(g_channel_type == CHANNEL_AWGN)
		{
			for(int i=g_shortening_start;i<=g_shortening_end;i++)
				g_received_LR[i-1] = exp(30.0);
		}
		else if(g_channel_type == CHANNEL_BEC)
		{
			for(int i=g_shortening_start;i<=g_shortening_end;i++)
				g_recv_codeword_hard[i-1] = g_codeword[i-1];
		}
	}
	else if(g_shortening==2)
	{
		if(g_channel_type == CHANNEL_AWGN)
		{
			for(int j=0;j<g_SC_CODE_L;j++)
			{
				if((j+1)%g_shortening_etha==0)
				{
					for(int i=0;i<g_shortening_num;i++)
					{
						int idx = j*g_SC_CODE_c*g_SC_CODE_z+i;
						g_received_LR[idx-1] = exp(30.0);
					}
				}
			}
		}
		else if(g_channel_type == CHANNEL_BEC)
		{
			for(int j=0;j<g_SC_CODE_L;j++)
			{
				if((j+1)%g_shortening_etha==0)
				{
					for(int i=0;i<g_shortening_num;i++)
					{
						int idx = j*g_SC_CODE_c*g_SC_CODE_z+i;
						g_recv_codeword_hard[idx-1] = g_codeword[idx-1];
					}
				}
			}
		}
	}

	if(g_punctuation)
	{
		if((g_channel_type == CHANNEL_AWGN) || (g_channel_type == CHANNEL_BSC))
		{
			if(g_punctuation==1)
			{
				for(int i=g_punctuation_start;i<=g_punctuation_end;i++)
					g_received_LR[i-1] =1;
			}
			else if(g_punctuation==2)
			{
				int idx = 0;
				for(int j=0;j<g_SC_CODE_L;j++)
				{
					for(int i=g_punctuation_start;i<=g_punctuation_end;i++)
					{
						g_received_LR[idx+i-1] =1;
					}
					idx += g_SC_CODE_M[j];
				}
			}
		}
		else if(g_channel_type == CHANNEL_BEC)
		{
			if(g_punctuation==1)
			{
				for(int i=g_punctuation_start;i<=g_punctuation_end;i++)
				{
					g_recv_codeword_hard[i-1] = ERASE_MARK;
				}
			}
			else if(g_punctuation==2)
			{
				int idx = 0;
				for(int j=0;j<g_SC_CODE_L;j++)
				{
					for(int i=g_punctuation_start;i<=g_punctuation_end;i++)
					{
						g_recv_codeword_hard[idx+i-1] = ERASE_MARK;
					}
					idx += g_SC_CODE_M[j];
				}
			}
			else if(g_punctuation==3)
			{
				int idx = 0;
				for(int j=0;j<g_SC_CODE_L;j++)
				{
					for(int i=g_punctuation_start;i<=g_punctuation_end;i++)
					{
						g_recv_codeword_hard[idx+i-1] = ERASE_MARK;
					}
					idx += g_SC_CODE_M[j];
				}
				for(int i=g_punctuation_start2;i<=g_punctuation_end2;i++)
				{
					g_recv_codeword_hard[i-1] = ERASE_MARK;
				}
			}
			else if(g_punctuation==4)
			{
				int idx = 0;
				for(int j=0;j<g_SC_CODE_L;j++)
				{
					for(int i=g_punctuation_start;i<=g_punctuation_end;i++)
					{
						g_recv_codeword_hard[idx+i-1] = ERASE_MARK;
					}
					idx += g_SC_CODE_M[j];
				}
				idx = g_SC_CODE_M[0]*g_SC_CODE_L;
				for(int j=0;j<g_SC_CODE_w-1;j++)
				{
					for(int i=g_punctuation_start1;i<=g_punctuation_end1;i++)
					{
						g_recv_codeword_hard[idx+i-1] = ERASE_MARK;
					}
					idx += g_SC_CODE_M[j];
				}
				for(int i=g_punctuation_start2;i<=g_punctuation_end2;i++)
				{
					g_recv_codeword_hard[i-1] = ERASE_MARK;
				}
			}
			else if(g_punctuation==5)
			{
				int idx_start = 0;
				int idx_end = g_SC_CODE_M[0]-1;
				for(int t=0;t<g_SC_CODE_L;t++)
				{
					int step = floor((double)(idx_end-idx_start+1)/g_punctuation_set[t]);

					if(g_punctuation_set[t]<0)
					{
						for(int i=0;i<=g_punctuation_set[t];i++)
							g_recv_codeword_hard[idx_start+i*step] = g_codeword[idx_start+i*step];
					}
					else if(g_punctuation_set[t]>0)
					{
						for(int i=0;i<=g_punctuation_set[t];i++)
							g_recv_codeword_hard[idx_start+i*step] = ERASE_MARK;
					}
					idx_start += g_SC_CODE_M[t];
					if(t<g_SC_CODE_L-1)
					{
						idx_end +=g_SC_CODE_M[t+1];
					}
				}
			}

		}
	}

	//===================================
	// Quantization
	//===================================
	if((g_decoder_type == DECODER_MSA) || 
		(g_decoder_type == DECODER_MSA_OFFSET))
	{
		// 유한한 bit를 사용하는 경우
		if(g_precision > 0)
		{
			for(int i = 0 ; i < g_CODE_N ; i++)
			{
				g_received_qLLR[i]	= Cal_MSA_Q(g_received_LLR[i], 0);
				g_L_Q[i]			= g_received_qLLR[i];
			}
		}
	}
	else if(g_decoder_type == DECODER_MSA_QUASI_UNIFORM)
	{
		for(int i = 0 ; i < g_CODE_N ; i++)
		{
			g_received_qLLR[i]	= Cal_MSA_Q(g_received_LLR[i], 1);
		}
	}

	return cnt;
}



//===============================================================================
//
//===============================================================================
int	LDPC_Decode(int * bIsCodeword)
{
	int iter = 0;
	if(g_save_type ==	SAVE_POSITION_BER && g_decoder_type == DECODER_BP)
	{
		iter = Run_Belief_Propagation_Decoder_SAVE(H, g_received_LR, g_bit_stream_trans, g_chks, bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_z);
	}
	else if(g_decoder_type == DECODER_BP)
	{
		iter = Run_Belief_Propagation_Decoder(H, g_received_LR, g_bit_stream_trans, g_chks, bIsCodeword);
	}
	else if(g_decoder_type == DECODER_GALLAGER_A)
	{
		iter = Run_Gallager_Decoder(H, g_recv_codeword_hard, g_bit_stream_trans, g_chks, 0, bIsCodeword);
	}
	else if(g_decoder_type == DECODER_GALLAGER_B_1)
	{
		iter = Run_Gallager_Decoder(H, g_recv_codeword_hard, g_bit_stream_trans, g_chks, 1, bIsCodeword);
	}
	else if(g_decoder_type == DECODER_GALLAGER_B_2)
	{
		iter = Run_Gallager_Decoder(H, g_recv_codeword_hard, g_bit_stream_trans, g_chks, 2, bIsCodeword);
	}
	else if((g_decoder_type == DECODER_MSA) || (g_decoder_type == DECODER_MSA_OFFSET) || (g_decoder_type == DECODER_MSA_QUASI_UNIFORM))
	{
		if(g_precision > 0)
			iter = Run_MSA_Decoder(H, g_received_qLLR, g_L_Q, g_bit_stream_trans, g_chks, bIsCodeword, fp_wrong_check, g_fp_W_V, g_wrong_variable_idx);
		else
			iter = Run_MSA_Decoder_INF(H, g_received_LLR, g_L, g_bit_stream_trans, g_chks, bIsCodeword);
	}
	else if(g_decoder_type == DECODER_PIPELINE)
	{
		iter = Run_Pipeline_Decoder(g_SC_CODE_b, g_SC_CODE_c,g_SC_CODE_w, g_SC_CODE_L, H, g_received_LR, g_bit_stream_trans, g_chks, bIsCodeword);
	}
	else if(g_decoder_type == DECODER_SW)
	{
		iter = Run_SW_Decoder(H,g_received_LR,g_bit_stream_trans,g_chks, bIsCodeword, g_position_BER, g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	else if(g_decoder_type == DECODER_SWW)
	{
		iter = Run_SWW_Decoder(H,g_received_LR,g_bit_stream_trans,g_chks, bIsCodeword, g_position_BER, g_code_type, g_SC_CODE_b ,g_SC_CODE_c ,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_z,g_SC_CODE_WIN,g_SC_CODE_cr);
	}
	else if(g_decoder_type == DECODER_SW_SKU)
	{
//		iter = Run_SW_Decoder_SKU(H,g_received_LR,g_bit_stream_trans,g_chks, bIsCodeword, g_position_BER, g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	else if(g_decoder_type == DECODER_BEC)
	{
		iter = Run_BEC_Decoder(H, g_recv_codeword_hard, g_bit_stream_trans, g_chks, bIsCodeword,g_target_VN);
	}
	else if(g_decoder_type == DECODER_BEC_SAVE)
	{
		iter = Run_BEC_Decoder_SAVE(H, g_recv_codeword_hard, g_bit_stream_trans, g_chks, bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_M);
	}
	else if(g_decoder_type == DECODER_BEC_TARGET)
	{
		iter = Run_BEC_Decoder_TARGET(H, g_recv_codeword_hard, g_bit_stream_trans, g_chks, bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_M,g_target_VN);
	}
	else if(g_decoder_type == DECODER_BEC_SW_SAVE)
	{
		if(g_targeting)
			iter = Run_BEC_SW_Decoder_TARGET(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks, bIsCodeword,g_position_BER, g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
		else
			iter = Run_BEC_SW_Decoder_SAVE(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks, bIsCodeword,g_position_BER, g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);

	}
	else if(g_decoder_type == DECODER_BEC_SW)
	{
		iter = Run_BEC_SW_Decoder(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks, bIsCodeword, g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	else if(g_decoder_type == DECODER_BEC_SW_OC)
	{
		iter = Run_BEC_SW_Decoder_OC(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks,bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M[0],g_SC_CODE_Mc[0],g_SC_CODE_etha);
	}
	else if(g_decoder_type == DECODER_BEC_SW_TWO)
	{
		iter = Run_BEC_SW_Decoder_Two(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks,bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	else if(g_decoder_type == DECODER_BEC_SW_TWO_CROSS)
	{
		iter = Run_BEC_SW_Decoder_Two_Cross(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks,bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	else if(g_decoder_type == DECODER_BEC_SW_TWO_INDI)
	{
		iter = Run_BEC_SW_Decoder_Two_Indi(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks,bIsCodeword,g_position_BER,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	else if(g_decoder_type == DECODER_BEC_SW_STEP)
	{
		iter = Run_BEC_SW_Decoder_Step(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks,bIsCodeword,g_position_BER,g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc,g_SC_CODE_etha);
	}
	else if(g_decoder_type == DECODER_BEC_SW_RA)
	{
		iter = Run_BEC_SW_Decoder_RA(H,g_recv_codeword_hard,g_bit_stream_trans,g_chks,bIsCodeword,g_position_BER,g_code_type,g_SC_CODE_L,g_SC_CODE_w,g_SC_CODE_WIN,g_SC_CODE_M,g_SC_CODE_Mc);
	}
	return iter;
}

//===============================================================================
//
//===============================================================================
int LDPC_Decode_2ND_Stage(int * bIsCodeword)
{
	int iter = 0;

	return iter;
}

//===============================================================================
// decoding 결과 bit error를 확인한다
//===============================================================================
int	LDPC_BIT_Check()
{


	int idx;
	int cnt = 0;
	int len;
	int offset;

	if(g_bSystematic == TRUE)		
	{
		len = g_CODE_K;
		for(int i = 0 ; i < len ; i++)
		{
			if(g_codeword[i] != g_bit_stream_trans[i])
				cnt++;

		}
	}
	else
	{
		for(int i = g_target_VN[0]-1 ; i < g_target_VN[1] ; i++)
		{
			if(g_codeword[i] != g_bit_stream_trans[i])
				cnt++;
		}
	}



	return cnt;
}

//===============================================================================
// channel을 통과한 codeword에 bit error가 몇개가 생겼는지를 체크
//===============================================================================
int	LDPC_Raw_Error_Check()
{
	int cnt = 0;
	int temp;
	int len;
	int offset;




	if(g_bSystematic == TRUE)	len = g_CODE_K;
	else						len = g_CODE_N;

	for(int i = 0 ; i < len ; i++)
	{
		if(g_channel_type == CHANNEL_AWGN)	// soft decision
		{
			if(g_recv_codeword_soft[i] >= 0)		temp = 0;
			else									temp = 1;
			if(g_codeword[i] != temp)				cnt++;
		}
		else if(g_channel_type == CHANNEL_BEC)
		{
			if(g_recv_codeword_hard[i] == ERASE_MARK)	 cnt++;
		}
		else	// hard decision
		{
			if(g_recv_codeword_hard[i] >= 0)		temp = 0;
			else									temp = 1;
			if(g_codeword[i]  != temp)				cnt++;
		}
	}




	return cnt;
}

//===============================================================================
//
//===============================================================================
int	LDPC_Codeword_Check(int raw)
{
	int r = FALSE;
	static char *	temp = NULL;
	static int		len  = 0;
	if(raw == TRUE)	//decode 이전
	{
		if(temp == NULL)
		{
			len = g_CODE_N;
			temp = (char*) calloc(len, sizeof(char));
		}

		if(g_channel_type == CHANNEL_AWGN)
		{
			for(int i = 0 ; i < g_CODE_N ; i++)
			{
				if(g_received_LLR[i] >= 0)		temp[i] = 0;
				else							temp[i] = 1;
			}
		}
		else
		{
			for(int i = 0 ; i < g_CODE_N ; i++)
			{
				if(g_recv_codeword_hard[i] >= 0)	temp[i] = 0;
				else								temp[i] = 1;
			}
		}

		if(check(H, temp, g_chks) == 0)
			r = TRUE;
	}
	else		//decode 후
	{
		if(check(H, g_bit_stream_trans, g_chks) == 0)	
			r = TRUE;
	}

	return r;
}


//===============================================================================
// 
//===============================================================================
void Print_word_state(long long index, FILE * fp)
{
	int idx_min;
	int idx_init;
	int offset;
	int N = mod2sparse_cols(H);
	int M = mod2sparse_rows(H);

	idx_min		= save_state_num;

	fprintf(fp, "< %lld >\n", index);
	fprintf(fp, "(iter\t:\t# error\t# UC\t/\terror loc)\n");
	for(int i = 0; i < tot_save_state_num - 1 ; i++)
	{
		fprintf(fp, "%d\t:\t%d \t%d\t/\t", i, g_err_num[i], g_unsat_check_num[i]);

		for(int j = 0 ; j < N ; j++)
		{
			if(g_word_state[i][j] != 0)
				fprintf(fp, "%d ", j);
		}
		fprintf(fp, "\t/\t");
		for(int j = 0 ; j < M ; j++)
		{
			if(g_check_state[i][j] != 0)
				fprintf(fp, "%d ", j);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");
}

//===============================================================================
// 
//===============================================================================
void Print_supp(FILE * fp, int * decoded)
{
	int N = mod2sparse_cols(H);
	mod2entry *e;

	fprintf(fp, "\nSupp(e)\n");
	for(int j = 0; j < N; j++)
	{
		if(decoded[j] != 0)
		{
			fprintf(fp, "%d - { ", j);
			for(e = mod2sparse_first_in_col(H,j);!mod2sparse_at_end(e);e = mod2sparse_next_in_col(e))
			{
				fprintf(fp,"%d ", e->row);
			}
			fprintf(fp, "}\n");
		}
	}
}


void Read_Wrong_Variable_Node(FILE * fp)
{
	int v = -1;
	int index;
	char temp[255];

	for(int i = 0 ; i < g_CODE_N ; i++)
	{
		g_wrong_variable_idx[i] = 0;
	}

	fscanf(fp, "%d", &v);
	//printf("A : %d ", v);

	memset(temp, 0, 255);
	fscanf(fp, "%s", temp);
	//printf("%s\n", temp);

	if(v > 0)
	{
		//printf("B :");
		for(int i = 0 ; i < v ; i++)
		{
			fscanf(fp, "%d", &index);
			g_wrong_variable_idx[index] = 1;
			//printf("%d ", index);

		}
		//printf("\n");
	}

	g_wrong_variable_num = v;
}

