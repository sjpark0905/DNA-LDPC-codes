/* RAND.C - Random number generation module. */

#include "rand.h"

static int					initialized = 0;		/* Has module been initialized? */
static VSLStreamStatePtr	state		= NULL;		/* Pointer to current state */
static VSLStreamStatePtr	save_state	= NULL;		/* Pointer to save state */


/* SET CURRENT STATE ACCORDING TO SEED. */
int rand_seed( int seed)
{ 
	initialized = 1;
	return	vslNewStream(&state, VSL_BRNG_MT2203, seed*37);
}

/* GENERATE UNIFORMLY FROM [0,1]. */
int rand_uniform (double *result, int number)
{
	if(!initialized)	rand_seed(0);
	return vdRngUniform(0,state,number,result,0.0,1.0);
}

int rand_gaussian (double *result, int number, double std_dev)
{
	if(!initialized)	rand_seed(0);
	return vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,state,number,result,0.0,std_dev);
}

int rand_uniformBit (unsigned int *result, int number)
{
	if(!initialized)	rand_seed(0);
	return viRngUniformBits(0,state,number,result);
}

int rand_save(char *filename)
{
	if(!initialized)	rand_seed(0);
	return vslSaveStreamF(state,filename);
}

int rand_load(char *filename)
{
	initialized = 1;
	return vslLoadStreamF(&state,filename);
}

int rand_copy()
{
	if(!initialized)	rand_seed(0);
	return vslCopyStream(&save_state,state);
}

int rand_save_err(char *filename)
{
	if(!initialized)	rand_seed(0);
	return vslSaveStreamF(save_state,filename);
}

// save : 0 - state 을 삭제
// save : 1 - save_state 을 삭제
int rand_del(int save)
{
	int result = 0;
	if(save == 0)
	{
		if(state != NULL)
		{
			result = vslDeleteStream(&state);
			state = NULL;
		}
	}
	else
	{
		if(save_state != NULL)
		{
			result         =       vslDeleteStream(&save_state);
			save_state     =       NULL;
		}
	}
	return result;
}

int rand_int(int r)
{
	int result = 0;
	unsigned int a;
	int v = 0;
	if(!initialized)	rand_seed(0);
	result = viRngUniformBits(0, state, 1, &a);
	v = a % r;
	return v;
}