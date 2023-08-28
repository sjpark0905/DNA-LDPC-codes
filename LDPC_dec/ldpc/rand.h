#ifndef __RAND__H__
#define __RAND__H__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mkl_vsl.h"


/* RAND.H - Interface to random number generation procedures. */

int rand_seed (int);		/* Initialize current state structure by seed */


/* GENERATORS FOR VARIOUS DISTRIBUTIONS. */

int rand_uniform (double *result, int number);	/* Uniform from [0,1) */
int rand_gaussian (double *result, int number, double std_dev);	 /* Gaussian with mean zero and unit variance */
int rand_uniformBit (unsigned int *result, int number);	 /* Uniform integer */
int rand_save(char *filename);
int rand_load(char *filename);
int rand_copy();
int rand_save_err(char *filename);
int rand_del(int save);
int rand_int(int r);

#endif