/* CHECK.C - Compute parity checks and other stats on decodings. */

/* Copyright (c) 2001 by Radford M. Neal 
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mod2sparse.h"
#include "check.h"



/* COMPUTE PARITY CHECKS.  Returns the number of parity checks violated by
   dblk.  The results of all the parity checks are stored in pchk. */

int check
( 
	mod2sparse *H,	/* Parity check matrix */
	char *dblk,		/* Guess for codeword */
	char *pchk		/* Place to store parity checks */
)
{
	int M, i, c;

	M = mod2sparse_rows(H);
	mod2sparse_mulvec(H, dblk, pchk);

	c = 0;	// unsatisfied check num
	for (i = 0; i < M; i++) 
	{ 
		c += pchk[i];
	}

	return c;
}

int check_bound
( 
	mod2sparse *H,	/* Parity check matrix */
	char *dblk,		/* Guess for codeword */
	char *pchk,		/* Place to store parity checks */
	int V_Start,
	int V_End,
	int C_Start,
	int C_End
)
{
	int M, i, c;

	M = mod2sparse_rows(H);
	mod2sparse_mulvec_bound(H, dblk, pchk,V_Start,V_End,C_Start,C_End);

	c = 0;	// unsatisfied check num
	for (i = C_Start; i < C_End; i++) 
	{ 
		c += pchk[i];
	}

	return c;
}