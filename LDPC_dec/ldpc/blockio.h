#ifndef __BLOCK_IO_H__
#define __BLOCK_IO_H__

/* BLOCKIO.H - Interface to block input/output routines. */

/* Copyright (c) 2000 by Radford M. Neal 
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

int  blockio_read  (FILE *, char *, int);
void blockio_write (FILE *, char *, int);

#endif