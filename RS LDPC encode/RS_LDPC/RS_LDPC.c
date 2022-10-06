/*
Name: NGUYEN VAN CUONG
ID: 2009126700
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int q, s, rho, gamma, H_pri;
int **gf_table;

void polynominal(int* poly,int s)
{
     int j;

     switch (s)
     {
          case 2: //p(x)=1+x+x^2
               poly[0]=1;
               poly[1]=1;
               break;
          case 3: //p(x)=1+x+x^3
               poly[0]=1;
               poly[1]=1;
               poly[2]=0;
               break;
          case 4: //p(x)=1+x+x^4
               poly[0]=1;
               poly[1]=1;
               poly[2]=0;
               poly[3]=0;
               break;
          case 5: //p(x)=1+x^2+x^5
               poly[0]=1;
               poly[1]=0;
               poly[2]=1;
               poly[3]=0;
               poly[4]=0;
               break;
          case 6: //p(x)=1+x+x^6
               poly[0]=1;
               poly[1]=1;
               poly[2]=0;
               poly[3]=0;
               poly[4]=0;
               poly[5]=0;
               break;
          case 7: //p(x)=1+x^3+x^7
               poly[0]=1;
               poly[1]=0;
               poly[2]=0;
               poly[3]=1;
               poly[4]=0;
               poly[5]=0;
               poly[6]=0;
               break;
          case 8: //p(x)=1+x^2+x^3+x^4+x^8
               poly[0]=1;
               poly[1]=0;
               poly[2]=1;
               poly[3]=1;
               poly[4]=1;
               poly[5]=0;
               poly[6]=0;
               poly[7]=0;
               break;
          case 9: //p(x)=1+x^4+x^9
               poly[0]=1;
               poly[1]=0;
               poly[2]=0;
               poly[3]=0;
               poly[4]=1;
               poly[5]=0;
               poly[6]=0;
               poly[7]=0;
               poly[8]=0;
               break;
         case 10: //p(x)=1+x^3+x^10
               poly[0]=1;
               poly[1]=0;
               poly[2]=0;
               poly[3]=1;
               poly[4]=0;
               poly[5]=0;
               poly[6]=0;
               poly[7]=0;
               poly[8]=0;
               poly[9]=0;
               break;
     }
}

void make_table()
{
	int i,j;
	int* poly;

	//init polinominal corresponding to s
	poly = (int*)calloc(s,sizeof(int));
	polynominal(poly,s);
	//initial elemtent 1
	gf_table[0][0]=1;
	for (j=1;j<s;j++)
        gf_table[0][j]=0;
    // ---------end of initial element 0 and 1
	for (i=1;i<q-1;i++)
	{ //from alpha^i-1 --> alpha^i
           gf_table[i][0]=0;
           for (j=1;j<s;j++)
               gf_table[i][j]=gf_table[i-1][j-1];
           if (gf_table[i-1][s-1]==1)
              for (j=0;j<s;j++)
                  gf_table[i][j]^= poly[j];
    }
}

int gf_add(int i, int j)
{
	// return k where alpha^k = alpha^i + alpha^j
	// if i=j, return -1
	// Be careful when i or j equals -1
	int k,t;
    int level;
    int* result;
    int ok;

    level=q-1;
    if (i==-1)
         return j%level;
    else if (j==-1)
         return i%level;
    else {
      //limiting in range of [0,level]
      i%=level;
      j%=level;
      if (i==j)
         return -1;
      result = (int*)calloc(s,sizeof(int));
      //find alpha^i+alpha^j
      for (k=0;k<s;k++)
	     result[k]=gf_table[i][k] ^ gf_table[j][k];
	  //find k whereas  alpha^k = alpha^i + alpha^j
      for (t=0;t<=level;t++)
      {
          for (k=0;k<s;k++)
          {
              ok=1;
              if (result[k]!=gf_table[t][k])
              {
                 ok=0;
                 break;
              }
          }
          if (ok==1)
             return t;
      }
   }
}

int gf_mult(int i, int j)
{
	// return k where alpha^k = alpha^i * alpha^j
	// Be careful when i or j equals -1
	int level;

	level=q-1;
	if (i==-1 || j==-1)
       return -1;
    else
       return (i+j)%level;
}

//function poly_mult will make multiply (alpha^b+x)(alpha^a0+alpha^a1*x+alpha^a2*x^2+....+alpha^an*x^n)
void poly_mult_xb(int b, int* poly, int n)
{
   int i;

   poly[n+1]=poly[n];

   for (i=n;i>0;i--)
      poly[i]=gf_add(poly[i-1],gf_mult(b,poly[i]));

   poly[0]=gf_mult(b,poly[0]);
}

void make_gen_poly(int* gen_poly)
{
	// store the coefficients of the generator polynomial in "gen_poly"
	// t are known
	int i;
	//init genator=x+alpha^1
	gen_poly[0]=1;
	gen_poly[1]=0;
	//Generate generator (x+alpha^1)...(x+alpha^(2*t))
	for (i=1;i<rho-2;i++)
	    poly_mult_xb(1+i,gen_poly,i);
}


void encode(int* gen_row1, int* gen_row2, int** codeword)
{
     int i,j,k;

	 for (i=-1; i<q-1; i++)
	 {
		 for (j=-1; j<q-1; j++)
		 {
			 for (k=0; k<rho; k++)
			 {
				 codeword[(i+1)*q + (j+1)][k] = gf_add(gf_mult(i,gen_row1[k]),gf_mult(j,gen_row2[k]));
			 }

		 }
	 }
}



int main(int argc, char *argv[])
{


	int *gen_poly;
	int *gen_row1;
	int *gen_row2;

	int **codeword;
	int **Cb;
	int *coset_num;

	int **H;
	int M, N;

	FILE *fp;
	char filename[255];
	int i,j,k,l;				
	int seed=1;
	int cnt;
	int selected;
	

	// input parameters

	  if (argc == 6)
	  {
		    s = atoi(argv[1]);
			rho = atoi(argv[2]);
			gamma = atoi(argv[3]);
			sprintf(filename,"%s",argv[4]);
			H_pri = atoi(argv[5]);
	  }
	  else
	  {
		    printf("Incorrect input arguments. Check again.\n");
		    printf("[Usage] ./RS_LDPC [s] [rho] [gamma] [out_filename]\n");
		    exit(0);
	  }


	// set paramters

	q = (int)pow(2.0,s);

	M = gamma*q;
	N = rho*q;


	// memory allocation

	gf_table = (int**)calloc((q-1)*s,sizeof(int));
	for(i=0; i<q-1; i++)
	{
		gf_table[i] = (int*)calloc(s,sizeof(int));
	}

	gen_poly = (int*)calloc(rho-1,sizeof(int));
	gen_row1 = (int*)calloc(rho,sizeof(int));
	gen_row2 = (int*)calloc(rho,sizeof(int));
	codeword = (int**)calloc(q*q*rho,sizeof(int));
	Cb = (int**)calloc(M*rho,sizeof(int));			//gamma개의  Cb 전체
	coset_num = (int*)calloc(q*q,sizeof(int));		//아마도 circulant matrix, 이걸 rho*gamma 크기로 circulant array.
	H = (int**)calloc(M*N,sizeof(int));

	for (i=0; i<q*q; i++)
	{
		codeword[i] = (int*)calloc(rho,sizeof(int));
	}

	for (i=0; i<M; i++)
	{
		Cb[i] = (int*)calloc(rho,sizeof(int));
	}

	for (i=0; i<M; i++)
	{
		H[i] = (int*)calloc(N,sizeof(int));
	}


	// initialize
	for (i=0; i<q*q; i++)
	{
		coset_num[i] = -1;
	}

	for (i=0; i<M; i++)
	{
		for (j=0; j<N; j++)
		{
			H[i][j] = 0;
		}
	}


	// make a generator polynomial

	make_table();
	make_gen_poly(gen_poly);
	
	// make two rows of the generator matrix

	for (i=0; i<rho-1; i++)
	{
		gen_row1[i] = gen_poly[i];
		gen_row2[i+1] = gen_poly[i];
	}
	gen_row1[rho-1] = -1;
	gen_row2[0] = -1;


	// make all codewords

	encode(gen_row1,gen_row2,codeword);

	// make the first coset Cb^(1) of size q X rho
	
	for (i=0; i<q*q; i++)		// select a codeword of weight rho
	{
		cnt=0;

		for (j=0; j<rho; j++)
		{
			if (codeword[i][j] != -1) cnt++;		//전체 codeword 중에서 weight가 rho인 codeword만 찾는 중
		}

		if (cnt==rho)
		{
			selected = i;
			break;									//찾으면 i를 저장한 후 break해서 바로 나옴
		}
	}

	for (i=0; i<q; i++)			// store the first coset Cb^(1)
	{
		for (j=0; j<rho; j++)
		{
			Cb[i][j] = gf_mult(i-1,codeword[selected][j]);		//gf_mult 에서 input이 -1이면 all-zero vector
		}
	}

	for (i=0; i<q; i++)			// check if a codeword is contained by the first coset Cb^(1)
	{
		for (j=0; j<q*q; j++)
		{
			cnt=0;
			for (k=0; k<rho; k++)
			{
				if (Cb[i][k]==codeword[j][k]) cnt++;		//Cb[i]에 속한 codeword를 모두 찾아서, codeword row 번호의 coset_num =0로 설정
			}
			if (cnt==rho)
			{
				coset_num[j] = 0;
				break;
			}
		}
	}

	// make the other cosets

	for (i=1; i<gamma; i++)		// select a codeword not belonging to previous cosets
	{
		for (j=0; j<q*q; j++)
		{
			if (coset_num[j] == -1)
			{
				selected = j;
				break;
			}
		}

		for (j=0; j<q; j++)			// store the (i+1)-th coset Cb^(i+1)
		{
			for (k=0; k<rho; k++)
			{
				Cb[j+i*q][k] = gf_add(Cb[j][k],codeword[selected][k]);
			}
		}

		for (j=0; j<q; j++)			// check if a codeword is contained by the (i+1)-th coset Cb^(i+1)
		{
			for (k=0; k<q*q; k++)
			{
				cnt=0;
				for (l=0; l<rho; l++)
				{
					if (Cb[j+i*q][l]==codeword[k][l]) cnt++;
				}
				if (cnt==rho)
				{
					coset_num[k] = i;
					break;
				}
			}
		}
	}


	// make H

	for (i=0; i<M; i++)
	{
		for (j=0; j<rho; j++)
		{
			H[i][j*q+Cb[i][j]+1] = 1;
		}
	}



	// make alists

	fp = fopen(filename,"w");

	fprintf(fp,"%d %d\n",M,N);
	fprintf(fp,"%d %d\n",rho,gamma);
	
	for (i=0; i<M; i++)
	{
		cnt=0;
		for (j=0; j<N; j++)
		{
			cnt+=H[i][j];
		}
		fprintf(fp,"%d ",cnt);
	}
	fprintf(fp,"\n");

	for (j=0; j<N; j++)
	{
		cnt=0;
		for (i=0; i<M; i++)
		{
			cnt+=H[i][j];
		}
		fprintf(fp,"%d ",cnt);
	}
	fprintf(fp,"\n");

	for (i=0; i<M; i++)
	{
		for (j=0; j<N; j++)
		{
			if (H[i][j] == 1) fprintf(fp,"%d ",j+1);
		}
		fprintf(fp,"\n");
	}
	
	for (j=0; j<N; j++)
	{
		for (i=0; i<M; i++)
		{
			if (H[i][j] == 1) fprintf(fp,"%d ",i+1);
		}
		fprintf(fp,"\n");
	}
	
	fclose(fp);


	// print out the results
	if (H_pri == 0)
	{
		
		printf("N = %d\nM = %d\n", N, M);
		printf("Generator polynomial\n");
		for (i = 0; i < rho - 1; i++)
		{
			printf("%d ", gen_poly[i]);
		}
		printf("\n");

		printf("Coset\n");
		for (i = 0; i < q * q; i++)
		{
			printf("%d ", coset_num[i]);
		}
		printf("\n");
		
	}
	
	else if (H_pri == 1)
	{
		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				printf("%d ", H[i][j]);
			}
			printf("\n");
		}
	}
	else
	{

		printf("%d %d", M, N);

	}

	/*
	// print  H

	fp = fopen(filename,"w");
	 fprintf(fp,"a=gf([");
	 for (i=0; i<M; i++)
	 {
		 for (j=0; j<N; j++)
		 {
			 fprintf(fp," %d",H[i][j]);
		 }
		 fprintf(fp,";");
	 }
	 fprintf(fp,"]); r_%d_%d_%d=rank(a);\n",s,rho,gamma);
	 fclose(fp);
	 */


	// free memory

	for(i=0; i<q-1; i++)
	{
		free(gf_table[i]);
	}

	for (i=0; i<q*q; i++)
	{
		free(codeword[i]);
	}

	for (i=0; i<M; i++)
	{
		free(Cb[i]);
	}

	for (i=0; i<M; i++)
	{
		free(H[i]);
	}

	
	free(gen_poly);
	free(gen_row1);
	free(gen_row2);
	free(coset_num);

	free(codeword);
	free(gf_table);
	free(Cb);
	free(H);

	
	
	return 0;
}
