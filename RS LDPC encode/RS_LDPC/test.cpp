#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int i, j, cnt;
char filename[30];
int M = 3;
int N = 7;
int H[3][7] = {
	{1, 0, 1, 0, 1, 0, 1},
	{0, 1, 1, 0, 0, 1, 1},
	{0, 0, 0, 1, 1, 1, 1},
};




int main()
{

FILE* fp;

sprintf(filename, "parity.alist");
fp = fopen(filename, "w");

fprintf(fp, "%d %d\n", M, N);

for (i = 0; i < M; i++)
{
	cnt = 0;
	for (j = 0; j < N; j++)
	{
		cnt += H[i][j];
	}
	fprintf(fp, "%d ", cnt);
}
fprintf(fp, "\n");

for (j = 0; j < N; j++)
{
	cnt = 0;
	for (i = 0; i < M; i++)
	{
		cnt += H[i][j];
	}
	fprintf(fp, "%d ", cnt);
}
fprintf(fp, "\n");

for (i = 0; i < M; i++)
{
	for (j = 0; j < N; j++)
	{
		if (H[i][j] == 1) fprintf(fp, "%d ", j + 1);
	}
	fprintf(fp, "\n");
}

for (j = 0; j < N; j++)
{
	for (i = 0; i < M; i++)
	{
		if (H[i][j] == 1) fprintf(fp, "%d ", i + 1);
	}
	fprintf(fp, "\n");
}

fclose(fp);
}