#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
#include <math.h>
#include <stdlib.h>

/*********************************************************************

Eighth Implementation Challenge

This file contains functions to generate, in TSPLIB format, random
clustered Euclidean TSP instances with coordinates in the range
0..MAXCOORD-1.

A portable random number generator, based on an algorithm of Knuth,
to generate the coordinates.

To use this program, run

	portcgen ncities seed > filename

where "ncities" is the number of cities and "seed" is an integer seed
for the generator.  The output goes to stdout, so you should redirect 
it to a file.
 
This code is based on previous work by David Johnson, Jon Bentley, and 
others.


						Lyle McGeoch
						lam@cs.amherst.edu

						June 24, 1999
 *********************************************************************/

#define MAXN 1000000
#define MAXCOORD 1000000
#define PRANDMAX 1000000000
#define CLUSTERFACTOR 100
#define SCALEFACTOR 1.0

int center[MAXN+1][2];
int a,b;
int arr[55];
void sprand(int);
double lprand(void);
double normal(void);

int main(int argc, char **argv)
{
	int N;
	int c;
	int i,j;
	int x, y;
	int seed;
	int nbase;
	double scale;

	if (argc < 2) {
		printf("Usage: portcgen num.of.cities seed > filename\n");
		exit(1);
		}
	N = atoi(argv[1]);
	seed = atoi(argv[2]);

	/* initialize random number generator */

	sprand(seed);

	nbase = N/CLUSTERFACTOR;
	scale = SCALEFACTOR/sqrt((double)N);

	printf ("NAME : portcgen-%d-%d\n", N, seed);
	printf ("COMMENT : portcgen N=%d, seed=%d\n", N, seed);
	printf ("TYPE : TSP\n");
	printf ("DIMENSION : %d\n", N);
	printf ("EDGE_WEIGHT_TYPE : EUC_2D\n");
	printf ("NODE_COORD_SECTION\n");

	for (i=1;i<=nbase;i++) 
		for (j=0;j<=1;j++)
	    		center[i][j] = (int) (lprand()/PRANDMAX*MAXCOORD);

	for (i=1;i<=N;i++) {
		c = (int) (lprand()/PRANDMAX*nbase) + 1;
		x = center[c][0] + (int) (normal()*scale*MAXCOORD);
		y = center[c][1] + (int) (normal()*scale*MAXCOORD);
		printf("%d %d %d\n",i,x,y);
	}
	return 0;
}

double normal(void)	/* Algorithm 3.4.1.P, p. 117, Knuth v. 2 */
{
	static int	goodstill = 0;
	static double	nextvar;
	double		s, t, v1, v2;

	if (goodstill) {
		goodstill = 0;
		return nextvar;
	}
	else {
		goodstill = 1;
		do {
			v1 = 2*lprand()/PRANDMAX - 1.0;
			v2 = 2*lprand()/PRANDMAX - 1.0;
			s = v1*v1 + v2*v2;
		} while (s >= 1.0);
		t = sqrt((-2.0 * log(s)) / s);
		nextvar = v1 * t;	/* Knuth's x1 */
		return v2 * t;		/* Knuth's x2 */
	}
}

void sprand (seed)
int seed;
{	int i,ii;
	int last,next;
	arr[0] = last = seed;
	next = 1;
	for (i=1;i<55;i++) {
		ii = (21*i)%55;
		arr[ii] = next;
		next = last - next;
		if (next < 0) next += PRANDMAX;
		last = arr[ii];
	}
	a = 0;
	b = 24;
	for (i=0;i<165;i++) last = lprand();
}

double lprand(void)
{	long t;
	if (a-- == 0) a = 54;
	if (b-- == 0) b = 54;
	t = arr[a] - arr[b];
	if (t<0) t+= PRANDMAX;
	arr[a] = t;
	return (double)t;
}

