#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
#include <stdlib.h>

/*********************************************************************

Eighth Implementation Challenge

This file contains functions to generate, in TSPLIB format,
Euclidean TSP instances with integer x and y coordinates 
uniformly and independently chosen in the range 0..MAXCOORD-1.

A portable random number generator, based on an algorithm of Knuth,
to generate the coordinates.

To use this program, run

	portgen ncities seed > filename

where "ncities" is the number of cities and "seed" is an integer seed
for the generator.  The output goes to stdout, so you should redirect 
it to a file.
 
This code is based on previous work by David Johnson, Jon Bentley, and 
others.


						Lyle McGeoch
						lam@cs.amherst.edu

						June 24, 1999
 *********************************************************************/

#define MAXCOORD 1000000
#define PRANDMAX 1000000000

int a,b;
int arr[55];
void sprand(int);
int lprand(void);

int main(int argc, char **argv)
{
	int factor = PRANDMAX/MAXCOORD;
	int N;
	int i;
	int x, y;
	int seed;

	if (argc < 2) {
		printf("Usage: portgen num.of.cities seed > filename\n");
		exit(1);
		}
	N = atoi(argv[1]);
	seed = atoi(argv[2]);

	/* initialize random number generator */

	sprand(seed);

	printf ("NAME : portgen-%d-%d\n", N, seed);
	printf ("COMMENT : portgen N=%d, seed=%d\n", N, seed);
	printf ("TYPE : TSP\n");
	printf ("DIMENSION : %d\n", N);
	printf ("EDGE_WEIGHT_TYPE : EUC_2D\n");
	printf ("NODE_COORD_SECTION\n");

	for (i=1;i<=N;i++) {
	    x = lprand()/factor;
	    y = lprand()/factor;
	    printf("%d %d %d\n", i, x, y);
	}
	return 0;
}

void sprand (int seed)
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

int lprand(void)
{	long t;
	if (a-- == 0) a = 54;
	if (b-- == 0) b = 54;
	t = arr[a] - arr[b];
	if (t<0) t+= PRANDMAX;
	arr[a] = t;
	return t;
}
