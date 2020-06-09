#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
#include <stdlib.h>

#define MAXCOORD 1000000

/*********************************************************************

Eighth Implementation Challenge

This file contains functions to generate, in TSPLIB format,
random distance matrices with distances in the range 0..MAXCOORD-1.
It uses a portable random number generator, based on an algorithm of Knuth.

To use this program, run

	portmgen ncities seed > filename

where "ncities" is the number of cities and "seed" is an integer seed
for the generator.  The output goes to stdout, so you should redirect 
it to a file.
 
The distance between two cities is calculated by applying a hash function,
based on the seed, to the indices of the cities.  Participants can run
their code on random distance matrices, without explicitly storing the
matrices, by using the techniques shown here to compute distances 
"on-the-fly."

This code is based on previous work by David Johnson, Jon Bentley, and 
others.


						Lyle McGeoch
						lam@cs.amherst.edu

						June 24, 1999
 *********************************************************************/

int a,b;
int arr[55];
int *key;
int param;
double factor;

void initializeKeys(int N, int seed) {

  int i;

  factor = MAXCOORD/2147483648.;
  param = 104*seed + 1;
  key = (int *) calloc (N+1, sizeof(int));

  for (i=1; i<=N; ++i) key[i] = 0x12345672*i + 1;
}

int distance (int i, int j) {
  
  int x, y, z;

  if (i == j) return 0;

  i = key[i];
  j = key[j];

  x = i&j;
  y = i|j;
  z = param;

  x *= z;
  y *= x;
  z *= y;

  z ^= param;

  x *= z;
  y *= x;
  z *= y;

  x = ((i+j)^z) & 0x7fffffff;
  return (int)(x*factor);
}

int main(int argc, char **argv)
{
	int N;
	int i,j;
	int seed;

	if (argc < 2) {
		printf("Usage: portmgen num.of.cities seed > filename\n");
		exit(1);
		}
	N = atoi(argv[1]);
	seed = atoi(argv[2]);

	initializeKeys(N, seed);

	printf ("NAME : portmgen-%d-%d\n", N, seed);
	printf ("COMMENT : portmgen N=%d, seed=%d\n", N, seed);
	printf ("TYPE : TSP\n");
	printf ("DIMENSION : %d\n", N);
	printf ("EDGE_WEIGHT_TYPE : EXPLICIT\n");
	printf ("EDGE_WEIGHT_FORMAT : UPPER_DIAG_ROW\n");
	printf ("EDGE_WEIGHT_SECTION\n");

	for (i=1;i<=N;i++) {
	  for (j=i; j<=N; ++j)
	    printf ("%d\n", distance(i,j));
	}
	return 0;
}
