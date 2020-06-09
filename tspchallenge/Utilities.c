#include <stdio.h>
#include <stdlib.h>

#include "Utilities.h"

void fatal (char *s) {
  fprintf (stderr, "%s\n", s);
  exit(1);
}

void **createMatrix (int i, int j, int s) {

  int k;
  char **result = (char **) (malloc(i*sizeof(result) + i*j*s));
  char *p = (char *)result + (i*sizeof(result));

  for (k=0; k<i; ++k) {
    result[k] = p;
    p += j*sizeof(s);
  }

  return (void **)result;
}

#define PRANDMAX 1000000000
int a,b;
int arr[55];

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

