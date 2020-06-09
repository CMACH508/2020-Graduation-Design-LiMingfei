#include <stdio.h>
#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"

static FILE *tourfp;           /* file pointer for output tour */

/*********************************************************************

Eighth Implementation Challenge

This program measures the length of a TSP tour.  It uses a reader (in
Instance.c) to read TSP instances in the TSPLIB format.

To use this program, run

	length filename tour_filename

where "filename" is a TSP instance in some supported TSPLIB format (see 
Instance.c for details), and tour_filename is a tour, a permuted list of the
numbers 1 thru N.

						Lyle McGeoch
						lam@cs.amherst.edu

						June 24, 1999

 *********************************************************************/


static void checktour(void) {       /* also saves tour in file */

  int *present;
  int i, first, prev, curr;
  double cost;
  char msg[200];

  cost = 0;

  present  = (int *) calloc(ncities+1, sizeof(int));
  for (i=1;i<=ncities;i++) present[i] = 0;

  if (!fscanf (tourfp, "%d", &first))
    fatal ("error reading tour file");

  if (first < 1 || first > ncities) {
    sprintf (msg, "City number (%d) is out of correct range", first);
    fatal (msg);
  }

  present[first] = 1;
  prev = first;

  for (i=2;i<=ncities;i++) {

    if (!fscanf (tourfp, "%d", &curr))
      fatal ("error reading tour file");

    if (curr < 1 || curr > ncities) {
      sprintf (msg, "City number (%d) is out of correct range", curr);
      fatal (msg);
    }

    if (present[curr]) {
      sprintf (msg, "City %d appears twice in tour", curr);
      fatal (msg);
    }

    present[curr] = 1;

    cost += distance (prev, curr);
    prev = curr;
  }

  cost += distance (prev, first);

  printf ("Tour length is %0.0f\n", cost);
}

/*******************************************************************/

int main (int argc, char **argv, char **envp) {
  
  if (argc < 3) fatal ("Usage:  length filename tour_filename");
  
  readInstance(argv[1]);

  tourfp = fopen (argv[2], "r");
  if (tourfp == 0)
    fatal ("Can't open tour file");

  checktour();

  return 0;
}


