#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>

#ifdef TIMING
#include <sys/times.h>
#endif

#include <limits.h>

#include "Instance.h"
#include "Utilities.h"

/*********************************************************************

Eighth Implementation Challenge

This file contains the main program for obtaining greedy (multifragment)
tours.  It uses a reader (in Instance.c) to read TSP instances in the 
TSPLIB format.

To use this program, run

	greedy filename [output_tour_filename]

where "filename" is a TSP instance in some supported TSPLIB format (see 
Instance.c for details).  The final tour, a permutation of the city numbers,
is stored in "output_tour_filename," if this parameter is given.  The
filename chosen can not begin with a digit.

To do timing of multiple iterations of the algorithm, you can also run

        greedy filename number_of_iterations

For Euclidean instances, k-d trees are used to find nearest neighbors.  For
matrix instances, a simple linear-time search is used whenever a nearest 
neighbor is needed.  (The latter isn't very efficient, but the processing
time for a matrix seems to be dominated by the time required to read the
matrix.)

The functions in this file fall into four classes.  In the order in
which they appear in this file, they are:

 	Functions for k-d trees
	Function to find the nearest neighbor in a matrix instance
	Functions implementing a heap
	Top-level functions (including main()) for computing greedy tours.

Distances between cities are assumed to be less than MAXINT.  Ints are
assumed to be at least 32 bits.

Random numbers are needed in the construction of the k-d trees.  They
are obtained from a portable random number generator.

This code is based on previous work by David Johnson, Jon Bentley, and 
others.

						Lyle McGeoch
						lam@cs.amherst.edu

						June 24, 1999
 *********************************************************************/

/************************ variables used in kd-trees *****************/

#define tloson(i) (2*(i))
#define thison(i) (2*(i)+1)

#define MAXINT (0x7fffffff)

typedef int xypair[2];
typedef xypair *EuclidGraph;

static EuclidGraph G      = 0;

static int	*tp       = 0;      	/* permutation vector */
static int	troot;          	/* root of k-d tree */
static int	tlevs;                  /* # levels in k-d tree: 0 .. tlevs-1; 
					   tlevs<MAXLEVS */
static int	*tbucket  = 0;	        /* Nodes. tbucket==0 => internal */
static int	*tdim     = 0;	        /* cut dimension */
static int	*tval     = 0;	        /* cut value */
static int	*tempty   = 0;	        /* 1 if entire subtree deleted */
static int	*tfirst   = 0;          /* first point in bucket */
static int	*tlast    = 0;          /* last point in bucket */
static int	*home     = 0;

static int	*degree   = 0;
static int	*tail     = 0;
static int	*present  = 0;   /* is this vertex temporarily deleted? */
static int	*nb1      = 0;
static int      *nb2      = 0;

static int	nntarget;	/* pt to find nn of */
static int	nndist;
static int      kdcutoff;
static int	nnbest;
static int	delct;
static int	nncount;

/*******************************************************************/
/*                                                                 */
/*           Functions implementing k-d trees                      */
/*                                                                 */
/*******************************************************************/

static EuclidGraph makeGraph (int n) {

  int i;
  
  if (isEuclidean) {
    G = (xypair *)calloc(n+1, sizeof(xypair));
    for (i=1; i<=n; ++i) {
      G[i][0] = (int) (xx[i] + 0.5);
      G[i][1] = (int) (yy[i] + 0.5);
    }
    return G;
  }
  else {
    fatal ("Can't run k-d trees on non-Euclidean instances");
    return 0;
  }
}

/*******************************************************************/

static int spread(int l, int u, int d) {
  int i;
  register int minv,maxv,t;

  minv = MAXINT; maxv = -MAXINT;

  for (i = l; i <= u; i++) {
    t = G[tp[i]][d];
    if      (t < minv) minv = t; 
    else if (t > maxv) maxv = t; 
  }
  return (maxv-minv);
}
  
/*******************************************************************/

static void swap(int i, int j) {
  int t;
  t = tp[i]; tp[i] = tp[j]; tp[j] = t;
}

/*******************************************************************/

static void partition(int l, int u, int k, int d) {
  int m, i, t, t1;
  if (l < u) {
    m = l + lprand() % (u-l+1);
    swap(l, m);
    t = G[tp[l]] [d];
    m = l;
    for (i = l+1; i <= u; i++)
      if (G[tp[i]] [d] < t) {
				/* swap(++m, i); */
	t1 = tp[++m]; tp[m] = tp[i]; tp[i] = t1;
      }
    swap(l, m);
    if      (m < k)	partition(m+1, u, k, d);
    else if (m > k)	partition(l, m-1, k, d);
  }
}

/*******************************************************************/

static void rbuild(int l, int u, int p, int level)
{
  int m, maxdim, j, maxval, ts;
  if (level >= tlevs) { /* bucket */
    tbucket[p] = 1;
    tempty[p] = 0;
    tfirst[p] = l;
    tlast[p] = u;
    for (j=l;j<=u;j++) {
      home[tp[j]] = p;
    }
  } else { /* internal node */
    tbucket[p] = 0;
    tempty[p] = 0;
    maxval = -MAXINT;
    maxdim = 0;            /* to keep compiler happy */
    for (j = 0; j < 2; j++)
      if ((ts=spread(l, u, j)) > maxval) {
	maxdim = j; maxval = ts;
      }
    m = (l+u)/2;
    partition(l, u, m, maxdim);
    tdim[p] = maxdim;
    tval[p] = G[tp[m]] [maxdim];
    rbuild(l,   m, tloson(p), level+1);
    rbuild(m+1, u, thison(p), level+1);
  }
}
  
/*******************************************************************/

static void Allocate(int size) {
  
  tp       = (int *) calloc(size+1, sizeof(int));
  home     = (int *) calloc(size+1, sizeof(int));

  tbucket  = (int *) calloc((1<<(tlevs+1))+1, sizeof(int));
  tdim     = (int *) calloc((1<<(tlevs+1))+1, sizeof(int));
  tval     = (int *) calloc((1<<(tlevs+1))+1, sizeof(int));
  tempty   = (int *) calloc((1<<(tlevs+1))+1, sizeof(int));
  tfirst   = (int *) calloc((1<<(tlevs+1))+1, sizeof(int));
  tlast    = (int *) calloc((1<<(tlevs+1))+1, sizeof(int));
}

/*******************************************************************/

static void buildtree(int n) {
  int i;
  tlevs = (int) (log((float) n)/log(2.0) - 2);
  Allocate(n);
  for (i = 1; i <= n; i++) tp[i] = i;
  troot = 1;
  rbuild(1, n, 1, 0);
}

/*******************************************************************/

static void rnn(int p) {
  int i,d,v;
  int thisdist;
  int ti;

  if (tempty[p]) return;

  if (tbucket[p]) {
    for (i=tfirst[p]; i<=tlast[p]; i++) {
      ti = tp[i];
      if (present[ti]) {
	thisdist = distance(ti,nntarget);
	if (thisdist < nndist) {
	  kdcutoff = nndist = thisdist;
	  nnbest = ti;
	  if (notIntegral && kdcutoff < MAXINT) kdcutoff++;
	}
      }
    }
  }
  else {
    int td;

    d = tdim[p]; v=tval[p];
    td = G[nntarget][d];
    if (td < v) {
      rnn(tloson(p));
      if (kdcutoff > v - td) rnn(thison(p));
    }
    else {
      rnn(thison(p));
      if (td - v < kdcutoff) rnn(tloson(p));
    }
  }
}
  
/*******************************************************************/

static void assert(int cond, char *s) {
  if (!cond) fatal(s);
}

/*******************************************************************/

static void deletept (int target) {
  register int i,p;
  int j;
  
  delct = 0;
  j = -1;
  p = home[target];
  for (i = tfirst[p]; i <= tlast[p]; i++)
    if (tp[i] == target) j = i;
  if (j != -1) {
    delct++;
    swap(j, tlast[p]);
    tlast[p]--;
    if (tfirst[p] > tlast[p]) tempty[p] = 1;
  }
  while (p>0) {
    p = p>>1;
    if (tempty[tloson(p)] && tempty[thison(p)])
      tempty[p] = 1;
    else break;
  }
  assert(delct == 1, "invalid del count");
}

/*******************************************************************/
/*                                                                 */
/*           Functions to find neighbors in matrix instances       */
/*                                                                 */
/*******************************************************************/


static int matnn (int q) {

  int bestdist, bestcity;
  int i, j;

  bestdist = MAXINT;
  bestcity = 0;

  for (i=1; i <= ncities; ++i) {
    if (present[i]) {
      if ((j=distMat[q][i]) < bestdist) {
	bestdist = j;
	bestcity = i;
      }
    }
  }
  nndist = bestdist;
  return bestcity;
}

/*******************************************************************/
/*                                                                 */
/*           Functions implementing heap                           */
/*                                                                 */
/*******************************************************************/

typedef struct {
  int here;
  int there;
  int d; } city;

#define hswap(heap,i,j) (swaptemp=heap[i], heap[i]=heap[j], heap[j]=swaptemp)

static int heapsize;
static city swaptemp;
static city *heap = 0;

/*******************************************************************/

static void sift_down (int from)
{
  /* Given that the subtrees of the heap rooted at the children
     of "from" have the heap property, establish that property on
     the subtree rooted at "from" */
  
  register int i,j,k,to;
  register city *heapptr;
  
  heapptr = heap;
  to = heapsize;
  i = from;
  j = i<<1;
  k = j+1;
  
  while (j<=to) {
    if (j == to) {
      if (heapptr[i].d > heapptr[j].d) hswap(heapptr,i,j);
      return;
    }
    if (heapptr[i].d < heapptr[j].d) {
      if (heapptr[i].d < heapptr[k].d) return;
      else {
	hswap(heapptr,i,k);
	i = k;
      }
    }
    else if (heapptr[j].d < heapptr[k].d) {
      hswap (heapptr,i,j);
      i = j;
    }
    else {
      hswap (heapptr,i,k);
      i = k;
    }
    j = i<<1;
    k = j+1;
  }
}

/*******************************************************************/

static void sift_up (int from)
    
{   
  register int i,j;
  
  j = from;
  i = j>>1;
  while (i>=1) {
    if (heap[j].d < heap[i].d)
      hswap(heap,i,j);
    else break;
    j = i;
    i = j>>1;
  }
}

/*******************************************************************/

static void Insert (int i, int j, int length) {

  /*  if (i < 1 || i > n) printf ("weird x in Insert\n");
      if (j < 1 || j > n) printf ("weird y %d in Insert\n", j);   */

  heap[++heapsize].here = i;
  heap[heapsize].there = j;
  heap[heapsize].d = length;
  sift_up (heapsize);
}

/*******************************************************************/

static city Deletemin() {
  city result = heap[1];
  heap[1] = heap[heapsize--];
  sift_down (1);
  return result;
}

/*******************************************************************/
/*                                                                 */
/*           Functions to compute greedy tour                      */
/*                                                                 */
/*******************************************************************/

static int nn(int q)	/* return nn to G[q] */
{

  nncount++;
  present[q] = 0;

  if (isEuclidean) {
    nntarget = q;
    kdcutoff = nndist = MAXINT;
    nnbest = 0;
    rnn(troot);
  }
  else {
    nnbest = matnn(q);
  }

  present[q] = 1;
  return(nnbest);
}

/******************** timing code ******** *************************/

#ifdef TIMING
	static struct tms buffer;
	static int utimer1;

	static double run_time(void) {
	  double result = 0;
	  times(&buffer);
	  result = (buffer.tms_utime-utimer1)/(double)CLK_TCK;
	  return result;
	}

	static void start_time(void) {
	  times(&buffer);
	  utimer1 = buffer.tms_utime;
	}
#endif

/*******************************************************************/

static void addedge(int x, int y) {
  degree[x]++; degree[y]++;

  /*  if (x < 1 || x > n || y < 1 || y > n) printf ("weird x or y\n");
  if (degree[x] > 2) printf ("degree of %d exceeds 2", x);
  if (degree[y] > 2) printf ("degree of %d exceeds 2", y);

  if (degree[x] < 0) printf ("degree of %d is negative", x);
  if (degree[y] < 0) printf ("degree of %d is negative", y);     */

  if (degree[x] == 1) nb1[x] = y;
  else nb2[x] = y;
  if (degree[y] == 1) nb1[y] = x;
  else nb2[y] = x;
}

/*******************************************************************/
static FILE *tourfp;           /* file pointer for output tour */
/*******************************************************************/

static void checktour(int n) {       /* also saves tour in file */
  int i,z1,z2,z3;

  z3 = 0;                /* to keep compiler happy */

  for (i=1;i<=n;i++) present[i] = 0;
  z1 = 1;

  if (tourfp) fprintf (tourfp, "%d\n", 1);

  z2 = nb1[1];
  present[z2] = 1;
  for (i=2;i<=n;i++) {

    if (tourfp) fprintf (tourfp, "%d\n", z2);

    z3 = nb1[z2];
    if (z3==z1) z3 = nb2[z2];
    ++present[z3];
    z1 = z2; z2 = z3;
  }
  if (z3 != 1) printf ("Unclosed cycle\n");
  for (i=1;i<=n;i++) {
    if (present[i] < 1) printf("Missing %d\n",i);
    if (present[i] > 1) printf("%d Visits to %d\n",present[i],i);
  }
}

/*
static void checkN (int j, int i) {

  if (nb1[j] == i) return;
  if (degree[j] == 2 && nb2[j] == i) return;
  printf ("Neighbor problem\n");
  exit(2);
}

static void checkAll(int counter) {

  int i;

  for (i=1; i<=ncities; ++i) {

    if (degree[i] != 1) {
      if (tail[i] != -1) break;
      if (degree[i] == 2) {
	checkN (nb1[i], i);
	checkN (nb2[i], i);
      }
    }
    else {
      if (tail[tail[i]] != i) break;
      checkN (nb1[i], i);
    }

  }
  if (i <= n) {
    printf ("Error with counter %d and degree %d\n", counter, degree[i]);
    exit(1);
  }
}
*/

/*******************************************************************/

static void doGreedy(void) {

  int i;
  int newnn, loopctr, tailx,taily;
  double cost = 0;
  int x, y;

  degree   = (int *) calloc(ncities+1, sizeof(int));
  tail     = (int *) calloc(ncities+1, sizeof(int));
  present  = (int *) calloc(ncities+1, sizeof(int));
  nb1	   = (int *) calloc(ncities+1, sizeof(int));
  nb2	   = (int *) calloc(ncities+1, sizeof(int));

  heap	   = (city *) calloc(ncities+1, sizeof(city));
  
  heapsize = 0;

  if (isEuclidean) {
    G = makeGraph(ncities);
    printf("greedy on Euclidean instance %s\n", filename);
#ifdef TIMING
    start_time();
#endif  
    buildtree(ncities);
    
    printf("%d buckets in kd-tree\n",1<<tlevs);
  }
  else {
    printf("greedy on matrix instance %s\n", filename);
#ifdef TIMING
    start_time();
#endif
  }

  for (i = 1; i <= ncities; i++) present[i] = 1;
  
  for (i = 1; i <= ncities; i++) {
    degree[i] = 0;
    newnn = nn(i);
    tail[i] = -1;
    Insert(i,newnn,nndist);
  }
  
  for (loopctr = 1; loopctr < ncities; loopctr++) {
    while (1) {
      city c = Deletemin();
      x = c.here;
      y = c.there;

      if (degree[x] >= 2) {
	continue;
      }
      if (degree[y] < 2 && y != tail[x])
	break;
      if (tail[x] >= 0) present[tail[x]] = 0;
      newnn = nn(x);
      if (tail[x] >= 0) present[tail[x]] = 1;
      Insert(x,newnn,nndist);
    }
    /* x, y is a new edge in tour */
    addedge(x,y);
    cost += distance(x,y);
    assert(x != tail[y], "bad edge");

    if (isEuclidean) {          /* fix kd tree */
      if (degree[y] >= 2) {
	deletept(y);
      }
      if (degree[x] >= 2) {
	deletept(x);
      }
    }
    else {                      /* update present */
      if (degree[y] >= 2)
	present[y] = 0;
      if (degree[x] >= 2)
	present[x] = 0;
    }

    tailx = tail[x];
    taily = tail[y];
    if (tailx == -1 && taily == -1) {
      tail[x] = y;
      tail[y] = x;
    } else if (tailx == -1 && taily > -1) {
      tail[y] = -1;
      tail[x] = taily;
      tail[taily] = x;
    } else if (tailx > -1 && taily == -1) {
      tail[x] = -1;
      tail[y] = tailx;
      tail[tailx] = y;
    } else /* (tailx > -1 && taily > -1) */ {
      tail[x] = tail[y] = -1;
      tail[tailx] = taily;
      tail[taily] = tailx;
    }
    if (degree[x] == 1) {
      if (tail[x] >= 0) present[tail[x]] = 0;
      newnn = nn(x);
      if (tail[x] >= 0) present[tail[x]] = 1;
      Insert(x,newnn,nndist);
    }
    /*    checkAll(loopctr);  */
  }
  
  for (x=1;x<=ncities;x++) {
    if (tail[x] >= 0) {
      y = tail[x];
      addedge(x,y);
      cost += distance(x,y);
      if (degree[x] != 2) printf("Error on last edge\n");
      if (degree[y] != 2) printf("Error on last edge\n");
      break;
    }
  }
  
  checktour(ncities);
  printf("nncount = %d\n",nncount);
    
  printf("\nTour cost = %0.0f  ", cost);
#ifdef TIMING  
  printf("Tour time = %4.2f", run_time());
#endif

  printf("\n\n");
}

static void myfree(void *p) {
  if (p != 0) free(p);
}

static void freeMemory() {

  /* free up dynamically allocated memory */

  myfree(G);
  myfree(tp);
  myfree(tbucket);
  myfree(tdim);
  myfree(tval);
  myfree(tempty);
  myfree(tfirst);
  myfree(tlast);
  myfree(home);
  myfree(degree);
  myfree(tail);
  myfree(present);
  myfree(nb1);
  myfree(nb2);
  myfree(heap);
}

  
/*******************************************************************/

int main (int argc, char **argv, char **envp) {

  int iterations;
  int i;

  setbuf(stdout,NULL);

  if (argc < 2) {
    fprintf (stderr, 
	     "Usage:  greedy filename [trials] or\n");
    fatal (  "        greedy filename [tour_output_filename]");
  }
  
#ifdef TIMING
  start_time();
#endif

  readInstance(argv[1]);

#ifdef TIMING
  printf ("Time to read instance = %4.2f\n\n", run_time());
#endif

  if (argc >= 3 && !isdigit(argv[2][0])) {
    tourfp = fopen (argv[2], "w");
    if (tourfp == 0) 
      fprintf (stderr,
	       "Error opening output tour file.  Tour won't be saved.\n");
  }
  else tourfp = 0;

  if (argc >= 3 && isdigit(argv[2][0])) {
    iterations = atoi(argv[2]);
    fprintf (stderr, "Beginning %d iterations.\n", iterations);
  }
  else iterations = 1;

  sprand(12345678); /* randomness is used only for construction of kd-tree */
  doGreedy();
  
  if (iterations > 1) {
    printf ("Doing more iterations...\n");
    freopen ("/dev/null", "w", stdout);   

    for (i=1; i<iterations; ++i) {
      freeMemory();
      sprand(12345678); 
      doGreedy();
    }
  }
  return 0;
}
