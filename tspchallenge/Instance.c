#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "Instance.h"
#include "Utilities.h"

/*********************************************************************

Eighth Implementation Challenge

This file contains functions to read TSP instances in some of the
formats used in TSPLIB.  All formats currently used in TSPLIB TSP
instances with more than 1000 cities are supported, including

	2-D Euclidean instances:  EDGE_WEIGHT_TYPE = EUC_2D

	2-D Euclidean instances with ceilings used in distance
		calculations:  EDGE_WEIGHT_TYPE = CEIL_2D

	Upper triangular distance matrices (including the diagonal):
		EDGE_WEIGHT_TYPE = EXPLICIT, 
		EDGE_WEIGHT_FORMAT = UPPER_DIAG_ROW


This code is based on previous work by David Johnson, Jon Bentley, and 
others.

						Lyle McGeoch
						lam@cs.amherst.edu

						June 24, 1999
 *********************************************************************/

char *filename;
int  ncities;
char *instanceName;

int  (*distance)(int, int);  /* pointer to function returning distance */
  
int  isEuclidean;      
int  notIntegral;      /* x,y coordinates include non-integral values */

double  *xx, *yy;      /* used only for Euclidean instances */
int  **distMat;        /* used only for matrix instances */

/*************************************************************************/

static FILE *fp;

#define BUFSIZE 512
typedef char mybuf[BUFSIZE];

static int linenum;
static mybuf buffer;         /* used for reading input lines */
static mybuf msg;            /* used for error messages */
static mybuf edgeWeightType;
static mybuf edgeWeightFormat;
static mybuf name;

/*************************************************************************/

static int roundDistance (int a, int b) {
  double xd = xx[a] - xx[b];
  double yd = yy[a] - yy[b];
  double r  = sqrt(xd*xd + yd*yd) + 0.5;
  return (int)r;
}

/*************************************************************************/

static int ceilDistance (int a, int b) {
  double xd = xx[a] - xx[b];
  double yd = yy[a] - yy[b];
  double r  = ceil(sqrt(xd*xd + yd*yd));
  return (int)r;
}

/*************************************************************************/

static int matDistance (int a, int b) {
  return distMat[a][b];
}

/*************************************************************************/

static int isWhiteSpace(char c) {
  return (c==' ') || (c=='\t');
}

/*************************************************************************/

static char *getline (void) {
  size_t length;
  size_t i;

  char *result = fgets (buffer, BUFSIZE, fp);
  if (result == 0) return result;

  length = strlen (result);
  if (length > 0 && result[length-1] == '\n') {
    result[--length] = 0;
  }
  else fatal ("Buffer size in getline() isn't large enough");
  
  for (i=0; i<length; ++i) 
    if (!isWhiteSpace(result[i])) return result;
  return getline();
}

/*************************************************************************/

static char *match (char *source, char *prefix) {

  size_t prefixLength = strlen (prefix);
  
  if (strncmp(source, prefix, prefixLength) != 0) return 0;
  source += prefixLength;
  while (isWhiteSpace(*source)) source++;
  if (*source == ':') source++;
  while (isWhiteSpace(*source)) source++;
  
  return source;
}

/*************************************************************************/

static int getInt (char *s) {

  int i;
  int result = sscanf (s, "%d", &i);
  
  if (result != 1) {
    sprintf (msg, "Error reading integer on line %d", linenum);
    fatal(msg);
  }
  return i;
}

/*************************************************************************/

static void readEdgeWeightSection(void) {

  int val;
  int i,j;

  if (match (edgeWeightType, "EXPLICIT")) {
    if (match (edgeWeightFormat, "UPPER_DIAG_ROW")) {
      
      distMat = (int **) createMatrix(ncities+1,ncities+1,sizeof(int));
      
      for (i=1; i<=ncities; ++i) 
	for (j=i; j<=ncities; ++j) {
	  fscanf (fp, "%d", &val);
	  distMat[i][j] = distMat[j][i] = val;
	}
      isEuclidean = 0;
      distance = matDistance;
      return;
    }
  }

  sprintf (msg, "Program doesn't support edge weights with type %s and format %s",
	   edgeWeightType, edgeWeightFormat);
  fatal(msg);
}


/*************************************************************************/

static int readInts(FILE *fp) {

  int d1, d2;
  int i,j;
  int result;

  for (i=1; i<=ncities; ++i) {

    result = fscanf (fp, "%d %d %d", &j, &d1, &d2);
    if (result != 3) return 0;

    if ((j<1) || (j >ncities))
      fatal ("City number is out of bounds");
    xx[j] = d1;
    yy[j] = d2;
  }

  return 1;
}
    
/*************************************************************************/

static void readDoubles(FILE *fp) {

  double d1, d2;
  int i,j;

  for (i=1; i<=ncities; ++i) {

    fscanf (fp, "%d %lf %lf", &j, &d1, &d2);

    if ((j<1) || (j >ncities))
      fatal ("City number is out of bounds");
    xx[j] = d1;
    yy[j] = d2;
    if (d1 != (int)d1) notIntegral=1;
    if (d2 != (int)d2) notIntegral=1;
  }
}
    
/*************************************************************************/

static void readNodeCoordSection(void) {
  
  int filePosition;

  if (match (edgeWeightType, "EUC_2D") ||
      match (edgeWeightType, "CEIL_2D")) {
    xx = (double *) calloc (ncities+1, sizeof(double));
    yy = (double *) calloc (ncities+1, sizeof(double));

    notIntegral = 0;

    filePosition = ftell(fp);

    if (!readInts(fp)) {
      fseek (fp, filePosition, SEEK_SET);
      readDoubles(fp);
    }

    isEuclidean = 1;
    if (match (edgeWeightType, "CEIL_2D")) 
      distance = ceilDistance;
    else
      distance = roundDistance;
    return;
  }
  else {
    
    sprintf (msg, "Program doesn't support node coordinates with type %s",
	     edgeWeightType);
    fatal(msg);
    return;
  }
}

/*************************************************************************/

void readInstance (char *s) {

  char *line;
  char *trailer;

  mybuf type;
  mybuf edgeDataFormat;
  mybuf nodeCoordType;
  mybuf displayDataType;

  instanceName = name;
  filename = s;
  fp = fopen (filename, "r");
  if (fp == NULL) fatal ("Couldn't open input file.\n");
  linenum = 0;

  ncities = 0;
  strcpy (instanceName, "NONE");
  strcpy (type, "NONE");
  strcpy (edgeWeightType, "NONE");
  strcpy (edgeDataFormat, "NONE");
  strcpy (nodeCoordType, "NONE");
  strcpy (displayDataType, "NONE");

  while ((line = getline())) {
    linenum++;
    
    if ((trailer = match (line, "NAME")))
      strcpy (instanceName, trailer);
    else if ((trailer = match (line, "TYPE")))
      strcpy (type, trailer);
    else if ((trailer = match (line, "COMMENT"))) {}
    else if ((trailer = match (line, "DIMENSION")))
      ncities = getInt (trailer);
    else if ((trailer = match (line, "EDGE_WEIGHT_TYPE")))
      strcpy (edgeWeightType, trailer);
    else if ((trailer = match (line, "EDGE_WEIGHT_FORMAT")))
      strcpy (edgeWeightFormat, trailer);
    else if ((trailer = match (line, "EDGE_DATA_FORMAT")))
      strcpy (edgeDataFormat, trailer);
    else if ((trailer = match (line, "NODE_COORD_TYPE")))
      strcpy (nodeCoordType, trailer);
    else if ((trailer = match (line, "DISPLAY_DATA_TYPE")))
      strcpy (displayDataType, trailer);
    else if ((trailer = match (line, "EOF"))) break;
    else if ((trailer = match (line, "EDGE_WEIGHT_SECTION")))
      readEdgeWeightSection();
    else if ((trailer = match (line, "NODE_COORD_SECTION")))
      readNodeCoordSection();
    else {
      sprintf (msg, "Unexpected input line: %s", line);
      fatal(msg);
    }
  }
  if (match (type, "TSP") == 0) {
    sprintf (msg, "This message has type %s, not TSP", type);
    fatal (msg);
  }

  if (distance == 0)
    fatal ("No matrix or coordinate information in input file");
}
