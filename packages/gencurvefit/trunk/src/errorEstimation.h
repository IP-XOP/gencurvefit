

/*
 functions contained in errorEstimation.c
 */
static int ludcmp(double**,int,int*,double*);
static int lubksb(double **a, int n, int *indx, double b[]);
static int matrixInversion(double **a, int N);
int getCovarianceMatrix(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP);
int updatePartialDerivative(double**, GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP);
int partialDerivative(double**, int, GenCurveFitRuntimeParamsPtr, GenCurveFitInternalsPtr,int);
int updateAlpha(double**,double**, GenCurveFitInternalsPtr goiP);
int calculateAlphaElement(int row, int col, double **alpha, double **derivativeMatrix, GenCurveFitInternalsPtr goiP);
int packAlphaSymmetric(double** alpha,GenCurveFitInternalsPtr);