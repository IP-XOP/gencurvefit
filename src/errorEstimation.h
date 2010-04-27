

/*
 functions contained in errorEstimation.c
 */
static int matrixInversion(double **a, int N, double *detA);
int getCovarianceMatrix(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP);
int updatePartialDerivative(double**, GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP);
int partialDerivative(double**, int, GenCurveFitRuntimeParamsPtr, GenCurveFitInternalsPtr,int);
int updateAlpha(double**,double**, GenCurveFitInternalsPtr goiP);
int calculateAlphaElement(int row, int col, double **alpha, double **derivativeMatrix, GenCurveFitInternalsPtr goiP);
int packAlphaSymmetric(double** alpha,GenCurveFitInternalsPtr);
static int choldc (double **a, int N, double *p);
static void cholsl(const double **a, int N, const double *p, double *b, double *x);
