// SVN date:    $Date$
// SVN author:  $Author$
// SVN rev.:    $Revision$
// SVN URL:     $HeadURL$
// SVN ID:      $Id$

/*
 *  errorEstimation.c
 *  GenCurvefit
 * 
 *	This code works out the covariance matrix, and therefore the fit errors for the genetic curvefit.  It does so by a matrix method
 *	THis matrix method does not affect the fit parameters in any way, but it is a gradient technique.
 *
 *  Created by andrew on 24/09/07.
 *  Copyright 2007 __Andrew Nelson and The Australian Nuclear Science and Technology Organisation__. All rights reserved.
 *
 */
#include "XOPStandardHeaders.h"
#include "GenCurveFit.h"
#include "errorEstimation.h"

#define TINY 1.0e-20

double factorial(double num){
	int ii;
	double result = 0;
	
	if( num<70){
		result = 1;
		for(ii = 1 ; ii < num+1 ; ii+=1){
			result *= (double)ii;
		}
	} else {
		result = sqrt(2 * 3.14159 * num) * pow(num/2.71828,num);
	}
	return result;
}

int getCovarianceMatrix(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	int err;
	double **derivativeMatrix = NULL;
	double hessianDeterminant = 0;
	int ii,jj;
	//	double temp;
	err = 0;
	
	if(goiP->covarianceMatrix == NULL){
		err = UNSPECIFIED_ERROR;
		goto done;
	}
	
	derivativeMatrix = (double**) malloc2d(goiP->numvarparams,goiP->unMaskedPoints,sizeof(double));
	if(derivativeMatrix == NULL){	err = NOMEM;	goto done;}
	
	if(err = updatePartialDerivative(derivativeMatrix, p, goiP)) goto done;
	
	if(err = updateAlpha(goiP->covarianceMatrix, derivativeMatrix, goiP)) goto done;

	if(err = matrixInversion(goiP->covarianceMatrix, goiP->numvarparams, &hessianDeterminant)) goto done;

	goiP->V_logBayes = exp(-0.5 * (*(goiP->chi2Array)) / (double)(goiP->unMaskedPoints - goiP->numvarparams));// * pow(4*3.14159,(double) goiP->numvarparams) ;//* factorial((double)goiP->numvarparams);
	goiP->V_logBayes = goiP->V_logBayes / (sqrt(hessianDeterminant));
	//	for(ii=0; ii < goiP->numvarparams ; ii+=1){
	//		temp = fabs(*(goiP->limits + *(goiP->varparams+ii) + goiP->totalnumparams)-*(goiP->limits + *(goiP->varparams+ii)));
	//		temp = temp / (0.5 * fabs(*(goiP->limits + *(goiP->varparams+ii) + goiP->totalnumparams)+(*(goiP->limits + *(goiP->varparams+ii)))));
	//		goiP->V_logBayes = goiP->V_logBayes / temp;
	//	}
	goiP->V_logBayes = log(goiP->V_logBayes);
		
	if(!p->WFlagEncountered)
		for(ii=0; ii< goiP->numvarparams ; ii++)
			for(jj=0 ; jj<goiP->numvarparams ; jj+=1)
				goiP->covarianceMatrix[ii][jj] *= *(goiP->chi2Array)/(goiP->unMaskedPoints - goiP->numvarparams);
		
done:
	if(derivativeMatrix != NULL)
		free(derivativeMatrix);
	
	return err;
}

/** Calculates the lower left elements for <code>this.alpha</code>. */
int updateAlpha(double **alpha, double **derivativeMatrix, GenCurveFitInternalsPtr goiP) {
	int err = 0, ii,jj;
	for (ii = 0; ii < goiP->numvarparams; ii++) {
		for (jj = 0; jj < ii+1 ; jj++) {
			if(err = calculateAlphaElement(ii, jj, alpha, derivativeMatrix, goiP)) return err;
		}
	}
	if(err = packAlphaSymmetric(alpha,goiP)) return err;
	return err;
}

int calculateAlphaElement(int row, int col, double **alpha, double **derivativeMatrix, GenCurveFitInternalsPtr goiP) {
	int err = 0;
	int ii;
	double result = 0;
	double num = 0;
	
	for (ii = 0; ii < goiP->unMaskedPoints ; ii++) {
		num = derivativeMatrix[row][ii]	* derivativeMatrix[col][ii];
		switch(goiP->weighttype){
			case -1:
				break;
			case 0:
				num *= ((*(goiP->dataSig+ii))*(*(goiP->dataSig+ii)));
				break;
			case 1:
				num /= ((*(goiP->dataSig+ii))*(*(goiP->dataSig+ii)));
				break;
		}
		result += num;
	}
	
	alpha[row][col] = result;
	return err;
}

/** packs the upper right elements of the alpha matrix, because the alpha matrix should be symmetrical*/
int packAlphaSymmetric(double** alpha,GenCurveFitInternalsPtr goiP){  
	int err = 0,ii,jj;
	
	for(ii=0 ; ii<goiP->numvarparams ; ii++){
		for(jj = goiP->numvarparams-1 ; jj > ii ; jj--){
			alpha[ii][jj] = alpha[jj][ii];
		}
	}
	return err;
}

int updatePartialDerivative(double **derivativeMatrix, GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	int err = 0;
	int ii;
	for(ii=0 ; ii< goiP->numvarparams ; ii++){
		if(err = partialDerivative(derivativeMatrix,ii,p,goiP,*(goiP->varparams+ii))) return err;
	} 
	return err;
}

int partialDerivative(double** derivativeMatrix,int derivativeMatrixRow, GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP,int parameterIndex){
	int err = 0;
	double param, diff;
	long indices[MAX_DIMENSIONS];
	double value[2];
	int jj;
	
	indices[0] = parameterIndex;
	if(err = MDGetNumericWavePointValue(goiP->GenCurveFitCoefs,indices,value)) return err;
	
	param = value[0];
	diff = 1.e-6*param;
	value[0] = param + diff;
	if(err = MDSetNumericWavePointValue(goiP->GenCurveFitCoefs,indices,value)) return err;
	
	if(err = calcModel(&goiP->fi,goiP->GenCurveFitCoefs,goiP->dataCalc,*(derivativeMatrix+derivativeMatrixRow),goiP->xcalc,goiP->independentVariable,goiP->numVarMD,goiP->isAAO,goiP->sp))
		return err;
	
	value[0] = param - diff;	
	if(err = MDSetNumericWavePointValue(goiP->GenCurveFitCoefs,indices,value)) return err;
	
	if(err = calcModel(&goiP->fi,goiP->GenCurveFitCoefs,goiP->dataCalc,goiP->dataTemp,goiP->xcalc,goiP->independentVariable,goiP->numVarMD,goiP->isAAO,goiP->sp))
		return err;
	
	for(jj=0; jj<goiP->unMaskedPoints ; jj++){
		derivativeMatrix[derivativeMatrixRow][jj] = (derivativeMatrix[derivativeMatrixRow][jj] - *(goiP->dataTemp+jj)) / (2*diff);
	}
	
	value[0] = param;	
	if(err = MDSetNumericWavePointValue(goiP->GenCurveFitCoefs,indices,value)) return err;
	
	return err;
}

/*
int matrixInversion(double **a, int N, double *detA){
	int err=0;
	int i,j;
	int *indx = NULL;
	double *col = NULL;
	double **tempA = NULL;
	double d;
	*detA = 1;

	
	
	indx = (int*)malloc(sizeof(int)*N);
	if(indx == NULL){
		err = NOMEM;
		goto done;
	}
	
	tempA = (double**)malloc2d(N,N,sizeof(double));
	if(tempA == NULL){
		err = NOMEM;
		goto done;
	}
	memcpy(tempA, a, sizeof(double) * N * N);
	
	col = (double*)malloc(sizeof(double)*N);
	if(col == NULL){
		err = NOMEM;
		goto done;
	}
	
	//perform the cholesky decomposition, purely to get the determinant
	if(err = choldc(tempA, N, col)) goto done;
	//the determinant of the original matrix is the square of the products of the elements in the cholesky diagonal
	for(i = 0 ; i < N ; i+=1)
		*detA *= tempA[i][i] * tempA[i][i];
	memcpy(tempA, a, sizeof(double) * N * N);

	
	if(err = ludcmp(tempA,N,indx,&d)) goto done;
	
	for(j=0 ; j<N ; j++){
		for(i=0 ; i<N ; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(tempA,N,indx,col);
		for(i=0 ; i<N ; i++){
			d = col[i];
			a[i][j] = col[i];
		};
	}
done:
	if(col!=NULL)
		free(col);
	if(indx!=NULL)
		free(indx);
	if(tempA!=NULL)
		free(tempA);
	return err;
}
*/


int matrixInversion(double **a, int N, double *detA){
	int err=0;
	int i,j;
	double *x = NULL;
	double *b = NULL;
	double *p = NULL;
	double **tempA = NULL;
	*detA = 1;
	
	x = (double*)malloc(sizeof(double) * N);
	if(x == NULL){
		err = NOMEM;
		goto done;
	}
	
	tempA = (double**)malloc2d(N,N,sizeof(double));
	if(tempA == NULL){
		err = NOMEM;
		goto done;
	}
	memcpy(tempA, a, sizeof(double) * N * N);
	
	p = (double*)malloc(sizeof(double) * N);
	if(p == NULL){
		err = NOMEM;
		goto done;
	}
	memset(p, 0, sizeof(double) * N);
	
	b = (double*)malloc(sizeof(double) * N);
	if(b == NULL){
		err = NOMEM;
		goto done;
	}
	memset(b, 0, sizeof(double) * N);
	
	//perform the cholesky decomposition
	if(err = choldc(tempA, N, p)) goto done;
	
	//now do the back substitution
	for(j = 0 ; j < N ; j++){
		memset(b, 0, sizeof(double) * N);
		b[j] = 1.0;
		memset(x, 0, sizeof(double) * N);
		cholsl(tempA, N, p, b, x);
		
		for(i = 0 ; i < N ; i++)
			a[i][j] = x[i];
	}
	
	//make the covariance matrix symmetric
	for(i=0 ; i < N ; i++)
		for(j = N-1 ; j > i ; j--)
			a[i][j] = a[j][i];
		
	
	//the determinant of the original matrix is the square of the products of the elements in the cholesky diagonal
	for(i = 0 ; i < N ; i+=1)
		*detA *= tempA[i][i] * tempA[i][i];
	
done:
	if(p)
		free(p);
	if(b)
		free(b);
	if(x)
		free(x);
	if(tempA)
		free(tempA);
	return err;
}


//Cholesky Decomposition
static int choldc (double **a, int N, double *p){
	int err = 0;	
	int ii, jj, kk;
	double sum = 0;
	
	
	for(ii = 0; ii < N ; ii++){
		for(jj = ii ; jj < N ; jj++){
			for(sum = a[ii][jj] , kk = ii - 1 ; kk >= 0; kk--) sum -= a[ii][kk] * a[jj][kk];
			if(ii == jj){
				if(sum <= 0.0)
					return UNSPECIFIED_ERROR;
				p[ii] = sqrt(sum);
			} else a[jj][ii] = sum/p[ii];
		}
	}
	
	
done:
	return err;
}

//Cholesky back substitution
static void cholsl(double **a, int N, const double *p, double *b, double *x){
	int ii, kk;
	double sum = 0;
	for(ii = 0 ; ii < N ; ii++){
		for(sum = b[ii], kk = ii-1 ; kk >=0 ; kk-- ) sum -=a[ii][kk] * x[kk];
		x[ii] = sum/p[ii];
	}
	for(ii = N-1 ; ii >= 0 ; ii--){
		for(sum = x[ii], kk = ii+1 ; kk < N ; kk++ ) sum -=a[kk][ii] * x[kk];
		x[ii] = sum/p[ii];
	}
}

static int ludcmp(double **a, int n, int *indx, double *d){
	int i, imax, j, k, err = 0;
	double big, dum, sum, temp;
	double *vv = NULL;
	
	vv = (double*)malloc(sizeof(double)*n);
	if(vv == NULL){
		err = 0;
		goto done;
	}
	
	*d = 1.0;
	
	for(i=0 ; i<n ; i++){
		big = 0.0;
		for(j=0 ; j<n ; j++)
			if((temp = fabs(a[i][j])) > big) big = temp;
		if(big == 0.0){
			err = UNSPECIFIED_ERROR;
			goto done;
		}
		vv[i] = 1.0/big;	
	} 
	
	for(j=0 ; j<n ; j++){
		for(i=0 ; i<j ; i++){
			sum = a[i][j];
			for (k=0 ; k<i ; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for(i=j ; i<n ; i++){
			sum = a[i][j];
			for(k=0; k<j ; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if( (dum=vv[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if(j != imax){
			for(k=0 ; k<n ; k++){
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if(a[j][j] == 0.0) a[j][j] = TINY;
		
		if(j != n-1){
			dum = 1.0/(a[j][j]);
			for(i=j+1; i<n ; i++) a[i][j] *=dum;
		}
	}
	
	
done:
	if(vv != NULL)
		free(vv);
	
	return err;
	
}

static void lubksb(double **a, int n, int *indx, double b[]){
	int i, ii=0, ip, j;
	double sum;
	
	for(i=1 ; i<=n ; i++){
		ip = indx[i-1];
		sum = b[ip];
		b[ip] = b[i-1];
		if(ii)
			for(j=ii ; j<=i-1 ; j++) sum -= a[i-1][j-1]*b[j-1];
		else if (sum) ii=i;
		b[i-1] = sum;
	}
	for(i=n ; i>=1 ; i--){
		sum = b[i-1];
		for(j=i+1 ; j<=n ; j++) sum -= a[i-1][j-1]*b[j-1];
		b[i-1] = sum/a[i-1][i-1];
	}	
}

