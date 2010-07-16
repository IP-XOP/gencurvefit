// SVN date:    $Date$
// SVN author:  $Author$
// SVN rev.:    $Revision$
// SVN URL:     $HeadURL$
// SVN ID:      $Id$

/*	GenCurveFit.c -- An XOP for curvefitting via Differential Evolution.
 See:
 Wormington, et. al., "Characterisation of structures from X-ray Scattering
 Data using Genetic Algorithms", Phil. Trans. R. Soc. Lond. A (1999) 357, 2827-2848
 
 And
 The Motofit packages: http://motofit.sourceforge.net/
 Nelson, "Co-refinement of multiple-contrast neutron/X-ray reflectivity data using Motofit",
 J. Appl. Cryst. (2006). 39,273-276.
 
 @copyright: Andrew Nelson and the Australian Nuclear Science and Technology Organisation 2007.
 
 */

/*
 put all the standard headers in GenCurveFit.h, including those from wavemetrics
 */
#include "XOPStandardHeaders.h"
#include "GenCurveFit.h"
#include "updateXOP.h"
#include "errorEstimation.h"
#include "dSFMT.h"


//gTheWindow is a window created to show the latest position of the fit
XOP_WINDOW_REF gTheWindow = NULL;

//a variable which sees if the fit is to be aborted.
int Abort_the_fit = 0;

/*
 ExecuteGenCurveFit performs the genetic curvefitting routines
 returns 0 if no error
 returns errorcode otherwise
 */
int
ExecuteGenCurveFit(GenCurveFitRuntimeParamsPtr p)
{
	/* the function that is called by IGOR. Here's where we farm the work out to different places.*/
	/*
	 GenCurveFitInternals carries the internal data structures for doing the fitting.
	 */
	GenCurveFitInternals goi;
	//initialise the structure to be zero
	
	/*
	 err carries the errors for all the operations
	 err2 carries the error code if we still want to return err, but we need to finish off
	 something else first.	
	 */
	int err = 0, err2=0;
	double value[2];
	long indices[MAX_DIMENSIONS];
	//variables listed below are purely for outputs sake.
	char varname[MAX_OBJ_NAME+1];
	long dimensionSizes[MAX_DIMENSIONS];
	double t1,t2;
	long lt1 = 0;
	char note[200], note_buffer1[MAX_WAVE_NAME+1], note_buffer2[MAX_WAVE_NAME+1], cmd[MAXCMDLEN+1];
	int output,ii,jj,isDisplayed, quiet, generateCovariance;
	uint32_t seed;
	
	//initialise all the internal data structures to NULL
	memset(&goi, 0, sizeof(goi));
	
	if( igorVersion < 503 )
		return REQUIRES_IGOR_500;
	
	//reset the abort condition
	Abort_the_fit = 0;
	
	strncpy(varname, "V_Fiterror", MAX_OBJ_NAME);
	if(FetchNumVar(varname, &t1, &t2)!=-1){
		if(!err)
			lt1 = 0;
		else 
			lt1 = 1;
		
		err = 0;
		if(err2 = SetIgorIntVar(varname, lt1, 1)){err = err2;goto done;};
	}	
	
	/*
	 Genetic Optimisation uses a lot of random numbers, we are seeding the generator here.
	 If you want to seed the generator we can 
	 */
	if(p->SEEDFlagEncountered)
		seed =  (uint32_t) p->SEEDFlag_seed;
	else 
		seed = (uint32_t) time(NULL);
	randomInteger(0, seed, 1);
	randomDouble(0, 0, seed , 1);
	
	/*
	 checkInput checks the input that IGOR sends to the XOP.  If everything is correct it returns 0, 
	 else it returns an error.  Errors can be caused by, e.g. different x and y wave lengths, etc.
	 */
	if(err = checkInput(p, &goi))
		goto done;
	
	/*
	 init_GenCurveFitInternals sets up the internal data arrays in the GenCurveFitInternals structure.  It holds
	 copies of the datapoints being fitted, arrays of the fitting structures, etc.  It calls malloc
	 several times to set aside memory for all this.  If this procedure works without a hitch then the function
	 returns 0.
	 */
	if(err = init_GenCurveFitInternals(p, &goi))
		goto done;
	
	/*
	 optimiseloop does the Differential Evolution, according to Storn and Price.  When this returns 0, then 
	 you have the best fit in the GenCurveFitInternals structure.  Otherwise it returns an error code.  If the user aborts
	 then the FIT_ABORTED error code is returned, but it is still possible to retrieve the best solution so far
	 */
	err = optimiseloop(&goi, p);
	
	/*
	 if there are no errors, or if the user aborted, then return the best fit.
	 If the data is displayed in the top graph append the fitcurve to the top graph
	 */
	if(err == 0 || err == FIT_ABORTED){
		if(err2 = ReturnFit( &goi,  p))
		{err = err2;goto done;}
		
		//return an error wave
		//make an error wave
		dimensionSizes[0] = goi.totalnumparams;
		dimensionSizes[1] = 0;
		if(err2 = MDMakeWave(&goi.W_sigma, "W_sigma",goi.cDF,dimensionSizes,NT_FP64, 1))
		{err = err2;goto done;}
		
		//set the error wave to zero
		for(ii=0; ii<goi.totalnumparams; ii+=1){
			indices[0] = ii;
			value[0] = 0;
			if(err2 = MDSetNumericWavePointValue(goi.W_sigma,indices,value))
			{err = err2;goto done;};					 
		}
		
		//do the errors
		generateCovariance = 0;
		if(p->MATFlagEncountered){
			generateCovariance = 1;
			if(p->MATFlagParamsSet[0] && (int) p->MATFlag_mat == 0)
				generateCovariance = 0;
		}
				
		goi.covarianceMatrix = (double**)malloc2d(goi.numvarparams, goi.numvarparams, sizeof(double));
		if(goi.covarianceMatrix == NULL)
		{ err = NOMEM; goto done;}
		
		if(generateCovariance && (!(err2 = getCovarianceMatrix(p, &goi)))){
			//set the error wave
			for(ii = 0; ii < goi.numvarparams ; ii += 1){
				indices[0] = *(goi.varparams+ii);
				value[0] = sqrt(goi.covarianceMatrix[ii][ii]);
				if(err2 = MDSetNumericWavePointValue(goi.W_sigma, indices, value))
					{err = err2;goto done;};
			}					 
			WaveHandleModified(goi.W_sigma);
			
			if(generateCovariance){
				//make the covariance matrix
				dimensionSizes[0] = goi.totalnumparams;
				dimensionSizes[1] = goi.totalnumparams;
				dimensionSizes[2] = 0;
				if(err2 = MDMakeWave(&goi.M_covariance, "M_Covar", goi.cDF, dimensionSizes, NT_FP64, 1))
					{err = err2; goto done;};
			
				for(ii = 0; ii < goi.totalnumparams; ii+=1){
					for(jj = 0 ; jj < goi.totalnumparams; jj += 1){
						indices[0] = ii;
						indices[1] = jj;
						value[0] = 0;
						if(err2 = MDSetNumericWavePointValue(goi.M_covariance, indices, value))
							{err = err2;goto done;};  
					}
				}
				
				for(ii = 0; ii < goi.numvarparams ; ii += 1){
					for(jj = 0; jj < goi.numvarparams ; jj += 1){
						indices[0] = *(goi.varparams + ii);
						indices[1] = *(goi.varparams + jj);
						value[0] = goi.covarianceMatrix[ii][jj];
						if(err2 = MDSetNumericWavePointValue(goi.M_covariance, indices, value))
							{err = err2;goto done;};
					}
				}
				
				WaveHandleModified(goi.M_covariance);
			}
		}
		
		
		if(err2 = isWaveDisplayed(p->dataWave.waveH, &isDisplayed))
		{err = err2;goto done;};
		if(isDisplayed && goi.numVarMD == 1){
			if(err2 = isWaveDisplayed(goi.OUT_data, &isDisplayed))
				{err = err2 ; goto done;};
			
			if(!isDisplayed){
				strncpy(cmd, "appendtograph/w=$(winname(0,1)) ", MAXCMDLEN);
				WaveName(goi.OUT_data,&note_buffer1[0]);
				strncat(cmd, &note_buffer1[0], MAXCMDLEN - strlen(note_buffer1));
				
				if(p->DFlagEncountered && p->XFlagEncountered){
					WaveName(p->XFlag_xx,&note_buffer2[0]);
					strncat(cmd, " vs ", MAXCMDLEN - strlen(cmd) - strlen(" vs "));
					strncat(cmd, note_buffer2, MAXCMDLEN - strlen(cmd) - strlen(note_buffer2) );
				}
				if(err2 = XOPSilentCommand(&cmd[0]))
					{err = err2; goto done;}
			}
		}
	}
	
	/*this section sets history and global variables
	 V_Chisq
	 V_Fiterror
	 If there are no errors returned above, then V_Chisq is produced in the current datafolder.
	 V_fiterror is for if there is a problem with the fitting.  
	 If this global variable exists in IGOR at runtime, and there is no XOPerror then V_fiterror = 0.  
	 If there is an XOPerror, bit 0 of V_fiterror is set. However, the XOP will return gracefully to IGOR,
	 allowing the user that called it to detect this and carry on.
	 In normal situations V_fiterror may not exist, so the GenCurveFit XOP will return an error message.
	 */
	
	if(!err){
		SetOperationNumVar("V_Chisq", *(goi.chi2Array));
		SetOperationNumVar("V_fitIters", (int)(goi.V_numfititers));
		SetOperationNumVar("V_npnts", (int)(goi.unMaskedPoints));
		SetOperationNumVar("V_nterms", (int)WavePoints(p->coefs));
		SetOperationNumVar("V_nheld", (int)(WavePoints(p->coefs) - goi.numvarparams));
		SetOperationNumVar("V_logBayes", goi.V_logBayes);
	}		
	
	strcpy(varname, "V_Fiterror");
	if(FetchNumVar(varname, &t1, &t2)!=-1){
		if(!err){
			lt1 = 0;
		} else {
			lt1 = 1;
		}
		err = 0;
		if(err2 = SetIgorIntVar(varname, lt1, 1))
		{err = err2;goto done;}
	}
	/*
	 This section puts a copy of the fit parameters into the history area, unless one sets quiet mode.
	 */
	quiet = 0;
	if(p->QFlagEncountered){
		quiet = 1;
		if(p->QFlagParamsSet[0] && (int) p->QFlag_quiet == 0)
			quiet = 0;
	}
	
	if(!quiet && (!err || err == FIT_ABORTED) && lt1==0 ){		
		if(!err)
			{output = snprintf(note, 199, "_______________________________\rGenetic Optimisation Successful\r");XOPNotice(note);}
		if(err == FIT_ABORTED)
			{output = snprintf(note, 199, "_______________________________\rGenetic Optimisation ABORTED\r");XOPNotice(note);}
		WaveName(p->dataWave.waveH, note_buffer1);
		
		output = snprintf(note, 
						  199,
						  "Fitting: %s to %s\r",
						  note_buffer1,
						  goi.fi.name);
		
		XOPNotice(note);
		
		output = snprintf(note,
						  199,
						  "V_fitIters = %li; V_Chisq = %g; V_npnts= %li; V_nterms= %li; V_nheld= %li; V_logBayes = %g\r",
						  goi.V_numfititers,
						  *(goi.chi2Array),
						  goi.unMaskedPoints,
						  WavePoints(p->coefs),
						  WavePoints(p->coefs) - goi.numvarparams,
						  goi.V_logBayes);
		
		XOPNotice(note);
		
		for(ii = 0; ii < WavePoints(p->coefs) ; ii += 1){
			indices[0] = ii;
			indices[1] = 0;
			if(err = MDGetNumericWavePointValue(goi.W_sigma, indices, value))
				goto done;
			
			output = snprintf(note,
							  199,
							  "\tw[%d]\t=\t%g   +/-   %g\r",
							  ii,
							  *(goi.gen_bestfitsofar+ii),
							  value[0]);
			
			XOPNotice(note);
		}
		output = snprintf(note, 199, "_______________________________\r");
		XOPNotice(note);
	}
	/*
	 freeAllocMem frees all the internal data structures which have had memory allocated to them.
	 this is ultra essential for no memory leaks.
	 */
done:
	
	freeAllocMem(&goi);
	return err;
}


int
//this will return if a wave contains NaN or INF.
checkNanInf(waveHndl wav){
	//this check examines to see if there are any NaN/Infs in a wave
	//this is really important if you want to calculate Chi2.
	int err = 0;
	long ii;
	double *dp = NULL;
	long points;
	if(wav == NULL)
		return NON_EXISTENT_WAVE;
	points = WavePoints(wav);
	
	dp = (double*) malloc(sizeof(double) * WavePoints(wav));
	if(dp == NULL)
		return NOMEM;
	
	if(err = MDGetDPDataFromNumericWave(wav,dp))
		goto done;
	
	for(ii = 0 ; ii < points ; ii += 1)
		if((err = IsNaN64(dp + ii)) || (err = IsINF64(dp + ii))) break;
	
	if(err > 0)
		err = INPUT_WAVES_CONTAINS_NANINF;
done:
	if(dp != NULL)
		free(dp);
	return err;
}


/*
 checkZeros sees how many zero points are in the wave
 returns 0 if no error
 returns errorcode otherwise
 */
int
checkZeros(waveHndl wavH,long* numzeros){
	int result = 0;
	long numBytes;
	double* dp = NULL;
	long points,ii;
	points = WavePoints(wavH);
	
	numBytes = points * sizeof(double); // Bytes needed for copy
	dp = (double*)malloc(numBytes);
	
	if (dp==NULL)
		return NOMEM;
	
	if (result = MDGetDPDataFromNumericWave(wavH, dp)) { // Get copy.
		free(dp);
		return result;
	}
	
	for(ii=0 ; ii<points ; ii+=1){
		if(*(dp+ii)==0)
			*numzeros+=1;
	}
	if(dp != NULL)
		free(dp);
	return result;
}

/*
 calcMaxLikelihood calculates Maximum Likelihood as a cost function for the fit
 dataObs		-	array containing the experimental data
 dataTemp		-	array containing theoretical model
 ysig		-	array containing the weighting
 len			-	the number of fit points
 *chi2		-	the calculated chi2 value
 weighttype	-	0 if ysig is weightarray, 1 if ysig is standard deviation
 returns 0 if no error
 returns errorcode otherwise
 */
static int
calcMaxLikelihood(double* dataObs, double* dataTemp, double* dataSig, long len, double* chi2, int weighttype){
	
	int err = 0;
	long ii;
	double abserr=0;
	if(dataObs == NULL || dataTemp == NULL || dataSig  == NULL)
		return UNSPECIFIED_ERROR;
	
	*chi2 = 0;
	
	for(ii=0 ; ii<len ; ii+=1){
		if(*(dataObs+ii) ==0 || *(dataTemp+ii)==0){
			abserr = *(dataTemp+ii) - *(dataObs+ii);
		} else {
			abserr = *(dataTemp+ii)-*(dataObs+ii)+(*(dataObs+ii)* log((*(dataObs+ii)/(*(dataTemp+ii)))));
		}
		*chi2 += abserr;
	}
	*chi2 *= 2;
	return err;
}


/*
 calcChi2 calculates chi2 for the fit
 dataObs		-	array containing the experimental data
 dataTemp		-	array containing theoretical model
 ysig		-	array containing the weighting
 len			-	the number of fit points
 *chi2		-	the calculated chi2 value
 weighttype	-	0 if ysig is weightarray, 1 if ysig is standard deviation
 returns 0 if no error
 returns errorcode otherwise
 */
static int
calcChi2(const double* dataObs, const double * dataTemp, const double* dataSig, long len, double* chi2, int weighttype){
	
	int err = 0;
	long ii;
	double abserr=0;
	if(dataObs == NULL || dataTemp == NULL || dataSig  == NULL)
		return UNSPECIFIED_ERROR;
	
	*chi2 = 0;
	
	switch(weighttype){
		case -1:
			for(ii=0 ; ii<len ; ii+=1){
				abserr = (*(dataObs+ii)-*(dataTemp+ii));
				*chi2 += abserr*abserr;
			}
			break;
		case 0:
			for(ii=0 ; ii<len ; ii+=1){
				abserr = (*(dataObs+ii)-*(dataTemp+ii));
				abserr *= (*(dataSig+ii));
				*chi2 += abserr*abserr;
			}		
			break;
		case 1:
			for(ii=0 ; ii<len ; ii+=1){
				abserr = (*(dataObs+ii)-*(dataTemp+ii));
				abserr /= (*(dataSig+ii));
				*chi2 += abserr*abserr;
			}
			break;
	}
	return err;
}


/*
 calcRobust calculates absolute errors for the fit
 dataObs		-	array containing the experimental data
 dataTemp		-	array containing theoretical model
 ysig		-	array containing the weighting
 len			-	the number of fit points
 *chi2		-	the calculated chi2 value
 weighttype	-	0 if ysig is weightarray, 1 if ysig is standard deviation
 returns 0 if no error
 returns errorcode otherwise
 */
static int
calcRobust(const double* dataObs, const double* dataTemp, const double* dataSig, long len, double* chi2, const int weighttype){
	
	int err = 0;
	long ii;
	double abserr=0;
	if(dataObs == NULL || dataTemp == NULL || dataSig  == NULL)
		return UNSPECIFIED_ERROR;
	
	*chi2 = 0;
	
	for(ii=0 ; ii<len ; ii+=1){
		abserr = (*(dataObs+ii)-*(dataTemp+ii));
		switch(weighttype){
			case 0:
				abserr *= (*(dataSig+ii));
				break;
			case 1:
				abserr /= (*(dataSig+ii));
				break;
		}
		*chi2 += fabs(abserr);
	}
	return err;
}

/*
 calcUserCostFunc calculates absolute errors for the fit
 minf		- user specified cost function
 yobs		-	wave to be sent to user cost function, contains observed data
 dataObs		-	array containing the original data
 dataCalc		-	wave containing the calculated model
 sobs			-	wave to be sent to user cost function, contains weight values
 dataSig		-	array containing the original weight values
 unMaskedPoints -	the number of data points being fitted
 gen_coefscopy  -	wave containing the parameters for that guess
 *chi2		-	the calculated chi2 value
 returns 0 if no error
 returns errorcode otherwise
 */
int calcUserCostFunc(FunctionInfo minf, 
					 waveHndl yobs,
					 const double *dataObs,
					 waveHndl dataCalc,
					 waveHndl sobs,
					 const double *dataSig,
					 long unMaskedPoints,
					 waveHndl GenCoefsCopy,
					 double *costVal){
	int err = 0;
	costFunc userCostFunc;
	long dataOffset;
	long numparams;
	double ret;

	numparams = WavePoints(GenCoefsCopy);
	userCostFunc.coefs = GenCoefsCopy;
	userCostFunc.yobs = yobs;
	userCostFunc.ycalc = dataCalc;
	userCostFunc.sobs = sobs;
	
	//copy the original data into the yobs wave created for the purpose
	//some sneaky users probably try to change it.
	if (err = MDAccessNumericWaveData(yobs, kMDWaveAccessMode0, &dataOffset)) 
		goto done; 
	memcpy((*yobs) + dataOffset, dataObs, unMaskedPoints*sizeof(double));
	
	//copy the original data into the yobs wave created for the purpose
	//some sneaky users probably try to change it.
	if (err = MDAccessNumericWaveData(sobs, kMDWaveAccessMode0, &dataOffset)) 
		goto done; 
	memcpy((*sobs) + dataOffset, dataSig, unMaskedPoints*sizeof(double));

	if(err = CallFunction(&minf, (void*) &userCostFunc, &ret))
		goto done;

	if(WavePoints(yobs) != unMaskedPoints || WavePoints(sobs) != unMaskedPoints || WavePoints(GenCoefsCopy)!= numparams){
		err = COSTFUNC_WAVES_CHANGED;
		goto done;
	}
	*costVal = ret;
	
	if(IsNaN64(&ret) || IsINF64(&ret)){
		err = COSTFUNC_DOESNT_RETURN_NUMBER;
		goto done;
	}

	
done:
	return err;
}



/*
 init_GenCurveFitInternals initialises the GenCurveFitInternals structure
 returns 0 if no error
 returns errorcode otherwise
 */
static int
init_GenCurveFitInternals(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	int err = 0, minpos,maxpos;
	long len;
	char *holdstr = NULL;
	long ii,jj,kk;
	char comparison[2];
	long dimensionSizes[MAX_DIMENSIONS+1],numdimensions;
	double value[2];
	long indices[MAX_DIMENSIONS];
	double bot,top;
	double chi2;
	long timeOutTicks=0;
	char xwavename[MAX_WAVE_NAME+1];
	char datawavename[MAX_WAVE_NAME+1];
	char reswavename[MAX_WAVE_NAME+1];
	char datawavestring[MAX_WAVE_NAME+1];
	char cmd[MAXCMDLEN];
	char letter[3];
	int toDisplay=0;
	int outPutType=NT_FP64;
	double temp1=0;
	long temp;
	waveStats wavStats;
	
	//do we want dynamic updates?
	goiP->noupdate = 0;
	if(p->NFlagEncountered){
		goiP->noupdate = 1;
		if(p->NFlagParamsSet[0] && (int) p->NFlag_noupdate == 0)
			goiP->noupdate = 0;
	}
	
	//initialise the chi2value
	goiP->chi2 = -1;
	
	//initialise an array to hold the parameters that are going to be varied.
	goiP->varparams = (int*) malloc(goiP->numvarparams*sizeof(int));
	//the total number of vectors in the population
	goiP->totalpopsize = goiP->numvarparams * (int)p->KFlag_popsize;
	
	if(goiP->varparams == NULL){
		err = NOMEM;
		goto done;
	}
	
	/*
	 Following section works out which parameters are going to vary
	 */
	len = GetHandleSize(p->holdstring);
	holdstr = (char*)malloc((len+1)*sizeof(char));
	if(holdstr == NULL){
		err = NOMEM;
		goto done;
	}
	comparison[1] = '\0';
	if(err = GetCStringFromHandle(p->holdstring,holdstr,len))
		goto done;
	jj=0;
	for(ii = 0L; ii<len ; ii++){
		comparison[0] = *(holdstr+ii);
		if(atoi(comparison)==0){
			*(goiP->varparams+jj)=ii;
			jj+=1;
		}
	}
	
	/*
	 goiP->temp is a utility array the same size as the input data
	 this needs to be specified at the top of the function
	 */
	goiP->temp = (double*)malloc(WavePoints(p->dataWave.waveH) * sizeof(double));
	if (goiP->temp == NULL){
		err = NOMEM;
		goto done;
	}
	
	/* get a full copy of the datawave */
	goiP->dataObsFull = (double*)malloc(WavePoints(p->dataWave.waveH) *sizeof(double));
	if (goiP->dataObsFull == NULL){
		err = NOMEM;
		goto done;
	}
	//get a copy of all the datawave.  This is so we can fill dataobs
	if (err = MDGetDPDataFromNumericWave(p->dataWave.waveH, goiP->dataObsFull)) // Get copy.
		goto done;
	
	/*
	 goiP->mask contains an array copy of the mask wave
	 */
	goiP->mask = (double*)malloc(WavePoints(p->dataWave.waveH) * sizeof(double));
	if (goiP->mask == NULL){
		err = NOMEM;
		goto done;
	}
	
	/*
	 if there is a range specified then we're fitting a subset of ywave
	 */
	if(p->dataWave.rangeSpecified){
		//the range was specified in point terms
		if(!p->dataWave.isPoint){
			err = SUBRANGE_SPECIFIED_ASX;
			goto done;
		} else {
			goiP->startPoint = (long)roundDouble(p->dataWave.startCoord);
			goiP->endPoint = (long)roundDouble(p->dataWave.endCoord);
			if(goiP->startPoint>goiP->endPoint){
				temp = goiP->startPoint;
				goiP->startPoint = goiP->endPoint;
				goiP->endPoint = temp;
			}
			if(goiP->startPoint<0)
				goiP->startPoint = 0;
			if(goiP->endPoint>WavePoints(p->dataWave.waveH)-1)
				goiP->endPoint=WavePoints(p->dataWave.waveH)-1;
		}
	} else {
		//if there is no range specified then we'll use the entire range
		goiP->startPoint = 0;
		goiP->endPoint = WavePoints(p->dataWave.waveH)-1;
	}
	/*
	 use the mask wave and the waverange specified to work out the points we need to fit
	 */
	goiP->unMaskedPoints = WavePoints(p->dataWave.waveH);
	if(p->MFlagEncountered){
		if(err = MDGetDPDataFromNumericWave(p->MFlag_maskwave, goiP->mask)) // Get copy.
			goto done;
	} else {
		//there was no mask wave specfied, use unit weighting
		for(ii = 0 ; ii< WavePoints(p->dataWave.waveH) ; ii+=1)
			*(goiP->mask + ii) = 1;
	}
	/* 
	 set up the mask array
	 need to correct goiP->unMaskedPoints as we go along, which specifies how many unmasked points there will be in the fit 
	 */
	for(ii=0;ii<WavePoints(p->dataWave.waveH);ii+=1){
		temp1 = *(goiP->mask+ii);
		if(*(goiP->mask+ii)==0 || IsNaN64(goiP->mask+ii) || ii < goiP->startPoint || ii>goiP->endPoint || IsNaN64(goiP->dataObsFull+ii) || IsINF64(goiP->dataObsFull+ii)){
			goiP->unMaskedPoints-=1;
			*(goiP->mask+ii)=0;
		}
	}
	
	//you can't fit the data if there's no fit points to be used.
	if(goiP->unMaskedPoints <1){
		err = INPUT_WAVES_NO_POINTS;
		goto done;
	}
	
	//now make the dataCalcwave and a copy of the original (unmasked) data+errors in the current datafolder
	dimensionSizes[0] = goiP->unMaskedPoints;
	dimensionSizes[1] = 0;
	dimensionSizes[2] = 0;
	if(err = MDMakeWave(&goiP->dataCalc,"GenCurveFit_dataCalc",goiP->cDF,dimensionSizes,NT_FP64, 1))
		goto done;
	if(err = MDMakeWave(&goiP->yobs,"GenCurveFit_yobs",goiP->cDF,dimensionSizes,NT_FP64, 1))
		goto done;
	if(err = MDMakeWave(&goiP->sobs,"GenCurveFit_sobs",goiP->cDF,dimensionSizes,NT_FP64, 1))
		goto done;
	
	if(goiP->isAAO){
		for(ii=0 ; ii<goiP->numVarMD ; ii+=1){
			sprintf(letter,"%li",ii);
			strcpy(xwavename,"GenCurveFit_xcalc");
			strcat(&xwavename[0],&letter[0]);
			if(err = MDMakeWave(&goiP->xcalc[ii],xwavename,goiP->cDF,dimensionSizes,NT_FP64, 1))
				goto done;
		}
	}
	
	////create a utility wave that will contains the x range of the original ywave
	dimensionSizes[0] = goiP->dataPoints;
	for(ii=0 ; ii<goiP->numVarMD ; ii+=1){
		strcpy(letter,"");
		sprintf(letter,"%li",ii);		
		strcpy(xwavename,"GenCurveFit_fullExtentOfData0");
		strcat(&xwavename[0],&letter[0]);
		if(err = MDMakeWave(&goiP->fullExtentOfData[ii],xwavename,goiP->cDF,dimensionSizes,NT_FP64, 1))
			goto done;
	}
	
	
	
	////make the temporary coefficients in the current datafolder
	dimensionSizes[0] = WavePoints(p->coefs);
	if(err = MDMakeWave(&goiP->GenCurveFitCoefs,"GenCurveFit_coefs",goiP->cDF,dimensionSizes,NT_FP64, 1))
		goto done;
	
	//initialise population vector
	goiP->gen_populationvector = (double**)malloc2d(goiP->totalpopsize,goiP->numvarparams,sizeof(double));
	if(goiP->gen_populationvector == NULL){
		err = NOMEM;
		goto done;
	}
	//initialise Chi2array
	goiP->chi2Array = (double*)malloc(goiP->totalpopsize * sizeof(double));
	if(goiP->chi2Array == NULL){
		err = NOMEM;
		goto done;
	}
	//initialise the trial vector
	goiP->gen_trial = (double*)malloc(goiP->numvarparams*sizeof(double));
	if(goiP->gen_trial == NULL){
		err = NOMEM;
		goto done;
	}
	//initialise the bprime vector
	goiP->gen_bprime = (double*)malloc(goiP->numvarparams*sizeof(double));
	if(goiP->gen_bprime == NULL){
		err = NOMEM;
		goto done;
	}
	//initialise the pvector
	goiP->gen_pvector = (double*)malloc(goiP->numvarparams*sizeof(double));
	if(goiP->gen_pvector == NULL){
		err = NOMEM;
		goto done;
	}
	//initialise space for a full array copy of the coefficients
	goiP->gen_coefsCopy = (double*)malloc(WavePoints(p->coefs)*sizeof(double));
	if(goiP->gen_coefsCopy == NULL){
		err = NOMEM;
		goto done;
	}
	//put the coefficients into the temporary space
	if(err = MDGetDPDataFromNumericWave(p->coefs, goiP->gen_coefsCopy))
		goto done;
	
	//initialise space for the best fit so far
	goiP->gen_bestfitsofar = (double*)malloc(WavePoints(p->coefs)*sizeof(double));
	if(goiP->gen_bestfitsofar == NULL){
		err = NOMEM;
		goto done;
	}
	//put the starting coefficients into the best fit so far
	if(err = MDGetDPDataFromNumericWave(p->coefs, goiP->gen_bestfitsofar))
		goto done;
	
	//initialise space for an array containing the unmasked fitpoint ydata 
	goiP->dataObs = (double*)malloc(goiP->unMaskedPoints * sizeof(double));
	if (goiP->dataObs == NULL){
		err = NOMEM;
		goto done;
	}
	
	//now fill up the dataObs array, if a point isn't being masked
	jj=0;
	for(ii=0;ii<WavePoints(p->dataWave.waveH);ii+=1){
		if(!(*(goiP->mask+ii) == 0 || IsNaN64(goiP->mask+ii))){
			*(goiP->dataObs+jj) = *(goiP->dataObsFull+ii);
			jj+=1;
		}
	}
	//copy those unmasked points into a dedicated wave
	if(err = MDStoreDPDataInNumericWave(goiP->yobs, goiP->dataObs))
	   goto done;
	
	//initialise array space for putting the calculated model in 
	goiP->dataTemp = (double*)malloc(goiP->unMaskedPoints * sizeof(double));
	if (goiP->dataTemp == NULL){
		err = NOMEM;
	    goto done;
	}
	//initialise array space for putting the limits in
	goiP->limits = (double*)malloc(WavePoints(p->limitswave) * sizeof(double));
	if (goiP->limits == NULL){
		err = NOMEM;
		goto done;
	}
	//put the limits in the dedicated limits array
	if (err = MDGetDPDataFromNumericWave(p->limitswave, goiP->limits))// Get copy.
		goto done;
	
	//initialise space for the weighting data
	goiP->dataSig = (double*)malloc(goiP->unMaskedPoints * sizeof(double));
	if (goiP->dataSig == NULL){
		err = NOMEM;
		goto done;
	}
	/*
	 this section initialises the weightwave, except for those that are masked
	 if there is no weightwave specified then set the weight wave to unity
	 */
	if(p->WFlagEncountered && p->WFlag_weighttype){
		if (err = MDGetDPDataFromNumericWave(p->WFlag_weighttype, goiP->temp)) // Get copy.
			goto done;
		jj=0;
		for(ii=0;ii<WavePoints(p->dataWave.waveH);ii+=1){
			if(!(*(goiP->mask+ii) == 0 || IsNaN64(goiP->mask+ii))){
				*(goiP->dataSig+jj) = *(goiP->temp+ii);
				jj+=1;
			}
		}
	} else {
		for(ii = 0 ; ii< goiP->unMaskedPoints ; ii+=1){
			*(goiP->dataSig+ii)=1;
		}
	}
	//copy those unmasked weighting points into a dedicated wave
	if(err = MDStoreDPDataInNumericWave(goiP->sobs, goiP->dataSig))
	   goto done;
	
	//initialise array space for x values
	goiP->independentVariable = (double*)malloc(goiP->unMaskedPoints*goiP->numVarMD*sizeof(double));
	if (goiP->independentVariable == NULL){
		err = NOMEM;
		goto done;
	}
	goiP->allIndependentVariable = (double*)malloc(goiP->dataPoints*goiP->numVarMD*sizeof(double));
	if (goiP->allIndependentVariable == NULL){
		err = NOMEM;
		goto done;
	}
	
	/*
	 if there was an xwave specified then fill up the temporary x array
	 with the unmasked points from the x-wave, otherwise use waveform scaling of ydata.
	 Either way, fill goiP->fullExtentOfData with the total range of x-values.
	 */
	if(p->XFlagEncountered){
		if(err = MDGetWaveDimensions(p->XFlag_xx, &numdimensions,dimensionSizes)) 
			goto done;
		if(goiP->numVarMD > 1 && numdimensions == 1){
			if (err = MDGetDPDataFromNumericWave(p->XFlag_xx, goiP->allIndependentVariable))// Get copy.
				goto done;
			for(ii=1 ; ii < goiP->numVarMD ; ii+=1){
				if (err = MDGetDPDataFromNumericWave(p->XFlagWaveH[ii-1], goiP->allIndependentVariable+(ii*goiP->dataPoints)))// Get copy.
					goto done;
			}
		} else {
			if (err = MDGetDPDataFromNumericWave(p->XFlag_xx, goiP->allIndependentVariable))// Get copy.
				goto done;
		}
		jj=0;
		for(ii=0 ; ii<goiP->dataPoints ; ii+=1){
			if(!(*(goiP->mask+ii) == 0 || IsNaN64(goiP->mask+ii))){
				for(kk=0 ; kk < goiP->numVarMD ; kk+=1){
					*(goiP->independentVariable+(kk*goiP->unMaskedPoints)+jj) = *(goiP->allIndependentVariable+(kk*goiP->dataPoints)+ii);
				}
				jj+=1;
			}
		}
	} else {	
		//by now the program should've aborted if you haven't specified xwaves and you
		//are fitting multivariate data
		jj=0;
		for(ii=0 ; ii<WavePoints(p->dataWave.waveH);ii+=1){
			*(goiP->allIndependentVariable+ii) = goiP->ystart + ((double)ii)*goiP->ydelta;
			if(!(*(goiP->mask+ii) == 0 || IsNaN64(goiP->mask+ii))){
				*(goiP->independentVariable+jj) = goiP->ystart + ((double)ii)*goiP->ydelta;
				jj+=1;
			}
		}
	}
	for(ii=0; ii<goiP->numVarMD ; ii+=1){
		if(goiP->fullExtentOfData[ii] == NULL){
			err = NOWAV;
			goto done;
		}
		if (err = MDStoreDPDataInNumericWave(goiP->fullExtentOfData[ii],goiP->allIndependentVariable+(ii*goiP->dataPoints)))//put the full extent of x vals into the utilitywave
			goto done;
	}
	
	//store the x array in an x wave used to calculate the theoretical model, but only if you are fitting all at once functions
	//creating these waves is necessary for all-at-once fits.
	if(goiP->isAAO){
		for(ii=0; ii<goiP->numVarMD ; ii+=1){
			if (err = MDStoreDPDataInNumericWave(goiP->xcalc[ii],goiP->independentVariable+(ii*goiP->unMaskedPoints)))//put the full extent of x vals into the utilitywave
				goto done;
		}
	}
	
	/*	setup output
	 if the Dflag is set, then that wave needs to be the same length as the original ywave.  If its not an error
	 will have already been returned.  If the flag is set, then no Xwave will be produced, as you can use the original X-wave.
	 if there is no flag set, then we will return the fit in new waves called coef_ywavename, fit_ywavename
	 BUT BEWARE, YOU MAY BE USING WAVESCALING.
	 */
	WaveName(p->dataWave.waveH,datawavestring);
	strcpy(datawavename,"fit_");
	strcpy(xwavename,"fitx_");
	strcpy(reswavename,"res_");
	for(ii=0;ii<MAX_WAVE_NAME-4;ii+=1){
		letter[0] = datawavestring[ii];
		datawavename[ii+4] =letter[0];
		xwavename[ii+5] =letter[0];
		reswavename[ii+4] = letter[0];
	}
	dimensionSizes[1] = 0;
	dimensionSizes[2] = 0;
	
	if(!p->DFlagEncountered){		//if there is no destination wave specified
		if(goiP->numVarMD == 1){	//and if the data is 1D
			dimensionSizes[0] = (long)p->LFlag_destLen;
			outPutType = WaveType(p->dataWave.waveH);	//to save memory use datawaves type
			if(err = MDMakeWave(&goiP->OUT_data,datawavename,goiP->cDF,dimensionSizes,outPutType, 1))
				goto done;
			minpos = findmin(goiP->allIndependentVariable,goiP->dataPoints);
			maxpos = findmax(goiP->allIndependentVariable,goiP->dataPoints);
			temp1 = *(goiP->allIndependentVariable+maxpos)-*(goiP->allIndependentVariable+minpos);
			temp1 /= floor(p->LFlag_destLen)-1;
			
			if(err = MDSetWaveScaling(goiP->OUT_data, ROWS, &temp1,goiP->allIndependentVariable+minpos))
				goto done;
			if(err = MDMakeWave(&goiP->tempWaveHndl_OUTx,xwavename,goiP->cDF,dimensionSizes,NT_FP64, 1))
				goto done;
			goiP->OUT_x[0] = goiP->tempWaveHndl_OUTx;
			
			for(ii=0 ; ii<(long)p->LFlag_destLen ; ii+=1){
				indices[0] = ii;
				value[0] = *(goiP->allIndependentVariable+minpos)+((double)ii)*temp1;
				if(err = MDSetNumericWavePointValue(goiP->OUT_x[0], indices, value))
					goto done;
			}
		} else {		//no destination wave and data is MD
			dimensionSizes[0] = goiP->dataPoints;
			if(err = MDMakeWave(&goiP->OUT_data,datawavename,goiP->cDF,dimensionSizes,NT_FP64, 1))
				goto done;
			for(ii=0;  ii<MAX_MDFIT_SIZE ; ii+=1) 
				goiP->OUT_x[ii] = goiP->fullExtentOfData[ii];
		}	
	} else {		//destination wave specified
		for(ii=0;  ii<MAX_MDFIT_SIZE ; ii+=1)
			goiP->OUT_x[ii] = goiP->fullExtentOfData[ii];
		goiP->OUT_data = p->DFlag_outputwave;	
	}
	
	if(p->RFlagEncountered){
		if(p->RFlag_resid != NULL){
			goiP->OUT_res = p->RFlag_resid;
		} else {
			dimensionSizes[1] = 0;
			dimensionSizes[0] = WavePoints(p->dataWave.waveH);
			if(err = MDMakeWave(&goiP->OUT_res,reswavename,goiP->cDF,dimensionSizes,NT_FP64, 1))
				goto done;
		}
	}
	
	//if you are doing updates, then append fit output to the topgraph (if data is shown there)
	if(!goiP->noupdate){
		//window for displaying output
		gTheWindow = CreateXOPWindow();
		if(gTheWindow == NULL){
			err = UNSPECIFIED_ERROR;
			return err;
		}
		
		if(goiP->numVarMD == 1){
			if(err = isWaveDisplayed(p->dataWave.waveH, &toDisplay))
				goto done;
			if(toDisplay){
				if(err = isWaveDisplayed(goiP->OUT_data, &toDisplay))
					goto done;
				if(!toDisplay){
					strcpy(cmd,"appendtograph/w=$(winname(0,1)) ");
					WaveName(goiP->OUT_data, &datawavename[0]);
					strcat(cmd, &datawavename[0]);
					
					if(p->DFlagEncountered){
						if(p->XFlagEncountered){
							WaveName(p->XFlag_xx,&xwavename[0]);
						} else {
							WaveName(goiP->OUT_x[0], &xwavename[0]);
						}
						strcat(cmd," vs ");
						strcat(cmd,&xwavename[0]);
					}
					if(err = XOPSilentCommand(&cmd[0]))
						goto done;
				}
			}
		}
	}
	
	//initialise population vector guesses, from between the limits
	for(ii=0; ii<goiP->totalpopsize ; ii+=1){
		for(jj=0 ; jj<goiP->numvarparams ; jj+=1){
			bot = *(goiP->limits+*(goiP->varparams+jj));
			top = *(goiP->limits+*(goiP->varparams+jj)+WavePoints(p->coefs));
			goiP->gen_populationvector[ii][jj] = randomDouble(bot,top, 0, 0);
		}
	}
	
	//now deal with all the options flags
	if(p->OPTFlagEncountered){
		//bit 0 of p->opt is set, then the initial guesses are used for the fit.
		if((long)p->OPTFlag_opt & (long)pow(2,0)){
			for(jj=0 ; jj<goiP->numvarparams ; jj+=1){
				ii = *(goiP->varparams+jj);
				goiP->gen_populationvector[0][jj] = *(goiP->gen_coefsCopy+ii);
			}
		}
	}
	
	//initialise Chi2array
	for(ii=0; ii<goiP->totalpopsize ; ii+=1){
		//perhaps the user wants to abort the fit straightaway, this GUI button does that.
		if(!goiP->noupdate){
			if (Abort_the_fit)
				return FIT_ABORTED;
		}
		//cmd-dot or abort button
		if(CheckAbort(timeOutTicks) == -1){
			err = FIT_ABORTED;
			goto done;
		}
		if(err = setPvectorFromPop(goiP, ii))
			goto done;
		if(err = insertVaryingParams(goiP, p))
			goto done;
		if(err = calcModel(&goiP->fi,
						   goiP->GenCurveFitCoefs,
						   goiP->dataCalc,
						   goiP->dataTemp,
						   goiP->xcalc,
						   goiP->independentVariable,
						   goiP->numVarMD,
						   goiP->isAAO,
						   goiP->sp))
			goto done;
		switch(goiP->METH){//which cost function
			case 0:
				if(err = calcChi2(goiP->dataObs,goiP->dataTemp,goiP->dataSig,goiP->unMaskedPoints,&chi2,goiP->weighttype))
					goto done;
				break;
			case 1:
				if(err = calcRobust(goiP->dataObs,goiP->dataTemp,goiP->dataSig,goiP->unMaskedPoints,&chi2,goiP->weighttype))
					goto done;
				break;
			case 2:
				if(err = calcUserCostFunc(goiP->minf, 
										  goiP->yobs,
										  goiP->dataObs,
										  goiP->dataCalc,
										  goiP->sobs,
										  goiP->dataSig,
										  goiP->unMaskedPoints,
										  goiP->GenCurveFitCoefs,
										  &chi2))
					goto done;
				break;
			default:
				break;
		}
		*(goiP->chi2Array+ii)= chi2;
	}
	//find best chi2 and put that into number 0 pos.
	wavStats = getWaveStats(goiP->chi2Array,goiP->totalpopsize,0);
	
	swapChi2values(goiP,0,wavStats.V_minloc);
	if(err = swapPopVector(goiP,goiP->totalpopsize,0,wavStats.V_minloc))
		goto done;
	
	if(err = setPvectorFromPop(goiP, 0))
		goto done;
	if(err = insertVaryingParams(goiP, p))
		goto done;
	memcpy(goiP->gen_bestfitsofar, goiP->gen_coefsCopy, sizeof(double) * goiP->totalnumparams);
	
done:
	if(holdstr != NULL)
		free(holdstr);
	
	return err;
}


/*
 freeAllocMem frees all the temporary arrays in the GenCurveFitInternals structure
 */
static void
freeAllocMem(GenCurveFitInternalsPtr goiP){
	int err=0,ii=0;
	
	waveHndl exists = NULL;
	if(goiP->temp!=NULL)
		free(goiP->temp);
	if(goiP->chi2Array!=NULL)
		free(goiP->chi2Array);
	if(goiP->gen_populationvector!=NULL)
		free(goiP->gen_populationvector);
	if(goiP->gen_coefsCopy!=NULL)
		free(goiP->gen_coefsCopy);
	if(goiP->gen_bestfitsofar!=NULL)
		free(goiP->gen_bestfitsofar);
	if(goiP->gen_bprime!=NULL)
		free(goiP->gen_bprime);
	if(goiP->gen_trial!=NULL)
		free(goiP->gen_trial);
	if(goiP->limits!=NULL)
		free(goiP->limits);
	if(goiP->mask!=NULL)
		free(goiP->mask);
	if(goiP->varparams!=NULL)
		free(goiP->varparams);
	if(goiP->independentVariable!=NULL)
		free(goiP->independentVariable);
	if(goiP->allIndependentVariable!=NULL)
		free(goiP->allIndependentVariable);
	if(goiP->dataObs!=NULL)
		free(goiP->dataObs);
	if(goiP->dataObsFull)
		free(goiP->dataObsFull);
	if(goiP->dataSig!=NULL)
		free(goiP->dataSig);
	if(goiP->dataTemp!=NULL)
		free(goiP->dataTemp);
	if(goiP->gen_pvector!=NULL)
		free(goiP->gen_pvector);
	if(goiP->covarianceMatrix!=NULL)
		free(goiP->covarianceMatrix);
	
	exists = FetchWaveFromDataFolder(goiP->cDF,"GenCurveFit_coefs");
	if(exists != NULL)
		err= KillWave(exists);	
	exists = FetchWaveFromDataFolder(goiP->cDF,"GenCurveFit_dataCalc");
	if(exists != NULL)
		err = KillWave(goiP->dataCalc);
	exists = FetchWaveFromDataFolder(goiP->cDF,"GenCurveFit_yobs");
	if(exists != NULL)
		err = KillWave(goiP->yobs);
	exists = FetchWaveFromDataFolder(goiP->cDF,"GenCurveFit_sobs");
	if(exists != NULL)
		err = KillWave(goiP->sobs);
	
	for(ii=0 ; ii<goiP->numVarMD ; ii+=1){
		if(goiP->xcalc[ii] != NULL)
			err = KillWave(goiP->xcalc[ii]);
	}
	if(goiP->tempWaveHndl_OUTx != NULL)
		err = KillWave(goiP->tempWaveHndl_OUTx);
	
	for(ii=0 ; ii<goiP->numVarMD ; ii+=1){
		if(goiP->fullExtentOfData[ii] != NULL)
			err = KillWave(goiP->fullExtentOfData[ii]);
	}
	//gTheWindow is a window created to show the latest position of the fit
	if (gTheWindow != NULL) {
		DestroyXOPWindow(gTheWindow);
		gTheWindow = NULL;
	}
}

static void
checkLimits(GenCurveFitInternalsPtr goiP,GenCurveFitRuntimeParamsPtr p){
	int ii;
	for(ii=0 ; ii<goiP->numvarparams ; ii+=1){
		if(*(goiP->gen_trial+ii) < *(goiP->limits+*(goiP->varparams+ii)) || *(goiP->gen_trial+ii) > *(goiP->limits+*(goiP->varparams+ii)+WavePoints(p->coefs)))
			*(goiP->gen_trial+ii) = randomDouble(*(goiP->limits+*(goiP->varparams+ii)),*(goiP->limits+*(goiP->varparams+ii)+WavePoints(p->coefs)), 0, 0);
	}
}


/*
 randomInteger returns an integer between 0 and upper EXclusive
 i.e. you will never get upper returned.
 */
static int
randomInteger (int upper, uint32_t seed, short initialise){
	static dsfmt_t generator;
	static int dispensed = 0;
	static double randomNumbers[1000];
	double randomNumber = 0.0;
	int val = 0;
	
	if(initialise){
		dsfmt_init_gen_rand(&generator, seed);
		dsfmt_fill_array_open_close(&generator, randomNumbers, 1000);
		dispensed = 0;
		return 0;
	}
	if(dispensed == 1000){
		dsfmt_fill_array_open_close(&generator, randomNumbers, 1000);
		dispensed = 0;
	}
	
	randomNumber = (double)upper * randomNumbers[dispensed];
	val = (int) randomNumber;
	dispensed += 1;
	
//	while (upper <= (val = rand() / (RAND_MAX/upper)));
	return val;

}

/*
 randomDouble returns a double value between lower <= x <= upper OR [lower,upper]
 */
static double
randomDouble(double lower, double upper, uint32_t seed, short initialise){
	static dsfmt_t generator;
	static int dispensed = 0;
	static double randomNumbers[1000];
	double randomNumber = 0.0;
	
	if(initialise){
		dsfmt_init_gen_rand(&generator, seed);
		dsfmt_fill_array_open_open(&generator, randomNumbers, 1000);
		dispensed = 0;
		return 0.0;
	}
	if(dispensed == 1000){
		dsfmt_fill_array_open_open(&generator, randomNumbers, 1000);
		dispensed = 0;
	}
	
	randomNumber = lower + randomNumbers[dispensed]*(upper-lower);
	dispensed += 1;
	return randomNumber;
//	return lower + rand()/(((double)RAND_MAX + 1)/(upper-lower));
}

/*
 insertVaryingParams inserts the current pvector into an array copy of the coefficients,
 then into a temporary wave
 returns 0 if no error
 returns errorcode otherwise
 */
static int
insertVaryingParams(GenCurveFitInternalsPtr goiP, GenCurveFitRuntimeParamsPtr p){
	int err=0,ii;
	
	for(ii=0 ; ii< goiP->numvarparams; ii+=1)
		*(goiP->gen_coefsCopy+*(goiP->varparams+ii)) =  *(goiP->gen_pvector+ii);
	
	err = MDStoreDPDataInNumericWave(goiP->GenCurveFitCoefs, goiP->gen_coefsCopy);
	return err;
}

/*
 extractVaryingParams copies the entire fit parameters from a wave into a temporary array
 it then extracts the varying parameters from that array and puts them into the pvector
 returns 0 if no error
 returns errorcode otherwise
 */
static int 
extractVaryingParams(GenCurveFitInternalsPtr goiP, GenCurveFitRuntimeParamsPtr p){
	int err=0,ii;
	
	if(err = MDGetDPDataFromNumericWave(goiP->GenCurveFitCoefs,goiP->gen_coefsCopy))
		return err;
	for(ii=0 ; ii< goiP->numvarparams; ii+=1){
		*(goiP->gen_pvector+ii) = *(goiP->gen_coefsCopy+*(goiP->varparams+ii));
	}
	return err;
}

/*
 calcModel calculates the theoretical curve for the model, given the coefficients
 fip			-	the function
 coefs		-	the coefficients to use in calculation
 output		-	where to put the theoretical fit (if using all-at-once function)
 outputPtr	-	where to put the theoretical fit (if you are using normal fit function)
 xx			-	wave containing the x values (if using all-at-once function)
 xpnts		-	array containing the x values (if you are using normal fit function)
 ndims		-	the dimensionality of the fit (i.e. how many independent variables there are
 isAAO		-	is the fit all-at-once?
 sp			-	if its a structure fit then this will hold the users structure
 returns 0 if no error
 returns errorcode otherwise
 */
int
calcModel(FunctionInfo *fip, waveHndl coefs, waveHndl output, double* outputPtr, waveHndl xx[MAX_MDFIT_SIZE], double* xpnts, int ndims,int isAAO,	fitfuncStruct* sp){
	int err = 0, ii,jj;
	int requiredParameterTypes[MAX_MDFIT_SIZE+2];
//	int badParameterNumber;
	allFitFunc allParameters;
	fitFunc parameters;
	
	//experimenting with faster wave access
	long dataOffset; 
	
	double result;
	long numfitpoints = WavePoints(output);
	
	// check if all the input isn't NULL
	if(coefs == NULL || output==NULL){
		err = UNSPECIFIED_ERROR;
		goto done;
	}
	if(fip == NULL){
		err = FITFUNC_NOT_SPECIFIED;	
		goto done;
	}
	if(outputPtr == NULL){
		err = UNSPECIFIED_ERROR;
		goto done;
	}
	
	switch(isAAO){
		case 0:
			if(xpnts == NULL){
				err = UNSPECIFIED_ERROR;
				goto done;
			}
			parameters.waveH = coefs;
			requiredParameterTypes[0] = WAVE_TYPE;
/*
			for(ii=0 ; ii<ndims ; ii+=1)
				requiredParameterTypes[ii+1] = NT_FP64;
			if (err = CheckFunctionForm(fip, ndims+1 , requiredParameterTypes,&badParameterNumber, NT_FP64))
				goto done;
*/			
			for(ii=0 ; ii<numfitpoints ; ii+=1){
				for(jj=0 ; jj<ndims ; jj+=1)
					parameters.x[jj] = *(xpnts + (jj * numfitpoints) + ii);

				// call the users fit function and put the result in the output array
				if (err = CallFunction(fip, (void*) &parameters, &result))
					goto done;
				*(outputPtr + ii) = result;
			}
			
			// copy the output array into an output wave
			//			if(err = MDStoreDPDataInNumericWave(output,outputPtr))
			//				goto done;
			
			if (err = MDAccessNumericWaveData(output, kMDWaveAccessMode0, &dataOffset)) 
				goto done; 
			memcpy((*output) + dataOffset , outputPtr, numfitpoints * sizeof(double));
			
			break;
		case 1:
			allParameters.waveC = coefs;
			allParameters.waveY = output;
			requiredParameterTypes[0] = WAVE_TYPE;
			requiredParameterTypes[1] = WAVE_TYPE;
			
			for(ii=0 ; ii<ndims ; ii+=1){
				requiredParameterTypes[ii+2] = WAVE_TYPE;
				if(xx[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				allParameters.waveX[ii] = xx[ii];
			}
//			if (err = CheckFunctionForm(fip, ndims + 2, requiredParameterTypes,&badParameterNumber, -1))
//				goto done;
			// call the users fit function and put the result in the output wave
			if (err = CallFunction(fip, (void*)&allParameters, &result))
				goto done;
			// the user may have changed the number of points in the output wave
			if(output == NULL || WavePoints(output) != numfitpoints){
				err = USER_CHANGED_FITWAVE;
				goto done;
			}
			
			// get the output wave and put it into the output array
			//			if(err = MDGetDPDataFromNumericWave(output,outputPtr))
			//				goto done;			
			if (err = MDAccessNumericWaveData(output, kMDWaveAccessMode0, &dataOffset)) 
				goto done; 
			memcpy(outputPtr , (*output) + dataOffset, numfitpoints * sizeof(double));
			break;
		case 2:
			if(sp == NULL){
				err = NULL_STRUCTURE;
				goto done;
			}
			(*sp).w = coefs;
			(*sp).yy = output;
			
			requiredParameterTypes[0] = FV_STRUCT_TYPE | FV_REF_TYPE;
			
			for(ii=0 ; ii<ndims ; ii+=1){
				if(xx[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				((*sp).xx[ii]) = (xx[ii]);
			}
			
//			if (err = CheckFunctionForm(fip, 1, requiredParameterTypes,&badParameterNumber, -1))
//				goto done;
			// call the users fit function and put the result in the output wave
			if (err = CallFunction(fip, (fitfuncStruct*)&sp, &result))
				goto done;
			//don't want any dangling references to waves
			sp->yy = NULL;
			sp->w = NULL;
			memset(sp->xx, 0, sizeof(sp->xx));
			
			// the user may have changed the number of points in the output wave
			if(output == NULL || WavePoints(output) != numfitpoints){
				err = USER_CHANGED_FITWAVE;
				goto done;
			}
			// get the output wave and put it into the output array
			//			if(err = MDGetDPDataFromNumericWave(output,outputPtr))
			//				goto done;
			if (err = MDAccessNumericWaveData(output, kMDWaveAccessMode0, &dataOffset)) 
				goto done; 
			memcpy(outputPtr , (*output) + dataOffset, numfitpoints * sizeof(double));
			break;
		default:
			err = UNSPECIFIED_ERROR;
			goto done;
			break;
	}
	
	// check that the fitfunction didn't return any NaN or INF
	if(err = checkNanInf(output)){
		err = FITFUNC_RETURNED_NANINF;
		goto done;
	}
done:		
	return err;
}


/*
 calcModelXY calculates the theoretical curve for the model, used for returning the results of the fit to igor
 fip			-	the function
 coefs		-	the coefficients to use in calculation
 output		-	where to put the theoretical curve
 xx			-	wave containing the x values
 ndims		-	the dimensionality of the fit (i.e. how many independent variables there are
 isAAO		-	is the fit all-at-once?
 sp			-	the structure for a structure fit
 returns 0 if no error
 returns errorcode otherwise
 
 it's needed because the output that is created isn't necessarily on the same X spacing as the fit
 */
static int
calcModelXY(FunctionInfo* fip, waveHndl coefs, waveHndl output, waveHndl xx[MAX_MDFIT_SIZE], int ndims,int isAAO, fitfuncStruct *sp){
	int err = 0, ii,jj;
	int requiredParameterTypes[MAX_MDFIT_SIZE+2];
	allFitFunc allParameters;
	fitFunc parameters;
	
	double result;
	double *tempX = NULL;
	double *tempY = NULL;
	long numfitpoints = WavePoints(output);
	
	// check if all the input isn't NULL
	if(coefs == NULL || output==NULL){
		err = UNSPECIFIED_ERROR;
		goto done;
	}
	if(fip == NULL){
		err = FITFUNC_NOT_SPECIFIED;	
		goto done;
	}
	
	switch(isAAO){
		case 0:
			tempX = (double*)malloc(ndims*numfitpoints*sizeof(double));
			if(tempX == NULL){
				err = NOMEM;
				goto done;
			}
			tempY = (double*)malloc(numfitpoints*sizeof(double));
			if(tempY == NULL){
				err = NOMEM;
				goto done;
			}
			for(ii=0 ; ii<ndims ; ii+=1){
				if(xx[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				if(err = MDGetDPDataFromNumericWave(xx[ii],tempX+(numfitpoints*ii)))
					goto done;
			}
			
			parameters.waveH = coefs;
			requiredParameterTypes[0] = WAVE_TYPE;
			for(ii=0 ; ii<ndims ; ii+=1)
				requiredParameterTypes[ii+1] = NT_FP64;
			
			for(ii=0 ; ii<numfitpoints ; ii+=1){
				for(jj=0 ; jj<ndims ; jj+=1){
					parameters.x[jj] = *(tempX+(jj*numfitpoints)+ ii);
				}
				// call the users fit function and put the result in the output array
				if (err = CallFunction(fip, (void*)&parameters, &result))
					goto done;
				*(tempY+ii) = result;
			}
			// copy the output array into an output wave
			if(err = MDStoreDPDataInNumericWave(output,tempY))
				goto done;
			break;
		case 1:
			allParameters.waveC = coefs;
			allParameters.waveY = output;
			requiredParameterTypes[0] = WAVE_TYPE;
			requiredParameterTypes[1] = WAVE_TYPE;
			
			for(ii=0 ; ii<ndims ; ii+=1){
				requiredParameterTypes[ii+2] = WAVE_TYPE;
				if(xx[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				allParameters.waveX[ii] = xx[ii];
			}
			// call the users fit function and put the result in the output wave
			if (err = CallFunction(fip, (void*)&allParameters, &result))
				goto done;
			
			// the user may have changed the number of points in the output wave
			if(output == NULL || WavePoints(output) != numfitpoints){
				err = USER_CHANGED_FITWAVE;
				goto done;
			}
			break;
		case 2:
			if(sp == NULL){
				err = NULL_STRUCTURE;
				goto done;
			}
			(*sp).w = coefs;
			(*sp).yy = output;
			
			requiredParameterTypes[0] = FV_STRUCT_TYPE | FV_REF_TYPE;
			
			for(ii=0 ; ii<ndims ; ii+=1){
				if(xx[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				((*sp).xx[ii]) = (xx[ii]);
			}
			
			// call the users fit function and put the result in the output wave
			if (err = CallFunction(fip, (fitfuncStruct*)&sp	, &result))
				goto done;
			//don't want any dangling references to waves
			sp->yy = NULL;
			sp->w = NULL;
			memset(sp->xx, 0, sizeof(sp->xx));
			
			// the user may have changed the number of points in the output wave
			if(output == NULL || WavePoints(output) != numfitpoints){
				err = USER_CHANGED_FITWAVE;
				goto done;
			}
			break;
		default:
			err = UNSPECIFIED_ERROR;
			goto done;
			break;
	}
	
	// check that the fitfunction didn't return any NaN or INF
	if(err = checkNanInf(output)){
		err = FITFUNC_RETURNED_NANINF;
		goto done;
	}
	
done:
	
	if(tempX != NULL)
		free(tempX);
	if(tempY != NULL)
		free(tempY);
	
	return err;
}



/*
 subtractTwoWaves subtracts wav2 from wav1
 returns 0 if no error
 returns errorcode otherwise
 */
static int
subtractTwoWaves(waveHndl wav1, waveHndl wav2){
	int err = 0;
	long ii;
	double *temp1 = NULL,*temp2 = NULL;
	double val=0,val2=0,val3=0;
	//check if the wave references are NULL
	if(wav1 == NULL)
		return NON_EXISTENT_WAVE;
	if(wav2 == NULL)
		return NON_EXISTENT_WAVE;
	if(WavePoints(wav1)!=WavePoints(wav2))
		return WAVE_LENGTH_MISMATCH;
	
	// we have to create temporary arrays to hold the wave data
	if((temp1 = (double*)malloc(sizeof(double)*WavePoints(wav1))) ==  NULL ){
		err = NOMEM;
		goto done;
	}
	if((temp2 = (double*)malloc(sizeof(double)*WavePoints(wav1))) ==  NULL ){
		err = NOMEM;
		goto done;
	}
	// get the data from the waves and put it into the temporary arrays
	if(err = MDGetDPDataFromNumericWave(wav1,temp1))
		goto done;
	if(err = MDGetDPDataFromNumericWave(wav2,temp2))
		goto done;
	// do the subtraction
	for(ii=0;ii<WavePoints(wav1);ii+=1){
		val = *(temp1+ii);
		val2 = *(temp2+ii);
		val3 = val - val2;
		*(temp1+ii) = val - val2;
	}
	// store the subtraction in wav1
	if(err = MDStoreDPDataInNumericWave(wav1,temp1))
		goto done;
	
	WaveHandleModified(wav1);
done:
	if(temp1!=NULL)
		free(temp1);
	if(temp2!=NULL)
		free(temp2);
	
	return err;
}
/*
 multiplies wav1 by scalar
 returns 0 if no error
 returns errorcode otherwise
 */
static int
scalarMultiply(waveHndl wav1, double scalar){
	int err = 0;
	long ii;
	double *temp1 = NULL;
	double val=0;
	//check if the wave references are NULL
	if(wav1 == NULL)
		return NON_EXISTENT_WAVE;
	
	// we have to create temporary arrays to hold the wave data
	if((temp1 = (double*)malloc(sizeof(double)*WavePoints(wav1))) ==  NULL ){
		err = NOMEM;
		goto done;
	}
	
	if(err = MDGetDPDataFromNumericWave(wav1,temp1))
		goto done;
	// do the subtraction
	for(ii=0;ii<WavePoints(wav1);ii+=1){
		val = *(temp1+ii)*scalar;
		*(temp1+ii) = val;
	}
	// store the subtraction in wav1
	if(err = MDStoreDPDataInNumericWave(wav1,temp1))
		goto done;
	
	WaveHandleModified(wav1);
	
done:
	if(temp1!=NULL)
		free(temp1);
	
	return err;
}

/*
 swapPopVector swaps the i vector with index j in the populationvector
 returns 0 if no error
 returns errorcode otherwise
 */
static int
swapPopVector(GenCurveFitInternalsPtr goiP,int popsize, int i, int j){
	double *tempparams = NULL;
	if(i<0 || j<0 || i>popsize-1 || j>popsize-1){
		return UNSPECIFIED_ERROR;
	} else {
		//do swap with pointers
		tempparams = *(goiP->gen_populationvector+j);
		*(goiP->gen_populationvector+j) = *(goiP->gen_populationvector+i);
		*(goiP->gen_populationvector+i) = tempparams;
		return 0;
	}
}
/*
 setPvector sets the pvector with a double array, checking to make sure the sizes are right
 returns 0 if no error
 returns errorcode otherwise
 */
static int 
setPvector(GenCurveFitInternalsPtr goiP,double* vector, int vectorsize){
	//goiP->gen_pvector[] = vector[p]
	memcpy(goiP->gen_pvector,vector,vectorsize*sizeof(double));
	return 0;
}
/*
 setPvectorFromPop sets the pvector from index vector from the population vector
 returns 0 if no error
 returns errorcode otherwise
 */

static int 
setPvectorFromPop(GenCurveFitInternalsPtr goiP, int vector){
	//goiP->gen_pvector[] = goiP->gen_populationvector[vector][p]
	memcpy(goiP->gen_pvector, *(goiP->gen_populationvector+vector), goiP->numvarparams*sizeof(double));
	return 0;
}

/*
 setPopVectorFromPVector sets the populationvector with index replace, with a double array
 returns 0 if no error
 returns errorcode otherwise
 */
static int 
setPopVectorFromPVector(GenCurveFitInternalsPtr goiP,double* vector, int vectorsize, int replace){
	//goiP->gen_populationvector[replace][] = vector[q]
	memcpy(*(goiP->gen_populationvector+replace), vector, vectorsize*sizeof(double));
	return 0;
}

/*
 swapChi2values swaps two values (i,j) in the goiP->chi2array 
 */
static void
swapChi2values(GenCurveFitInternalsPtr goiP, int i, int j){
	double temp = *(goiP->chi2Array+i);
	*(goiP->chi2Array+i) = *(goiP->chi2Array+j);
	*(goiP->chi2Array+j) = temp;
}
/*
 findmin finds the minimum value in a pointer array
 returns minimum position.
 */
static int
findmin(double* sort, int sortsize){
	int ii = 0 , minpos = 0;
	double minval = *(sort+ii);
	for(ii=0 ; ii<sortsize ; ii+=1){
		if(*(sort+ii) < minval){
			minval = *(sort+ii);
			minpos = ii;
		}
	}
	return minpos;
}
/*
 findMax finds the minimum value in a pointer array
 returns max position.
 */
static int
findmax(double* sort, int sortsize){
	int ii = 0 , maxpos = 0;
	double maxval = *(sort+ii);
	for(ii=0 ; ii<sortsize ; ii+=1){
		if(*(sort+ii) > maxval){
			maxval = *(sort+ii);
			maxpos = ii;
		}
	}
	return maxpos;
}


/*
 createTrialVector makes a mutated vector.  It fills the trialVector from the current pvector and from bPrime,
 in modulo.
 bPrime is created from two random population vectors and the best fit vector.
 */
static void 
createTrialVector(GenCurveFitInternalsPtr goiP, GenCurveFitRuntimeParamsPtr p, int currentpvector){
	void (*theStrategy)(GenCurveFitInternalsPtr, int);
	
	switch(goiP->STGY){
		case 0:
			theStrategy = Best1Bin;
			break;
		case 1:
			theStrategy = Best1Exp;
			break;
		case 2:
			theStrategy = Rand1Exp;
		case 3:
			theStrategy = RandToBest1Exp;
			break;
		case 4:
			theStrategy = Best2Exp;
			break;
		case 5:
			theStrategy = Rand2Exp;
			break;
		case 6:
			theStrategy = RandToBest1Bin;
			break;
		case 7:
			theStrategy = Best2Bin;
			break;
		case 8:
			theStrategy = Rand2Bin;
			break;
		default:
			theStrategy = Best1Bin;
			break;
	}
	
	theStrategy(goiP, currentpvector);
}

void SelectSamples(int popsize, int candidate,int *r1,int *r2, int *r3,int *r4,int *r5){
	if (r1){
		do
			*r1 = randomInteger(popsize, 0, 0);	
		while (*r1 == candidate);
	}
	
	if (r2)	{
		do
			*r2 = randomInteger(popsize, 0, 0);
		while ((*r2 == candidate) || (*r2 == *r1));
	}
	
	if (r3){
		do
			*r3 = randomInteger(popsize, 0, 0);
		while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
	}
	
	if (r4){
		do
			*r4 = randomInteger(popsize, 0, 0);
		while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
	}
	
	if (r5){
		do
			*r5 = randomInteger(popsize, 0, 0);
		while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
			   || (*r5 == *r2) || (*r5 == *r1));
	}
	
	return;
}


void Best1Bin(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2;
	int n, i;
	
	SelectSamples(goiP->totalpopsize, candidate, &r1, &r2, NULL, NULL, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));
	
	for (i=0; i < goiP->numvarparams; i++) {
		if ((randomDouble(0, 1, 0, 0) < goiP->recomb) || (i == (goiP->numvarparams - 1)))
			*(goiP->gen_trial + n) = goiP->gen_populationvector[0][n]
									+ goiP->k_m * (goiP->gen_populationvector[r1][n] - goiP->gen_populationvector[r2][n]);

		n = (n + 1) % goiP->numvarparams;
	}
	return;
}

void Best1Exp(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2;
	int n, i;
	
	SelectSamples(goiP->totalpopsize, candidate, &r1, &r2, NULL, NULL, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);

	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));
	
	for (i=0 ; (randomDouble(0, 1, 0, 0) < goiP->recomb) && (i < goiP->numvarparams); i++){
		*(goiP->gen_trial + n) = goiP->gen_populationvector[0][n]
								+ goiP->k_m * (goiP->gen_populationvector[r1][n] - goiP->gen_populationvector[r2][n]);
		
		n = (n + 1) % goiP->numvarparams;
	}
	return;
}

void Rand1Exp(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2, r3;
	int n, i;
	
	SelectSamples(goiP->totalpopsize, candidate,&r1,&r2,&r3, NULL, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));

	for (i=0; (randomDouble(0, 1, 0, 0) < goiP->recomb) && (i < goiP->numvarparams); i++) {
		*(goiP->gen_trial + n) = goiP->gen_populationvector[r1][n]
		+ goiP->k_m * (goiP->gen_populationvector[r2][n] - goiP->gen_populationvector[r3][n]);

		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}

void RandToBest1Exp(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2;
	int n,  i;
	
	SelectSamples(goiP->numvarparams, candidate,&r1,&r2, NULL, NULL, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));
	
	for (i=0; (randomDouble(0, 1, 0, 0) < goiP->recomb) && (i < goiP->numvarparams); i++) {
		*(goiP->gen_trial + n) += goiP->k_m * (goiP->gen_populationvector[0][n] - *(goiP->gen_trial + n))
		+ goiP->k_m * (goiP->gen_populationvector[r1][n]
				   - goiP->gen_populationvector[r2][n]);
		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}

void Best2Exp(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2, r3, r4;
	int n, i;
	
	SelectSamples(goiP->numvarparams, candidate,&r1,&r2,&r3,&r4, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));

	for (i=0; (randomDouble(0, 1, 0, 0) < goiP->recomb) && (i < goiP->numvarparams); i++) {
		*(goiP->gen_trial + n) = goiP->gen_populationvector[0][n] +
		goiP->k_m * (goiP->gen_populationvector[r1][n]
				 + goiP->gen_populationvector[r2][n]
				 - goiP->gen_populationvector[r3][n]
				 - goiP->gen_populationvector[r4][n]);
		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}

void Rand2Exp(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2, r3, r4, r5;
	int n, i;
	
	SelectSamples(goiP->numvarparams, candidate,&r1,&r2,&r3,&r4,&r5);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));

	for (i=0; (randomDouble(0, 1, 0, 0) < goiP->recomb) && (i < goiP->numvarparams); i++) {
		*(goiP->gen_trial + n) = goiP->gen_populationvector[r1][n]
		+ goiP->k_m * (goiP->gen_populationvector[r2][n]
				   + goiP->gen_populationvector[r3][n]
				   - goiP->gen_populationvector[r4][n]
				   - goiP->gen_populationvector[r5][n]);
		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}

void RandToBest1Bin(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2;
	int n, i;
	
	SelectSamples(goiP->numvarparams, candidate,&r1,&r2, NULL, NULL, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));
	for (i=0; i < goiP->numvarparams; i++) 
	{
		if ((randomDouble(0, 1, 0, 0) < goiP->recomb) || (i  == (goiP->numvarparams - 1)))
			*(goiP->gen_trial + n) += goiP->k_m * (goiP->gen_populationvector[0][n] - *(goiP->gen_trial + n))
			+ goiP->k_m * (goiP->gen_populationvector[r1][n]
					   - goiP->gen_populationvector[r2][n]);
		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}

void Best2Bin(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2, r3, r4;
	int n, i;
	
	SelectSamples(goiP->numvarparams, candidate,&r1,&r2,&r3,&r4, NULL);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));
	for (i=0; i < goiP->numvarparams; i++) 
	{
		if ((randomDouble(0, 1, 0, 0) < goiP->recomb) || (i  == (goiP->numvarparams - 1)))
			*(goiP->gen_trial + n) = goiP->gen_populationvector[0][n]
			+ goiP->k_m * (goiP->gen_populationvector[r1][n]
					   + goiP->gen_populationvector[r2][n]
					   - goiP->gen_populationvector[r3][n]
					   - goiP->gen_populationvector[r4][n]);
		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}

void Rand2Bin(GenCurveFitInternalsPtr goiP, int candidate){
	int r1, r2, r3, r4, r5;
	int n, i;
	
	SelectSamples(goiP->numvarparams, candidate,&r1,&r2,&r3,&r4,&r5);
	n = randomInteger(goiP->numvarparams, 0, 0);
	
	memcpy(goiP->gen_trial, *(goiP->gen_populationvector + candidate), goiP->numvarparams * sizeof(double));
	for (i=0; i < goiP->numvarparams; i++) 
	{
		if ((randomDouble(0, 1, 0, 0) < goiP->recomb) || (i  == (goiP->numvarparams - 1)))
			*(goiP->gen_trial + n) = goiP->gen_populationvector[r1][n]
			+ goiP->k_m * (goiP->gen_populationvector[r2][n]
					   + goiP->gen_populationvector[r3][n]
					   - goiP->gen_populationvector[r4][n]
					   - goiP->gen_populationvector[r5][n]);
		n = (n + 1) % goiP->numvarparams;
	}
	
	return;
}


/*
 ensureConstraints takes the current trial vector and makes sure that all the individual 
 parameters lie inbetween the upper and lower limits.
 returns void
 */
static void
ensureConstraints(GenCurveFitInternalsPtr goiP,GenCurveFitRuntimeParamsPtr p){
	int ii;
	long points = WavePoints(p->coefs);
	
	for(ii=0 ; ii < goiP->numvarparams ; ii+=1){
		if(*(goiP->gen_trial+ii) <*(goiP->limits+*(goiP->varparams+ii)) || *(goiP->gen_trial+ii)>(*(goiP->limits+*(goiP->varparams+ii)+points))){
			*(goiP->gen_trial+ii) = randomDouble(*(goiP->limits+*(goiP->varparams+ii)),(*(goiP->limits+*(goiP->varparams+ii)+points)), 0, 0);
		}
	}
}

/*
 optimiseloop performs the optimisation.  It takes the initial population and mutates it until we find the best fit solution
 returns 0 if no errors
 returns errorcode otherwise.
 */
static int
optimiseloop(GenCurveFitInternalsPtr goiP, GenCurveFitRuntimeParamsPtr p){
	
	long ii,kk,mm;
	int err=0;
	int currentpvector;
	double chi2pvector,chi2trial;
	int acceptMoveGrudgingly = 0;
	waveStats wavStats;
	
	//an array for dumping the population at each iteration
	MemoryStruct dumpRecord;
	dumpRecord.memory = NULL;
	dumpRecord.size = 0;
	
	//Display the coefficients so far.
	if(!goiP->noupdate){
		DisplayWindowXOP1Message(gTheWindow, WavePoints(p->coefs), goiP->gen_coefsCopy, *(goiP->chi2Array), goiP->fi.name, goiP->V_numfititers, goiP->convergenceNumber);
	}
	
	if(p->DUMPFlagEncountered && (!p->TEMPFlagEncountered)){
		for(mm=0 ; mm < goiP->totalpopsize ; mm+=1){
			WriteMemoryCallback((goiP->gen_populationvector[mm]), sizeof(double), goiP->numvarparams, &dumpRecord);
			if(dumpRecord.memory == NULL){
				err = NOMEM;goto done;
			}
		}
	}
	
	// the user sets how many times through the entire population
	for(kk=1; kk<=p->KFlag_iterations ; kk+=1){
		goiP->V_numfititers = kk;
		
		if(!goiP->noupdate){
			if(err = setPvectorFromPop(goiP, 0))
				return err;
			if(err = insertVaryingParams(goiP, p))
				return err;
			DisplayWindowXOP1Message(gTheWindow, 
									 WavePoints(p->coefs),
									 goiP->gen_bestfitsofar,
									 *(goiP->chi2Array),
									 goiP->fi.name, 
									 goiP->V_numfititers,
									 goiP->convergenceNumber);
		}	
		//iterate over all the individual members of the population
		for(ii = 0 ; ii < goiP->totalpopsize ; ii += 1){
			// perhaps the user wants to abort the fit using gui panel?
			if(Abort_the_fit){
				err = FIT_ABORTED;
				goto done;
			}
			
			//cmd-dot or abort button
			if(CheckAbort(0)==-1){
				err = FIT_ABORTED;
				goto done;
			}
			
			currentpvector=ii;
			//now set up the trial vector using a wave from the populationvector and bprime
			//first set the pvector 
			// create a mutated trial vector from the best fit and two random population members
			createTrialVector(goiP, p, currentpvector);
			// make sure the trial vector parameters lie between the user defined limits
			ensureConstraints(goiP, p);
			
			chi2pvector = *(goiP->chi2Array + ii);
			/*
			 find out the chi2 value of the trial vector		
			 */
			if(err = setPvector(goiP, goiP->gen_trial, goiP->numvarparams))
				goto done;
			if(err = insertVaryingParams(goiP, p))
				goto done;
			if(err = calcModel(&goiP->fi,
							   goiP->GenCurveFitCoefs,
							   goiP->dataCalc,
							   goiP->dataTemp,
							   goiP->xcalc,
							   goiP->independentVariable,
							   goiP->numVarMD,
							   goiP->isAAO,
							   goiP->sp))
				goto done;
			switch(goiP->METH){//which cost function
				case 0:
					if(err = calcChi2(goiP->dataObs, goiP->dataTemp, goiP->dataSig, goiP->unMaskedPoints, &chi2trial,goiP->weighttype))
						goto done;
					break;
				case 1:
					if(err = calcRobust(goiP->dataObs, goiP->dataTemp, goiP->dataSig, goiP->unMaskedPoints, &chi2trial,goiP->weighttype))
						goto done;
					break;
				case 2:
					if(err = calcUserCostFunc(goiP->minf, 
											  goiP->yobs,
											  goiP->dataObs,
											  goiP->dataCalc,
											  goiP->sobs,
											  goiP->dataSig,
											  goiP->unMaskedPoints,
											  goiP->GenCurveFitCoefs,
											  &chi2trial))
					   goto done;
					break;
				default:
					break;
			}
			
			/*
			 if the /TEMP switch is specified then the user may want to accept the move, a la Monte Carlo.
			*/
			if(p->TEMPFlagEncountered && p->TEMPFlag_opt > 0){
				if( exp(-chi2trial/chi2pvector / p->TEMPFlag_opt) < randomDouble(0, 1, 0, 0))
					acceptMoveGrudgingly = 1;
				else
					acceptMoveGrudgingly = 0;
				if(p->DUMPFlagEncountered && ii > 0){
					WriteMemoryCallback(goiP->gen_pvector, sizeof(double), goiP->numvarparams, &dumpRecord);
					if(dumpRecord.memory == NULL){
						err = NOMEM; goto done;
					}
				}
			}
			
			/*
			 if the chi2 of the trial vector is less than the current populationvector then replace it
			 */
			if(chi2trial < chi2pvector || (acceptMoveGrudgingly && ii)){
				if(err = setPopVectorFromPVector(goiP, goiP->gen_pvector, goiP->numvarparams, currentpvector))
					goto done;
				*(goiP->chi2Array+ii) = chi2trial;
				/*
				 if chi2 of the trial vector is less than that of the best fit, then replace the best fit vector
				 */
				if(chi2trial < *(goiP->chi2Array)){		//if this trial vector is better than the current best then replace it
					if(err = setPopVectorFromPVector(goiP, goiP->gen_pvector, goiP->numvarparams,0))
						goto done;
					/*
					 work out the fractional decrease in chi2
					 */
					wavStats = getWaveStats(goiP->chi2Array, goiP->totalpopsize, 1);
					
					/*
					 put a copy of the best fit so far into an array.
					 */
					memcpy(goiP->gen_bestfitsofar, goiP->gen_coefsCopy, sizeof(double) * goiP->totalnumparams);
					
					/*
					 update the best chi2 if you've just found a better fit (but not yet reached termination)
					 */
					*(goiP->chi2Array) = chi2trial;
					
					/*
					 if you're in update mode then update fit curve and the coefficients
					 */
					if(!goiP->noupdate){
						//DisplayWindowXOP1Message calls code in updateXOP<x>.c
						//this gives a window that gives the user the current chi2 value
						//and the number of iterations.
						goiP->convergenceNumber = p->TOLFlag_tol / (wavStats.V_stdev/ wavStats.V_avg);
						DisplayWindowXOP1Message(gTheWindow, 
												 WavePoints(p->coefs), 
												 goiP->gen_bestfitsofar,
												 *(goiP->chi2Array), 
												 goiP->fi.name,
												 goiP->V_numfititers,
												 goiP->convergenceNumber
												 );
						
						if(err = ReturnFit(goiP, p))
							goto done;
					}
					
					/*
					 if the fractional decrease is less than the tolerance abort the fit.
					 */
					if( wavStats.V_stdev/wavStats.V_avg < p->TOLFlag_tol){	//if the fractional decrease is less and 0.5% stop.
						if(p->DUMPFlagEncountered){
							for(mm=0 ; mm < goiP->totalpopsize ; mm+=1){
								WriteMemoryCallback((goiP->gen_populationvector[mm]), sizeof(double), goiP->numvarparams, &dumpRecord);
								if(dumpRecord.memory == NULL){
									err = NOMEM; goto done;
								}
							}
						}
						goto done;
					}
				}
			}
		}
		if(p->DUMPFlagEncountered && (!p->TEMPFlagEncountered)){
			for(mm=0 ; mm < goiP->totalpopsize ; mm+=1){
				WriteMemoryCallback((goiP->gen_populationvector[mm]), sizeof(double), goiP->numvarparams, &dumpRecord);
				if(dumpRecord.memory == NULL){
					err = NOMEM; goto done;
				}
			}
		}
	}
	
done:
	if(dumpRecord.memory){
		err = dumpRecordToWave(goiP,&dumpRecord);
		free(dumpRecord.memory);
	}
	
	return err;
}

/*
 dumpRecordToWave puts the dumped population array into a wave, if the /DUMP flag was specified.
 */
int dumpRecordToWave(GenCurveFitInternalsPtr goiP,	MemoryStruct *dumpRecord){
	int err = 0;
	
	waveHndl dump;
	long dimensionSizes[MAX_DIMENSIONS+1]; // Array of dimension sizes 
	
	memset(dimensionSizes, 0, sizeof(dimensionSizes));
	dimensionSizes[0] = goiP->numvarparams;
	dimensionSizes[1] = dumpRecord->size / (sizeof(double) * goiP->numvarparams);
//	dimensionSizes[2] = dumpRecord->size/(sizeof(double)*goiP->totalpopsize*goiP->numvarparams);
	
	if(err = MDMakeWave(&dump,"M_gencurvefitpopdump",goiP->cDF,dimensionSizes, NT_FP64, 1))
		return err;
	
	if(err = MDStoreDPDataInNumericWave(dump, (double*)dumpRecord->memory))
		return err;
	
	return err;
}

/*
 ReturnFit updates the model fits and coefficients, then informs IGOR
 that we've changed the waves.
 returns 0 if no errors
 returns errorcode otherwise
 */
static int
ReturnFit(GenCurveFitInternalsPtr goiP, GenCurveFitRuntimeParamsPtr p){
	int err = 0;
	
	if(err = MDStoreDPDataInNumericWave(p->coefs, goiP->gen_bestfitsofar))
		return err;
	WaveHandleModified(p->coefs);
	
	switch(goiP->numVarMD){
		case 1:
			if(p->DFlagEncountered && p->XFlagEncountered){
				if(err = calcModelXY(&goiP->fi,p->coefs,goiP->OUT_data,goiP->fullExtentOfData,goiP->numVarMD, goiP->isAAO,goiP->sp))
					return err;
			} else {
				if(err = calcModelXY(&goiP->fi,p->coefs,goiP->OUT_data,goiP->OUT_x,goiP->numVarMD, goiP->isAAO,goiP->sp))
					return err;
			}
			break;
		default:
			if(err = calcModelXY(&goiP->fi,p->coefs,goiP->OUT_data,goiP->fullExtentOfData,goiP->numVarMD, goiP->isAAO,goiP->sp))
				return err;
			break;
	}
	
	WaveHandleModified(goiP->OUT_data);
	
	if(p->RFlagEncountered){
		if(err = calcModelXY(&goiP->fi,p->coefs,goiP->OUT_res,goiP->fullExtentOfData,goiP->numVarMD,goiP->isAAO,goiP->sp))
			return err;
		if(err = subtractTwoWaves(goiP->OUT_res,p->dataWave.waveH))
			return err;
		if(err = scalarMultiply(goiP->OUT_res, -1))
			return err;
		WaveHandleModified(goiP->OUT_res);
	}
	DoUpdate();
done:
	
	return err;
};


/*
 identicalWaves tests whether two waveHandles refer to the same wave.
 returns 0 if no errors
 returns errorcode otherwise
 if(wav1 == wav2) then isSame=1 
 */
int
identicalWaves(waveHndl wav1, waveHndl wav2, int* isSame){
	int err = 0;
	char wav1Name[MAX_WAVE_NAME+1];
	char wav2Name[MAX_WAVE_NAME+1];
	DataFolderHandle df1H,df2H;
	long df1,df2;
	if(wav1 == NULL || wav2 == NULL){
		return 0;
	}
	*isSame = 0;
	WaveName(wav1,wav1Name);
	WaveName(wav2,wav2Name);
	
	if(err = GetWavesDataFolder(wav1,&df1H))
		return err;
	if(err = GetWavesDataFolder(wav2,&df2H))
		return err;
	if(err= GetDataFolderIDNumber(df1H,&df1))
		return err;
	if(err= GetDataFolderIDNumber(df2H,&df2))
		return err;
	
	if(CmpStr(wav1Name,wav2Name)==0 && df1 == df2)
		*isSame = 1;
	
	return err;
}

/*
 isWaveDisplayed tests if wav is displayed in the top graph
 returns 0 if no error
 returns errorcode otherwise
 if(wav is displayed in top graph) then isDisplayed=1
 */
static int
isWaveDisplayed(waveHndl wav, int *isDisplayed){
	char cmd[MAXCMDLEN+1];
	char varName[MAX_OBJ_NAME+1];
	char gwaveName[MAX_WAVE_NAME+1];
	double re=-1,imag=-1;
	int err=0;
	*isDisplayed = 0;
	
	if(wav == NULL)
		return NON_EXISTENT_WAVE;
	
	strcpy(varName, "TEMPGenCurveFit_GLOBALVAR");
	re = -1;
	if(err = SetIgorFloatingVar(varName, &re, 1))
		return err;
	strcpy(cmd,"TEMPGenCurveFit_GLOBALVAR=whichlistitem(\"");
	
	WaveName(wav,gwaveName);
	strcat(cmd,gwaveName);
	strcat(cmd,"\",(tracenamelist(winname(0,1),\";\",1)))");
	if(err = XOPSilentCommand(&cmd[0]))
		return err;
	if(FetchNumVar(varName, &re, &imag)==-1)
		return EXPECTED_VARNAME;
	if(re != -1) 
		*isDisplayed = 1;
	strcpy(cmd,"Killvariables/z TEMPGenCurveFit_GLOBALVAR");
	if(err = XOPSilentCommand(&cmd[0]))
		return err;
	return 0;
}
/*
 arraySd returns the standard deviation of a pointer to an array of doubles
 */
static double
arraySD(double* data, long datasize){
	long ii=0;
	double sd = 0;
	double nx2=0,nx=0,mean=0;
	for(ii=0;ii<datasize;ii+=1){
		nx += (*(data+ii));
	}
	mean = nx/(double)datasize;
	for(ii=0;ii<datasize;ii+=1){
		nx2 += pow((*(data+ii)-mean),2);
	}
	sd = nx2/(double)datasize;
	sd = sqrt(sd);
	
	return sd;
}
/*
 arrayMean returns the arithmetic mean of a pointer array
 */
static double
arrayMean(double* data, long datasize){
	long ii=0;
	double mean = 0;
	double nx=0; 
	for(ii=0;ii<datasize;ii+=1){
		nx+= (*(data+ii));
	}
	mean = nx/(double)datasize;
	
	return mean;
}
static long
numInArray3SD(double* data, double sd, long datasize){
	double mean = arrayMean(data,datasize);
	
	long ii=0;
	long num=0;
	for(ii=0;ii<datasize;ii+=1){
		if(abs(*(data+ii)-mean)<(3*sd))
			num+=1;
	}
	return num;
}
static int
getRange (WaveRange range,long *startPoint,long *endPoint){
	int direction;
	int err = 0;
	*startPoint = (double)range.startCoord;
	*endPoint = (double)range.endCoord;
	direction = 1;
	if (range.rangeSpecified) {
		WaveRangeRec wr;
		MemClear(&wr,sizeof(WaveRangeRec));
		wr.x1 = range.startCoord;
		wr.x2 = range.endCoord;
		wr.rangeMode = 3;
		wr.isBracket = range.isPoint;
		wr.gotRange = 1;
		wr.waveHandle = range.waveH;
		wr.minPoints = 2;
		if (err = CalcWaveRange(&wr))
			return err;
		*startPoint = wr.p1;
		*endPoint = wr.p2;
		direction = wr.wasBackwards ? -1:1;
	}
	return err;
}

/*
 roundDouble returns a rounded value for val
 */
static double
roundDouble(double val){
	double retval;
	if(val>0){
		if(val-floor(val) < 0.5){
			retval = floor(val);
		} else {
			retval = ceil(val);
		}
	} else {
		if(val-floor(val) <= 0.5){
			retval = floor(val);
		} else {
			retval = ceil(val);
		}
	}
	return retval;
}


int
WindowMessage(void){
	long item0;										// Item from IO list.
	long message;
	
	message = GetXOPMessage();
	
	if ((message & XOPWINDOWCODE) == 0)
		return 0;
	
	item0 = GetXOPItem(0);
	
	switch (message) {
#ifdef _MACINTOSH_			// [
		case UPDATE:
		{
		}
			break;			
		case ACTIVATE:
		{
		}
			break;			
		case CLICK:
		{
		}
			break;
#endif						// _MACINTOSH_ ]			
		case CLOSE:
			HideAndDeactivateXOPWindow(gTheWindow);
			break;
			
		case NULLEVENT:				// This is not sent on Windows. Instead, similar processing is done in response to the WM_MOUSEMOVE message.
			ArrowCursor();
			break;
	}												// Ignore other window messages.
	return 1;
}

static waveStats getWaveStats(double *sort, long length,int moment){
	long ii=0;
	double minval = *sort, maxval = *sort;
	long minpos=0,maxpos=0;
	double nx2=0,nx=0;
	struct waveStats retval;
	
	switch(moment){
		case 0:
			for(ii=0;ii<length;ii+=1){
				if(*(sort+ii)>maxval){
					maxval = *(sort+ii);
					maxpos = ii;
				}
				if(*(sort+ii)<minval){
					minval = *(sort+ii);
					minpos = ii;
				}
			}
			retval.V_maxloc = maxpos;
			retval.V_minloc = minpos;
			break;
		case 1:
			for(ii=0;ii<length;ii+=1){
				nx += (*(sort+ii));
				nx2 += *(sort+ii)*(*(sort+ii));
				if(*(sort+ii)>maxval){
					maxval = *(sort+ii);
					maxpos = ii;
				}
				if(*(sort+ii)<minval){
					minval = *(sort+ii);
					minpos = ii;
				}
			}
			retval.V_maxloc = maxpos;
			retval.V_minloc = minpos;
			retval.V_avg = nx/(double)length;
			retval.V_stdev = sqrt((nx2/(double)length)-(retval.V_avg*retval.V_avg));
			break;
	}
	return retval;
}



