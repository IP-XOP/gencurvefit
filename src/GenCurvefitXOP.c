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
#include "GenCurveFitXOP.h"
#include "gencurvefit.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*
 a variable which sees if the fit is to be aborted.
 */
int Abort_the_fit = 0;

/*
 some memory utilities
 */

//create a platform independent routine for continuous reallocation of memory, appending data to it
void *myrealloc(void *src_ptr, size_t size)
{
    /* There might be a realloc() out there that doesn't like reallocing
	 NULL pointers, so we take care of it here */
    if(src_ptr)
		return realloc(src_ptr, size);
    else
		return malloc(size);
}

size_t
WriteMemoryCallback(void *ptr, size_t size, size_t nmemb, void *data)
{
    size_t realsize = size * nmemb;
    struct MemoryStruct *mem = (struct MemoryStruct *)data;
	
    mem->memory = (char *)myrealloc(mem->memory, mem->size + realsize + 1);
    if (mem->memory) {
		memcpy(&(mem->memory[mem->size]), ptr, realsize);
		mem->size += realsize;
		mem->memory[mem->size] = 0;
    }
    return realsize;
}

/*
 the fitfunction for the genetic optimisation
 */
int lgencurvefit_fitfunction(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long datapoints, unsigned int numDataDims){
	int err = 0;
	
	CountInt ii;
	unsigned int jj;
	int requiredParameterTypes[MAX_MDFIT_SIZE + 2];
	allFitFunc allParameters;
	fitFunc parameters;
	double result;
	double *dp;
    GenCurveFitInternals *goiP = (GenCurveFitInternals*) userdata;

	memset(&parameters, 0, sizeof(parameters));
	memset(&allParameters, 0, sizeof(allParameters));

	//cmd-dot or abort button
	if(CheckAbort(0) == -1){
		err = FIT_ABORTED;
		goto done;
	}

	/*
	 check if all the input isn't NULL
	 */
	if(goiP->coefs == NULL || goiP->dataCalc == NULL){
		err = UNSPECIFIED_ERROR;
		goto done;
	}
	
	/*
	 Place the parameters into the coefficient wave
	 */
    dp = WaveData(goiP->GenCurveFitCoefs);
    memcpy(dp, coefs, sizeof(double) * numcoefs);

//	if(err = MDStoreDPDataInNumericWave(goiP->GenCurveFitCoefs, (double*) coefs))
//		goto done;
			
	switch(goiP->isAAO){
		case 0:
			parameters.waveH = goiP->GenCurveFitCoefs;
			requiredParameterTypes[0] = WAVE_TYPE;
			
			for(ii = 0 ; ii < datapoints ; ii += 1){
				for(jj = 0 ; jj < numDataDims ; jj += 1)
					parameters.x[jj] = xdata[jj][ii];
				
				// call the users fit function and put the result in the output array
				if (err = CallFunction(&goiP->fi, (void*) &parameters, &result))
					goto done;
				model[ii] = result;
			}
						
			break;
		case 1:
			allParameters.waveC = goiP->GenCurveFitCoefs;
			allParameters.waveY = goiP->dataCalc;
			requiredParameterTypes[0] = WAVE_TYPE;
			requiredParameterTypes[1] = WAVE_TYPE;
			
			for(ii = 0 ; ii < numDataDims ; ii += 1){
				requiredParameterTypes[ii + 2] = WAVE_TYPE;
				if(goiP->xcalc[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				allParameters.waveX[ii] = goiP->xcalc[ii];
			}

			// call the users fit function and put the result in the output wave
			if (err = CallFunction(&goiP->fi, (void*)&allParameters, &result))
				goto done;

			// the user may have changed the number of points in the output wave
			if(goiP->dataCalc == NULL || WavePoints(goiP->dataCalc ) != datapoints){
				err = USER_CHANGED_FITWAVE;
				goto done;
			}
			dp = WaveData(goiP->dataCalc);
			memcpy(model, dp, sizeof(double) * datapoints);
			
			break;
		case 2:			
			if(goiP->sp == NULL){
				err = NULL_STRUCTURE;
				goto done;
			}
			goiP->sp->w = goiP->GenCurveFitCoefs;
			goiP->sp->yy = goiP->dataCalc;
			
			
			for(ii = 0 ; ii < numDataDims ; ii+=1){
				if(goiP->xcalc[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				goiP->sp->xx[ii] = goiP->xcalc[ii];
			}

			if (err = CallFunction(&goiP->fi, (fitfuncStruct*)&goiP->sp, &result))
				goto done;
			//don't want any dangling references to waves
			goiP->sp->yy = NULL;
			goiP->sp->w = NULL;
			memset(goiP->sp->xx, 0, sizeof(goiP->sp->xx));
			
			// the user may have changed the number of points in the output wave
			if(goiP->dataCalc == NULL || WavePoints(goiP->dataCalc) != datapoints){
				err = USER_CHANGED_FITWAVE;
				goto done;
			}
			
			dp = WaveData(goiP->dataCalc);
			memcpy(model, dp, sizeof(double) * datapoints);
			
			break;
		default:
			err = UNSPECIFIED_ERROR;
			goto done;
			break;
	}
	
	// check that the fitfunction didn't return any NaN or INF
	if(err = checkNanInfArray(model, datapoints)){
		err = FITFUNC_RETURNED_NANINF;
		goto done;
	}

done:		

	return err;
}

/*
 the cost function for the genetic optimisation
 */
double lgencurvefit_costfunction(void *userdata, const double *coefs, unsigned int numcoefs, const double *data, const double *model, const double *errors, long datapoints){
	double val;

	GenCurveFitInternals *goiP = (GenCurveFitInternals*) userdata;
	costFunc userCostFunc;
	double *dp;
	int err;
	
	switch(goiP->METH){//which cost function
		case 0:
			val = chisquared(userdata, coefs, numcoefs, data, model, errors, datapoints);			
			break;
		case 1:
			val = robust(userdata, coefs, numcoefs, data, model, errors, datapoints);
			break;
		case 2:
			userCostFunc.coefs = goiP->GenCurveFitCoefs;
			userCostFunc.yobs =	 goiP->yobs;
			userCostFunc.ycalc = goiP->dataCalc;
			userCostFunc.sobs =  goiP->sobs;
			
			//copy the original data into the yobs wave created for the purpose
			//some sneaky users probably try to change it.
			dp = WaveData(goiP->yobs);
			memcpy(dp, data, datapoints * sizeof(double));
			
			//copy the original data into the sobs wave created for the purpose
			//some sneaky users probably try to change it.
			dp = WaveData(goiP->sobs);
			memcpy(dp, errors, datapoints * sizeof(double));

			//copy the original data into the sobs wave created for the purpose
			//some sneaky users probably try to change it.
			dp = WaveData(goiP->dataCalc);
			memcpy(dp, model, datapoints * sizeof(double));
            
			if(err = CallFunction(&goiP->minf, (void*) &userCostFunc, &val))
				goto done;
			
			if(WavePoints(goiP->yobs) != datapoints || WavePoints(goiP->sobs) != datapoints || WavePoints(goiP->GenCurveFitCoefs)!= numcoefs){
				err = COSTFUNC_WAVES_CHANGED;
				goto done;
			}
			
			if(IsNaN64(&val) || IsINF64(&val)){
				err = COSTFUNC_DOESNT_RETURN_NUMBER;
				goto done;
			}
			
			break;
		default:
			val = chisquared(userdata, coefs, numcoefs, data, model, errors, datapoints);			
			break;
	}
	done:
			
	return val;
};


/*
 the update function for the genetic optimisation
 */
int lgencurvefit_updatefunction(void *userdata,
								const double *coefs,
								 unsigned int numcoefs,
								 unsigned int iterations,
								 double cost,
								 unsigned int updatetime,
								 double convergenceNumber,
								 const double** population,
								 const unsigned int *varparams,
								 unsigned int numvarparams,
								 unsigned int totalpopsize,
								 const double *costmap){
	int err = 0;
	double val = 0;
	
	GenCurveFitInternals *goiP = (GenCurveFitInternals*) userdata;
	updtFunc userUpdateFunc;
	
	goiP->V_numfititers = iterations;
	
	/*
	 Display the coefficients so far.
	 */

	// perhaps the user wants to abort the fit using gui panel?
	if(Abort_the_fit){
		err = FIT_ABORTED;
		goto done;
	}

	//We want a record of all the improvements in the fit
	if(goiP->dump && (updatetime == 2 || updatetime == 1)){
		WriteMemoryCallback((void*) coefs, sizeof(double), goiP->totalnumparams, &(goiP->dumpRecord));
		if(goiP->dumpRecord.memory == NULL){
			err = NOMEM; goto done;
		}
	}
	
	//store the best fit coefficients in the original coef wave and update the output data.
	if((updatetime == 1 && !goiP->noupdate) || updatetime == 16){
		goiP->cost = cost;

		if(err = MDStoreDPDataInNumericWave(goiP->OUT_coefs, (void*) coefs))
			return err;
		WaveHandleModified(goiP->OUT_coefs);
		
		if(err = calcModelXY(&goiP->fi, goiP->OUT_coefs, goiP->OUT_data, goiP->OUT_x, goiP->numVarMD, goiP->isAAO, goiP->sp))
			return err;
		
		WaveHandleModified(goiP->OUT_data);
		
		if(goiP->OUT_res){
			if(err = calcModelXY(&goiP->fi, goiP->OUT_coefs, goiP->OUT_res, goiP->fullExtentOfData, goiP->numVarMD, goiP->isAAO, goiP->sp))
				return err;
			if(err = subtractTwoWaves(goiP->OUT_res, goiP->originalYobs))
				return err;
			if(err = scalarMultiply(goiP->OUT_res, -1))
				return err;
			WaveHandleModified(goiP->OUT_res);
		}
		if(RunningInMainThread() && !goiP->noupdate)
		   DoUpdate();
	}
	
	//the user defined a user IGOR update function so call that.
	if(!goiP->noupdate && goiP->useIgorUpdateFunction && population && costmap){
		userUpdateFunc.currentbestfit = goiP->OUT_coefs;
		userUpdateFunc.population = goiP->M_population;
		userUpdateFunc.costmap = goiP->W_costmap;
		userUpdateFunc.updatetime = updatetime;

		memcpy(WaveData(goiP->M_population), &(population[0][0]), totalpopsize * numvarparams * sizeof(double));
		memcpy(WaveData(goiP->W_costmap), costmap, totalpopsize * sizeof(double) );

		WaveHandleModified(goiP->M_population);
		WaveHandleModified(goiP->W_costmap);

		if(err = CallFunction(&goiP->igorUpdateFunction, (void*) &userUpdateFunc, &val))
			goto done;
		if(val)
			err = (int)val;
	}
	
	
done:
	return err;
}

/* 
 function obtains the covariance matrix for the fit.
 */
int getGCovarianceMatrix(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	int err = 0;
	double hessianDeterminant = 0;
	
	if(err = getCovarianceMatrix(&goiP->covarianceMatrix,
								 &hessianDeterminant,
								 goiP,
								 lgencurvefit_fitfunction,
                                 lgencurvefit_costfunction,
								 goiP->coefs,
								 (int) goiP->totalnumparams,
								 goiP->holdvector,
								 goiP->dataObs,
								 goiP->dataSig,
								 (const double**) goiP->independentVariable,
								 (long) goiP->unMaskedPoints,
								 (int) goiP->numVarMD,
								 !p->WFlagEncountered))
		return err;
				
	goiP->V_logBayes = exp(-0.5 * (goiP->cost) / (double)(goiP->unMaskedPoints - goiP->numvarparams));// * pow(4*3.14159,(double) goiP->numvarparams) ;//* factorial((double)goiP->numvarparams);
	goiP->V_logBayes = goiP->V_logBayes / (sqrt(hessianDeterminant));
	goiP->V_logBayes = log(goiP->V_logBayes);

	return err;
}

/*
 converts a libgencurvefit error to an XOP error
 All libgencurvefit errors are negative, all XOP errors are positive
 */

int convertlgencurvefitErrorToXOP(int lgencurvefiterror){
	int err = lgencurvefiterror;
	
	if(lgencurvefiterror > -1)
		return err;
	
	switch (lgencurvefiterror) {
		case NO_MEMORY:
			err = NOMEM;
			break;
		case INCORRECT_LIMITS:
			err = LIMITS_INVALID;
			break;
		case HOLDVECTOR_COEFS_MISMATCH:
			err = HOLDSTRING_NOT_RIGHT_SIZE;
			break;
		case NO_VARYING_PARAMS:
			err = ALL_COEFS_BEING_HELD;
			break;
		case WRONG_NUMBER_OF_PARAMS:
			break;
		case COEFS_MUST_BE_WITHIN_LIMITS:
			err = LIMITS_INVALID;
			break;
		case PROBLEM_CALCULATING_COVARIANCE:
			err = SINGULAR_MATRIX;
			break;
		case NO_FIT_FUNCTION_SPECIFIED:
			err = FITFUNC_NOT_SPECIFIED;
			break;
		case NO_Y_ARRAY:
			break;
		case NO_X_ARRAY:
			break;
		case NO_E_ARRAY:
			break;
		case NO_COEFS_ARRAY:
			break;
		case NO_LIMITS_ARRAY:
			err = LIMITS_INVALID;
			break;
		case SINGULAR_MATRIX_ERROR:
			err = SINGULAR_MATRIX;
			break;
		default:
			break;
	}
	
	
	return err;
}


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
	
	/*
	 err carries the errors for all the operations
	 err2 carries the error code if we still want to return err, but we need to finish off
	 something else first.	
	 */
	int err = 0, err2 = 0;
	double value[2];
    int numDimensions;
	CountInt indices[MAX_DIMENSIONS + 1];
    int wtype = 0;
	
	//libgencurvefit options
	gencurvefitOptions gco;
	
	//variables listed below are purely for outputs sake.
	char varname[MAX_OBJ_NAME + 1];
	CountInt dimensionSizes[MAX_DIMENSIONS + 1];
	double t1,t2, chi2;
	double *wP;
    int updateStatus;
	long lt1 = 0;
	char note[200], note_buffer1[MAX_WAVE_NAME + 1], note_buffer2[MAX_WAVE_NAME + 1], cmd[MAXCMDLEN + 1];
	int output, ii, jj, isDisplayed, quiet;
	
	//initialise all the internal data structures to NULL
	memset(&goi, 0, sizeof(goi));
	
	//you have to use IGOR > 9.00
	if( igorVersion < 900 )
		return REQUIRES_IGOR_900;

    //reset the options for libgencurvefit
	memset(&gco, 0, sizeof(gencurvefitOptions));
	
	//reset the abort condition
	Abort_the_fit = 0;
	
	strncpy(varname, "V_Fiterror", MAX_OBJ_NAME);
	if(FetchNumVar(varname, &t1, &t2)!=-1){
		if(!err)
			lt1 = 0;
		else 
			lt1 = 1;
		
		err = 0;
		if(err2 = SetIgorIntVar(varname, lt1, 1)){
            err = err2;
            goto done;
            
        };
	}	
	
	/*
	 Genetic Optimisation uses a lot of random numbers, we are seeding the generator here.
	 If you want to seed the generator we can 
	 */
	if(p->SEEDFlagEncountered)
		gco.seed = (int) p->SEEDFlag_seed;
	else 
		gco.seed = -1;
	
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
	
	if(p->OPTFlagEncountered && (((long)p->OPTFlag_opt) & (long)pow(2, 0)))
		gco.useinitialguesses = 1;
    
    
    if(p->POPFlagEncountered && p->initial_popwave != NULL){
        //check how many points are in the wave
        wtype = WaveType(p->initial_popwave);
        if(wtype != NT_FP64){
            err = INITIAL_POPULATION_DP;
            goto done;
        }
        
        if(err = MDGetWaveDimensions(p->initial_popwave, &numDimensions, indices))
            goto done;
        if(numDimensions != 2 || indices[0] != goi.totalnumparams || indices[1] < 1){
            err = INCORRECT_INITIAL_POPULATION;
            goto done;
        }
        gco.initial_population_rows = indices[1];
        gco.initial_population = (double *)WaveData(p->initial_popwave);
    }
	
	if(p->DITHFlagEncountered && p->DITHFlagParamsSet[0] && p->DITHFlagParamsSet[1]){
        gco.dither[0] = p->dith1;
        gco.dither[1] = p->dith2;
    } else {
        gco.dither[0] = -1.;
        gco.dither[1] = -1.;
    }
    
	gco.updatefun = (updatefunction) &lgencurvefit_updatefunction;
	gco.k_m = goi.k_m;
	gco.recomb = goi.recomb;
	gco.strategy = goi.STGY;
	gco.tolerance = goi.tolerance;
	gco.popsizeMultiplier = goi.popsize;
	gco.iterations = goi.iterations;
	gco.updatefrequency = 31;
	
	/*
	 you want to polish the fit (at the end) using LevenbergMarquardt
	 */
	if(p->POLFlagEncountered)
		gco.polishUsingLM = 1;
	
	/*
	 if you want to do a Monte Carlo fit.
	 */
	if(p->MCFlagEncountered && p->WFlagEncountered)
		gco.monteCarlo = 1;
	
	//make an error wave before we start the fit, to make sure it's always there
	dimensionSizes[0] = goi.totalnumparams;
	dimensionSizes[1] = 0;
	if(err = MDMakeWave(&goi.W_sigma, "W_sigma", goi.cDF, dimensionSizes, NT_FP64, 1))
			goto done;
	for(ii = 0 ; ii < goi.totalnumparams ; ii+=1){
		indices[0] = ii;
		value[0] = 0;
		if(err = MDSetNumericWavePointValue(goi.W_sigma, indices, value))
			goto done;
	}
	WaveHandleModified(goi.W_sigma);
		
//	if(p->MATFlagEncountered && !(p->MATFlagParamsSet[0] && (int) p->MATFlag_mat == 0)){
    //make the covariance matrix
    dimensionSizes[0] = goi.totalnumparams;
    dimensionSizes[1] = goi.totalnumparams;
    dimensionSizes[2] = 0;
    if(err = MDMakeWave(&goi.M_covariance, "M_Covar", goi.cDF, dimensionSizes, NT_FP64, 1))
            goto done;
        
    wP = WaveData(goi.M_covariance);
        
    memset(wP, 0, sizeof(double) * goi.totalnumparams * goi.totalnumparams);
        
    WaveHandleModified(goi.M_covariance);
//	}
	
	/*
	 optimiseloop does the Differential Evolution, according to Storn and Price.  When this returns 0, then 
	 you have the best fit in the GenCurveFitInternals structure.  Otherwise it returns an error code.  If the user aborts
	 then the FIT_ABORTED error code is returned, but it is still possible to retrieve the best solution so far
	 */
	if(err = genetic_optimisation(lgencurvefit_fitfunction,
							   lgencurvefit_costfunction,
							   (unsigned int) goi.totalnumparams,
							   goi.coefs,
							   goi.holdvector,
							   (const double**) goi.limits,
							   (long) goi.unMaskedPoints,
							   goi.dataObs,
							   (const double**) goi.independentVariable,
							   goi.dataSig,
							   goi.numVarMD,
							   &chi2,
							   &gco,
							   &goi))
		goto done;
	
	//make sure the coefficients for the fit are updated, by calling the update function
	Abort_the_fit = 0;
	if(err = lgencurvefit_updatefunction(&goi,
										 goi.coefs,
										 (unsigned int) goi.totalnumparams,
										 goi.V_numfititers,
										 chi2,
										 16,
										 1,
										 NULL,
										 NULL,
										 0,
										 0,
										 NULL))
	   goto done;
	
	//output the dumprecord
	if(!err && goi.dump && goi.dumpRecord.memory)
		if(err2 = dumpRecordToWave(&goi, &goi.dumpRecord)){
			err = err2;
			goto done;
		}
	
	/*
	 if there are no errors, or if the user aborted, then return the best fit.
	 If the data is displayed in the top graph append the fitcurve to the top graph
	 */
	if((err == 0 || err == FIT_ABORTED) && !getGCovarianceMatrix(p, &goi) && goi.covarianceMatrix){
		//set the error wave
		jj = 0;
		for(ii = 0; ii < goi.totalnumparams ; ii += 1){
			indices[0] = ii;
			if(!goi.holdvector[ii]){
				value[0] = sqrt(goi.covarianceMatrix[goi.varParams[jj]][goi.varParams[jj]]);
				jj++;
			} else {
				value[0] = 0;
			}

			if(err = MDSetNumericWavePointValue(goi.W_sigma, indices, value))
				goto done;
			
		}
		WaveHandleModified(goi.W_sigma);
		
//		if(p->MATFlagEncountered && !(p->MATFlagParamsSet[0] && (int) p->MATFlag_mat == 0)){
        //make the covariance matrix
        dimensionSizes[0] = goi.totalnumparams;
        dimensionSizes[1] = goi.totalnumparams;
        dimensionSizes[2] = 0;
        wP = WaveData(goi.M_covariance);			
        memcpy(wP, *goi.covarianceMatrix, sizeof(double) * goi.totalnumparams * goi.totalnumparams);
        WaveHandleModified(goi.M_covariance);
//		}
	}
	//append the fit to the graph
	if((!err || err == FIT_ABORTED) && RunningInMainThread()){
		if(err = isWaveDisplayed(p->dataWave.waveH, &isDisplayed))
			goto done;
		
		if(isDisplayed && goi.numVarMD == 1){
			if(err = isWaveDisplayed(goi.OUT_data, &isDisplayed))
				goto done;
			
			if(!isDisplayed){
				strncpy(cmd, "appendtograph/w=$(winname(0,1)) ", MAXCMDLEN);
				WaveName(goi.OUT_data,&note_buffer1[0]);
				strncat(cmd, &note_buffer1[0], MAXCMDLEN - strlen(note_buffer1));
				
				if(p->DFlagEncountered && p->XFlagEncountered){
					WaveName(p->XFlag_xx,&note_buffer2[0]);
					strncat(cmd, " vs ", MAXCMDLEN - strlen(cmd) - strlen(" vs "));
					strncat(cmd, note_buffer2, MAXCMDLEN - strlen(cmd) - strlen(note_buffer2) );
				}
                PauseUpdate(&updateStatus);
				if(err = XOPSilentCommand(&cmd[0])){
                    ResumeUpdate(&updateStatus);
                    goto done;
                }
                ResumeUpdate(&updateStatus);
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
		SetOperationNumVar("V_Chisq", chi2);
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
						  "V_fitIters = %li; V_Chisq = %g; V_npnts= %lld; V_nterms= %lld; V_nheld= %lld; V_logBayes = %g\r",
						  goi.V_numfititers,
						  chi2,
						  (SInt64) goi.unMaskedPoints,
						  (SInt64) WavePoints(p->coefs),
						  (SInt64) WavePoints(p->coefs) - goi.numvarparams,
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
							  goi.coefs[ii],
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
	if(err)
		err = convertlgencurvefitErrorToXOP(err);
    
	freeAllocMem(&goi);
	return err;
}

int
//this will return if an array contains NaN or INF.
checkNanInfArray(double *array, CountInt datapoints){
	//this check examines to see if there are any NaN/Infs in a wave
	//this is really important if you want to calculate Chi2.
	int err = 0;
	CountInt ii;
	
	for(ii = 0 ; ii < datapoints ; ii += 1)
		if((err = IsNaN64(array + ii)) || (err = IsINF64(array + ii)))
		   break;
	
	if(err > 0)
		err = INPUT_WAVES_CONTAINS_NANINF;

	return err;
}

int
//this will return if a wave contains NaN or INF.
checkNanInf(waveHndl wav){
	//this check examines to see if there are any NaN/Infs in a wave
	//this is really important if you want to calculate Chi2.
	int err = 0;
	CountInt ii;
	double *dp = NULL;
	CountInt points;
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
checkZeros(waveHndl wavH, long* numzeros){
	int result = 0;
	size_t numBytes;
	double* dp = NULL;
	CountInt points, ii;
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
 init_GenCurveFitInternals initialises the GenCurveFitInternals structure
 returns 0 if no error
 returns errorcode otherwise
 */
static int
init_GenCurveFitInternals(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	int err = 0;
	char *holdstr = NULL;
	CountInt ii, jj, kk, minpos, maxpos;
	CountInt dimensionSizes[MAX_DIMENSIONS+1];
	int numdimensions;
	double value[2];
	CountInt indices[MAX_DIMENSIONS];
	char xwavename[MAX_WAVE_NAME + 1];
	char datawavename[MAX_WAVE_NAME + 1];
	char reswavename[MAX_WAVE_NAME + 1];
	char datawavestring[MAX_WAVE_NAME + 1];
	char cmd[MAXCMDLEN];
	char letter[3];
	int toDisplay = 0;
	int outPutType = NT_FP64;
	double *dp;
    int updateStatus;
	CountInt temp;
	double temp1;
	DataFolderHandle tempWavesDFH;
	
	waveHndl gcf_dataCalc, gcf_yobs, gcf_sobs, gcf_xcalc[MAX_MDFIT_SIZE], gcf_fullExtentOfData[MAX_MDFIT_SIZE], gcf_GenCurveFitCoefs;
	waveHndl gcf_tempWaveHndl_OUTx, gcf_W_costmap, gcf_M_population;
	
    tempWavesDFH = (DataFolderHandle) -1;

	//do we want dynamic updates?
	goiP->noupdate = 0;
	if(p->NFlagEncountered){
		goiP->noupdate = 1;
		if(p->NFlagParamsSet[0] && (int) p->NFlag_noupdate == 0)
			goiP->noupdate = 0;
	}
	
	//if you're not running in the main thread don't do updates
	if(!RunningInMainThread())
		goiP->noupdate = 1;

	
	/* get a full copy of the datawave */
	goiP->dataObsFull = (double*)malloc(WavePoints(p->dataWave.waveH) * sizeof(double));
	if (goiP->dataObsFull == NULL){
		err = NOMEM;
		goto done;
	}
	//get a copy of all the datawave.  This is so we can fill dataobs
	if (err = MDGetDPDataFromNumericWave(p->dataWave.waveH, goiP->dataObsFull)) // Get copy.
		goto done;
	
	goiP->originalYobs = p->dataWave.waveH;
	
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
			goiP->startPoint = (CountInt)roundDouble(p->dataWave.startCoord);
			goiP->endPoint = (CountInt)roundDouble(p->dataWave.endCoord);
			if(goiP->startPoint > goiP->endPoint){
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
		for(ii = 0 ; ii < WavePoints(p->dataWave.waveH) ; ii += 1)
			*(goiP->mask + ii) = 1;
	}
	/* 
	 set up the mask array
	 need to correct goiP->unMaskedPoints as we go along, which specifies how many unmasked points there will be in the fit 
	 */
	for(ii = 0 ; ii < WavePoints(p->dataWave.waveH) ; ii += 1){
		temp1 = *(goiP->mask+ii);
		if(*(goiP->mask + ii) == 0 || IsNaN64(goiP->mask + ii) || ii < goiP->startPoint || ii>goiP->endPoint || IsNaN64(goiP->dataObsFull + ii) || IsINF64(goiP->dataObsFull + ii)){
			goiP->unMaskedPoints -= 1;
			*(goiP->mask + ii)=0;
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
	
	
	if(err = MDMakeWave(&gcf_dataCalc,"GenCurveFit_dataCalc", tempWavesDFH, dimensionSizes, NT_FP64, 1))
		goto done;
	if(err = MDMakeWave(&gcf_yobs, "GenCurveFit_yobs", tempWavesDFH, dimensionSizes, NT_FP64, 1))
		goto done;
	if(err = MDMakeWave(&gcf_sobs, "GenCurveFit_sobs", tempWavesDFH, dimensionSizes,NT_FP64, 1))
		goto done;
	
	if(err = HoldWave(gcf_dataCalc))
		goto done;
	if(err = HoldWave(gcf_yobs))
		goto done;
	if(err = HoldWave(gcf_sobs))
		goto done;
	goiP->dataCalc= gcf_dataCalc;
	goiP->yobs = gcf_yobs;
	goiP->sobs = gcf_sobs;
	
	if(goiP->isAAO){
		for(ii=0 ; ii<goiP->numVarMD ; ii+=1){
			sprintf(letter,"%li",ii);
			strcpy(xwavename,"GenCurveFit_xcalc");
			strcat(&xwavename[0],&letter[0]);
			if(err = MDMakeWave(&gcf_xcalc[ii], xwavename, tempWavesDFH, dimensionSizes,NT_FP64, 1))
				goto done;
			if(err = HoldWave(gcf_xcalc[ii]))
				goto done;
			goiP->xcalc[ii] = gcf_xcalc[ii];
		}
	}
	
	/*
	 create a utility wave that will contains the x range of the original ywave
	 */
	dimensionSizes[0] = goiP->dataPoints;
	for(ii=0 ; ii<goiP->numVarMD ; ii+=1){
		strcpy(letter, "");
		sprintf(letter, "%li", ii);		
		strcpy(xwavename, "GenCurveFit_fullExtentOfData0");
		strcat(&xwavename[0], &letter[0]);
		if(err = MDMakeWave(&gcf_fullExtentOfData[ii], xwavename, tempWavesDFH, dimensionSizes, NT_FP64, 1))
			goto done;
		if(err = HoldWave(gcf_fullExtentOfData[ii]))
			goto done;
		goiP->fullExtentOfData[ii] = gcf_fullExtentOfData[ii];
	}
	
	
	/*
	 make the temporary coefficients in the current datafolder
	 */
	dimensionSizes[0] = WavePoints(p->coefs);
	if(err = MDMakeWave(&gcf_GenCurveFitCoefs, "GenCurveFit_coefs", tempWavesDFH, dimensionSizes, NT_FP64, 1))
		goto done;
	goiP->OUT_coefs = p->coefs;
	if(err = HoldWave(gcf_GenCurveFitCoefs))
		goto done;
	goiP->GenCurveFitCoefs = gcf_GenCurveFitCoefs;
	
	/*
	 initialise space for an array containing the unmasked fitpoint ydata 
	 */
	goiP->dataObs = (double*)malloc(goiP->unMaskedPoints * sizeof(double));
	if (goiP->dataObs == NULL){
		err = NOMEM;
		goto done;
	}
	
	//now fill up the dataObs array, if a point isn't being masked
	jj = 0;
	for(ii = 0 ; ii < WavePoints(p->dataWave.waveH) ; ii += 1){
		if(!(*(goiP->mask + ii) == 0 || IsNaN64(goiP->mask + ii))){
			*(goiP->dataObs + jj) = *(goiP->dataObsFull + ii);
			jj += 1;
		}
	}
	//copy those unmasked points into a dedicated wave, only useful for user specified cost functions
	if(err = MDStoreDPDataInNumericWave(goiP->yobs, goiP->dataObs))
	   goto done;
	
	//a temporary array used in calculating parameter uncertainties.
	goiP->dataTemp = (double*) malloc(sizeof(double) * goiP->unMaskedPoints);
	if(!goiP->dataTemp){
		err = NOMEM;
		goto done;
	}
	
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
		jj = 0;
		for(ii = 0 ; ii < goiP->dataPoints ; ii += 1){
			if(!(goiP->mask[ii] == 0 || IsNaN64(goiP->mask + ii))){
				indices[0] = ii;
				if(err = MDGetNumericWavePointValue(p->WFlag_weighttype, indices, &temp1))
					goto done;
				
				if((int)goiP->weighttype == 0)
					temp1 = 1 / temp1;
				
				goiP->dataSig[jj] = temp1;
				jj++;
			}
		}
	} else {
		for(ii = 0 ; ii < goiP->unMaskedPoints ; ii++)
			goiP->dataSig[ii] = 1;
	}
	
	//copy it into the sobs wave (only useful for user specified cost functions)
	dp = WaveData(goiP->sobs);
	memcpy(dp, goiP->dataSig, sizeof(double) * goiP->unMaskedPoints);
	
	
	//initialise array space for x values
	goiP->independentVariable = (double**)malloc2d(goiP->numVarMD, (int) goiP->unMaskedPoints, sizeof(double));
	if (goiP->independentVariable == NULL){
		err = NOMEM;
		goto done;
	}
	
	goiP->allIndependentVariable = (double**)malloc2d(goiP->numVarMD, (int) goiP->dataPoints, sizeof(double));
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
		if(err = MDGetWaveDimensions(p->XFlag_xx, &numdimensions, dimensionSizes)) 
			goto done;
		if(goiP->numVarMD > 1 && numdimensions == 1){
			if (err = MDGetDPDataFromNumericWave(p->XFlag_xx, *goiP->allIndependentVariable))// Get copy.
				goto done;
			
			for(ii = 1 ; ii < goiP->numVarMD ; ii += 1)
				if (err = MDGetDPDataFromNumericWave(p->XFlagWaveH[ii-1], *(goiP->allIndependentVariable + ii)))// Get copy.
					goto done;

		} else
			if (err = MDGetDPDataFromNumericWave(p->XFlag_xx, *goiP->allIndependentVariable))// Get copy.
				goto done;
		
		jj = 0;
		for(ii = 0 ; ii < goiP->dataPoints ; ii += 1){
			if(!(goiP->mask[ii] == 0 || IsNaN64(goiP->mask + ii))){
				for(kk = 0 ; kk < goiP->numVarMD ; kk += 1)
					goiP->independentVariable[kk][jj] = goiP->allIndependentVariable[kk][ii];
				jj+=1;
			}
		}
	} else {	
		//by now the program should've aborted if you haven't specified xwaves and you
		//are fitting multivariate data
		jj = 0;
		for(ii = 0 ; ii < WavePoints(p->dataWave.waveH) ; ii += 1){
			goiP->allIndependentVariable[0][ii] = goiP->ystart + ((double)ii) * goiP->ydelta;
			if(!(goiP->mask[ii] == 0 || IsNaN64(goiP->mask + ii))){
				goiP->independentVariable[0][jj] = goiP->ystart + ((double)ii)*goiP->ydelta;
				jj+=1;
			}
		}
	}
	for(ii = 0; ii < goiP->numVarMD ; ii += 1){
		if(goiP->fullExtentOfData[ii] == NULL){
			err = NOWAV;
			goto done;
		}
		if (err = MDStoreDPDataInNumericWave(goiP->fullExtentOfData[ii], *(goiP->allIndependentVariable + ii)))//put the full extent of x vals into the utilitywave
			goto done;
	}
	
	//store the x array in an x wave used to calculate the theoretical model, but only if you are fitting all at once functions
	//creating these waves is necessary for all-at-once fits.
	if(goiP->isAAO){
		for(ii = 0; ii < goiP->numVarMD ; ii += 1){
			if (err = MDStoreDPDataInNumericWave(goiP->xcalc[ii], *(goiP->independentVariable + ii)))//put the full extent of x vals into the utilitywave
				goto done;
		}
	}
	
	/*	
	 setup output
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
			dimensionSizes[0] = (CountInt)p->LFlag_destLen;
			outPutType = WaveType(p->dataWave.waveH);	//to save memory use datawaves type
			if(err = MDMakeWave(&goiP->OUT_data, datawavename, goiP->cDF, dimensionSizes, outPutType, 1))
				goto done;
			minpos = findmin(*goiP->allIndependentVariable, goiP->dataPoints);
			maxpos = findmax(*goiP->allIndependentVariable, goiP->dataPoints);
			temp1 = (*goiP->allIndependentVariable)[maxpos] - (*goiP->allIndependentVariable)[minpos];
			temp1 /= floor(p->LFlag_destLen) - 1;
			
			if(err = MDSetWaveScaling(goiP->OUT_data, ROWS, &temp1, &(goiP->allIndependentVariable[0][minpos])))
				goto done;
			if(err = MDMakeWave(&gcf_tempWaveHndl_OUTx, xwavename, tempWavesDFH, dimensionSizes, NT_FP64, 1))
				goto done;
			if(err = HoldWave(gcf_tempWaveHndl_OUTx))
					goto done;
			goiP->tempWaveHndl_OUTx = gcf_tempWaveHndl_OUTx;
			
			goiP->OUT_x[0] = goiP->tempWaveHndl_OUTx;
			
			for(ii = 0 ; ii < (CountInt)p->LFlag_destLen ; ii += 1){
				indices[0] = ii;
				value[0] = goiP->allIndependentVariable[0][minpos]+((double)ii) * temp1;
				if(err = MDSetNumericWavePointValue(goiP->OUT_x[0], indices, value))
					goto done;
			}
		} else {		//destination wave and data is MD
			dimensionSizes[0] = goiP->dataPoints;
			if(err = MDMakeWave(&goiP->OUT_data, datawavename, goiP->cDF, dimensionSizes, NT_FP64, 1))
				goto done;
			for(ii = 0 ; ii < MAX_MDFIT_SIZE ; ii += 1) 
				goiP->OUT_x[ii] = goiP->fullExtentOfData[ii];
		}	
	} else {		//destination wave specified
		for(ii = 0; ii < MAX_MDFIT_SIZE ; ii += 1)
			goiP->OUT_x[ii] = goiP->fullExtentOfData[ii];
		goiP->OUT_data = p->DFlag_outputwave;	
	}
	
	if(p->RFlagEncountered){
		if(p->RFlag_resid != NULL){
			goiP->OUT_res = p->RFlag_resid;
		} else {
			dimensionSizes[1] = 0;
			dimensionSizes[0] = WavePoints(p->dataWave.waveH);
			if(err = MDMakeWave(&goiP->OUT_res, reswavename, goiP->cDF, dimensionSizes, NT_FP64, 1))
				goto done;
		}
	}
	
	//if you are doing updates, then append fit output to the topgraph (if data is shown there)
	if(!goiP->noupdate){
		if(goiP->numVarMD == 1 && RunningInMainThread()){
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
							WaveName(p->XFlag_xx, &xwavename[0]);
						} else {
							WaveName(goiP->OUT_x[0], &xwavename[0]);
						}
						if(!strcmp(xwavename, "_free_")){
							strcat(cmd," vs ");
							strcat(cmd, &xwavename[0]);
						}
					}
                    PauseUpdate(&updateStatus);
					if(err = XOPSilentCommand(&cmd[0])){
						ResumeUpdate(&updateStatus);
                        goto done;
                    }
                    ResumeUpdate(&updateStatus);
				}
			}
		}
		//user also specified that they wanted a userupdate function.
		if(goiP->useIgorUpdateFunction){
			dimensionSizes[1] = 0;
			dimensionSizes[0] = goiP->popsize * goiP->numvarparams;
			if(err = MDMakeWave(&gcf_W_costmap, "TEMP_costmap", tempWavesDFH, dimensionSizes, NT_FP64, 1))
				goto done;
			

			dimensionSizes[2] = 0;
			dimensionSizes[1] = goiP->popsize * goiP->numvarparams;
			dimensionSizes[0] =  goiP->numvarparams;
			if(err = MDMakeWave(&gcf_M_population, "TEMP_population", tempWavesDFH, dimensionSizes, NT_FP64, 1))
				goto done;
			
			if(err = HoldWave(gcf_W_costmap))
				goto done;
			if(err = HoldWave(gcf_M_population))
				goto done;
			goiP->W_costmap = gcf_W_costmap;
			goiP->M_population = gcf_M_population;
			
		}
	}
		

done:
	if(holdstr != NULL)
		free(holdstr);
	
	return err;
}


/*
 freeAllocMem frees all the temporary arrays in the GenCurveFitInternals structure
 */
static void
freeAllocMem(GenCurveFitInternalsPtr goiP) {
	int err = 0, ii = 0;

	if (goiP->varParams)
		free(goiP->varParams);
	if (goiP->dataTemp)
		free(goiP->dataTemp);
	if (goiP->holdvector)
		free(goiP->holdvector);
	if (goiP->coefs)
		free(goiP->coefs);
	if (goiP->limits != NULL)
		free(goiP->limits);
	if (goiP->mask != NULL)
		free(goiP->mask);
	if (goiP->independentVariable != NULL)
		free(goiP->independentVariable);
	if (goiP->allIndependentVariable != NULL)
		free(goiP->allIndependentVariable);
	if (goiP->dataObs != NULL)
		free(goiP->dataObs);
	if (goiP->dataObsFull)
		free(goiP->dataObsFull);
	if (goiP->dataSig != NULL)
		free(goiP->dataSig);
	if (goiP->covarianceMatrix != NULL)
		free(goiP->covarianceMatrix);
	if (goiP->dumpRecord.memory)
		free(goiP->dumpRecord.memory);

    err = ReleaseWave(&goiP->GenCurveFitCoefs);
    err = ReleaseWave(&goiP->dataCalc);
    err = ReleaseWave(&goiP->yobs);
    err = ReleaseWave(&goiP->sobs);
    err = ReleaseWave(&goiP->M_population);
    err = ReleaseWave(&goiP->W_costmap);
    for (ii = 0; ii < goiP->numVarMD; ii += 1)
        err = ReleaseWave(&goiP->xcalc[ii]);
    err = ReleaseWave(&goiP->tempWaveHndl_OUTx);
    for (ii = 0; ii < goiP->numVarMD; ii += 1)
        err = ReleaseWave(&goiP->fullExtentOfData[ii]);

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
	CountInt numfitpoints = WavePoints(output);
	
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
			tempX = (double*)malloc(ndims * numfitpoints * sizeof(double));
			if(tempX == NULL){
				err = NOMEM;
				goto done;
			}
			tempY = (double*)malloc(numfitpoints * sizeof(double));
			if(tempY == NULL){
				err = NOMEM;
				goto done;
			}
			for(ii = 0 ; ii < ndims ; ii += 1){
				if(xx[ii] == NULL){
					err = UNSPECIFIED_ERROR;
					goto done;
				}
				if(err = MDGetDPDataFromNumericWave(xx[ii], tempX+(numfitpoints*ii)))
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
	CountInt ii;
	double *temp1 = NULL, *temp2 = NULL;
	double val = 0, val2 = 0, val3 = 0;
	//check if the wave references are NULL
	if(wav1 == NULL)
		return NON_EXISTENT_WAVE;
	if(wav2 == NULL)
		return NON_EXISTENT_WAVE;
	if(WavePoints(wav1) != WavePoints(wav2))
		return WAVE_LENGTH_MISMATCH;
	
	// we have to create temporary arrays to hold the wave data
	if((temp1 = (double*)malloc(sizeof(double) * WavePoints(wav1))) ==  NULL ){
		err = NOMEM;
		goto done;
	}
	if((temp2 = (double*)malloc(sizeof(double) * WavePoints(wav1))) ==  NULL ){
		err = NOMEM;
		goto done;
	}
	// get the data from the waves and put it into the temporary arrays
	if(err = MDGetDPDataFromNumericWave(wav1, temp1))
		goto done;
	if(err = MDGetDPDataFromNumericWave(wav2, temp2))
		goto done;
	// do the subtraction
	for(ii = 0 ; ii < WavePoints(wav1) ; ii += 1){
		val = *(temp1+ii);
		val2 = *(temp2+ii);
		val3 = val - val2;
		*(temp1 + ii) = val - val2;
	}
	// store the subtraction in wav1
	if(err = MDStoreDPDataInNumericWave(wav1, temp1))
		goto done;
	
	WaveHandleModified(wav1);
done:
	if(temp1 != NULL)
		free(temp1);
	if(temp2 != NULL)
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
	CountInt ii;
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
	
	if(err = MDGetDPDataFromNumericWave(wav1, temp1))
		goto done;
	// do the multiplication
	for(ii = 0 ; ii < WavePoints(wav1) ; ii += 1){
		val = *(temp1 + ii) * scalar;
		*(temp1 + ii) = val;
	}
	// store the multiplication in wav1
	if(err = MDStoreDPDataInNumericWave(wav1, temp1))
		goto done;
	
	WaveHandleModified(wav1);
	
done:
	if(temp1 != NULL)
		free(temp1);
	
	return err;
}

/*
 findmin finds the minimum value in a pointer array
 returns minimum position.
 */
static CountInt
findmin(double* sort, CountInt sortsize){
	CountInt ii = 0 , minpos = 0;
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
static CountInt
findmax(double* sort, CountInt sortsize){
	CountInt ii = 0 , maxpos = 0;
	double maxval = sort[0];
	for(ii = 0 ; ii < sortsize ; ii += 1){
		if(sort[ii] > maxval){
			maxval = sort[ii];
			maxpos = ii;
		}
	}
	return maxpos;
}


/*
 dumpRecordToWave puts the dumped population array into a wave, if the /DUMP flag was specified.
 */
int dumpRecordToWave(GenCurveFitInternalsPtr goiP,	MemoryStruct *dumpRecord){
	int err = 0;
	
	waveHndl dump;
	CountInt dimensionSizes[MAX_DIMENSIONS + 1]; // Array of dimension sizes 
	
	memset(dimensionSizes, 0, sizeof(dimensionSizes));
	dimensionSizes[0] = goiP->totalnumparams;
	dimensionSizes[1] = dumpRecord->size / (sizeof(double) * goiP->totalnumparams);
	
	if(err = MDMakeWave(&dump,"M_gencurvefitpopdump", goiP->cDF, dimensionSizes, NT_FP64, 1))
		return err;
	
	if(err = MDStoreDPDataInNumericWave(dump, (double*)dumpRecord->memory))
		return err;
	
	return err;
}


/*
 identicalWaves tests whether two waveHandles refer to the same wave.
 returns 0 if no errors
 returns errorcode otherwise
 if(wav1 == wav2) then isSame=1 
 */
int
identicalWaves(waveHndl wav1, waveHndl wav2){
	if(wav1 == wav2)
        return 1;
    else
        return 0;
}

/*
 isWaveDisplayed tests if wav is displayed in the top graph
 returns 0 if no error
 returns errorcode otherwise
 if(wav is displayed in top graph) then isDisplayed=1
 if it's not displayed isDisplayed = 0
 */
static int
isWaveDisplayed(waveHndl wav, int *isDisplayed){
	char cmd[MAXCMDLEN+1];
	char varName[MAX_OBJ_NAME+1];
	char gwaveName[MAX_WAVE_NAME+1];
	double re=-1,imag=-1;
	int err=0;
    int updateStatus = 0;
	*isDisplayed = 0;
	   
	//we'll say that it's not displayed anyway if you're in a threaded environment
	if(!RunningInMainThread())
		return 0;
	
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
    
    PauseUpdate(&updateStatus);
	if(err = XOPSilentCommand(&cmd[0])){
        ResumeUpdate(&updateStatus);
        return err;
    }
    
	if(FetchNumVar(varName, &re, &imag)==-1)
		return EXPECTED_VARNAME;
	if(re != -1) 
		*isDisplayed = 1;
	strcpy(cmd,"Killvariables/z TEMPGenCurveFit_GLOBALVAR");

	if(err = XOPSilentCommand(&cmd[0])){
        ResumeUpdate(&updateStatus);
		return err;
    }
    ResumeUpdate(&updateStatus);

	return 0;
}

static int
getRange (WaveRange range, CountInt *startPoint, CountInt *endPoint){
	int direction;
	int err = 0;
	*startPoint = (CountInt)range.startCoord;
	*endPoint = (CountInt)range.endCoord;
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

static waveStats getWaveStats(double *sort, CountInt length, int moment){
	CountInt ii=0;
	double minval = *sort, maxval = *sort;
	CountInt minpos=0,maxpos=0;
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



