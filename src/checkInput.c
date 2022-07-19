// SVN date:    $Date$
// SVN author:  $Author$
// SVN rev.:    $Revision$
// SVN URL:     $HeadURL$
// SVN ID:      $Id$

/*
 *  checkInput.c
 *  GenCurvefit
 *
 *  Created by andrew on 25/04/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "XOPStandardHeaders.h"
#include "GenCurveFitXOP.h"
#include "gencurvefit.h"

/*
 this function checks the input from all the parameters IGOR gives it.
 returns 0 if error
 returns error code otherwise
 */
int checkInput(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	int numdimensions;
	CountInt indices[MAX_DIMENSIONS + 1];
	int  err =0;
	int badParameterNumber;
	long numzeros = 0;
	char *holdstr = NULL;
	int requiredParameterTypes[MAX_MDFIT_SIZE + 2];
	int METH=0, ii=0, jj;
	double value[2];
		
	/*
	 get the current datafolder
	 */
	if(err = GetCurrentDataFolder(&goiP->cDF))
		goto done;
	
	/*
	 start analysing the fitfunction
	 */
	goiP->isAAO = 0;
	goiP->numVarMD = -1;
	
	if (p->fitfunEncountered) {
        //couldn't get function information
		if(err = GetFunctionInfo(p->fitfun, &goiP->fi))
			goto done;

		if(!goiP->fi.isThreadSafe && !RunningInMainThread()){
			err = ONLY_THREADSAFE;
			goto done;
		}
		
		goiP->functionname = p->fitfun;

		// function is not proper fitfunc
		if(goiP->fi.totalNumParameters < 1 || goiP->fi.totalNumParameters > MAX_MDFIT_SIZE){
			err = INVALID_FIT_FUNC;
			goto done;
		}
		
		switch(goiP->fi.totalNumParameters){
			case 1:
				if(p->sp){
					goiP->sp = p->sp;
					
					if(p->sp == NULL){
						err = NULL_STRUCTURE;
						goto done;
					}					
					if (p->sp->version != kfitfuncStructVersion) {
						err = INCOMPATIBLE_STRUCT_VERSION; 
						goto done; 
					} 
				} else {
					err = NEED_STRC;
					goto done;
				}
				goiP->numVarMD = p->sp->numVarMD;
				goiP->sp->numVarMD = goiP->numVarMD;
				
				requiredParameterTypes[0] = FV_STRUCT_TYPE | FV_REF_TYPE;
				if(err = CheckFunctionForm(&goiP->fi,
                                           goiP->fi.totalNumParameters,
                                           requiredParameterTypes,
                                           &badParameterNumber,
                                           NT_FP64)){
					err = INVALID_FIT_FUNC;
					goto done;
				}
				goiP->isAAO = 2;
				break;
			default:		//either normal fitfunc or all at once
							//first argument is always a wave containing the coefficients.
				requiredParameterTypes[0] = WAVE_TYPE;
				
				for(ii = 1 ; ii <goiP->fi.totalNumParameters ; ii+=1)
					requiredParameterTypes[ii] = NT_FP64;
				
				goiP->numVarMD = goiP->fi.totalNumParameters-1;
				goiP->isAAO = 0;
				err = CheckFunctionForm(&goiP->fi,
                                        goiP->fi.totalNumParameters,
                                        requiredParameterTypes,
                                        &badParameterNumber,
                                        NT_FP64);
				
				if(err){ //it may be all-at-once
					for(ii = 0 ; ii <goiP->fi.totalNumParameters ; ii+=1)
						requiredParameterTypes[ii] = WAVE_TYPE;
					goiP->numVarMD = goiP->fi.totalNumParameters-2;
					goiP->isAAO = 1;
					if(err = CheckFunctionForm(&goiP->fi,
                                               goiP->fi.totalNumParameters,
                                               requiredParameterTypes,
                                               &badParameterNumber,
                                               NT_FP64)){
						err = INVALID_FIT_FUNC;					
						goto done;
					}
				}
					//fit function always has to return a number, not complex or text, even if its all-at-once.
					if(goiP->fi.returnType != NT_FP64){
						err = FITFUNC_DOESNT_RETURN_NUMBER;
						goto done;
					}
					break;
		}
	} else {
        
		//no fit function was specified.
		err = FITFUNC_NOT_SPECIFIED;
		goto done;
	}
	
	if(p->STGYFlagEncountered){
		goiP->STGY = (int) p->stgy;
		if(goiP->STGY < 0 || goiP->STGY > 8)
			goiP->STGY = 0;
	}
	
	if(p->LFlagEncountered){
		if(IsNaN64(&p->destLen) || IsINF64(&p->destLen) || p->destLen<1){
			err = BAD_FLAG_NUM;
			goto done;
		}
	} else {
		p->destLen = 200;
	}
	
	if (p->dataWaveEncountered) {
		//if y wave doesn't exist go no further
		if(p->dataWave.waveH == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		
		//the ywave has to be SP or DP.
		if(!((WaveType(p->dataWave.waveH) == NT_FP64) || (WaveType(p->dataWave.waveH) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}

		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->dataWave.waveH, &numdimensions,indices)) 
			goto done;
		
		goiP->dataPoints = indices[0];
		//if the dataWave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		
		//if there are no points to fit, then you can't do anything.
		if(goiP->dataPoints == 0){
			err = INPUT_WAVES_NO_POINTS;
			goto done;
		}
		////we're not going to do the fit if there are any NaN or INFS in the input.
		//if(err = checkNanInf(p->dataWave.waveH)){
		//	err = INPUT_WAVES_CONTAINS_NANINF;
		//	goto done;
		//} 
	}
	
	if (p->coefsEncountered) {
		//if ceofs wave doesn't exist go no further
		if(p->coefs == NIL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		//if the coefsWave isn't Double precision
		if(!((WaveType(p->coefs) == NT_FP64) || (WaveType(p->coefs) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}
		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->coefs, &numdimensions,indices))
			goto done;
		goiP->totalnumparams = indices[0];
		
		//if the coefswave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		//if there are no coefficients, then you can't do anything.
		if(indices[0] == 0){
			err =COEF_HAS_NO_POINTS;
			goto done;
		}
		//all the parameters have to be usable numbers.
		if(err = checkNanInf(p->coefs)){
			err = COEF_HAS_NANINF; 
			goto done;
		}
				
		//put the coefficients into the temporary space
		goiP->coefs = (double*) malloc (sizeof(double) * WavePoints(p->coefs));
		if(goiP->coefs == NULL){
			err = NOMEM;
			goto done;
		}
		if(err = MDGetDPDataFromNumericWave(p->coefs, goiP->coefs))
			goto done;
	}
	
	//get the ywave scaling just in case the xwave isn't specified.
	if (err = MDGetWaveScaling(p->dataWave.waveH, ROWS, &goiP->ydelta, &goiP->ystart)) // Get X scaling
		goto done;
	
	/*
	 check if the independent variables are specified.
	 */
	if (p->XFlagEncountered) {
		if(p->xx == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		if(goiP->numVarMD > 1){
			if(err = MDGetWaveDimensions(p->xx, &numdimensions,indices))
				goto done;
			switch(numdimensions){
				case 1:
					ii=0;
					while(p->XFlagParamsSet[ii] != 0){
						ii+=1;
					}
						if(ii != goiP->numVarMD){
							err = SPARSE_INDEPENDENT_VARIABLE;
							goto done;
						}
						
						if(p->XFlagParamsSet[goiP->numVarMD-1] == 0){
							err =  SPARSE_INDEPENDENT_VARIABLE;
							goto done;
						}
						for(ii=0 ; ii < goiP->numVarMD-1 ; ii+=1){
							if(p->XFlagWaveH[ii] == NULL){
								err = NON_EXISTENT_WAVE;
								goto done;
							}
							//if the xwave isn't Double precision
							if(!((WaveType(p->XFlagWaveH[ii]) == NT_FP64) || (WaveType(p->XFlagWaveH[ii]) == NT_FP32))){
								err = REQUIRES_SP_OR_DP_WAVE;
								goto done;
							}
							//check how many points are in the wave
							if(err = MDGetWaveDimensions(p->XFlagWaveH[ii], &numdimensions,indices)) 
								goto done;
							//if the xwave isn't 1D then we can't fit it.
							if(numdimensions>1){
								err = INPUT_WAVES_NOT_1D;
								goto done;
							}
							//if it isn't the same size as the datawave abort.
							if(indices[0] != goiP->dataPoints){
								err = WAVES_NOT_SAME_LENGTH;
								goto done;
							}
							//if the xwave contains NaN or INF, then stop.
							if(err = checkNanInf(p->XFlagWaveH[ii])){
								err = INPUT_WAVES_CONTAINS_NANINF; 
								goto done;
							}
						}
						break;
				case 2:
					if(p->XFlagParamsSet[1] != 0){
						err = SPARSE_INDEPENDENT_VARIABLE;
						goto done;
					}
					if(indices[0] != goiP->dataPoints){
						err = WAVES_NOT_SAME_LENGTH;
						goto done;
					}
					if(indices[1] != goiP->numVarMD){
						err = SPARSE_INDEPENDENT_VARIABLE;
						goto done;
					}
					break;
				default:
					err = SPARSE_INDEPENDENT_VARIABLE;
					goto done;
					break;
			}
		} else {	//we are fitting 1D data
			if(goiP->numVarMD != 1){
				err = SPARSE_INDEPENDENT_VARIABLE;
				goto done;
			}
			//if the xwave isn't Double precision
			if(!((WaveType(p->xx) == NT_FP64) || (WaveType(p->xx) == NT_FP32))){
				err = REQUIRES_SP_OR_DP_WAVE;
				goto done;
			}
			//check how many points are in the wave
			if(err = MDGetWaveDimensions(p->xx, &numdimensions,indices))
				goto done;
			//if the ywave isn't 1D then we can't fit it.
			if(numdimensions>1){
				err = INPUT_WAVES_NOT_1D;
				goto done;
			}
			//if it isn't the same size as the ywave abort.
			if(indices[0] != goiP->dataPoints){
				err = WAVES_NOT_SAME_LENGTH;
				goto done;
			}
			//if the xwave contains NaN or INF, then stop.
			if(err = checkNanInf(p->xx)){
				err = INPUT_WAVES_CONTAINS_NANINF; 
				goto done;
			}
		}
	} else {
		//we're just going to work with the y wave scaling
		//if the scaling of the y wave is buggered up
		if(goiP->numVarMD != 1){
			err = SPARSE_INDEPENDENT_VARIABLE;
			goto done;
		}
		if(IsNaN64(&goiP->ystart) || IsNaN64(&goiP->ydelta) || IsINF64(&goiP->ystart) || IsINF64(&goiP->ydelta)){
			err = SCALING_OF_YWAVE;
			goto done;
		}
	}
	// if the weightwave contains standard deviations then type=1 otherwise
	// the wave contains weighting values, goiP->weighttype = 0.
	// if no weight wave is specified then goiP->weighttype = -1.
	if(p->IFlagEncountered){
		if(p->IFlagParamsSet[0])
			goiP->weighttype = (int) p->iflag;
	}  else {
		goiP->weighttype = 0;
	}
	
	/*
	 was there a weight (sd) wave specified?
	 */
	if (p->weighttype) {
		if(p->weighttype == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		//if the weightwave isn't Double precision
		if(!((WaveType(p->weighttype) == NT_FP64) || (WaveType(p->weighttype) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}
		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->weighttype, &numdimensions,indices))
			goto done;
		//if the weight wave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		//if it isn't the same size as the ywave abort.
		if(indices[0] != goiP->dataPoints){
			err = WAVES_NOT_SAME_LENGTH;
			goto done;
		}
		//check the weightwave for NaN/INF
		if(err = checkNanInf(p->weighttype)){
			err = INPUT_WAVES_CONTAINS_NANINF;
			goto done;
		}
		//check if there are any zeros in the weightwave
		//this is because you will get a divide by zero error if chi2 uses the value
		//as a denominator
		if(err = checkZeros(p->weighttype, &numzeros))
			goto done;
		if(goiP->weighttype == 1 && numzeros>0){
			err = STANDARD_DEV_IS_ZERO;
			goto done;
		}
	} else {
		goiP->weighttype = -1;
	}
	
	/*
	 was there a mask wave specified?  Set to 0 or NaN to mask points from a fit.
	 */
	if (p->MFlagEncountered) {
		if(p->maskwave == NIL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		//if the maskwave isn't Double precision
		if(!((WaveType(p->maskwave) == NT_FP64) || (WaveType(p->maskwave) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}
		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->maskwave, &numdimensions,indices))
			goto done;
		//if the weight wave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		//if it isn't the same size as the ywave abort.
		if(indices[0] != goiP->dataPoints){
			err = WAVES_NOT_SAME_LENGTH;
			goto done;
		}
	}
	
	/*
	 check if we are producing residuals.
	 */
	if (p->RFlagEncountered) {
		if(p->resid != NULL){
			if(!((WaveType(p->resid) == NT_FP64) || (WaveType(p->resid) == NT_FP32))){
				err = REQUIRES_SP_OR_DP_WAVE;
				goto done;
			}
			//check how many points are in the wave
			if(err = MDGetWaveDimensions(p->resid, &numdimensions,indices)) 
				goto done;
			//if the ywave isn't 1D then we can't fit it.
			if(numdimensions>1){
				err = INPUT_WAVES_NOT_1D;
				goto done;
			}
			//if it isn't the same size as the ywave abort.
			if(indices[0] != goiP->dataPoints){
				err = WAVES_NOT_SAME_LENGTH;
				goto done;
			}
		}
	}
	
	/*
	 these parameters control how the differential evolution operates.
	 */
	if (p->KFlagEncountered) {
		//can't have duff genetic optimisation input
		if(IsNaN64(&p->iterations) || IsINF64(&p->iterations) || p->iterations<1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
		if(IsNaN64(&p->popsize) || IsINF64(&p->popsize) || p->popsize <1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
		if(IsNaN64(&p->km) || IsINF64(&p->km) || p->km<=0 || p->km > 1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
		if(IsNaN64(&p->recomb) || IsINF64(&p->recomb) || p->recomb<=0 || p->recomb>1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
	} else {
		p->iterations = 100.;
		p->popsize = 20.;
		p->km = 0.7;
		p->recomb = 0.5;
	}
	goiP->recomb = p->recomb;
	goiP->k_m = p->km;
	goiP->iterations = (long) p->iterations;
	goiP->popsize = (long) p->popsize;
	

	/*
	 the cost function for minimisation is now specified.  
	 */
	if(p->METHFlagEncountered){
		METH = (int)p->method;
		switch(METH){
			case 0:		//this will be Chi2
				goiP->METH = METH;
				break;
			case 1:		//robust fitting
				goiP->METH = METH;
				break;
			default:
				err = INCORRECT_COST_FUNCTION;
				goto done;
				break;
		}
	} else {
		//default is least squares
		goiP->METH = 0;
	}
	
	/*
	 the user specifies a costfunction
	 */
//	if (p->MINFFlagEncountered) {
//		//couldn't get function information
//		if(err = GetFunctionInfo(p->minfun, &goiP->minf))
//			goto done;
//
//			// function is not proper costfunc
//		if(goiP->minf.totalNumParameters != 4){
//			err = INVALID_COST_FUNCTION;
//			goto done;
//		}
//
//		if(!RunningInMainThread() && (!goiP->minf.isThreadSafe)){
//			err = ONLY_THREADSAFE;
//			goto done;
//		}
//
//		//fit function always has to return a number, not complex or text, even if its all-at-once.
//		if(goiP->minf.returnType != NT_FP64){
//			err = COSTFUNC_DOESNT_RETURN_NUMBER;
//			goto done;
//		}
//		requiredParameterTypes[0] = WAVE_TYPE;
//		requiredParameterTypes[1] = WAVE_TYPE;
//		requiredParameterTypes[2] = WAVE_TYPE;
//		requiredParameterTypes[3] = WAVE_TYPE;
//		if(err = CheckFunctionForm(&goiP->minf, 4, requiredParameterTypes, &badParameterNumber, NT_FP64)){
//			err = INVALID_COST_FUNCTION;
//			goto done;
//		}
//		goiP->METH = 2;
//	}
	
	/*
	 a holdstring is used to work out which parameters are being fitted.  i.e. "0001000111".
	0=fit, 1=hold
	 the holdvector gets passed to the genetic optimisation function
	 */
	goiP->holdvector = (unsigned int*) malloc(goiP->totalnumparams * sizeof(unsigned int));
	if(goiP->holdvector == NULL){
		err = NOMEM;
		goto done;
	}
	for(ii = 0 ; ii < goiP->totalnumparams ; ii++)
		goiP->holdvector[ii] = 1;
	
	//user can specify a holdwave (leads to smaller command lines)
	if(p->HOLDFlagEncountered && p->holdwav != NULL){
		memset(indices, 0, MAX_DIMENSIONS + 1);
		if(WavePoints(p->holdwav) != goiP->totalnumparams){
			err = HOLDSTRING_NOT_RIGHT_SIZE;
			goto done;
		}
		for(ii = 0 ; ii < goiP->totalnumparams ; ii++){
			indices[0] = ii;
			if(err = MDGetNumericWavePointValue(p->holdwav, indices, value))
				goto done;
			if(value[0] == 0 || IsNaN64(value) || IsINF64(value)){
				goiP->holdvector[ii] = 0;
				goiP->numvarparams += 1;
			}
		}
	} else if (!p->holdstring) {
		//please specify holdstring
		err = HOLDSTRING_NOT_SPECIFIED;
		goto done;
	} else {
		//if we have a holdstring we want to use it.
		if(p->holdstring !=NULL)
		{
			char comparison[2];
			BCInt len;
			long ii;
			int val;
			
			memset(comparison, 0, sizeof(char) * 2);
			
			len = WMGetHandleSize(p->holdstring);
			//if specified the holdstring should be the same length as the coefficient wave
			if(len != goiP->totalnumparams){
				err = HOLDSTRING_NOT_RIGHT_SIZE;
				goto done;
			}
			
			holdstr = (char*)malloc((len + 1) * sizeof(char));
			if(holdstr == NULL){
				err = NOMEM;
				goto done;
			}

			/*
			 get the holdstring from the operation handle
			 */
			if(err = GetCStringFromHandle(p->holdstring, holdstr, (int) len)){
				goto done;
			}
			goiP->numvarparams = 0;
			
			/*
			 we have to check that the holdstring is simply 0 or 1's.
			 */
			for(ii = 0L ; ii < len ; ii++){
				//if its not a digit its not a correct holdstring
				if(!isdigit(*(holdstr + ii))){
					err = HOLDSTRING_INVALID;
					goto done;
				}
				comparison[0] = holdstr[ii];
				val = atoi(comparison);
				// you may have an invalid holdstring
				if(!(val == 0 || val == 1)){
					err = HOLDSTRING_INVALID;
					goto done;
				} else {
					// if the holdstring = '0' then you want to vary that parameter
					if(val == 0){
						goiP->numvarparams +=1;
						goiP->holdvector[ii] = 0;
					}
				}
			}

		} else {
			/*
			 please specify holdstring
			 */
			err = HOLDSTRING_NOT_SPECIFIED;
			goto done;
		}
	}
	
	/*
	 if all the parameters are being held then go no further.
	 */
	if(goiP->numvarparams == 0){
		err = ALL_COEFS_BEING_HELD;
		goto done;
	}
	
	goiP->varParams = malloc(sizeof(unsigned int) * goiP->numvarparams);
	if(!goiP->varParams){
		err = NOMEM;
		goto done;
	}
	
	jj = 0;
	for(ii = 0 ; ii < goiP->totalnumparams ; ii++){
		if(goiP->holdvector[ii] == 0){
			goiP->varParams[jj] = ii;
			jj++;
		}
	}
	
	/*
	 the fractional tolerance for stopping the fit.
	 */
	if (p->TOLFlagEncountered) {
		if(IsNaN64(&p->tol)
		   || IsINF64(&p->tol)
		   || p->tol < 0){
			err = STOPPING_TOL_INVALID;
			goto done;
		}
		goiP->tolerance = p->tol;
	} else
		goiP->tolerance = 0.0001;

	/*
	 the user specifies an update function
	 */
//	if (p->UPDTFlagEncountered) {
//		//couldn't get function information
//		if(err = GetFunctionInfo(p->igorUpdateFunc, &goiP->igorUpdateFunction))
//			goto done;
//
//		if(!RunningInMainThread() && (!goiP->igorUpdateFunction.isThreadSafe)){
//			err = ONLY_THREADSAFE;
//			goto done;
//		}
//
//		// function is not proper fitfunc
//		if(goiP->igorUpdateFunction.totalNumParameters != 4){
//			err = INVALID_UPDATE_FUNCTION;
//			goto done;
//		}
//
//		//update function has to return a number
//		if(goiP->igorUpdateFunction.returnType != NT_FP64){
//			err = UPDTFUNC_DOESNT_RETURN_NUMBER;
//			goto done;
//		}
//		requiredParameterTypes[0] = WAVE_TYPE;
//		requiredParameterTypes[1] = WAVE_TYPE;
//		requiredParameterTypes[2] = WAVE_TYPE;
//		requiredParameterTypes[3] = NT_FP64;
//
//		if(err = CheckFunctionForm(&goiP->igorUpdateFunction, 4, requiredParameterTypes, &badParameterNumber, NT_FP64)){
//			err = INVALID_UPDATE_FUNCTION;
//			goto done;
//		}
//		goiP->useIgorUpdateFunction = 1;
//	}
	
	/*
	 Checkout the limits wave
	 */
	if (p->limitswaveEncountered) {		
		//if the limitswave handle doesn't exist.
		if(p->limitswave == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		//if the ywave isn't Double precision
		if(!((WaveType(p->limitswave) == NT_FP64) || (WaveType(p->limitswave) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}
		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->limitswave, &numdimensions,indices)) 
			goto done;
		
		//we need an upper and lower boundary.
		if(numdimensions != 2){
			err = LIMITS_WRONG_DIMS;
			goto done;
		}
		//if it isn't the same size as input coefs abort.
		if(indices[0] != goiP->totalnumparams){
			err = LIMITS_WRONG_DIMS;
			goto done;
		}
		//if there any Nan/Inf's there will be a problem
		if(err = checkNanInf(p->limitswave)){
			err = INPUT_WAVES_CONTAINS_NANINF;
			goto done;
		}
		
		/*
		 now we need to check that the coefficients lie between the limits and
		 that the limits are sane
		 */
		goiP->limits = (double**) malloc2d(2, (int) goiP->totalnumparams, sizeof(double));
		if(goiP->limits == NULL){
			err = NOMEM;
			goto done;
		}
		if(err = MDGetDPDataFromNumericWave(p->limitswave, *(goiP->limits)))
		   goto done;
		   
	} else {
		/*
		 if you don't have a limits wave you can't do anything
		 */
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	
	//DFlag will be the output
	if(p->DFlagEncountered){
		if(p->outputwave == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		// the output wave has to be the same size as the input fit wave
		if(WavePoints(p->outputwave) != WavePoints(p->dataWave.waveH)){
			err = OUTPUT_WAVE_WRONG_SIZE;
			goto done;
		}
		// the input waves and output wave shouldn't be the same
		if(identicalWaves(p->outputwave, p->coefs)){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(identicalWaves(p->outputwave, p->dataWave.waveH)){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(identicalWaves(p->outputwave, p->xx)){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(identicalWaves(p->outputwave, p->limitswave)){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(identicalWaves(p->outputwave, p->maskwave)){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(identicalWaves(p->outputwave, p->weighttype)){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
	}
	
	if(p->DUMPFlagEncountered)
		goiP->dump = 1;
	
done:
	if(holdstr)
		free(holdstr);

	return err;
}


