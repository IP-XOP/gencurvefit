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
#include "GenCurvefit.h"


/*
 this function checks the input from all the parameters IGOR gives it.
 returns 0 if error
 returns error code otherwise
 */
int checkInput(GenCurveFitRuntimeParamsPtr p, GenCurveFitInternalsPtr goiP){
	long numdimensions;
	long dimsize[MAX_DIMENSIONS+1];
	long  coefsdimsize = 0;
	int  err =0;
	int badParameterNumber;
	int sameWave;
	long numzeros = 0;
	char *comparison = NULL;
	double* dpL = NULL;
	double* dpC = NULL;
	char *holdstr = NULL;
	int requiredParameterTypes[MAX_MDFIT_SIZE+2];
	int METH=0, ii=0;
	
	//get the current datafolder
	if(err = GetCurrentDataFolder(&goiP->cDF)){
		goto done;
	}
	
	//start analysing the fitfunction
	goiP->isAAO = 0;
	goiP->numVarMD = -1;
	
	if (p->fitfunEncountered) {
		//couldn't get function information
		if(err = GetFunctionInfo(p->fitfun, &goiP->fi))
			goto done;
		
		// function is not proper fitfunc
		if(goiP->fi.totalNumParameters < 1 || goiP->fi.totalNumParameters > MAX_MDFIT_SIZE){
			err = INVALID_FIT_FUNC;
			goto done;
		}
		
		switch(goiP->fi.totalNumParameters){
			case 1:
				if(p->STRCFlag_sp){
					goiP->sp = p->STRCFlag_sp;
					
					if(p->STRCFlag_sp == NULL){
						err = NULL_STRUCTURE;
						goto done;
					}					
					if (p->STRCFlag_sp->version != kfitfuncStructVersion) { 
						err = INCOMPATIBLE_STRUCT_VERSION; 
						goto done; 
					} 
				} else {
					err = NEED_STRC;
					goto done;
				}
				goiP->numVarMD = p->STRCFlag_sp->numVarMD;
				goiP->sp->numVarMD = goiP->numVarMD;
				
				requiredParameterTypes[0] = FV_STRUCT_TYPE | FV_REF_TYPE;
				if(err = CheckFunctionForm(&goiP->fi,goiP->fi.totalNumParameters,requiredParameterTypes,&badParameterNumber,-1)){
					err = INVALID_FIT_FUNC;
					goto done;
				}
					goiP->isAAO = 2;
				break;
			default:		//either normal fitfunc or all at once
							//first argument is always a wave containing the coefficients.
				requiredParameterTypes[0] = WAVE_TYPE;
				
				for(ii = 1 ; ii <goiP->fi.totalNumParameters ; ii+=1){
					requiredParameterTypes[ii] = NT_FP64;
				}
					goiP->numVarMD = goiP->fi.totalNumParameters-1;
				goiP->isAAO = 0;
				err = CheckFunctionForm(&goiP->fi,goiP->fi.totalNumParameters,requiredParameterTypes,&badParameterNumber,-1);
				
				if(err){ //it may be all-at-once
					for(ii = 0 ; ii <goiP->fi.totalNumParameters ; ii+=1)
						requiredParameterTypes[ii] = WAVE_TYPE;
					goiP->numVarMD = goiP->fi.totalNumParameters-2;
					goiP->isAAO = 1;
					if(err = CheckFunctionForm(&goiP->fi,goiP->fi.totalNumParameters,requiredParameterTypes,&badParameterNumber,-1)){
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
	
	if(p->LFlagEncountered){
		if(IsNaN64(&p->LFlag_destLen) || IsINF64(&p->LFlag_destLen) || p->LFlag_destLen<1){
			err = BAD_FLAG_NUM;
			goto done;
		}
	} else {
		p->LFlag_destLen = 200;
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
		if(err = MDGetWaveDimensions(p->dataWave.waveH, &numdimensions,dimsize)) 
			goto done;
		
		goiP->dataPoints = dimsize[0];
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
		if(err = MDGetWaveDimensions(p->coefs, &numdimensions,dimsize))
			goto done;
		goiP->totalnumparams = dimsize[0];
		
		//if the coefswave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		//if there are no coefficients, then you can't do anything.
		if(dimsize[0] == 0){
			err =COEF_HAS_NO_POINTS;
			goto done;
		}
		//all the parameters have to be usable numbers.
		if(err = checkNanInf(p->coefs)){
			err = COEF_HAS_NANINF; 
			goto done;
		}
		coefsdimsize = dimsize[0];
	}
	
	//get the ywave scaling just in case the xwave isn't specified.
	if (err = MDGetWaveScaling(p->dataWave.waveH, ROWS, &goiP->ydelta, &goiP->ystart)) // Get X scaling
		goto done;
	
	//check if the independent variables are specified.
	if (p->XFlagEncountered) {
		if(p->XFlag_xx == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		if(goiP->numVarMD > 1){
			if(err = MDGetWaveDimensions(p->XFlag_xx, &numdimensions,dimsize)) 
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
							if(err = MDGetWaveDimensions(p->XFlagWaveH[ii], &numdimensions,dimsize)) 
								goto done;
							//if the xwave isn't 1D then we can't fit it.
							if(numdimensions>1){
								err = INPUT_WAVES_NOT_1D;
								goto done;
							}
							//if it isn't the same size as the datawave abort.
							if(dimsize[0] != goiP->dataPoints){
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
					if(dimsize[0] != goiP->dataPoints){
						err = err = WAVES_NOT_SAME_LENGTH;
						goto done;
					}
					if(dimsize[1] != goiP->numVarMD){
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
			if(!((WaveType(p->XFlag_xx) == NT_FP64) || (WaveType(p->XFlag_xx) == NT_FP32))){
				err = REQUIRES_SP_OR_DP_WAVE;
				goto done;
			}
			//check how many points are in the wave
			if(err = MDGetWaveDimensions(p->XFlag_xx, &numdimensions,dimsize)) 
				goto done;
			//if the ywave isn't 1D then we can't fit it.
			if(numdimensions>1){
				err = INPUT_WAVES_NOT_1D;
				goto done;
			}
			//if it isn't the same size as the ywave abort.
			if(dimsize[0] != goiP->dataPoints){
				err = WAVES_NOT_SAME_LENGTH;
				goto done;
			}
			//if the xwave contains NaN or INF, then stop.
			if(err = checkNanInf(p->XFlag_xx)){
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
			goiP->weighttype = p->IFlag_weighttype;
	}  else {
		goiP->weighttype = 0;
	}
	
	//was there a weight (sd) wave specified?
	if (p->WFlag_weighttype) {
		if(p->WFlag_weighttype == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		//if the weightwave isn't Double precision
		if(!((WaveType(p->WFlag_weighttype) == NT_FP64) || (WaveType(p->WFlag_weighttype) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}
		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->WFlag_weighttype, &numdimensions,dimsize)) 
			goto done;
		//if the weight wave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		//if it isn't the same size as the ywave abort.
		if(dimsize[0] != goiP->dataPoints){
			err = WAVES_NOT_SAME_LENGTH;
			goto done;
		}
		//check the weightwave for NaN/INF
		if(err = checkNanInf(p->WFlag_weighttype)){
			err = INPUT_WAVES_CONTAINS_NANINF;
			goto done;
		}
		//check if there are any zeros in the weightwave
		//this is because you will get a divide by zero error if chi2 uses the value
		//as a denominator
		if(err = checkZeros(p->WFlag_weighttype, &numzeros))
			goto done;
		if(goiP->weighttype == 1 && numzeros>0){
			err = STANDARD_DEV_IS_ZERO;
			goto done;
		}
	} else {
		goiP->weighttype = -1;
	}
	
	//was there a mask wave specified?  Set to 0 or NaN to mask points from a fit.
	if (p->MFlagEncountered) {
		if(p->MFlag_maskwave == NIL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		//if the maskwave isn't Double precision
		if(!((WaveType(p->MFlag_maskwave) == NT_FP64) || (WaveType(p->MFlag_maskwave) == NT_FP32))){
			err = REQUIRES_SP_OR_DP_WAVE;
			goto done;
		}
		//check how many points are in the wave
		if(err = MDGetWaveDimensions(p->MFlag_maskwave, &numdimensions,dimsize)) 
			goto done;
		//if the weight wave isn't 1D then we can't fit it.
		if(numdimensions>1){
			err = INPUT_WAVES_NOT_1D;
			goto done;
		}
		//if it isn't the same size as the ywave abort.
		if(dimsize[0] != goiP->dataPoints){
			err = WAVES_NOT_SAME_LENGTH;
			goto done;
		}
	}
	
	//check if we are producing residuals.
	if (p->RFlagEncountered) {
		if(p->RFlag_resid != NULL){
			if(!((WaveType(p->RFlag_resid) == NT_FP64) || (WaveType(p->RFlag_resid) == NT_FP32))){
				err = REQUIRES_SP_OR_DP_WAVE;
				goto done;
			}
			//check how many points are in the wave
			if(err = MDGetWaveDimensions(p->RFlag_resid, &numdimensions,dimsize)) 
				goto done;
			//if the ywave isn't 1D then we can't fit it.
			if(numdimensions>1){
				err = INPUT_WAVES_NOT_1D;
				goto done;
			}
			//if it isn't the same size as the ywave abort.
			if(dimsize[0] != goiP->dataPoints){
				err = WAVES_NOT_SAME_LENGTH;
				goto done;
			}
		}
	}
	
	//these parameters control how the differential evolution operates.
	if (p->KFlagEncountered) {
		//can't have duff genetic optimisation input
		if(IsNaN64(&p->KFlag_iterations) || IsINF64(&p->KFlag_iterations) || p->KFlag_iterations<1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
		if(IsNaN64(&p->KFlag_popsize) || IsINF64(&p->KFlag_popsize) || p->KFlag_popsize <1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
		if(IsNaN64(&p->KFlag_km) || IsINF64(&p->KFlag_km) || p->KFlag_km<=0 || p->KFlag_km > 1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
		if(IsNaN64(&p->KFlag_recomb) || IsINF64(&p->KFlag_recomb) || p->KFlag_recomb<=0 || p->KFlag_recomb>1){
			err = GenCurveFit_PARS_INCORRECT;
			goto done;
		}
	} else {
		p->KFlag_iterations = 100.;
		p->KFlag_popsize = 20.;
		p->KFlag_km = 0.7;
		p->KFlag_recomb = 0.5;
	}
	
	//the cost function for minimisation is now specified.  
	if(p->METHFlagEncountered){
		METH = (int)p->METHFlag_method;
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
	
	// the user specifies a costfunction
	if (p->MINFFlagEncountered) {
		//couldn't get function information
		if(err = GetFunctionInfo(p->MINFFlag_minfun, &goiP->minf))
			goto done;
			
			// function is not proper fitfunc
		if(goiP->minf.totalNumParameters != 4){
			err = INVALID_COST_FUNCTION;
			goto done;
		}
			
			//fit function always has to return a number, not complex or text, even if its all-at-once.
		if(goiP->minf.returnType != NT_FP64){
			err = COSTFUNC_DOESNT_RETURN_NUMBER;
			goto done;
		}
		requiredParameterTypes[0] = WAVE_TYPE;
		requiredParameterTypes[1] = WAVE_TYPE;
		requiredParameterTypes[2] = WAVE_TYPE;
		requiredParameterTypes[3] = WAVE_TYPE;
		if(err = CheckFunctionForm(&goiP->minf, 4, requiredParameterTypes, &badParameterNumber, NT_FP64)){
			err = INVALID_COST_FUNCTION;
			goto done;
		}
		goiP->METH = 2;
	}
	
	//a holdstring is used to work out which parameters are being fitted.  i.e. "0001000111".
	//0=fit, 1=hold
	if (!p->holdstring) {
		//please specify holdstring
		err = HOLDSTRING_NOT_SPECIFIED;
		goto done;
	} else {
		//if we have a holdstring we want to use it.
		if(p->holdstring !=NULL)
		{
			long len;
			long ii;
			int val;
			
			len = GetHandleSize(p->holdstring);
			//if specified the holdstring should be the same length as the coefficient wave
			if(len != coefsdimsize){
				err = HOLDSTRING_NOT_RIGHT_SIZE;
				goto done;
			}
			
			holdstr = (char*)malloc((len+1)*sizeof(char));
			if(holdstr == NULL){
				err = NOMEM;
				goto done;
			}
			comparison = (char*)malloc(2*sizeof(char));
			if(comparison == NULL){
				err = NOMEM;
				goto done;
			}
			
			*(comparison+1) = '\0';
			//get the holdstring from the operation handle
			if(err = GetCStringFromHandle(p->holdstring,holdstr,len)){
				goto done;
			}
			goiP->numvarparams = 0;
			//we have to check that the holdstring is simply 0 or 1's.
			for(ii = 0L; ii<len ; ii++){
				//if its not a digit its not a correct holdstring
				if(!isdigit(*(holdstr+ii))){
					err = HOLDSTRING_INVALID;
					goto done;
				}
				*comparison = *(holdstr+ii);
				val = atoi(comparison);
				// you may have an invalid holdstring
				if(!(val == 0 || val==1)){
					err = HOLDSTRING_INVALID;
					goto done;
				} else {
					// if the holdstring = '0' then you want to vary that parameter
					if(val == 0)
						goiP->numvarparams +=1;
				}
			}
			//if all the parameters are being held then go no further.
			if(goiP->numvarparams == 0){
				err = ALL_COEFS_BEING_HELD;
				goto done;
			}
		} else {
			//please specify holdstring
			err = HOLDSTRING_NOT_SPECIFIED;
			goto done;
		}
	}
	
	//the fractional tolerance for stopping the fit.
	if (p->TOLFlagEncountered) {
		if(IsNaN64(&p->TOLFlag_tol) || IsINF64(&p->TOLFlag_tol) || p->TOLFlag_tol<0){
			err = STOPPING_TOL_INVALID;
			goto done;
		}
	} else {
		p->TOLFlag_tol = 0.0001;
	}
	
	// Main parameters.
	if (p->limitswaveEncountered) {
		long len;
		long ii;
		int val;
		long numBytesL,numBytesC;
		
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
		if(err = MDGetWaveDimensions(p->limitswave, &numdimensions,dimsize)) 
			goto done;
		//we need an upper and lower boundary.
		if(numdimensions != 2){
			err = LIMITS_WRONG_DIMS;
			goto done;
		}
		//if it isn't the same size as input coefs abort.
		if(dimsize[0] != coefsdimsize){
			err = LIMITS_WRONG_DIMS;
			goto done;
		}
		//if there any Nan/Inf's there will be a problem
		if(err = checkNanInf(p->limitswave)){
			err = INPUT_WAVES_CONTAINS_NANINF;
			goto done;
		}
		
		//now we need to check that the coefficients lie between the limits and
		//that the limits are sane
		
		numBytesL = WavePoints(p->limitswave) * sizeof(double); // Bytes needed for copy
		numBytesC = WavePoints(p->coefs) * sizeof(double);
		
		dpL = (double*)malloc(numBytesL);
		dpC = (double*)malloc(numBytesC);
		if (dpC==NULL || dpL == NULL){
			err = NOMEM;
			goto done;
		}
		if (err = MDGetDPDataFromNumericWave(p->coefs, dpC)) { // Get copy.
			goto done;
		}
		if (err = MDGetDPDataFromNumericWave(p->limitswave, dpL)) // Get copy.
			goto done;
		//get the holdstring
		len = GetHandleSize(p->holdstring);
		
		*(comparison+1) = '\0';
		if(err = GetCStringFromHandle(p->holdstring,holdstr,len)){
			goto done;
		}
		for(ii = 0L; ii<len ; ii++){
			*comparison = *(holdstr+ii);
			val = atoi(comparison);
			if(val==0){
				//lowerlim must be < upperlim and param must be in between.
				if(p->OPTFlagEncountered && (((long)p->OPTFlag_opt) & (long)pow(2,0))){
					if(*(dpC+ii)<*(dpL + ii) || *(dpC+ii) > *(dpL + ii+len) || *(dpL + ii+len) <*(dpL + ii)){
						err = LIMITS_INVALID;
						goto done;
					}
				} else {//it doesn't need to be inbetween because we generate our own values
					if(*(dpL + ii+len) <*(dpL + ii)){
						err = LIMITS_INVALID;
						goto done;
					}
				}
			}
		}
	} else {
		//if you don't have a limits wave you can't do anything
		err = NON_EXISTENT_WAVE;
		goto done;
	}
	
	//DFlag will be the output
	if(p->DFlagEncountered){
		if(p->DFlag_outputwave == NULL){
			err = NON_EXISTENT_WAVE;
			goto done;
		}
		// the output wave has to be the same size as the input fit wave
		if(WavePoints(p->DFlag_outputwave) != WavePoints(p->dataWave.waveH)){
			err = OUTPUT_WAVE_WRONG_SIZE;
			goto done;
		}
		// the input waves and output wave shouldn't be the same
		if(err = identicalWaves(p->DFlag_outputwave,p->coefs,&sameWave)) goto done;
		if(sameWave == 1){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(err = identicalWaves(p->DFlag_outputwave,p->dataWave.waveH,&sameWave))	goto done;
		if(sameWave == 1){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(err = identicalWaves(p->DFlag_outputwave,p->XFlag_xx,&sameWave))	goto done;
		if(sameWave == 1){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(err = identicalWaves(p->DFlag_outputwave,p->limitswave,&sameWave))	goto done;
		if(sameWave == 1){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(err = identicalWaves(p->DFlag_outputwave,p->MFlag_maskwave,&sameWave)) goto done;
		if(sameWave == 1){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
		if(err = identicalWaves(p->DFlag_outputwave,p->WFlag_weighttype,&sameWave)) goto done;
		if(sameWave == 1){
			err = OUTPUT_WAVE_OVERWRITING_INPUT;
			goto done;
		}
	}
	// free all memory allocated during this function
done:
		if(comparison != NULL)
			free(comparison);
	if(holdstr != NULL)
		free(holdstr);
	if(dpC !=NULL)
		free(dpC);
	if(dpL !=NULL)
		free(dpL);
	return err;
}


