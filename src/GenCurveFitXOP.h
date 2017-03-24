// SVN date:    $Date: 2010-07-27 19:25:00 +1000 (Tue, 27 Jul 2010) $
// SVN author:  $Author: andyfaff $
// SVN rev.:    $Revision: 1372 $
// SVN URL:     $HeadURL: svn://svn.igorexchange.com/packages/gencurvefit/trunk/src/GenCurveFit.h $
// SVN ID:      $Id: GenCurveFit.h 1372 2010-07-27 09:25:00Z andyfaff $

/*
	
GenCurvefit.c -- An XOP for curvefitting via Differential Evolution.
@copyright: Andrew Nelson and the Australian Nuclear Science and Technology Organisation 2007.

*/
#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include <time.h>
#include <stdlib.h>

#ifdef WINIGOR
#define snprintf sprintf_s
#endif

//maximum dimension of fit
#define MAX_MDFIT_SIZE 50
#define kfitfuncStructVersion 1000 

// Custom error codes
#define REQUIRES_IGOR_610 1 + FIRST_XOP_ERR
#define NON_EXISTENT_WAVE 2 + FIRST_XOP_ERR
#define REQUIRES_SP_OR_DP_WAVE 3 + FIRST_XOP_ERR
#define WAVES_NOT_SAME_LENGTH 4 + FIRST_XOP_ERR
#define INPUT_WAVES_NOT_1D 5 + FIRST_XOP_ERR
#define INPUT_WAVES_NO_POINTS 6 + FIRST_XOP_ERR
#define INPUT_WAVES_CONTAINS_NANINF 7 + FIRST_XOP_ERR
#define COEF_HAS_NO_POINTS 8 + FIRST_XOP_ERR
#define COEF_HAS_NANINF 9 + FIRST_XOP_ERR
#define SCALING_OF_YWAVE 10 + FIRST_XOP_ERR
#define GenCurveFit_PARS_INCORRECT 11 + FIRST_XOP_ERR
#define HOLDSTRING_NOT_SPECIFIED 12 + FIRST_XOP_ERR
#define HOLDSTRING_INVALID 13 + FIRST_XOP_ERR
#define STOPPING_TOL_INVALID 14 + FIRST_XOP_ERR
#define INVALID_FIT_FUNC 15 + FIRST_XOP_ERR
#define FITFUNC_DOESNT_RETURN_NUMBER 16 + FIRST_XOP_ERR
#define SPARSE_INDEPENDENT_VARIABLE 17 + FIRST_XOP_ERR
#define INVALID_AFITFUNC_INPUT 18 + FIRST_XOP_ERR
#define FITFUNC_NOT_SPECIFIED 19 + FIRST_XOP_ERR
#define LIMITS_WRONG_DIMS 20 + FIRST_XOP_ERR
#define LIMITS_INVALID 21 + FIRST_XOP_ERR
#define FITFUNC_RETURNED_NANINF 22 + FIRST_XOP_ERR
#define UNSPECIFIED_ERROR 23 + FIRST_XOP_ERR
#define HOLDSTRING_NOT_RIGHT_SIZE 24 + FIRST_XOP_ERR
#define ALL_COEFS_BEING_HELD 25 + FIRST_XOP_ERR
#define FIT_ABORTED 26 + FIRST_XOP_ERR
#define OUTPUT_WAVE_WRONG_SIZE 27 + FIRST_XOP_ERR
#define OUTPUT_WAVE_OVERWRITING_INPUT 28 + FIRST_XOP_ERR
#define STANDARD_DEV_IS_ZERO 29 + FIRST_XOP_ERR
#define USER_CHANGED_FITWAVE 30 + FIRST_XOP_ERR
#define INCORRECT_COST_FUNCTION 31 + FIRST_XOP_ERR
#define SUBRANGE_SPECIFIED_ASX 32 + FIRST_XOP_ERR
#define NULL_STRUCTURE 33 + FIRST_XOP_ERR
#define NEED_STRC 34 + FIRST_XOP_ERR
#define INVALID_COST_FUNCTION 35 + FIRST_XOP_ERR
#define COSTFUNC_DOESNT_RETURN_NUMBER 36 + FIRST_XOP_ERR
#define COSTFUNC_WAVES_CHANGED 37 + FIRST_XOP_ERR
#define INVALID_UPDATE_FUNCTION 38 + FIRST_XOP_ERR
#define UPDTFUNC_DOESNT_RETURN_NUMBER 39 + FIRST_XOP_ERR
#define INCORRECT_INITIAL_POPULATION 40 + FIRST_XOP_ERR
#define INITIAL_POPULATION_DP 41 + FIRST_XOP_ERR

/*
Structure fitfuncStruct   
Wave w
wave y
wave x[50]
 
int16 numVarMD
wave ffsWaves[50]
wave ffsTextWaves[10]
variable ffsvar[5]
string ffsstr[5]
nvar ffsnvars[5]
svar ffssvars[5]
funcref allatoncefitfunction ffsfuncrefs[10]
uint32 ffsversion    // Structure version. 
EndStructure 
*/
#pragma pack(2)
struct fitfuncStruct { 
 waveHndl w;
 waveHndl yy;
 waveHndl xx[MAX_MDFIT_SIZE];
 short numVarMD;

 waveHndl otherNumWaves[50];
 waveHndl otherTextWaves[10];
 
 double var[5];
 Handle str[5];
 
 NVARRec nvars[5];
 SVARRec svars[5];
 void* funcRef[10];
 
 UInt32 version;     // Structure version.

}; 
typedef struct fitfuncStruct fitfuncStruct; 
typedef struct fitfuncStruct* fitfuncStructPtr; 

struct GenCurveFitRuntimeParams {
	// Flag parameters.
    
    // Parameters for /MC flag group.
	int MCFlagEncountered;
	// There are no fields for this group because it has no parameters.
    
	// Parameters for /HOLD flag group.
	int HOLDFlagEncountered;
	waveHndl holdwav;
	int HOLDFlagParamsSet[1];
    
	// Parameters for /POL flag group.
	int POLFlagEncountered;
	// There are no fields for this group because it has no parameters.
    
	// Parameters for /STGY flag group.
	int STGYFlagEncountered;
	double stgy;
	int STGYFlagParamsSet[1];
    
	// Parameters for /MINF flag group.
	int MINFFlagEncountered;
	char minfun[MAX_OBJ_NAME + 1];
	int MINFFlagParamsSet[1];
    
	// Parameters for /DITH flag group.
	int DITHFlagEncountered;
	double dith1;
	double dith2;
	int DITHFlagParamsSet[2];
    
	// Parameters for /UPDT flag group.
	int UPDTFlagEncountered;
	char UPDTFlag_igorUpdateFunc[MAX_OBJ_NAME + 1];
	int UPDTFlagParamsSet[1];
	
	// Parameters for /DUMP flag group.
	int DUMPFlagEncountered;
	// There are no fields for this group because it has no parameters.
	
	// Parameters for /STRC flag group.
	int STRCFlagEncountered;
	fitfuncStruct* STRCFlag_sp;
	int STRCFlagParamsSet[1];

	// Parameters for /OPT flag group.
	int OPTFlagEncountered;
	double OPTFlag_opt;
	int OPTFlagParamsSet[1];

	// Parameters for /MAT flag group.
	int MATFlagEncountered;
	double MATFlag_mat;
	int MATFlagParamsSet[1];
	// There are no fields for this group because it has no parameters.

	// Parameters for /Q flag group.
	int QFlagEncountered;
	double QFlag_quiet;
	int QFlagParamsSet[1];

	// Parameters for /N flag group.
	int NFlagEncountered;
	double NFlag_noupdate;
	int NFlagParamsSet[1];
	
	// There are no fields for this group because it has no parameters.

	// Parameters for /SEED flag group.
	int SEEDFlagEncountered;
	double SEEDFlag_seed;
	int SEEDFlagParamsSet[1];

	// Parameters for /L flag group.
	int LFlagEncountered;
	double LFlag_destLen;
	int LFlagParamsSet[1];

	// Parameters for /R flag group.
	int RFlagEncountered;
	waveHndl RFlag_resid;					// Optional parameter.
	int RFlagParamsSet[1];

	// Parameters for /METH flag group.
	int METHFlagEncountered;
	double METHFlag_method;
	int METHFlagParamsSet[1];

	// Parameters for /X flag group.
	int XFlagEncountered;
	waveHndl XFlag_xx;
	waveHndl XFlagWaveH[49];				// Optional parameter.
	int XFlagParamsSet[50];

	// Parameters for /D flag group.
	int DFlagEncountered;
	waveHndl DFlag_outputwave;
	int DFlagParamsSet[1];

	// Parameters for /W flag group.
	int WFlagEncountered;
	waveHndl WFlag_weighttype;
	int WFlagParamsSet[1];

	// Parameters for /I flag group.
	int IFlagEncountered;
	double IFlag_weighttype;
	int IFlagParamsSet[1];

	// Parameters for /M flag group.
	int MFlagEncountered;
	waveHndl MFlag_maskwave;
	int MFlagParamsSet[1];

	// Parameters for /K flag group.
	int KFlagEncountered;
	double KFlag_iterations;
	double KFlag_popsize;
	double KFlag_km;
	double KFlag_recomb;
	int KFlagParamsSet[4];

	// Parameters for /TOL flag group.
	int TOLFlagEncountered;
	double TOLFlag_tol;
	int TOLFlagParamsSet[1];

    // Parameters for /HOLD flag group.
    int POPFlagEncountered;
    waveHndl initial_popwave;
    int POPFlagParamsSet[1];
    
	// Main parameters.

	// Parameters for simple main group #0.
	int fitfunEncountered;
	char fitfun[MAX_OBJ_NAME+1];
	int fitfunParamsSet[1];

	// Parameters for simple main group #1.
	int dataWaveEncountered;
	WaveRange dataWave;
	int dataWaveParamsSet[1];

	// Parameters for simple main group #2.
	int coefsEncountered;
	waveHndl coefs;
	int coefsParamsSet[1];

	// Parameters for simple main group #3.
	int holdstringEncountered;
	Handle holdstring;
	int holdstringParamsSet[1];

	// Parameters for simple main group #4.
	int limitswaveEncountered;
	waveHndl limitswave;
	int limitswaveParamsSet[1];

	// These are postamble fields that Igor sets.
	int calledFromFunction;					// 1 if called from a user function, 0 otherwise.
	int calledFromMacro;					// 1 if called from a macro, 0 otherwise.
	UserFunctionThreadInfoPtr tp;	

};
typedef struct GenCurveFitRuntimeParams GenCurveFitRuntimeParams;
typedef struct GenCurveFitRuntimeParams* GenCurveFitRuntimeParamsPtr;
#pragma pack()

struct MemoryStruct {
	char *memory;
	size_t size;
};
typedef struct MemoryStruct MemoryStruct;

//this structure contains all the internal memory arrays necessary for the fit to proceed.
struct GenCurveFitInternals{
	//how many parameters are being varied
	int numvarparams;
	//total number of parameters
	CountInt totalnumparams;
	//which parameters are varying;
	unsigned int *varParams;
	//how many dimensions you are trying to fit
	int numVarMD;
	
	//the holdvector for which parameters you want to fix
	unsigned int *holdvector;

	double recomb;
	double k_m;
	double tolerance;
	long popsize;
	long iterations;
	
	//lowest cost
	double cost;
	
	//strategy for minimisation
	int STGY;

	//cost function indicator for minimisation
	int METH;
	FunctionInfo minf;

	//number of fititerations done
	long V_numfititers;
	//do you want to perform dynamic updates?
	int noupdate;
	
	//logarithm of the Bayes Posterior probability
	double V_logBayes;
	
	//the coefficients being fitted
	double *coefs;
	//a full copy of the y data being fitted
	double *dataObsFull;
	//the number of dataPoints;
	CountInt dataPoints;
	//a temporary array the same length as the unmasked dataset;
	double *dataTemp;
	//a copy of the y data being fitted, without the masked points.
	double *dataObs;
	//which weight type?
	int weighttype;
	//the corresponding errors for each of those points.
	double *dataSig;
	//corresponding independent variable points for each of the masked data points.
	double **independentVariable;
	//corresponding independent variable points for each of the unmasked data points.
	double **allIndependentVariable;
	//a copy of the limits that are being used.
	double **limits;
	//an array specifying which y points are not included.
	double *mask;
	//if a range is specified for the ywave needs start and endpoints
	CountInt startPoint;
	CountInt endPoint;
	//how many points are being fitted.
	CountInt unMaskedPoints;
	//the ywave scaling, in case no x-wave is specified.
	double ystart,ydelta;
	//the function being fitted.
	FunctionInfo fi;
	char *functionname;

	//is it a structure fit?
	fitfuncStruct* sp;

	//is it all at once (normal fit func=0,AAO=1,Structfit=2)?
	int isAAO;
	//an estimated covariance Matrix
	double **covarianceMatrix;
	//a handle for the covariance matrix
	waveHndl M_covariance;
	//a handle for the error wave
	waveHndl W_sigma;
	
	//an array for dumping the population at each iteration
	MemoryStruct dumpRecord;

	//a user specified update function
	int useIgorUpdateFunction;
	FunctionInfo igorUpdateFunction;
	//a wave containing the current population
	waveHndl M_population;
	//a wave containing the current chi2 values for each of the population
	waveHndl W_costmap;
	
	//Did you want to dump the population?
	int dump;
	
	//a number that says if the fit is converging, only of use for tracking progress
	//not used anywhere else
	double convergenceNumber;

	//the current datafolder needs to be stored, so we have a place to put temporary waves.
	//dataCalc, xcalc,GenCurveFitCoefs are temporary waves created so that we can call a function.
	
	//a handle to the current data folder
	DataFolderHandle cDF;
	//a handle to the calculated data. This is filled in calcModel.  This will have the same size as yobs and sobs.
	waveHndl dataCalc;
	//Wavehandle to the original data
	waveHndl originalYobs;
	//wavehandles to the unmasked y_obs points, with errors, that you are fitting to.  This may be just a copy of the original data
	waveHndl yobs;	//GenCurveFit_yobs
	waveHndl sobs;
	
	// a handle to the xwave used to calculate the model (excluding masked points).
	// Only made for an all at once FF.
	waveHndl xcalc[MAX_MDFIT_SIZE]; 
	// the coefficients used to calculate the model
	waveHndl GenCurveFitCoefs;
	// the full range of xpoints being used (including masked points).
	waveHndl fullExtentOfData[MAX_MDFIT_SIZE];
	// a temporary wave handle.  This will be made if /D is specified and if the independent data is univariate
	waveHndl tempWaveHndl_OUTx;
	
	//Wave Handles for the output, i.e. the fit waves.
	//the output y wave
	waveHndl OUT_data;
	//the output xwave
	waveHndl OUT_x[MAX_MDFIT_SIZE];	//these aren't actual waves, but handles to prexisting waves.
	//the output residual wave
	waveHndl OUT_res;
	//the output coefficients
	waveHndl OUT_coefs;
};
typedef struct GenCurveFitInternals GenCurveFitInternals;
typedef struct GenCurveFitInternals* GenCurveFitInternalsPtr;

/*
	A structure to hold statistics of a wave
*/
struct waveStats {
	double V_avg;
	double V_stdev;
	CountInt V_maxloc;
	CountInt V_minloc;
};
typedef struct waveStats waveStats;
typedef struct waveStats* waveStatsPtr;

#pragma pack(2)
struct fitFunc { // Used to pass parameters to the function.
	waveHndl waveH; // For the first function parameter.
	double x[MAX_MDFIT_SIZE];
};
typedef struct fitFunc fitFunc;
typedef struct fitFunc* fitFuncPtr;

struct allFitFunc { // Used to pass parameters to the function.
	waveHndl waveC; // For the coefficients.
	waveHndl waveY;	// for filling up by the function
	waveHndl waveX[MAX_MDFIT_SIZE];	// supplies independent values for function
};

typedef struct allFitFunc allFitFunc;
typedef struct allFitFunc* allFitFuncPtr;

struct costFunc { // Used to pass parameters to the function.
	waveHndl coefs;
	waveHndl yobs;
	waveHndl ycalc;
	waveHndl sobs;
};
typedef struct costFunc costFunc;
typedef struct costFunc* costFuncPtr;

struct updtFunc { // Used to pass parameters to the function.
	waveHndl currentbestfit;
	waveHndl population;
	waveHndl costmap;
	double updatetime;
};
typedef struct updtFunc updtFunc;
typedef struct updtFunc* updtFuncPtr;

#pragma pack()

/*
	Functions contained in Gencurvefit.c
*/
int lgencurvefit_fitfunction(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long datapoints, unsigned int numDataDims);
int ExecuteGenCurveFit(GenCurveFitRuntimeParamsPtr p);
int checkInput(GenCurveFitRuntimeParamsPtr, GenCurveFitInternalsPtr);
int checkNanInf(waveHndl);
int checkNanInfArray(double *array, CountInt datapoints);
int checkZeros(waveHndl , long* );
static void freeAllocMem(GenCurveFitInternalsPtr goiP);
static int calcModelXY(FunctionInfo*, waveHndl , waveHndl , waveHndl[MAX_MDFIT_SIZE] , int ,int ,fitfuncStruct*);
static CountInt findmin(double* , CountInt );
static CountInt findmax(double* , CountInt );
static int init_GenCurveFitInternals(GenCurveFitRuntimeParamsPtr, GenCurveFitInternalsPtr);
int identicalWaves(waveHndl , waveHndl);
static int subtractTwoWaves(waveHndl, waveHndl   );
static int scalarMultiply(waveHndl wav1, double scalar);
static int isWaveDisplayed(waveHndl, int *);
static int getRange (WaveRange ,CountInt *,CountInt *);
static double roundDouble(double);
static waveStats getWaveStats(double*,CountInt,int);
int dumpRecordToWave(GenCurveFitInternalsPtr goiP,	MemoryStruct *dumpRecord);

/* Prototypes */
HOST_IMPORT int XOPMain(IORecHandle ioRecHandle);