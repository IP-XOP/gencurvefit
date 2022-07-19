/*
 *  main.c
 *  GenCurvefit
 */
#include "XOPStandardHeaders.h"
#include "GenCurveFitXOP.h"

static int
RegisterGenCurveFit(void){
    const char* cmdTemplate = "GenCurveFit/MC/HOLD=wave:holdwav/POL/STGY=number:stgy/MINF=name:minfun/DITH={number:dith1, number:dith2}/UPDT=name:igorUpdateFunc/DUMP/STRC=structure:{sp,fitfuncStruct}/opt=number:opt /MAT[=number:mat] /q[=number:quiet] /n[=number:noupdate] /SEED=number:seed /L=number:destLen /R[=Wave:resid] /meth=number:method /X={Wave:xx[,Wave[49]]}  /D=wave:outputwave /W=wave:weighttype /I=[number:iflag] /M=wave:maskwave /k={number:iterations, number:popsize, number:km, number:recomb}/TOL=number:tol/POP=wave:initial_popwave name:fitfun, waveRange:dataWave, wave:coefs, string:holdstring, wave:limitswave";
    char const * runtimeNumVarList = "V_Chisq;V_fitIters;V_npnts;V_nterms;V_nheld;V_logBayes";
    char const * runtimeStrVarList = "";

	// NOTE: If you change this template, you must change the GenCurveFitRuntimeParams structure as well.
	return RegisterOperation(cmdTemplate, runtimeNumVarList, runtimeStrVarList, sizeof(GenCurveFitRuntimeParams), (void*)ExecuteGenCurveFit, kOperationIsThreadSafe);
}

static int
RegisterOperations(void)		// Register any operations with Igor.
{
	int result;
	
	// Register XOP1 operation.
	if (result = RegisterGenCurveFit())
		return result;
	
	// There are no more operations added by this XOP.
	
	return 0;
}


/*	XOPEntry()
This is the entry point from the host application to the XOP for all
messages after the INIT message.
*/
static void
XOPEntry(void)
{	
	XOPIORecResult result = 0;
		
	switch (GetXOPMessage()) {
		// We don't need to handle any messages for this XOP.
	}
	SetXOPResult(result);
}

/*	main(ioRecHandle)

This is the initial entry point at which the host application calls XOP.
The message sent by the host must be INIT.

main does any necessary initialization and then sets the XOPEntry field of the
ioRecHandle to the address to be called for future messages.
*/
HOST_IMPORT int
XOPMain(IORecHandle ioRecHandle){
	int result;
	
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.

	if (result = RegisterOperations()) {
		SetXOPResult(result);
		return EXIT_FAILURE;
	}
    
    if (igorVersion < 900){
        SetXOPResult(IGOR_OBSOLETE);
        return EXIT_FAILURE;
    } else {
        SetXOPResult(0L);
        return EXIT_SUCCESS;
    }
}
