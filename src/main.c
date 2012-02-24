// SVN date:    $Date$
// SVN author:  $Author$
// SVN rev.:    $Revision$
// SVN URL:     $HeadURL$
// SVN ID:      $Id$

/*
 *  main.c
 *  GenCurvefit
 */
#include "XOPStandardHeaders.h"
#include "GenCurveFitXOP.h"
#include "updateXOP.h"

static int
RegisterGenCurveFit(void){
	char* cmdTemplate;
	char* runtimeNumVarList;
	char* runtimeStrVarList;

	// NOTE: If you change this template, you must change the GenCurveFitRuntimeParams structure as well.
	cmdTemplate = "GenCurveFit/MC/HOLD=wave:holdwav/POL/STGY=number:stgy/TEMP=number:temp/MINF=name:minfun/UPDT=name:igorUpdateFunc/DUMP/STRC=structure:{sp,fitfuncStruct}/opt=number:opt /MAT[=number:mat] /q[=number:quiet] /n[=number:noupdate] /SEED=number:seed /L=number:destLen /R[=Wave:resid] /meth=number:method /X={Wave:xx[,Wave[49]]}  /D=wave:outputwave /W=wave:weighttype /I=[number:weighttype] /M=wave:maskwave /k={number:iterations, number:popsize, number:km, number:recomb}/TOL=number:tol name:fitfun, waveRange:dataWave, wave:coefs, string:holdstring, wave:limitswave";
	runtimeNumVarList = "V_Chisq;V_fitIters;V_npnts;V_nterms;V_nheld;V_logBayes";
	runtimeStrVarList = "";
	return RegisterOperation(cmdTemplate, runtimeNumVarList, runtimeStrVarList, sizeof(GenCurveFitRuntimeParams), (void*)ExecuteGenCurveFit, 0);
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
	long result = 0;
	
	if (WindowMessage())							// Handle all messages related to XOP window.
		return;
	
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
HOST_IMPORT void
main(IORecHandle ioRecHandle){
	int result;
	
	XOPInit(ioRecHandle);							// Do standard XOP initialization.
	
	SetXOPEntry(XOPEntry);							// Set entry point for future calls.
	
	//CreateXOP Window Class for the update window.  This is not needed for Mac usage.
#ifdef _WINDOWS_
	{
		if (result = CreateXOPWindowClass()) {
			SetXOPResult(result);
			return;
		}
	}
#endif
	
	if (result = RegisterOperations()) {
		SetXOPResult(result);
		return;
	}
	
	SetXOPResult(0);
}
