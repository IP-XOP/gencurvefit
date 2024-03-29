/*
GenCurveFit.c -- An XOP for curvefitting via Differential Evolution.
@copyright: Andrew Nelson and the Australian Nuclear Science and Technology Organisation 2007.
*/
#include "XOPStandardHeaders.r"

resource 'vers' (1) {						// XOP version info.
	0x01, 0x20, release, 0x00, 0,			// version bytes and country integer.
	"1.",
	"1, � andyfaff@gmail.com all rights reserved."
};

resource 'vers' (2) {						// Igor version info.
	0x08, 0x00, release, 0x00, 0,			// Version bytes and country integer.
	"8.00",
	"(for Igor 8.00 or later)"
};

resource 'STR#' (1100) {					/* custom error messages */
	{
	//[1]
	"GenCurvefit requires Igor Pro 9.0 or later.",
	//[2]
	"Non-Existent Wave.",
	//[3]
	"Requires Double or Single Precision Wave.",
	//[4]
	"Waves not the same length.",
	//[5]
	"Input data must be 1D.",
	//[6]
	"Input data has no points.",
	//[7]
	"Input data contains NaN/INF.",
	//[8]
	"Coefficient Wave has no points.",
	//[9]
	"Coefficient Wave contains NaN/INF.",
	//[10]
	"Scaling of Y wave is NaN/INF.",
	//[11]
	"/K --> genetic optimisation parameters not specified correctly.",
	//[12]
	"Please specify holdstring.",
	//[13]
	"Holdstring should only contain 0 and 1 (fit,don't fit).",
	//[14]
	"Stopping tolerance must be 0<tol and not NaN/INF .",
	//[15]
	"Invalid FitFunction.",
	//[16]
	"Fitfunction doesn't return a number.",
	//[17]
	"Multivariate fitfunction supplied, but wrong independent variable information supplied.",
	//[18]
	"Parameters for fitfunc must be f(Wave,Wave,Wave).",
	//[19]
	"Fit function wasn't specified.",
	//[20]
	"Limits Wave must have dimensions NUMROWS=numpnts(coefs), NUMCOLS=2.",
	//[21]
	"Limits and/or coefficient waves not specified correctly.\n lowerlim<coef<upperlim.",
	//[22]
	"FitFunction returned NaN/INF for one or more values.",
	//[23]
	"Unspecified Internal Gencurvefit Error.",
	//[24]
	"Holdstring should be same length as coefficient wave.",
	//[25]
	"All parameters are being held.",
	//[26]
	"Genetic Optimisation Aborted.",
	//[27]
	"Output wave requires same number of points as ywave (/R/D).",
	//[28]
	"You are trying to overwrite an Input/Output wave with an output wave (/R/D).",
	//[29]
	"you are trying to weight with standard deviations, but one of them is zero - DIVIDE BY ZERO.",
	//[30]
	"User redimensioned/killed Internal Gencurvefit fitwave.",
	//[31]
	"Wrong Index for cost function.",
	//[32]
	"Subrange limited to point specification (you specified xscaling) - remember to sort wave.",
	//[33]
	"Null structure specified in GenCurvefit.",
	//[34]
	"GenCurvefit requires a proper structure in the fitfunction.",
	//[35]
	"The costfunction should be of the form cost(params, y_obs, y_calc, s_obs)",
	//[36]
	"The cost function should return a double precision number",
	//[37]
	"Do not change the length of the waves in the cost function",
	//[38]
	"The update function must be of the form updt(wave, wave, wave, updatetime)",
	//[39]
	"The update function should return a number",
    //[40]
    "If you specify an initial population is should have dimensions [numcoefs][N], where N > 0",
    //[41]
    "Initial population needs to be a double precision wave",
	}
};

resource 'STR#' (1101) {					// Misc strings for XOP.
	{
		"-1",								// This item is no longer supported by the Carbon XOP Toolkit.
		"No Menu Item",						// This item is no longer supported by the Carbon XOP Toolkit.
		"GenCurveFit Help",					// Name of XOP's help file.
	}
};

// No menu item

resource 'XOPI' (1100) {
	XOP_VERSION,							// XOP protocol version.
	DEV_SYS_CODE,							// Development system information.
    XOP_FEATURE_FLAGS,                        // Tells Igor about XOP features
    XOPI_RESERVED,                            // Reserved - must be zero.
    XOP_TOOLKIT_VERSION,                    // XOP Toolkit version.
};

resource 'XOPC' (1100) {
	{
		"GenCurveFit",								// Name of operation.
		XOPOp+UtilOP+compilableOp+threadSafeOp,			// Operation's category.
	}
};


