

/*
 Functions contained in updateXOP<x>.c
 */
void DrawXOPWindow(XOP_WINDOW_REF w);
void DisplayWindowXOP1Message(XOP_WINDOW_REF w, int numcoefs, const double* coefs, double chi2,char* fitfunc,long fititers, unsigned int updatetime, double convergenceNumber);
XOP_WINDOW_REF CreateXOPWindow(void);
void DestroyXOPWindow(XOP_WINDOW_REF w);

#ifdef _WINDOWS_			// [
int CreateXOPWindowClass(void);
#endif	
