

/*
 Functions contained in updateXOP<x>.c
 */
void DrawXOPWindow(IgorWindowRef w);
void DisplayWindowXOP1Message(IgorWindowRef w, int numcoefs, const double* coefs, double chi2,char* fitfunc,long fititers, unsigned int updatetime, double convergenceNumber);
IgorWindowRef CreateXOPWindow(void);
void DestroyXOPWindow(IgorWindowRef w);

#ifdef _WINDOWS_			// [
int CreateXOPWindowClass(void);
#endif	
