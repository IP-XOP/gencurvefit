// SVN date:    $Date$
// SVN author:  $Author$
// SVN rev.:    $Revision$
// SVN URL:     $HeadURL$
// SVN ID:      $Id$

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h
#include <Carbon/Carbon.h>

//constants for NIB HiViews.
#define kMyHIViewSignature 'pViw'// the HiView window signature
#define kMyHIViewFieldID 123	//the HiView window ID
#define kMyHIViewButtonSignature 'pBut'  //the HiView Button signature
#define kMyHIViewButtonID 124   //the HiView Button ID


#define MSG_LEN 254

void DestroyXOPWindow(XOP_WINDOW_REF w);

struct displayData{
	WindowPtr theWindow;
	int numcoefs;
	const double *coefs;
	double chi2;
	const char *fitfunc;
	long fititers;
};
typedef struct displayData displayData;

//a pointer to hold the data to display
static displayData *theDisplayData = NULL;

//the HiView for the text container in the progress window.
static HIViewRef theHIView = NULL;

void DisplayWindowXOP1Message(WindowPtr theWindow,int numcoefs, const double* coefs, double chi2, const char* fitfunc,long fititers)
{	
	OSStatus err;
	theDisplayData = (displayData*)malloc(sizeof(displayData));
	if(!theDisplayData)
		return;
	
	theDisplayData->theWindow = theWindow;
	theDisplayData->numcoefs = numcoefs;
	theDisplayData->coefs = coefs;
	theDisplayData->chi2 = chi2;
	theDisplayData->fitfunc = fitfunc;
	theDisplayData->fititers = fititers;
	
	err = HIViewSetNeedsDisplay (theHIView, true);	
}

OSStatus MyDrawEventHandler (EventHandlerCallRef myHandler, EventRef event, void *userData){
    OSStatus status = noErr;
    CGContextRef myContext;
    HIRect      bounds;
	
	extern displayData *theDisplayData;
	float w, h;
	
	int err2, posx, posy, fsize = 12;
	char message[MSG_LEN + 1];
	int vertspace = 20, space = 105, ii, jj;	

    status = GetEventParameter (event,
								kEventParamCGContextRef,
								typeCGContextRef,
								NULL,
								sizeof (CGContextRef),
								NULL,
								&myContext);
	
    require_noerr(status, CantGetGraphicsContext);
	
    HIViewGetBounds ((HIViewRef) userData, &bounds);
    require_noerr(status, CantGetBoundingRectangle);
	
	w = bounds.size.width;
	h = bounds.size.height;
	
	CGContextTranslateCTM (myContext, 0, bounds.size.height);
	CGContextScaleCTM (myContext, 1.0, -1.0);
	
    CGContextSelectFont (myContext, 
						 "Helvetica-Bold",
						 fsize,
						 kCGEncodingMacRoman);
    CGContextSetCharacterSpacing (myContext, 2); // 4
    CGContextSetTextDrawingMode (myContext, kCGTextFill); // 5
    CGContextSetRGBFillColor (myContext, 0, 0, 0, 1); // 6
	CGContextSetRGBStrokeColor (myContext, 0, 0, 0, 1); // 7
	
	if(theDisplayData){
		posx = bounds.origin.x;
		posy = bounds.size.height - fsize;
		
		strncpy(message, "Fitting to: ", MSG_LEN);
		strncat(message, theDisplayData->fitfunc, MSG_LEN - strlen(message));
		CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
		
		//print out the number of iterations and the Chi2 value
		posx = bounds.origin.x;
		posy -= vertspace;
		
		err2 = snprintf(message, MSG_LEN, "Iterations: %li", theDisplayData->fititers);
		CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 

		posy -= vertspace;
		err2 = snprintf(message, MSG_LEN, "Chi2: %-6.4g", theDisplayData->chi2);
		CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 

		posy -= vertspace;
		strncpy(message, "Coefficients", MSG_LEN);
		CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 

		//now print out the first 25 coefficients.  In order to get grid like behaviour there's
		//a lot of moving about.
		for(jj=0 ; (jj<(int)ceil((double)(theDisplayData->numcoefs)/5) && jj<6) ; jj+=1){
			posy -= vertspace;
			for(ii=0 ; (ii<5 && ii < ((theDisplayData->numcoefs)-(jj*5))) ; ii+=1){
				posx = bounds.origin.x + ii*space;
				err2 = snprintf(message, MSG_LEN, "%-6.4g", *(theDisplayData->coefs + ii + (jj*5)));
				CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
			}
		}
		
	}

CantGetGraphicsContext:	
CantGetBoundingRectangle:
    return status;
}

OSStatus MyButtonHandler (EventHandlerCallRef myHandler, EventRef event, void *userData){
	OSStatus err = 0;
	extern int Abort_the_fit;
	Abort_the_fit = 1;
	return err;
}


WindowPtr CreateXOPWindow(void){
	WindowPtr theWindow = NULL;
	extern HIViewRef theHIView;	//drawing window
	IBNibRef nibRef;
    OSStatus err;

	//HiViewRef for the "Abort" push button
	ControlRef theHIViewButton;

	static const EventTypeSpec myHIViewSpec[] = {kEventClassControl, kEventControlDraw };
	static const HIViewID myHIViewID = { kMyHIViewSignature, kMyHIViewFieldID };
	static const EventTypeSpec myHIViewButtonSpec[] = {kEventClassControl, kEventControlHit };
	static const ControlID myHIViewButtonID = { kMyHIViewButtonSignature, kMyHIViewButtonID };
		
	
	// Create a Nib reference passing the name of the nib file (without the .nib extension)
    // CreateNibReference only searches into the application bundle.
	CFBundleRef xop_bundle = CFBundleGetBundleWithIdentifier(CFSTR("gencurvefit"));
	err= CreateNibReferenceWithCFBundle ( xop_bundle, CFSTR("main"), &nibRef);
    require_noerr( err, done );	

	err = CreateWindowFromNib(nibRef, CFSTR("main"), &theWindow);
    require_noerr(err, done );	

	// We don't need the nib reference anymore.
    DisposeNibReference(nibRef);
	
	//get a HiView for the window we create
	HIViewFindByID (HIViewGetRoot(theWindow), myHIViewID, &theHIView);
	
    err = InstallEventHandler (GetControlEventTarget (theHIView),
							   NewEventHandlerUPP (MyDrawEventHandler),
							   GetEventTypeCount (myHIViewSpec),
							   myHIViewSpec,
							   (void *) theHIView,
							   NULL);
	require_noerr(err, done );	
	
	// Get the button control
	err = GetControlByID( theWindow, &myHIViewButtonID, &theHIViewButton );
	//err = HIViewFindByID(HIViewGetRoot(theWindow),myHIViewButtonID, &theHIViewButton);
	require_noerr(err, done );								   	

	err = InstallControlEventHandler( theHIViewButton,
									 NewEventHandlerUPP( MyButtonHandler ), 
									 GetEventTypeCount (myHIViewButtonSpec),
									 myHIViewButtonSpec,
									 (void*) theHIViewButton,
									 NULL);

    require_noerr(err, done );

	//display the window, even if there is nothing in it.
	ShowWindow(theWindow);
done:
	return theWindow;
}

//cleans up all the XOP windows.
void DestroyXOPWindow(WindowPtr theWindow) {
	CFRelease(theWindow);
	theWindow = NULL;
}

