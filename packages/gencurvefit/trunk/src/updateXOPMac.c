// SVN date:    $Date$
// SVN author:  $Author$
// SVN rev.:    $Revision$
// SVN URL:     $HeadURL$
// SVN ID:      $Id$

#include "XOPStandardHeaders.h"			// Include ANSI headers, Mac headers, IgorXOP.h, XOP.h and XOPSupport.h

//Carbon headers required for NIB/HiView windows
#include <Carbon/Carbon.h>

//constants for NIB HiViews.
//these fields/signatures should be the same as those given to the 
//UI items in Interface Builder
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
	double convergenceNumber;
};
typedef struct displayData displayData;

//a global variable to hold the data to display
//this is necessary because we can't send it to the event handler directly.
static displayData theDisplayData;

void DisplayWindowXOP1Message(WindowPtr theWindow,int numcoefs, const double* coefs, double chi2, const char* fitfunc,long fititers, double convergenceNumber)
{	
	OSStatus err;
	HIViewRef theHIView;	//drawing window
	static const HIViewID myHIViewID = { kMyHIViewSignature, kMyHIViewFieldID };

	//stick all the data to be drawn in a structure
	memset(&theDisplayData, 0, sizeof(theDisplayData));
	
	theDisplayData.theWindow = theWindow;
	theDisplayData.numcoefs = numcoefs;
	theDisplayData.coefs = coefs;
	theDisplayData.chi2 = chi2;
	theDisplayData.fitfunc = fitfunc;
	theDisplayData.fititers = fititers;
	theDisplayData.convergenceNumber = convergenceNumber;

	//get a HiView for the window we created
	HIViewFindByID (HIViewGetRoot(theWindow), myHIViewID, &theHIView);

	//so we can then force the HIView to re-draw
	err = HIViewSetNeedsDisplay (theHIView, true);
	
}

//an event handler for redraw events on the window
OSStatus MyDrawEventHandler (EventHandlerCallRef myHandler, EventRef event, void *userData){
    OSStatus status = noErr;
    CGContextRef myContext;
    HIRect      bounds;
	
	//the data to be displayed
	extern displayData theDisplayData;
	float w, h;
	int err2, posx, posy, fsize = 12;
	char message[MSG_LEN + 1];
	//spacings to layout the text
	int vertspace = 20, space = 105, ii, jj;	

	//get the graphics context for drawing on.
    status = GetEventParameter (event,
								kEventParamCGContextRef,
								typeCGContextRef,
								NULL,
								sizeof (CGContextRef),
								NULL,
								&myContext);
	
    require_noerr(status, CantGetGraphicsContext);
	
	//the bounds of the rectangle for drawing on
    HIViewGetBounds ((HIViewRef) userData, &bounds);
    require_noerr(status, CantGetBoundingRectangle);
	
	//height and width of the HIView window
	w = bounds.size.width;
	h = bounds.size.height;
	
	//the HIView coordinate system is reversed for y,
	//compared to the normal Quartz2D graphics system
	//this code reverses it.
	CGContextTranslateCTM (myContext, 0, bounds.size.height);
	CGContextScaleCTM (myContext, 1.0, -1.0);
	
	//set up the text to draw: font, colour, etc.
    CGContextSelectFont (myContext, 
						 "Helvetica-Bold",
						 fsize,
						 kCGEncodingMacRoman);
    CGContextSetCharacterSpacing (myContext, 2); // 4
    CGContextSetTextDrawingMode (myContext, kCGTextFill); // 5
    CGContextSetRGBFillColor (myContext, 0, 0, 0, 1); // 6
	CGContextSetRGBStrokeColor (myContext, 0, 0, 0, 1); // 7
	
	//we want to start drawing at the top of the box.
	posx = bounds.origin.x;
	posy = bounds.size.height - fsize;
	
	//now print out the data we want to display in the box.
	
	if(theDisplayData.fitfunc){
		//can't use strncpy/strncat on strings that don't exist.
		strncpy(message, "Fitting to: ", MSG_LEN);
		strncat(message, theDisplayData.fitfunc, MSG_LEN - strlen(message));
		CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
	}
	//print out the number of iterations and the Chi2 value
	posx = bounds.origin.x;
	posy -= vertspace;
	
	err2 = snprintf(message, MSG_LEN, "Iterations: %li", theDisplayData.fititers);
	CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
	
	posy -= vertspace;
	err2 = snprintf(message, MSG_LEN, "Chi2: %-6.4g", theDisplayData.chi2);
	CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
	
	posx += 2*space;
	err2 = snprintf(message, MSG_LEN, "Convergence (fit stops when >1): %-6.4g", theDisplayData.convergenceNumber);
	CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
	
	posx = bounds.origin.x;
	posy -= vertspace;
	strncpy(message, "Coefficients", MSG_LEN);
	CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
	
	//now print out the first 25 coefficients.  In order to get grid like behaviour there's
	//a lot of moving about.
	for(jj=0 ; theDisplayData.coefs && (jj<(int)ceil((double)(theDisplayData.numcoefs)/5) && jj<6) ; jj+=1){
		posy -= vertspace;
		for(ii=0 ; (ii<5 && ii < ((theDisplayData.numcoefs)-(jj*5))) ; ii+=1){
			posx = bounds.origin.x + ii*space;
			err2 = snprintf(message, MSG_LEN, "%-6.4g", *(theDisplayData.coefs + ii + (jj*5)));
			CGContextShowTextAtPoint (myContext, posx, posy, message, strlen(message) ); 
		}
	}
	
CantGetGraphicsContext:	
CantGetBoundingRectangle:
    return status;
}

//This event handles all the button presses for the abort button in the window.
//All it does is set a global variable which signals the main fitting routine to stop.
OSStatus MyButtonHandler (EventHandlerCallRef myHandler, EventRef event, void *userData){
	OSStatus err = 0;
	extern int Abort_the_fit;
	Abort_the_fit = 1;
	return err;
}

//create the XOP window from the NIB file.
//note that the Info.plist file of the XOP has a key called "Bundle Identifier".
//This key must be the same as the string in the CFBundleGetBundleWithIdentifier function call,
//otherwise the NIB doesn't load.
WindowPtr CreateXOPWindow(void){
	WindowPtr theWindow = NULL;
	IBNibRef nibRef;
    OSStatus err;
	HIViewRef theHIView;
	
	//HiViewRef for the "Abort" push button
	ControlRef theHIViewButton;

	static const EventTypeSpec myHIViewSpec[] = {kEventClassControl, kEventControlDraw };
	static const HIViewID myHIViewID = { kMyHIViewSignature, kMyHIViewFieldID };
	static const EventTypeSpec myHIViewButtonSpec[] = {kEventClassControl, kEventControlHit };
	static const ControlID myHIViewButtonID = { kMyHIViewButtonSignature, kMyHIViewButtonID };
	
	// Create a Nib reference passing the name of the nib file (without the .nib extension)
    // CreateNibReference only searches into the application bundle.
	CFBundleRef xop_bundle = CFBundleGetBundleWithIdentifier(CFSTR("gencurvefit"));
	//the window in the NIB file must be the same as the 2nd argument ("main" in this case).
	err= CreateNibReferenceWithCFBundle ( xop_bundle, CFSTR("main"), &nibRef);
    require_noerr( err, done );	

	//create a window from the NIB reference.  
	//the window in the NIB file must be the same as the 2nd argument ("main" in this case).
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

