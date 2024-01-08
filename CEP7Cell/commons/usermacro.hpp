#ifndef USERMACRO_HPP
#define USERMACRO_HPP

// File for defining the global macro.
// Do not put any includes in here.

// Default current time and future time
#define T0 0.0
#define T1 0.1

// Set this value to use solver instead of analytical
#define SOLVER

// Set this value if you want combine the geometry and scalar into one VTK file
#define VTK_COMBINED

// Define cell model. 
// Uncomment the models that you want to use.
// Comment the others.

#define CRN1998
//#define HHUXLEY1952
//#define ORUDY2011_STATIC
//#define ORUDY2017_DYNAMIC
//#define TN2004ENDO
//#define TN2004EPI
//#define TN2004M
//#define TN2006ENDO
//#define TN2006EPI
//#define TN2006M


#endif
