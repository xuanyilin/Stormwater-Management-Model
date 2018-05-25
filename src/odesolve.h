//-----------------------------------------------------------------------------
//  odesolve.h
//
//  Header file for ODE solver contained in odesolve.c
//
//-----------------------------------------------------------------------------

#ifndef ODESOLVE_H
#define ODESOLVE_H


//-----------------------------------------------------------------------------
//    Local declarations
//-----------------------------------------------------------------------------
typedef struct
{
    int      nmax;      // max. number of equations
    double*  y;         // dependent variable
    double*  yscal;     // scaling factors
    double*  yerr;      // integration errors
    double*  ytemp;     // temporary values of y
    double*  dydx;      // derivatives of y
    double*  ak;        // derivatives at intermediate points
} TOdesolveShared;

// functions that open, close, and use the ODE solver
int  odesolve_open(SWMM_Project *sp, int n);

void odesolve_close(SWMM_Project *sp);

int  odesolve_integrate(SWMM_Project *sp, double ystart[], int n, double x1,
        double x2, double eps, double h1,
        void (*derivs)(SWMM_Project*, double, double*, double*));


#endif //ODESOLVE_H
