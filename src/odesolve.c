//-----------------------------------------------------------------------------
//   odesolve.c
//
//   Fifth-order Runge-Kutta integration with adaptive step size control
//   based on code from Numerical Recipes in C (Cambridge University
//   Press, 1992).
//
//   Date:     11/15/06
//   Author:   L. Rossman
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include "headers.h"
#include "odesolve.h"


#define MAXSTP 10000
#define REALLY_TINY   1.0e-30
#define SAFETY 0.9
#define PGROW  -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4    // = (5/SAFETY)^(1/PGROW)


// function that integrates over an error-controlled stepsize
int rkqs(SWMM_Project *sp,double* x, int n, double htry, double eps, double* hdid,
        double* hnext, void (*derivs)(SWMM_Project*, double, double*, double*));

// function that performs the Runge-Kutta integration step
void rkck(SWMM_Project *sp, double x, int n, double h,
        void (*derivs)(SWMM_Project*, double, double*, double*));


//-----------------------------------------------------------------------------
//    open the ODE solver to solve system of n equations
//    (return 1 if successful, 0 if not)
//-----------------------------------------------------------------------------
int odesolve_open(SWMM_Project *sp, int n)
{
    TOdesolveShared *odslv = &sp->OdesolveShared;

    odslv->nmax  = 0;
    odslv->y     = (double *) calloc(n, sizeof(double));
    odslv->yscal = (double *) calloc(n, sizeof(double));
    odslv->dydx  = (double *) calloc(n, sizeof(double));
    odslv->yerr  = (double *) calloc(n, sizeof(double));
    odslv->ytemp = (double *) calloc(n, sizeof(double));
    odslv->ak    = (double *) calloc(5*n, sizeof(double));

    if ( !odslv->y || !odslv->yscal || !odslv->dydx ||
            !odslv->yerr || !odslv->ytemp || !odslv->ak )
        return 0;

    odslv->nmax = n;
    return 1;
}


//-----------------------------------------------------------------------------
//    close the ODE solver
//-----------------------------------------------------------------------------
void odesolve_close(SWMM_Project *sp)
{
    TOdesolveShared *odslv = &sp->OdesolveShared;

    if ( odslv->y ) free(odslv->y);
    odslv->y = NULL;

    if ( odslv->yscal ) free(odslv->yscal);
    odslv->yscal = NULL;

    if ( odslv->dydx ) free(odslv->dydx);
    odslv->dydx = NULL;

    if ( odslv->yerr ) free(odslv->yerr);
    odslv->yerr = NULL;

    if ( odslv->ytemp ) free(odslv->ytemp);
    odslv->ytemp = NULL;

    if ( odslv->ak ) free(odslv->ak);
    odslv->ak = NULL;

    odslv->nmax = 0;
}


int odesolve_integrate(SWMM_Project *sp, double ystart[], int n, double x1,
        double x2, double eps, double h1,
        void (*derivs)(SWMM_Project*, double, double*, double*))
//---------------------------------------------------------------
//   Driver function for Runge-Kutta integration with adaptive
//   stepsize control. Integrates starting n values in ystart[]
//   from x1 to x2 with accuracy eps. h1 is the initial stepsize
//   guess and derivs is a user-supplied function that computes
//   derivatives dy/dx of y. On completion, ystart[] contains the
//   new values of y at the end of the integration interval.
//---------------------------------------------------------------
{
    int    i, errcode, nstp;
    double hdid, hnext;
    double x = x1;
    double h = h1;

    TOdesolveShared *odslv = &sp->OdesolveShared;

    if (odslv->nmax < n) return 1;
    for (i=0; i<n; i++) odslv->y[i] = ystart[i];
    for (nstp=1; nstp<=MAXSTP; nstp++)
    {
        derivs(sp,x,odslv->y,odslv->dydx);
        for (i=0; i<n; i++)
            odslv->yscal[i] = fabs(odslv->y[i]) + fabs(odslv->dydx[i]*h) + REALLY_TINY;
        if ((x+h-x2)*(x+h-x1) > 0.0) h = x2 - x;
        errcode = rkqs(sp, &x,n,h,eps,&hdid,&hnext,derivs);
        if (errcode) break;
        if ((x-x2)*(x2-x1) >= 0.0)
        {
            for (i=0; i<n; i++) ystart[i] = odslv->y[i];
            return 0;
        }
        if (fabs(hnext) <= 0.0) return 2;
        h = hnext;
    }
    return 3;
}


int rkqs(SWMM_Project *sp, double* x, int n, double htry, double eps, double* hdid,
         double* hnext, void (*derivs)(SWMM_Project*, double, double*, double*))
//---------------------------------------------------------------
//   Fifth-order Runge-Kutta integration step with monitoring of
//   local truncation error to assure accuracy and adjust stepsize.
//   Inputs are current value of x, trial step size (htry), and
//   accuracy (eps). Outputs are stepsize taken (hdid) and estimated
//   next stepsize (hnext). Also updated are the values of y[].
//---------------------------------------------------------------
{
    int i;
    double err, errmax, h, htemp, xnew, xold = *x;

    TOdesolveShared *odslv = &sp->OdesolveShared;

    // --- set initial stepsize
    h = htry;
    for (;;)
    {
        // --- take a Runge-Kutta-Cash-Karp step
        rkck(sp, xold, n, h, derivs);

        // --- compute scaled maximum error
        errmax = 0.0;
        for (i=0; i<n; i++)
        {
            err = fabs(odslv->yerr[i]/odslv->yscal[i]);
            if (err > errmax) errmax = err;
        }
        errmax /= eps;

        // --- error too large; reduce stepsize & repeat
        if (errmax > 1.0)
        {
            htemp = SAFETY*h*pow(errmax,PSHRNK);
            if (h >= 0)
            {
                if (htemp > 0.1*h) h = htemp;
                else h = 0.1*h;
            }
            else
            {
                if (htemp < 0.1*h) h = htemp;
                else h = 0.1*h;
            }
            xnew = xold + h;
            if (xnew == xold) return 2;
            continue;
        }

        // --- step succeeded; compute size of next step
        else
        {
            if (errmax > ERRCON) *hnext = SAFETY*h*pow(errmax,PGROW);
            else *hnext = 5.0*h;
            *x += (*hdid=h);
            for (i=0; i<n; i++) odslv->y[i] = odslv->ytemp[i];
            return 0;
        }
    }
}


void rkck(SWMM_Project *sp, double x, int n, double h,
        void (*derivs)(SWMM_Project*, double, double*, double*))
//----------------------------------------------------------------------
//   Uses the Runge-Kutta-Cash-Karp method to advance y[] at x
//   over stepsize h.
//----------------------------------------------------------------------
{
    TOdesolveShared *odslv = &sp->OdesolveShared;

    const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
           b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42= -0.9, b43=1.2,
           b51= -11.0/54.0, b52=2.5, b53= -70.0/27.0, b54=35.0/27.0,
           b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0,
           b64=44275.0/110592.0, b65=253.0/4096.0, c1=37.0/378.0,
           c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
           dc5= -277.0/14336.0;
    const double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
           dc4=c4-13525.0/55296.0, dc6=c6-0.25;
    int i;
    double *ak2 = (odslv->ak);
    double *ak3 = ((odslv->ak)+(n));
    double *ak4 = ((odslv->ak)+(2*n));
    double *ak5 = ((odslv->ak)+(3*n));
    double *ak6 = ((odslv->ak)+(4*n));

    for (i=0; i<n; i++)
        odslv->ytemp[i] = odslv->y[i] + b21*h*odslv->dydx[i];
    derivs(sp, x+a2*h,odslv->ytemp,ak2);

    for (i=0; i<n; i++)
        odslv->ytemp[i] = odslv->y[i] + h*(b31*odslv->dydx[i]+b32*ak2[i]);
    derivs(sp, x+a3*h,odslv->ytemp,ak3);

    for (i=0; i<n; i++)
        odslv->ytemp[i] = odslv->y[i] + h*(b41*odslv->dydx[i]+b42*ak2[i] + b43*ak3[i]);
    derivs(sp, x+a4*h,odslv->ytemp,ak4);

    for (i=0; i<n; i++)
        odslv->ytemp[i] = odslv->y[i] + h*(b51*odslv->dydx[i]+b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);
    derivs(sp, x+a5*h,odslv->ytemp,ak5);

    for (i=0; i<n; i++)
        odslv->ytemp[i] = odslv->y[i] + h*(b61*odslv->dydx[i]+b62*ak2[i] + b63*ak3[i] + b64*ak4[i]
                   + b65*ak5[i]);
    derivs(sp, x+a6*h,odslv->ytemp,ak6);

    for (i=0; i<n; i++)
        odslv->ytemp[i] = odslv->y[i] + h*(c1*odslv->dydx[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]);

    for (i=0; i<n; i++)
        odslv->yerr[i] = h*(dc1*odslv->dydx[i] +dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]);
}
