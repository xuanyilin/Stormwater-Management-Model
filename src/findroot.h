//-----------------------------------------------------------------------------
//   findroot.h
//
//   Header file for root finding method contained in findroot.c
//
//   Last modified on 11/19/13.
//-----------------------------------------------------------------------------

#ifndef FINDROOT_H
#define FINDROOT_H


int findroot_Newton(SWMM_Project *sp, double x1, double x2, double* rts, double xacc,
        void (*func) (SWMM_Project *sp, double x, double* f, double* df,
                void* p), void* p);

double findroot_Ridder(SWMM_Project *sp, double x1, double x2, double xacc,
	                   double (*func)(SWMM_Project *sp, double, void* p), void* p);


#endif //FINDROOT_H
