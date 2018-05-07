/*
 * treatment.h
 *
 *  Created on: May 7, 2018
 *      Author: mtryby
 */

#ifndef SRC_TREATMNT_H_
#define SRC_TREATMNT_H_


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    int     ErrCode;                // treatment error code
    int     J;                      // index of node being analyzed
    double  Dt;                     // curent time step (sec)
    double  Q;                      // node inflow (cfs)
    double  V;                      // node volume (ft3)
    double* R;                      // array of pollut. removals
    double* Cin;                    // node inflow concentrations
//static TTreatment* Treatment; // defined locally in treatmnt_treat()         //(5.1.008)
} TTreatmntShared;


#endif /* SRC_TREATMNT_H_ */
