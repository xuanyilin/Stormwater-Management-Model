/*
 * dynwave.h
 *
 *  Created on: May 2, 2018
 *      Author: mtryby
 */

#ifndef SRC_DYNWAVE_H_
#define SRC_DYNWAVE_H_


//-----------------------------------------------------------------------------
//  Data Structures
//-----------------------------------------------------------------------------
typedef struct
{
    char    converged;                 // TRUE if iterations for a node done
    double  newSurfArea;               // current surface area (ft2)
    double  oldSurfArea;               // previous surface area (ft2)
    double  sumdqdh;                   // sum of dqdh from adjoining links
    double  dYdT;                      // change in depth w.r.t. time (ft/sec)
} TXnode;

//-----------------------------------------------------------------------------
//  Shared Variables
//-----------------------------------------------------------------------------
typedef struct
{
    double  VariableStep;           // size of variable time step (sec)
    TXnode* Xnode;                  // extended nodal information

    double  Omega;                  // actual under-relaxation parameter
    int     Steps;                  // number of Picard iterations
} TDynwaveShared;


#endif /* SRC_DYNWAVE_H_ */
