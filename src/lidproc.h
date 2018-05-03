/*
 * lidproc.h
 *
 *  Created on: May 3, 2018
 *      Author: mtryby
 */

#ifndef SRC_LIDPROC_H_
#define SRC_LIDPROC_H_


//-----------------------------------------------------------------------------
//  Local Variables
//-----------------------------------------------------------------------------
typedef struct
{
    TLidUnit*  theLidUnit;     // ptr. to a subcatchment's LID unit
    TLidProc*  theLidProc;     // ptr. to a LID process

    double     Tstep;          // current time step (sec)
    //static double   Rainfall;       // current rainfall rate (ft/s)              //(5.1.008)
    double     EvapRate;       // evaporation rate (ft/s)
    double     MaxNativeInfil; // native soil infil. rate limit (ft/s)

    double     SurfaceInflow;  // precip. + runon to LID unit (ft/s)
    double     SurfaceInfil;   // infil. rate from surface layer (ft/s)
    double     SurfaceEvap;    // evap. rate from surface layer (ft/s)
    double     SurfaceOutflow; // outflow from surface layer (ft/s)
    double     SurfaceVolume;  // volume in surface storage (ft)

    double     PaveEvap;       // evap. from pavement layer (ft/s)          //(5.1.008)
    double     PavePerc;       // percolation from pavement layer (ft/s)    //(5.1.008)
    double     PaveVolume;     // volume stored in pavement layer  (ft)     //(5.1.008)

    double     SoilEvap;       // evap. from soil layer (ft/s)
    double     SoilPerc;       // percolation from soil layer (ft/s)
    double     SoilVolume;     // volume in soil/pavement storage (ft)

    double     StorageInflow;  // inflow rate to storage layer (ft/s)
    double     StorageExfil;   // exfil. rate from storage layer (ft/s)     //(5.1.011)
    double     StorageEvap;    // evap.rate from storage layer (ft/s)
    double     StorageDrain;   // underdrain flow rate layer (ft/s)
    double     StorageVolume;  // volume in storage layer (ft)

    double     Xold[MAX_LAYERS];  // previous moisture level in LID layers  //(5.1.008)
} TLidprocShared;


#endif /* SRC_LIDPROC_H_ */
