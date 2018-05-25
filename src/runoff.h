/*
 * runoff.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_RUNOFF_H_
#define SRC_RUNOFF_H_


typedef struct
{
    //-----------------------------------------------------------------------------
    //  Exportable variables
    //-----------------------------------------------------------------------------
    char    HasWetLids;  // TRUE if any LIDs are wet (used in lidproc.c)           //(5.1.008)
    double* OutflowLoad; // exported pollutant mass load (used in surfqual.c)

    //-----------------------------------------------------------------------------
    // Shared variables
    //-----------------------------------------------------------------------------
    char  IsRaining;                // TRUE if precip. falls on study area
    char  HasRunoff;                // TRUE if study area generates runoff
    char  HasSnow;                  // TRUE if any snow cover on study area
    int   Nsteps;                   // number of runoff time steps taken
    int   MaxSteps;                 // final number of runoff time steps
    long  MaxStepsPos;              // position in Runoff interface file
                                    //    where MaxSteps is saved
} TRunoffShared;


#endif /* SRC_RUNOFF_H_ */
