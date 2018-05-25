/*
 * swmm-internal.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_SWMM_INTERNAL_H_
#define SRC_SWMM_INTERNAL_H_


int  swmm_IsOpenFlag(SWMM_Project *sp);
int  swmm_IsStartedFlag(SWMM_Project *sp);

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    int  IsOpenFlag;           // TRUE if a project has been opened
    int  IsStartedFlag;        // TRUE if a simulation has been started
    int  SaveResultsFlag;      // TRUE if output to be saved to binary file
    int  ExceptionCount;       // number of exceptions handled
    int  DoRunoff;             // TRUE if runoff is computed
    int  DoRouting;            // TRUE if flow routing is computed
} TSwmmShared;


#endif /* SRC_SWMM_INTERNAL_H_ */
