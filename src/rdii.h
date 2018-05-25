/*
 * rdii.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_RDII_H_
#define SRC_RDII_H_


typedef struct                         // Data for a single unit hydrograph
{                                      // -------------------------------------
   double*   pastRain;                 // array of past rainfall values
   char*     pastMonth;                // month in which past rainfall occurred
   int       period;                   // current UH time period
   int       hasPastRain;              // true if > 0 past periods with rain
   int       maxPeriods;               // max. past rainfall periods
   long      drySeconds;               // time since last nonzero rainfall
   double    iaUsed;                   // initial abstraction used (in or mm)
}  TUHData;

typedef struct                         // Data for a unit hydrograph group
{                                      //---------------------------------
   int       isUsed;                   // true if UH group used by any nodes
   int       rainInterval;             // time interval for RDII processing (sec)
   double    area;                     // sewered area covered by UH's gage (ft2)
   double    rdii;                     // rdii flow (in rainfall units)
   DateTime  gageDate;                 // calendar date of rain gage period
   DateTime  lastDate;                 // date of last rdii computed
   TUHData   uh[3];                    // data for each unit hydrograph
}  TUHGroup;

//-----------------------------------------------------------------------------
// Shared Variables
//-----------------------------------------------------------------------------
typedef struct
{
    TUHGroup*  UHGroup;             // processing data for each UH group
    int        RdiiStep;            // RDII time step (sec)
    int        NumRdiiNodes;        // number of nodes w/ RDII data
    int*       RdiiNodeIndex;       // indexes of nodes w/ RDII data
    REAL4*     RdiiNodeFlow;        // inflows for nodes with RDII          //(5.1.003)
    int        RdiiFlowUnits;       // RDII flow units code
    DateTime   RdiiStartDate;       // start date of RDII inflow period
    DateTime   RdiiEndDate;         // end date of RDII inflow period
    double     TotalRainVol;        // total rainfall volume (ft3)
    double     TotalRdiiVol;        // total RDII volume (ft3)
    int        RdiiFileType;        // type (binary/text) of RDII file
} TRdiiShared;


#endif /* SRC_RDII_H_ */
