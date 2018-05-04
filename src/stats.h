/*
 * stats.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_STATS_H_
#define SRC_STATS_H_


#define MAX_STATS 5

//-----------------------
// SYSTEM-WIDE STATISTICS
//-----------------------
typedef struct
{
   double        minTimeStep;
   double        maxTimeStep;
   double        avgTimeStep;
   double        avgStepCount;
   double        steadyStateCount;
}  TSysStats;

//-------------------------
// MAXIMUM VALUE STATISTICS
//-------------------------
typedef struct
{
   int           objType;         // either NODE or LINK
   int           index;           // node or link index
   double        value;           // value of node or link statistic
}  TMaxStats;

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    TSysStats       SysStats;
    TMaxStats       MaxMassBalErrs[MAX_STATS];
    TMaxStats       MaxCourantCrit[MAX_STATS];
    TMaxStats       MaxFlowTurns[MAX_STATS];
    double          SysOutfallFlow;
} TStatsShared;


#endif /* SRC_STATS_H_ */
