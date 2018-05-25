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


//------------------------
// SUBCATCHMENT STATISTICS
//------------------------
typedef struct
{
    double       precip;
    double       runon;
    double       evap;
    double       infil;
    double       runoff;
    double       maxFlow;
    double*      surfaceBuildup;
}  TSubcatchStats;


//----------------
// NODE STATISTICS
//----------------
typedef struct
{
   double        avgDepth;
   double        maxDepth;
   DateTime      maxDepthDate;
   double        maxRptDepth;                                                  //(5.1.008)
   double        volFlooded;
   double        timeFlooded;
   double        timeSurcharged;
   double        timeCourantCritical;
   double        totLatFlow;
   double        maxLatFlow;
   double        maxInflow;
   double        maxOverflow;
   double        maxPondedVol;
   DateTime      maxInflowDate;
   DateTime      maxOverflowDate;
}  TNodeStats;


//-------------------
// STORAGE STATISTICS
//-------------------
typedef struct
{
   double        initVol;
   double        avgVol;
   double        maxVol;
   double        maxFlow;
   double        evapLosses;
   double        exfilLosses;
   DateTime      maxVolDate;
}  TStorageStats;


//-------------------
// OUTFALL STATISTICS
//-------------------
typedef struct
{
   double       avgFlow;
   double       maxFlow;
   double*      totalLoad;
   int          totalPeriods;
}  TOutfallStats;


//----------------
// PUMP STATISTICS
//----------------
typedef struct
{
   double       utilized;
   double       minFlow;
   double       avgFlow;
   double       maxFlow;
   double       volume;
   double       energy;
   double       offCurveLow;
   double       offCurveHigh;
   int          startUps;
   int          totalPeriods;
}  TPumpStats;


//----------------
// LINK STATISTICS
//----------------
typedef struct
{
   double        maxFlow;
   DateTime      maxFlowDate;
   double        maxVeloc;
   //DateTime      maxVelocDate;  //deprecated                                 //(5.1.008)
   double        maxDepth;
   double        timeNormalFlow;
   double        timeInletControl;
   double        timeSurcharged;
   double        timeFullUpstream;
   double        timeFullDnstream;
   double        timeFullFlow;
   double        timeCapacityLimited;
   double        timeInFlowClass[MAX_FLOW_CLASSES];
   double        timeCourantCritical;
   long          flowTurns;
   int           flowTurnSign;
}  TLinkStats;


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


//-----------------------------------------------------------------------------
//  Exportable variables (shared with statsrpt.c)
//-----------------------------------------------------------------------------
typedef struct
{
    TSubcatchStats*     SubcatchStats;
    TNodeStats*         NodeStats;
    TLinkStats*         LinkStats;
    TStorageStats*      StorageStats;
    TOutfallStats*      OutfallStats;
    TPumpStats*         PumpStats;
    double              MaxOutfallFlow;
    double              MaxRunoffFlow;
} TStatsExport;


#endif /* SRC_STATS_H_ */
