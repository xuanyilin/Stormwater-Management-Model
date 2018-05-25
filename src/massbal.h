/*
 * massbal.h
 *
 *  Created on: May 3, 2018
 *      Author: mtryby
 */

#ifndef SRC_MASSBAL_H_
#define SRC_MASSBAL_H_


//-------------------------------
// CUMULATIVE RUNOFF TOTALS
//-------------------------------
typedef struct
{                                 // All volume totals are in ft3.             //(5.1.008)
   double        rainfall;        // rainfall volume
   double        evap;            // evaporation loss
   double        infil;           // infiltration loss
   double        runoff;          // runoff volume
   double        drains;          // LID drains                                //(5.1.008)
   double        runon;           // runon from outfalls                       //(5.1.008)
   double        initStorage;     // inital surface storage
   double        finalStorage;    // final surface storage
   double        initSnowCover;   // initial snow cover
   double        finalSnowCover;  // final snow cover
   double        snowRemoved;     // snow removal
   double        pctError;        // continuity error (%)
}  TRunoffTotals;


//--------------------------
// CUMULATIVE LOADING TOTALS
//--------------------------
typedef struct
{                                 // All loading totals are in lbs.
   double        initLoad;        // initial loading
   double        buildup;         // loading added from buildup
   double        deposition;      // loading added from wet deposition
   double        sweeping;        // loading removed by street sweeping
   double        bmpRemoval;      // loading removed by BMPs
   double        infil;           // loading removed by infiltration
   double        runoff;          // loading removed by runoff
   double        finalLoad;       // final loading
   double        pctError;        // continuity error (%)
}  TLoadingTotals;


//------------------------------
// CUMULATIVE GROUNDWATER TOTALS
//------------------------------
typedef struct
{                                 // All GW flux totals are in feet.
   double        infil;           // surface infiltration
   double        upperEvap;       // upper zone evaporation loss
   double        lowerEvap;       // lower zone evaporation loss
   double        lowerPerc;       // percolation out of lower zone
   double        gwater;          // groundwater flow
   double        initStorage;     // initial groundwater storage
   double        finalStorage;    // final groundwater storage
   double        pctError;        // continuity error (%)
}  TGwaterTotals;


//----------------------------
// CUMULATIVE ROUTING TOTALS
//----------------------------
typedef struct
{                                  // All routing totals are in ft3.
   double        dwInflow;         // dry weather inflow
   double        wwInflow;         // wet weather inflow
   double        gwInflow;         // groundwater inflow
   double        iiInflow;         // RDII inflow
   double        exInflow;         // direct inflow
   double        flooding;         // internal flooding
   double        outflow;          // external outflow
   double        evapLoss;         // evaporation loss
   double        seepLoss;         // seepage loss
   double        reacted;          // reaction losses
   double        initStorage;      // initial storage volume
   double        finalStorage;     // final storage volume
   double        pctError;         // continuity error
}  TRoutingTotals;


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    TRunoffTotals    RunoffTotals;    // overall surface runoff continuity totals
    TLoadingTotals*  LoadingTotals;   // overall WQ washoff continuity totals
    TGwaterTotals    GwaterTotals;    // overall groundwater continuity totals
    TRoutingTotals   FlowTotals;      // overall routed flow continuity totals
    TRoutingTotals*  QualTotals;      // overall routed WQ continuity totals
    TRoutingTotals   StepFlowTotals;  // routed flow totals over time step
    TRoutingTotals   OldStepFlowTotals;
    TRoutingTotals*  StepQualTotals;  // routed WQ totals over time step
} TMassbalShared;


//-----------------------------------------------------------------------------
//  Exportable variables
//-----------------------------------------------------------------------------
typedef struct
{
    double*  NodeInflow;              // total inflow volume to each node (ft3)
    double*  NodeOutflow;             // total outflow volume from each node (ft3)
    double   TotalArea;               // total drainage area (ft2)
} TMassbalExport;


#endif /* SRC_MASSBAL_H_ */
