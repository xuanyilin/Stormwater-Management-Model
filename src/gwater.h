/*
 * gwater.h
 *
 *  Created on: May 2, 2018
 *      Author: mtryby
 */

#ifndef SRC_GWATER_H_
#define SRC_GWATER_H_


//-------------------
// AQUIFER OBJECT
//-------------------
typedef struct
{
    char*       ID;               // aquifer name
    double      porosity;         // soil porosity
    double      wiltingPoint;     // soil wilting point
    double      fieldCapacity;    // soil field capacity
    double      conductivity;     // soil hyd. conductivity (ft/sec)
    double      conductSlope;     // slope of conductivity v. moisture curve
    double      tensionSlope;     // slope of tension v. moisture curve
    double      upperEvapFrac;    // evaporation available in upper zone
    double      lowerEvapDepth;   // evap depth existing in lower zone (ft)
    double      lowerLossCoeff;   // coeff. for losses to deep GW (ft/sec)
    double      bottomElev;       // elevation of bottom of aquifer (ft)
    double      waterTableElev;   // initial water table elevation (ft)
    double      upperMoisture;    // initial moisture content of unsat. zone
    int         upperEvapPat;     // monthly upper evap. adjustment factors
} TAquifer;


////  Added to release 5.1.008.  ////                                          //(5.1.008)
//-----------------------
// GROUNDWATER STATISTICS
//-----------------------
typedef struct
{
    double       infil;           // total infiltration (ft)
    double       evap;            // total evaporation (ft)
    double       latFlow;         // total lateral outflow (ft)
    double       deepFlow;        // total flow to deep aquifer (ft)
    double       avgUpperMoist;   // avg. upper zone moisture
    double       finalUpperMoist; // final upper zone moisture
    double       avgWaterTable;   // avg. water table height (ft)
    double       finalWaterTable; // final water table height (ft)
    double       maxFlow;         // max. lateral outflow (cfs)
}  TGWaterStats;


//------------------------
// GROUNDWATER OBJECT
//------------------------
typedef struct
{
    int           aquifer;        // index of associated gw aquifer
    int           node;           // index of node receiving gw flow
    double        surfElev;       // elevation of ground surface (ft)
    double        a1, b1;         // ground water outflow coeff. & exponent
    double        a2, b2;         // surface water outflow coeff. & exponent
    double        a3;             // surf./ground water interaction coeff.
    double        fixedDepth;     // fixed surface water water depth (ft)
    double        nodeElev;       // elevation of receiving node invert (ft)
    double        bottomElev;     // bottom elevation of lower GW zone (ft)
    double        waterTableElev; // initial water table elevation (ft)
    double        upperMoisture;  // initial moisture content of unsat. zone
    //----------------------------
    double        theta;          // upper zone moisture content
    double        lowerDepth;     // depth of saturated zone (ft)
    double        oldFlow;        // gw outflow from previous time period (fps)  //(5.1.011)
    double        newFlow;        // gw outflow from current time period (fps)   //(5.1.011)
    double        evapLoss;       // evaporation loss rate (ft/sec)
    double        maxInfilVol;    // max. infil. upper zone can accept (ft)
    TGWaterStats  stats;          // gw statistics                             //(5.1.008)
} TGroundwater;


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
//  NOTE: all flux rates are in ft/sec, all depths are in ft.
typedef struct
{
    double    Area;            // subcatchment area (ft2)                   //(5.1.008)
    double    Infil;           // infiltration rate from surface
    double    MaxEvap;         // max. evaporation rate
    double    AvailEvap;       // available evaporation rate
    double    UpperEvap;       // evaporation rate from upper GW zone
    double    LowerEvap;       // evaporation rate from lower GW zone
    double    UpperPerc;       // percolation rate from upper to lower zone
    double    LowerLoss;       // loss rate from lower GW zone
    double    GWFlow;          // flow rate from lower zone to conveyance node
    double    MaxUpperPerc;    // upper limit on UpperPerc
    double    MaxGWFlowPos;    // upper limit on GWFlow when its positve
    double    MaxGWFlowNeg;    // upper limit on GWFlow when its negative
    double    FracPerv;        // fraction of surface that is pervious
    double    TotalDepth;      // total depth of GW aquifer
    double    Theta;           // moisture content of upper zone
    double    HydCon;          // unsaturated hydraulic conductivity (ft/s) //(5.1.010)
    double    Hgw;             // ht. of saturated zone
    double    Hstar;           // ht. from aquifer bottom to node invert
    double    Hsw;             // ht. from aquifer bottom to water surface
    double    Tstep;           // current time step (sec)
    TAquifer  A;               // aquifer being analyzed
    TGroundwater* GW;          // groundwater object being analyzed
    MathExpr* LatFlowExpr;     // user-supplied lateral GW flow expression  //(5.1.007)
    MathExpr* DeepFlowExpr;    // user-supplied deep GW flow expression     //(5.1.007)
} TGwaterShared;


#endif /* SRC_GWATER_H_ */
