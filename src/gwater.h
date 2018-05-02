/*
 * gwater.h
 *
 *  Created on: May 2, 2018
 *      Author: mtryby
 */

#ifndef SRC_GWATER_H_
#define SRC_GWATER_H_


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
