/*
 * massbal.h
 *
 *  Created on: May 3, 2018
 *      Author: mtryby
 */

#ifndef SRC_MASSBAL_H_
#define SRC_MASSBAL_H_


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


#endif /* SRC_MASSBAL_H_ */
