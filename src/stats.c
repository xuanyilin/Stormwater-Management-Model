//-----------------------------------------------------------------------------
//   stats.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             09/15/14   (Build 5.1.007)
//             03/19/15   (Build 5.1.008)
//             08/01/16   (Build 5.1.011)
//             03/14/17   (Build 5.1.012)
//   Author:   L. Rossman (EPA)
//             R. Dickinson (CDM)
//
//   Simulation statistics functions.
//
//   Build 5.1.007:
//   - Exfiltration losses added to storage node statistics.
//
//   Build 5.1.008:
//   - Support for updating groundwater statistics added.
//   - Support for updating maximum reported nodal depths added.
//   - OpenMP parallelization applied to updating node and link flow statistics.
//   - Updating of time that conduit is upstrm/dnstrm full was modified.
//
//   Build 5.1.011:
//   - Surcharging is now evaluated only under dynamic wave flow routing and
//     storage nodes cannot be classified as surcharged.
//
//   Build 5.1.012:
//   - Time step statistics now evaluated only in non-steady state periods.
//   - Check for full conduit flow now accounts for number of barrels.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#if defined(_OPENMP)
  #include <omp.h>                                                             //(5.1.008)
#endif
#include "headers.h"
#include "swmm5.h"

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
//extern double*         NodeInflow;     // defined in massbal.c
//extern double*         NodeOutflow;    // defined in massbal.c

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  stats_open                    (called from swmm_start in swmm5.c)
//  stats_close                   (called from swmm_end in swmm5.c)
//  stats_report                  (called from swmm_end in swmm5.c)
//  stats_updateSubcatchStats     (called from subcatch_getRunoff)
//  stats_updateGwaterStats       (called from gwater_getGroundwater)          //(5.1.008)
//  stats_updateFlowStats         (called from routing_execute)
//  stats_updateCriticalTimeCount (called from getVariableStep in dynwave.c)
//  stats_updateMaxNodeDepth      (called from output_saveNodeResults)         //(5.1.008)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void stats_updateNodeStats(SWMM_Project *sp, int node, double tStep,
        DateTime aDate);
static void stats_updateLinkStats(SWMM_Project *sp, int link, double tStep,
        DateTime aDate);
static void stats_findMaxStats(SWMM_Project *sp);
static void stats_updateMaxStats(TMaxStats maxStats[], int i, int j, double x);

//=============================================================================

int  stats_open(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: opens the simulation statistics system.
//
{
    int j, k, p;

    TStatsShared *stts = &sp->StatsShared;
    TStatsExport *sttsx = &sp->StatsExport;

    // --- set all pointers to NULL
    sttsx->NodeStats = NULL;
    sttsx->LinkStats = NULL;
    sttsx->StorageStats = NULL;
    sttsx->OutfallStats = NULL;
    sttsx->PumpStats = NULL;

    // --- allocate memory for & initialize subcatchment statistics
    sttsx->SubcatchStats = NULL;
    if ( sp->Nobjects[SUBCATCH] > 0 )
    {
        sttsx->SubcatchStats = (TSubcatchStats *) calloc(sp->Nobjects[SUBCATCH],
                                               sizeof(TSubcatchStats));
        if ( !sttsx->SubcatchStats )
        {
            report_writeErrorMsg(sp, ERR_MEMORY, "");
            return sp->ErrorCode;
        }
        for (j=0; j<sp->Nobjects[SUBCATCH]; j++)
        {
            sttsx->SubcatchStats[j].precip  = 0.0;
            sttsx->SubcatchStats[j].runon   = 0.0;
            sttsx->SubcatchStats[j].evap    = 0.0;
            sttsx->SubcatchStats[j].infil   = 0.0;
            sttsx->SubcatchStats[j].runoff  = 0.0;
            sttsx->SubcatchStats[j].maxFlow = 0.0;

            if ( sp->Nobjects[POLLUT] > 0 )
            {
                sttsx->SubcatchStats[j].surfaceBuildup =
                    (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
                if ( !sttsx->SubcatchStats[j].surfaceBuildup )
                {
                    report_writeErrorMsg(sp, ERR_MEMORY, "");
                    return sp->ErrorCode;
                }
                for ( p = 0; p < sp->Nobjects[POLLUT]; p++ )
                    sttsx->SubcatchStats[j].surfaceBuildup[p] = 0.0;
            }
            else sttsx->SubcatchStats[j].surfaceBuildup = NULL;
        }

////  Added to release 5.1.008.  ////                                          //(5.1.008)
////
        for (j=0; j<sp->Nobjects[SUBCATCH]; j++)
        {
            if ( sp->Subcatch[j].groundwater == NULL ) continue;
            sp->Subcatch[j].groundwater->stats.avgUpperMoist = 0.0;
            sp->Subcatch[j].groundwater->stats.avgWaterTable = 0.0;
            sp->Subcatch[j].groundwater->stats.infil = 0.0;
            sp->Subcatch[j].groundwater->stats.latFlow = 0.0;
            sp->Subcatch[j].groundwater->stats.deepFlow = 0.0;
            sp->Subcatch[j].groundwater->stats.evap = 0.0;
            sp->Subcatch[j].groundwater->stats.maxFlow = 0.0;
        }
////
    }

    // --- allocate memory for node & link stats
    if ( sp->Nobjects[LINK] > 0 )
    {
        sttsx->NodeStats = (TNodeStats *) calloc(sp->Nobjects[NODE], sizeof(TNodeStats));
        sttsx->LinkStats = (TLinkStats *) calloc(sp->Nobjects[LINK], sizeof(TLinkStats));
        if ( !sttsx->NodeStats || !sttsx->LinkStats )
        {
            report_writeErrorMsg(sp, ERR_MEMORY, "");
            return sp->ErrorCode;
        }
    }

    // --- initialize node stats
    if ( sttsx->NodeStats ) for ( j = 0; j < sp->Nobjects[NODE]; j++ )
    {
        sttsx->NodeStats[j].avgDepth = 0.0;
        sttsx->NodeStats[j].maxDepth = 0.0;
        sttsx->NodeStats[j].maxDepthDate = sp->StartDateTime;
        sttsx->NodeStats[j].maxRptDepth = 0.0;                                        //(5.1.008)
        sttsx->NodeStats[j].volFlooded = 0.0;
        sttsx->NodeStats[j].timeFlooded = 0.0;
        sttsx->NodeStats[j].timeSurcharged = 0.0;
        sttsx->NodeStats[j].timeCourantCritical = 0.0;
        sttsx->NodeStats[j].totLatFlow = 0.0;
        sttsx->NodeStats[j].maxLatFlow = 0.0;
        sttsx->NodeStats[j].maxInflow = 0.0;
        sttsx->NodeStats[j].maxOverflow = 0.0;
        sttsx->NodeStats[j].maxPondedVol = 0.0;
        sttsx->NodeStats[j].maxInflowDate = sp->StartDateTime;
        sttsx->NodeStats[j].maxOverflowDate = sp->StartDateTime;
    }

    // --- initialize link stats
    if ( sttsx->LinkStats ) for ( j = 0; j < sp->Nobjects[LINK]; j++ )
    {
        sttsx->LinkStats[j].maxFlow = 0.0;
        sttsx->LinkStats[j].maxVeloc = 0.0;
        sttsx->LinkStats[j].maxDepth = 0.0;
        sttsx->LinkStats[j].timeSurcharged = 0.0;
        sttsx->LinkStats[j].timeFullUpstream = 0.0;
        sttsx->LinkStats[j].timeFullDnstream = 0.0;
        sttsx->LinkStats[j].timeFullFlow = 0.0;
        sttsx->LinkStats[j].timeCapacityLimited = 0.0;
        sttsx->LinkStats[j].timeCourantCritical = 0.0;
        for (k=0; k<MAX_FLOW_CLASSES; k++)
            sttsx->LinkStats[j].timeInFlowClass[k] = 0.0;
        sttsx->LinkStats[j].flowTurns = 0;
        sttsx->LinkStats[j].flowTurnSign = 0;
    }

    // --- allocate memory for & initialize storage unit statistics
    if ( sp->Nnodes[STORAGE] > 0 )
    {
        sttsx->StorageStats = (TStorageStats *) calloc(sp->Nnodes[STORAGE],
                           sizeof(TStorageStats));
        if ( !sttsx->StorageStats )
        {
            report_writeErrorMsg(sp, ERR_MEMORY, "");
            return sp->ErrorCode;
        }
        else for ( k = 0; k < sp->Nobjects[NODE]; k++ )
        {
            if ( sp->Node[k].type != STORAGE ) continue;
            j = sp->Node[k].subIndex;
            sttsx->StorageStats[j].initVol = sp->Node[k].newVolume;
            sttsx->StorageStats[j].avgVol = 0.0;
            sttsx->StorageStats[j].maxVol = 0.0;
            sttsx->StorageStats[j].maxFlow = 0.0;
            sttsx->StorageStats[j].evapLosses = 0.0;
            sttsx->StorageStats[j].exfilLosses = 0.0;                                 //(5.1.007)
            sttsx->StorageStats[j].maxVolDate = sp->StartDateTime;
        }
    }

    // --- allocate memory for & initialize outfall statistics
    if ( sp->Nnodes[OUTFALL] > 0 )
    {
        sttsx->OutfallStats = (TOutfallStats *) calloc(sp->Nnodes[OUTFALL],
                           sizeof(TOutfallStats));
        if ( !sttsx->OutfallStats )
        {
            report_writeErrorMsg(sp, ERR_MEMORY, "");
            return sp->ErrorCode;
        }
        else for ( j = 0; j < sp->Nnodes[OUTFALL]; j++ )
        {
            sttsx->OutfallStats[j].avgFlow = 0.0;
            sttsx->OutfallStats[j].maxFlow = 0.0;
            sttsx->OutfallStats[j].totalPeriods = 0;
            if ( sp->Nobjects[POLLUT] > 0 )
            {
                sttsx->OutfallStats[j].totalLoad =
                    (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
                if ( !sttsx->OutfallStats[j].totalLoad )
                {
                    report_writeErrorMsg(sp, ERR_MEMORY, "");
                    return sp->ErrorCode;
                }
                for (k=0; k<sp->Nobjects[POLLUT]; k++)
                    sttsx->OutfallStats[j].totalLoad[k] = 0.0;
            }
            else sttsx->OutfallStats[j].totalLoad = NULL;
        }
    }

    // --- allocate memory & initialize pumping statistics
    if ( sp->Nlinks[PUMP] > 0 )
    { 
        sttsx->PumpStats = (TPumpStats *) calloc(sp->Nlinks[PUMP], sizeof(TPumpStats));
        if ( !sttsx->PumpStats )
        {
            report_writeErrorMsg(sp, ERR_MEMORY, "");
            return sp->ErrorCode;
        }
        else for ( j = 0; j < sp->Nlinks[PUMP]; j++ )
        {
            sttsx->PumpStats[j].utilized = 0.0;
            sttsx->PumpStats[j].minFlow  = 0.0;
            sttsx->PumpStats[j].avgFlow  = 0.0;
            sttsx->PumpStats[j].maxFlow  = 0.0;
            sttsx->PumpStats[j].volume   = 0.0;
            sttsx->PumpStats[j].energy   = 0.0;
            sttsx->PumpStats[j].startUps = 0;
            sttsx->PumpStats[j].offCurveLow = 0.0;
            sttsx->PumpStats[j].offCurveHigh = 0.0;
        } 
    } 

    // --- initialize system stats
    sttsx->MaxRunoffFlow = 0.0;
    sttsx->MaxOutfallFlow = 0.0;
    stts->SysStats.maxTimeStep = 0.0;
    stts->SysStats.minTimeStep = sp->RouteStep;
    stts->SysStats.avgTimeStep = 0.0;
    stts->SysStats.avgStepCount = 0.0;
    stts->SysStats.steadyStateCount = 0.0;
    return 0;
}

//=============================================================================

void  stats_close(SWMM_Project *sp)
//
//  Input:   none
//  Output:  
//  Purpose: closes the simulation statistics system.
//
{
    int j;

    TStatsExport *sttsx = &sp->StatsExport;

    if ( sttsx->SubcatchStats)
    {
        for ( j=0; j<sp->Nobjects[SUBCATCH]; j++ )
            FREE(sttsx->SubcatchStats[j].surfaceBuildup);
        FREE(sttsx->SubcatchStats);
    }
    FREE(sttsx->NodeStats);
    FREE(sttsx->LinkStats);
    FREE(sttsx->StorageStats);
    if ( sttsx->OutfallStats )
    {
        for ( j=0; j<sp->Nnodes[OUTFALL]; j++ )
            FREE(sttsx->OutfallStats[j].totalLoad);
        FREE(sttsx->OutfallStats);
    }
    FREE(sttsx->PumpStats);
}

//=============================================================================

void  stats_report(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reports simulation statistics.
//
{
    TStatsShared *stts = &sp->StatsShared;

    // --- report flow routing accuracy statistics
    if ( sp->Nobjects[LINK] > 0 && sp->RouteModel != NO_ROUTING )
    {
        stats_findMaxStats(sp);
        report_writeMaxStats(sp, stts->MaxMassBalErrs, stts->MaxCourantCrit, MAX_STATS);
        report_writeMaxFlowTurns(sp, stts->MaxFlowTurns, MAX_STATS);
        report_writeSysStats(sp, &stts->SysStats);
    }

    // --- report summary statistics
    statsrpt_writeReport(sp);
}

//=============================================================================

void   stats_updateSubcatchStats(SWMM_Project *sp, int j, double rainVol,
        double runonVol, double evapVol, double infilVol, double runoffVol,
        double runoff)
//
//  Input:   j = subcatchment index
//           rainVol   = rainfall + snowfall volume (ft3)
//           runonVol  = runon volume from other subcatchments (ft3)
//           evapVol   = evaporation volume (ft3)
//           infilVol  = infiltration volume (ft3)
//           runoffVol = runoff volume (ft3)
//           runoff    = runoff rate (cfs)
//  Output:  none
//  Purpose: updates totals of runoff components and the surface buildup
//           of pollutants for a specific subcatchment.
//
{
    int p;
    
    TStatsExport *sttsx = &sp->StatsExport;

    sttsx->SubcatchStats[j].precip += rainVol;
    sttsx->SubcatchStats[j].runon  += runonVol;
    sttsx->SubcatchStats[j].evap   += evapVol;
    sttsx->SubcatchStats[j].infil  += infilVol;
    sttsx->SubcatchStats[j].runoff += runoffVol;
    sttsx->SubcatchStats[j].maxFlow = MAX(sttsx->SubcatchStats[j].maxFlow, runoff);
    
    for ( p = 0; p < sp->Nobjects[POLLUT]; p++ )
    {
        sttsx->SubcatchStats[j].surfaceBuildup[p] = subcatch_getBuildup(sp, j, p);
    }
        
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

void  stats_updateGwaterStats(SWMM_Project *sp, int j, double infil, double evap,
        double latFlow, double deepFlow, double theta, double waterTable,
        double tStep)
{
    sp->Subcatch[j].groundwater->stats.infil += infil * tStep;
    sp->Subcatch[j].groundwater->stats.evap += evap * tStep;
    sp->Subcatch[j].groundwater->stats.latFlow += latFlow * tStep;
    sp->Subcatch[j].groundwater->stats.deepFlow += deepFlow * tStep;
    sp->Subcatch[j].groundwater->stats.avgUpperMoist += theta * tStep;
    sp->Subcatch[j].groundwater->stats.avgWaterTable += waterTable * tStep;
    sp->Subcatch[j].groundwater->stats.finalUpperMoist = theta;
    sp->Subcatch[j].groundwater->stats.finalWaterTable = waterTable;
    if ( fabs(latFlow) > fabs(sp->Subcatch[j].groundwater->stats.maxFlow) )
    {
        sp->Subcatch[j].groundwater->stats.maxFlow = latFlow;
    }
}

//=============================================================================

void  stats_updateMaxRunoff(SWMM_Project *sp)
//
//   Input:   none
//   Output:  updates global variable MaxRunoffFlow
//   Purpose: updates value of maximum system runoff rate.
//
{
    int j;
    double sysRunoff = 0.0;
    
    TStatsExport *sttsx = &sp->StatsExport;

    for (j=0; j<sp->Nobjects[SUBCATCH]; j++)
        sysRunoff += sp->Subcatch[j].newRunoff;
    sttsx->MaxRunoffFlow = MAX(sttsx->MaxRunoffFlow, sysRunoff);
}    

//=============================================================================

////  New function added for release 5.1.008.  ////                            //(5.1.008)

void   stats_updateMaxNodeDepth(SWMM_Project *sp, int j, double depth)
//
//   Input:   j = node index
//            depth = water depth at node at current reporting time (ft)
//   Output:  none
//   Purpose: updates a node's maximum depth recorded at reporting times.
//
{
    TStatsExport *sttsx = &sp->StatsExport;

    if ( sttsx->NodeStats != NULL )
        sttsx->NodeStats[j].maxRptDepth = MAX(sttsx->NodeStats[j].maxRptDepth, depth);
}

//=============================================================================

void   stats_updateFlowStats(SWMM_Project *sp, double tStep, DateTime aDate,
        int stepCount, int steadyState)
//
//  Input:   tStep = routing time step (sec)
//           aDate = current date/time
//           stepCount = # steps required to solve routing at current time period
//           steadyState = TRUE if steady flow conditions exist
//  Output:  none
//  Purpose: updates various flow routing statistics at current time period.
//
{
    int   j;

    TStatsShared *stts = &sp->StatsShared;
    TStatsExport *sttsx = &sp->StatsExport;

    // --- update stats only after reporting period begins
    if ( aDate < sp->ReportStart ) return;
    stts->SysOutfallFlow = 0.0;

    // --- update node & link stats
#pragma omp parallel num_threads(sp->NumThreads)                                   //(5.1.008)
{
    #pragma omp for                                                            //(5.1.008)
    for ( j=0; j<sp->Nobjects[NODE]; j++ )
        stats_updateNodeStats(sp, j, tStep, aDate);
    #pragma omp for                                                            //(5.1.008)
    for ( j=0; j<sp->Nobjects[LINK]; j++ )
        stats_updateLinkStats(sp, j, tStep, aDate);
}

////  Following code segment modified for release 5.1.012.  ////               //(5.1.012)

    // --- update count of times in steady state
    stts->SysStats.steadyStateCount += steadyState;

    // --- update time step stats if not in steady state
	if ( steadyState == FALSE )
	{
        // --- skip initial time step for min. value)
        if ( sp->OldRoutingTime > 0 )
        {
            stts->SysStats.minTimeStep = MIN(stts->SysStats.minTimeStep, tStep);
        }
        stts->SysStats.avgTimeStep += tStep;
        stts->SysStats.maxTimeStep = MAX(stts->SysStats.maxTimeStep, tStep);

        // --- update iteration step count stats
        stts->SysStats.avgStepCount += stepCount;
	}

////

    // --- update max. system outfall flow
	sttsx->MaxOutfallFlow = MAX(sttsx->MaxOutfallFlow, stts->SysOutfallFlow);
}

//=============================================================================
   
void stats_updateCriticalTimeCount(SWMM_Project *sp, int node, int link)
//
//  Input:   node = node index
//           link = link index
//  Output:  none
//  Purpose: updates count of times a node or link was time step-critical.
//
{
    TStatsExport *sttsx = &sp->StatsExport;

    if      ( node >= 0 )
        sttsx->NodeStats[node].timeCourantCritical += 1.0;
    else if ( link >= 0 )
        sttsx->LinkStats[link].timeCourantCritical += 1.0;
}

//=============================================================================

////  Function modified for release 5.1.008.  ////                             //(5.1.008)

void stats_updateNodeStats(SWMM_Project *sp, int j, double tStep, DateTime aDate)
//
//  Input:   j = node index
//           tStep = routing time step (sec)
//           aDate = current date/time
//  Output:  none
//  Purpose: updates flow statistics for a node.
//
{
    int    k, p;
    double newVolume = sp->Node[j].newVolume;
    double newDepth = sp->Node[j].newDepth;
    int    canPond = (sp->AllowPonding && sp->Node[j].pondedArea > 0.0);

    TStatsShared *stts = &sp->StatsShared;
    TStatsExport *sttsx = &sp->StatsExport;

    // --- update depth statistics
    sttsx->NodeStats[j].avgDepth += newDepth;
    if ( newDepth > sttsx->NodeStats[j].maxDepth )
    {
        sttsx->NodeStats[j].maxDepth = newDepth;
        sttsx->NodeStats[j].maxDepthDate = aDate;
    }
    
    // --- update flooding, ponding, and surcharge statistics
    if ( sp->Node[j].type != OUTFALL )
    {
        if ( newVolume > sp->Node[j].fullVolume || sp->Node[j].overflow > 0.0 )
        {
            sttsx->NodeStats[j].timeFlooded += tStep;
            sttsx->NodeStats[j].volFlooded += sp->Node[j].overflow * tStep;
            if ( canPond ) sttsx->NodeStats[j].maxPondedVol =
                MAX(sttsx->NodeStats[j].maxPondedVol,
                    (newVolume - sp->Node[j].fullVolume));
        }

        // --- for dynamic wave routing, classify a non-storage node as        //(5.1.011)
        //     surcharged if its water level exceeds its crown elev.           //(5.1.011)
        if ( sp->RouteModel == DW && sp->Node[j].type != STORAGE &&                    //(5.1.011)
             newDepth + sp->Node[j].invertElev + FUDGE >= sp->Node[j].crownElev )
        {
            sttsx->NodeStats[j].timeSurcharged += tStep;
        }
    }

    // --- update storage statistics
    if ( sp->Node[j].type == STORAGE )
    {
        k = sp->Node[j].subIndex;
        sttsx->StorageStats[k].avgVol += newVolume;
        sttsx->StorageStats[k].evapLosses +=
            sp->Storage[sp->Node[j].subIndex].evapLoss; 
        sttsx->StorageStats[k].exfilLosses +=
            sp->Storage[sp->Node[j].subIndex].exfilLoss; 

        newVolume = MIN(newVolume, sp->Node[j].fullVolume);
        if ( newVolume > sttsx->StorageStats[k].maxVol )
        {
            sttsx->StorageStats[k].maxVol = newVolume;
            sttsx->StorageStats[k].maxVolDate = aDate;
        }
        sttsx->StorageStats[k].maxFlow = MAX(sttsx->StorageStats[k].maxFlow, sp->Node[j].outflow);
    }

    // --- update outfall statistics
    if ( sp->Node[j].type == OUTFALL ) 
    {
        k = sp->Node[j].subIndex;
        if ( sp->Node[j].inflow >= MIN_RUNOFF_FLOW )
        {
            sttsx->OutfallStats[k].avgFlow += sp->Node[j].inflow;
            sttsx->OutfallStats[k].maxFlow = MAX(sttsx->OutfallStats[k].maxFlow, sp->Node[j].inflow);
            sttsx->OutfallStats[k].totalPeriods++;
        }
        for (p=0; p<sp->Nobjects[POLLUT]; p++)
        {
            sttsx->OutfallStats[k].totalLoad[p] += sp->Node[j].inflow *
                sp->Node[j].newQual[p] * tStep;
        }
        stts->SysOutfallFlow += sp->Node[j].inflow;
    }

    // --- update inflow statistics
    sttsx->NodeStats[j].totLatFlow += ( (sp->Node[j].oldLatFlow + sp->Node[j].newLatFlow) *
                                 0.5 * tStep );
    if ( fabs(sp->Node[j].newLatFlow) > fabs(sttsx->NodeStats[j].maxLatFlow) )
        sttsx->NodeStats[j].maxLatFlow = sp->Node[j].newLatFlow;
    if ( sp->Node[j].inflow > sttsx->NodeStats[j].maxInflow )
    {
        sttsx->NodeStats[j].maxInflow = sp->Node[j].inflow;
        sttsx->NodeStats[j].maxInflowDate = aDate;
    }

    // --- update overflow statistics
    if ( sp->Node[j].overflow > sttsx->NodeStats[j].maxOverflow )
    {
        sttsx->NodeStats[j].maxOverflow = sp->Node[j].overflow;
        sttsx->NodeStats[j].maxOverflowDate = aDate;
    }
}

//=============================================================================

void  stats_updateLinkStats(SWMM_Project *sp, int j, double tStep, DateTime aDate)
//
//  Input:   j = link index
//           tStep = routing time step (sec)
//           aDate = current date/time
//  Output:  none
//  Purpose: updates flow statistics for a link.
//
{
    int    k;
    double q, v;
    double dq;

    TStatsExport *sttsx = &sp->StatsExport;

    // --- update max. flow
    dq = sp->Link[j].newFlow - sp->Link[j].oldFlow;
    q = fabs(sp->Link[j].newFlow);
    if ( q > sttsx->LinkStats[j].maxFlow )
    {
        sttsx->LinkStats[j].maxFlow = q;
        sttsx->LinkStats[j].maxFlowDate = aDate;
    }

    // --- update max. velocity
    v = link_getVelocity(sp, j, q, sp->Link[j].newDepth);
    if ( v > sttsx->LinkStats[j].maxVeloc )
    {
        sttsx->LinkStats[j].maxVeloc = v;
        //LinkStats[j].maxVelocDate = aDate;                                   //(5.1.008)
    }

    // --- update max. depth
    if ( sp->Link[j].newDepth > sttsx->LinkStats[j].maxDepth )
    {
        sttsx->LinkStats[j].maxDepth = sp->Link[j].newDepth;
    }

    if ( sp->Link[j].type == PUMP )
    {
        if ( q >= sp->Link[j].qFull )
            sttsx->LinkStats[j].timeFullFlow += tStep;
        if ( q > MIN_RUNOFF_FLOW )
        {
            k = sp->Link[j].subIndex;
            sttsx->PumpStats[k].minFlow = MIN(sttsx->PumpStats[k].minFlow, q);
            sttsx->PumpStats[k].maxFlow = sttsx->LinkStats[j].maxFlow;
            sttsx->PumpStats[k].avgFlow += q;
            sttsx->PumpStats[k].volume += q*tStep;
            sttsx->PumpStats[k].utilized += tStep;
            sttsx->PumpStats[k].energy += link_getPower(sp, j)*tStep/3600.0;
            if ( sp->Link[j].flowClass == DN_DRY )
                sttsx->PumpStats[k].offCurveLow += tStep;
            if ( sp->Link[j].flowClass == UP_DRY )
                sttsx->PumpStats[k].offCurveHigh += tStep;
            if ( sp->Link[j].oldFlow < MIN_RUNOFF_FLOW )
                sttsx->PumpStats[k].startUps++;
            sttsx->PumpStats[k].totalPeriods++;
            sttsx->LinkStats[j].timeSurcharged += tStep;
            sttsx->LinkStats[j].timeFullUpstream += tStep;
            sttsx->LinkStats[j].timeFullDnstream += tStep;
        }
    }
    else if ( sp->Link[j].type == CONDUIT )
    {

        // --- update time under normal flow & inlet control 
        if ( sp->Link[j].normalFlow ) sttsx->LinkStats[j].timeNormalFlow += tStep;
        if ( sp->Link[j].inletControl ) sttsx->LinkStats[j].timeInletControl += tStep;
    
        // --- update flow classification distribution
        k = sp->Link[j].flowClass;
        if ( k >= 0 && k < MAX_FLOW_CLASSES )
        {
            ++sttsx->LinkStats[j].timeInFlowClass[k];
        }

        // --- update time conduit is full
        k = sp->Link[j].subIndex;
        if ( q >= sp->Link[j].qFull * (double)sp->Conduit[k].barrels )                 //(5.1.012)
            sttsx->LinkStats[j].timeFullFlow += tStep;
        if ( sp->Conduit[k].capacityLimited )
            sttsx->LinkStats[j].timeCapacityLimited += tStep;

////  Following section modified for release 5.1.008.  ////                    //(5.1.008)
////
        switch (sp->Conduit[k].fullState)
        {
        case ALL_FULL:
            sttsx->LinkStats[j].timeSurcharged += tStep;
            sttsx->LinkStats[j].timeFullUpstream += tStep;
            sttsx->LinkStats[j].timeFullDnstream += tStep;
            break;
        case UP_FULL:
            sttsx->LinkStats[j].timeFullUpstream += tStep;
            break;
        case DN_FULL:
            sttsx->LinkStats[j].timeFullDnstream += tStep;
        }
////
    }

    // --- update flow turn count
    k = sttsx->LinkStats[j].flowTurnSign;
    sttsx->LinkStats[j].flowTurnSign = SGN(dq);
    if ( fabs(dq) > 0.001 &&  k * sttsx->LinkStats[j].flowTurnSign < 0 )
        sttsx->LinkStats[j].flowTurns++;
}

//=============================================================================

void  stats_findMaxStats(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: finds nodes & links with highest mass balance errors
//           & highest times Courant time-step critical.
//
{
    int    j;
    double x;

    TStatsShared *stts = &sp->StatsShared;
    TStatsExport *sttsx = &sp->StatsExport;
    TMassbalExport *mssblx = &sp->MassbalExport;

    // --- initialize max. stats arrays
    for (j=0; j<MAX_STATS; j++)
    {
        stts->MaxMassBalErrs[j].objType = NODE;
        stts->MaxMassBalErrs[j].index   = -1;
        stts->MaxMassBalErrs[j].value   = -1.0;
        stts->MaxCourantCrit[j].index   = -1;
        stts->MaxCourantCrit[j].value   = -1.0;
        stts->MaxFlowTurns[j].index     = -1;
        stts->MaxFlowTurns[j].value     = -1.0;
    }

    // --- find links with most flow turns 
    if ( sp->StepCount > 2 )
    {
        for (j=0; j<sp->Nobjects[LINK]; j++)
        {
            x = 100.0 * sttsx->LinkStats[j].flowTurns / (2./3.*(sp->StepCount-2));
            stats_updateMaxStats(stts->MaxFlowTurns, LINK, j, x);
        }
    }

    // --- find nodes with largest mass balance errors
    for (j=0; j<sp->Nobjects[NODE]; j++)
    {
        // --- skip terminal nodes and nodes with negligible inflow
        if ( sp->Node[j].degree <= 0  ) continue;
        if ( mssblx->NodeInflow[j] <= 0.1 ) continue;

        // --- evaluate mass balance error
        //     (Note: NodeInflow & NodeOutflow include any initial and final
        //            stored volumes, respectively).
        if ( mssblx->NodeInflow[j]  > 0.0 )
            x = 1.0 - mssblx->NodeOutflow[j] / mssblx->NodeInflow[j];
        else if ( mssblx->NodeOutflow[j] > 0.0 ) x = -1.0;
        else                             x = 0.0;
        stats_updateMaxStats(stts->MaxMassBalErrs, NODE, j, 100.0*x);
    }

    // --- stop if not using a variable time step
    if ( sp->RouteModel != DW || sp->CourantFactor == 0.0 ) return;

    // --- find nodes most frequently Courant critical
    if ( sp->StepCount == 0 ) return;                                              //(5.1.008)
    for (j=0; j<sp->Nobjects[NODE]; j++)
    {
        x = sttsx->NodeStats[j].timeCourantCritical / sp->StepCount;
        stats_updateMaxStats(stts->MaxCourantCrit, NODE, j, 100.0*x);
    }

    // --- find links most frequently Courant critical
    for (j=0; j<sp->Nobjects[LINK]; j++)
    {
        x = sttsx->LinkStats[j].timeCourantCritical / sp->StepCount;
        stats_updateMaxStats(stts->MaxCourantCrit, LINK, j, 100.0*x);
    }
}

//=============================================================================

void  stats_updateMaxStats(TMaxStats maxStats[], int i, int j, double x)
//
//  Input:   maxStats[] = array of critical statistics values
//           i = object category (NODE or LINK)
//           j = object index
//           x = value of statistic for the object
//  Output:  none
//  Purpose: updates the collection of most critical statistics
//
{
    int   k;
    TMaxStats maxStats1, maxStats2;
    maxStats1.objType = i;
    maxStats1.index   = j;
    maxStats1.value   = x;
    for (k=0; k<MAX_STATS; k++)
    {
        if ( fabs(maxStats1.value) > fabs(maxStats[k].value) )
        {
            maxStats2 = maxStats[k];
            maxStats[k] = maxStats1;
            maxStats1 = maxStats2;
        }
    }
}

//=============================================================================
//
int stats_getNodeStat(SWMM_Project *sp, int index, TNodeStats *nodeStats)
//
// Input:    index
//           element = element to return
// Return:   value
// Purpose:  Gets a Node Stat for toolkitAPI
//
{
	int errorcode = 0;

    TStatsExport *sttsx = &sp->StatsExport;

	// Check if Open
	if (swmm_IsOpenFlag(sp) == FALSE)
	{
		errorcode = ERR_API_INPUTNOTOPEN;
	}

	// Check if Simulation is Running
	else if (swmm_IsStartedFlag(sp) == FALSE)
	{
		errorcode = ERR_API_SIM_NRUNNING;
	}

	// Check if object index is within bounds
	else if (index < 0 || index >= sp->Nobjects[NODE])
	{
		errorcode = ERR_API_OBJECT_INDEX;
	}

	else
	{
		memcpy(nodeStats, &sttsx->NodeStats[index], sizeof(TNodeStats));
	}
	return errorcode;
}

int stats_getStorageStat(SWMM_Project *sp, int index, TStorageStats *storageStats)
//
// Input:    subindex
//           element = element to return
// Return:   value
// Purpose:  Gets a Storage Stat for toolkitAPI
//
{
	int errorcode = 0;

    TStatsExport *sttsx = &sp->StatsExport;

	// Check if Open
	if (swmm_IsOpenFlag(sp) == FALSE)
	{
		errorcode = ERR_API_INPUTNOTOPEN;
	}

	// Check if Simulation is Running
	else if (swmm_IsStartedFlag(sp) == FALSE)
	{
		errorcode = ERR_API_SIM_NRUNNING;
	}

	// Check if object index is within bounds
	else if (index < 0 || index >= sp->Nobjects[NODE])
	{
		errorcode = ERR_API_OBJECT_INDEX;
	}

	// Check Node Type is storage
	else if (sp->Node[index].type != STORAGE)
	{
		errorcode = ERR_API_WRONG_TYPE;
	}

	else
	{
		// fetch sub index
		int k = sp->Node[index].subIndex;
		// Copy Structure
		memcpy(storageStats, &sttsx->StorageStats[k], sizeof(TStorageStats));
	}
	return errorcode;
}

int stats_getOutfallStat(SWMM_Project *sp, int index, TOutfallStats *outfallStats)
//
// Input:    subindex
//           element = element to return
// Return:   value
// Purpose:  Gets a Outfall Stat for toolkitAPI
//
{
	int errorcode = 0;
    int p;

    TStatsExport *sttsx = &sp->StatsExport;

	// Check if Open
	if (swmm_IsOpenFlag(sp) == FALSE)
	{
		errorcode = ERR_API_INPUTNOTOPEN;
	}

	// Check if Simulation is Running
	else if (swmm_IsStartedFlag(sp) == FALSE)
	{
		errorcode = ERR_API_SIM_NRUNNING;
	}

	// Check if object index is within bounds
	else if (index < 0 || index >= sp->Nobjects[NODE])
	{
		errorcode = ERR_API_OBJECT_INDEX;
	}

	// Check Node Type is outfall
	else if (sp->Node[index].type != OUTFALL)
	{
		errorcode = ERR_API_WRONG_TYPE;
	}

	else
	{
		// fetch sub index
		int k = sp->Node[index].subIndex;
		// Copy Structure
		memcpy(outfallStats, &sttsx->OutfallStats[k], sizeof(TOutfallStats));
		
		// Perform Deep Copy of Pollutants Results
        if (sp->Nobjects[POLLUT] > 0)
        {
            outfallStats->totalLoad =
                (double *)calloc(sp->Nobjects[POLLUT], sizeof(double));
            if (!outfallStats->totalLoad)
            {
                errorcode = ERR_MEMORY;
            }
            if (errorcode == 0)
            {
                for (p = 0; p < sp->Nobjects[POLLUT]; p++)
                    outfallStats->totalLoad[p] = sttsx->OutfallStats[k].totalLoad[p];
            }
        }
        else outfallStats->totalLoad = NULL;
    }
    return errorcode;
}

int stats_getLinkStat(SWMM_Project *sp, int index, TLinkStats *linkStats)
//
// Input:    index
//           element = element to return
// Return:   value
// Purpose:  Gets a Link Stat for toolkitAPI
//
{
	int errorcode = 0;

    TStatsExport *sttsx = &sp->StatsExport;

	// Check if Open
	if (swmm_IsOpenFlag(sp) == FALSE)
	{
		errorcode = ERR_API_INPUTNOTOPEN;
	}

	// Check if Simulation is Running
	else if (swmm_IsStartedFlag(sp) == FALSE)
	{
		errorcode = ERR_API_SIM_NRUNNING;
	}

	// Check if object index is within bounds
	else if (index < 0 || index >= sp->Nobjects[LINK])
	{
		errorcode = ERR_API_OBJECT_INDEX;
	}

	else
	{
		// Copy Structure
		memcpy(linkStats, &sttsx->LinkStats[index], sizeof(TLinkStats));
	}
	return errorcode;
}

int stats_getPumpStat(SWMM_Project *sp, int index, TPumpStats *pumpStats)
//
// Input:    subindex
//           element = element to return
// Return:   value
// Purpose:  Gets a Pump Stat for toolkitAPI
//
{
	int errorcode = 0;

    TStatsExport *sttsx = &sp->StatsExport;

	// Check if Open
	if (swmm_IsOpenFlag(sp) == FALSE)
	{
		errorcode = ERR_API_INPUTNOTOPEN;
	}

	// Check if Simulation is Running
	else if (swmm_IsStartedFlag(sp) == FALSE)
	{
		errorcode = ERR_API_SIM_NRUNNING;
	}

	// Check if object index is within bounds
	else if (index < 0 || index >= sp->Nobjects[LINK])
	{
		errorcode = ERR_API_OBJECT_INDEX;
	}
	
	// Check if pump
	else if (sp->Link[index].type != PUMP)
	{
		errorcode = ERR_API_WRONG_TYPE;
	}

	else
	{
		// fetch sub index
		int k = sp->Link[index].subIndex;
		// Copy Structure
		memcpy(pumpStats, &sttsx->PumpStats[k], sizeof(TPumpStats));
	}
	return errorcode;
}

int stats_getSubcatchStat(SWMM_Project *sp ,int index, TSubcatchStats *subcatchStats)
//
// Input:    index
//           element = element to return
// Return:   value
// Purpose:  Gets a Subcatchment Stat for toolkitAPI
//
{
	int errorcode = 0;
    int p;

    TStatsExport *sttsx = &sp->StatsExport;

	// Check if Open
	if (swmm_IsOpenFlag(sp) == FALSE)
	{
		errorcode = ERR_API_INPUTNOTOPEN;
	}

	// Check if Simulation is Running
	else if (swmm_IsStartedFlag(sp) == FALSE)
	{
		errorcode = ERR_API_SIM_NRUNNING;
	}

	// Check if object index is within bounds
	else if (index < 0 || index >= sp->Nobjects[SUBCATCH])
	{
		errorcode = ERR_API_OBJECT_INDEX;
	}

	else
	{
		// Copy Structure
		memcpy(subcatchStats, &sttsx->SubcatchStats[index], sizeof(TSubcatchStats));
        
        // Perform Deep Copy of Pollutant Buildup Results
        if (sp->Nobjects[POLLUT] > 0)
        {
            subcatchStats->surfaceBuildup =
                (double *)calloc(sp->Nobjects[POLLUT], sizeof(double));
            if (!subcatchStats->surfaceBuildup)
            {
                errorcode = ERR_MEMORY;
            }
            if (errorcode == 0)
            {
                for (p = 0; p < sp->Nobjects[POLLUT]; p++)
                    subcatchStats->surfaceBuildup[p] = sttsx->SubcatchStats[index].surfaceBuildup[p];
            }
        }
        else subcatchStats->surfaceBuildup = NULL;
	}
	return errorcode;
}

