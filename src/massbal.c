//-----------------------------------------------------------------------------
//   massbal.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/14  (Build 5.1.001)
//             09/15/14  (Build 5.1.007)
//             04/02/15  (Build 5.1.008)
//             08/05/15  (Build 5.1.010)
//             08/01/16  (Build 5.1.011)
//             03/14/17  (Build 5.1.012)
//   Author:   L. Rossman (EPA)
//             M. Tryby (EPA)
//
//   Mass balance functions
//
//   Build 5.1.007:
//   - Mass balances modified to to correctly handle negative external inflows.
//   - Volume from minimum surface area at nodes included in mass balances.
//
//   Build 5.1.008:
//   - massbal_updateRunoffTotals() modified.
//   - LID drain flows and returned outfall flows added to components of
//     runoff mass balance.
//   - Seepage pollutant loss added into mass balances.
//
//   Build 5.1.010:
//   - Remaining pollutant mass in "dry" elements now added to final storage.
//
//   Build 5.1.011:
//   - Final stored pollutant mass in links ignored for Steady Flow routing.
//
//   Build 5.1.012:
//   - Terminal storage nodes no longer treated as non-storage terminal
//     nodes are when updating total outflow volume.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"
#include "swmm5.h"

//-----------------------------------------------------------------------------
//  Constants   
//-----------------------------------------------------------------------------
static const double MAX_RUNOFF_BALANCE_ERR = 10.0;
static const double MAX_FLOW_BALANCE_ERR   = 10.0;

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  massbal_open                (called from swmm_start in swmm5.c)
//  massbal_close               (called from swmm_end in swmm5.c)
//  massbal_report              (called from swmm_end in swmm5.c)
//  massbal_updateRunoffTotals  (called from subcatch_getRunoff)
//  massbal_updateDrainTotals   (called from evalLidUnit in lid.c)             //(5.1.008)
//  massbal_updateLoadingTotals (called from subcatch_getBuildup)
//  massbal_updateGwaterTotals  (called from updateMassBal in gwater.c)
//  massbal_updateRoutingTotals (called from routing_execute)
//  massbal_initTimeStepTotals  (called from routing_execute)
//  massbal_addInflowFlow       (called from routing.c)
//  massbal_addInflowQual       (called from routing.c)
//  massbal_addOutflowFlow      (called from removeOutflows in routing.c)
//  massbal_addOutflowQual      (called from removeOutflows in routing.c)
//  massbal_addNodeLosses       (called from removeStorageLosses in routing.c)
//  massbal_addLinkLosses       (called from removeConduitLosses in routing.c)
//  massbal_addReactedMass      (called from qualrout.c & treatmnt.c)
//  massbal_addSeepageLoss      (called from routing.c)                        //(5.1.008)
//  massbal_addToFinalStorage   (called from qualrout.c)                       //(5.1.008)
//  massbal_getStepFlowError    (called from routing.c)

//-----------------------------------------------------------------------------
//  Local Functions   
//-----------------------------------------------------------------------------
double massbal_getBuildup(SWMM_Project *sp, int pollut);
double massbal_getStorage(SWMM_Project *sp, char isFinalStorage);
double massbal_getStoredMass(SWMM_Project *sp, int pollut);
double massbal_getLoadingError(SWMM_Project *sp);
double massbal_getGwaterError(SWMM_Project *sp);
double massbal_getQualError(SWMM_Project *sp);


//=============================================================================

int massbal_open(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: opens and initializes mass balance continuity checking.
//
{
    int j, n;

    TMassbalShared *mssbl = &sp->MassbalShared;
    TMassbalExport *mssblx = &sp->MassbalExport;

    // --- initialize global continuity errors
    sp->RunoffError = 0.0;
    sp->GwaterError = 0.0;
    sp->FlowError   = 0.0;
    sp->QualError   = 0.0;

    // --- initialize runoff totals
    mssbl->RunoffTotals.rainfall    = 0.0;
    mssbl->RunoffTotals.evap        = 0.0;
    mssbl->RunoffTotals.infil       = 0.0;
    mssbl->RunoffTotals.runoff      = 0.0;
    mssbl->RunoffTotals.runon       = 0.0;                                            //(5.1.008)
    mssbl->RunoffTotals.drains      = 0.0;                                            //(5.1.008)
    mssbl->RunoffTotals.snowRemoved = 0.0;
    mssbl->RunoffTotals.initStorage = 0.0;
    mssbl->RunoffTotals.initSnowCover = 0.0;
    mssblx->TotalArea = 0.0;
    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        mssbl->RunoffTotals.initStorage += subcatch_getStorage(sp, j);
        mssbl->RunoffTotals.initSnowCover += snow_getSnowCover(sp, j);
        mssblx->TotalArea += sp->Subcatch[j].area;
    }

    // --- initialize groundwater totals
    mssbl->GwaterTotals.infil        = 0.0;
    mssbl->GwaterTotals.upperEvap    = 0.0;
    mssbl->GwaterTotals.lowerEvap    = 0.0;
    mssbl->GwaterTotals.lowerPerc    = 0.0;
    mssbl->GwaterTotals.gwater       = 0.0;
    mssbl->GwaterTotals.initStorage  = 0.0;
    mssbl->GwaterTotals.finalStorage = 0.0;
    for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
    {
        mssbl->GwaterTotals.initStorage += gwater_getVolume(sp, j) * sp->Subcatch[j].area;
    }

    // --- initialize node flow & storage totals
    mssbl->FlowTotals.dwInflow = 0.0;
    mssbl->FlowTotals.wwInflow = 0.0;
    mssbl->FlowTotals.gwInflow = 0.0;
    mssbl->FlowTotals.iiInflow = 0.0;
    mssbl->FlowTotals.exInflow = 0.0;
    mssbl->FlowTotals.flooding = 0.0;
    mssbl->FlowTotals.outflow  = 0.0;
    mssbl->FlowTotals.evapLoss = 0.0;
    mssbl->FlowTotals.seepLoss = 0.0;
    mssbl->FlowTotals.reacted  = 0.0;
    mssbl->FlowTotals.initStorage = 0.0;
    for (j = 0; j < sp->Nobjects[NODE]; j++)
        mssbl->FlowTotals.initStorage += sp->Node[j].newVolume;
    for (j = 0; j < sp->Nobjects[LINK]; j++)
        mssbl->FlowTotals.initStorage += sp->Link[j].newVolume;
    mssbl->StepFlowTotals = mssbl->FlowTotals;

    // --- add contribution of minimum surface area (i.e., manhole area)
    //     to initial storage under dynamic wave routing
    if ( sp->RouteModel == DW )
    {
        for (j = 0; j < sp->Nobjects[NODE]; j++)
	{
            if ( sp->Node[j].type != STORAGE &&
                sp->Node[j].initDepth <= sp->Node[j].crownElev - sp->Node[j].invertElev )  //(5.1.007)
                mssbl->FlowTotals.initStorage += sp->Node[j].initDepth * sp->MinSurfArea;
	}
    }

    // --- initialize arrays to null
    mssbl->LoadingTotals = NULL;
    mssbl->QualTotals = NULL;
    mssbl->StepQualTotals = NULL;
    mssblx->NodeInflow = NULL;
    mssblx->NodeOutflow = NULL;

    // --- allocate memory for WQ washoff continuity totals
    n = sp->Nobjects[POLLUT];
    if ( n > 0 )
    {
        mssbl->LoadingTotals = (TLoadingTotals *) calloc(n, sizeof(TLoadingTotals));
        if ( mssbl->LoadingTotals == NULL )
        {
             report_writeErrorMsg(sp, ERR_MEMORY, "");
             return sp->ErrorCode;
        }
        for (j = 0; j < n; j++)
        {
            mssbl->LoadingTotals[j].initLoad      = massbal_getBuildup(sp, j);
            mssbl->LoadingTotals[j].buildup       = 0.0;
            mssbl->LoadingTotals[j].deposition    = 0.0;
            mssbl->LoadingTotals[j].sweeping      = 0.0;
            mssbl->LoadingTotals[j].infil         = 0.0;
            mssbl->LoadingTotals[j].bmpRemoval    = 0.0;
            mssbl->LoadingTotals[j].runoff        = 0.0;
            mssbl->LoadingTotals[j].finalLoad     = 0.0;
        }
    }

    // --- allocate memory for nodal WQ continuity totals
    if ( n > 0 )
    {
        mssbl->QualTotals = (TRoutingTotals *) calloc(n, sizeof(TRoutingTotals));
        mssbl->StepQualTotals = (TRoutingTotals *) calloc(n, sizeof(TRoutingTotals));
         if ( mssbl->QualTotals == NULL || mssbl->StepQualTotals == NULL )
         {
             report_writeErrorMsg(sp, ERR_MEMORY, "");
             return sp->ErrorCode;
         }
     }

    // --- initialize WQ totals
    for (j = 0; j < n; j++)
    {
        mssbl->QualTotals[j].dwInflow = 0.0;
        mssbl->QualTotals[j].wwInflow = 0.0;
        mssbl->QualTotals[j].gwInflow = 0.0;
        mssbl->QualTotals[j].exInflow = 0.0;
        mssbl->QualTotals[j].flooding = 0.0;
        mssbl->QualTotals[j].outflow  = 0.0;
        mssbl->QualTotals[j].evapLoss = 0.0;
        mssbl->QualTotals[j].seepLoss = 0.0;
        mssbl->QualTotals[j].reacted  = 0.0;
        mssbl->QualTotals[j].initStorage = massbal_getStoredMass(sp, j);
    }

    // --- initialize totals used over a single time step
    massbal_initTimeStepTotals(sp);

    // --- allocate memory for nodal flow continuity
    if ( sp->Nobjects[NODE] > 0 )
    {
        mssblx->NodeInflow = (double *) calloc(sp->Nobjects[NODE], sizeof(double));
        if ( mssblx->NodeInflow == NULL )
        {
             report_writeErrorMsg(sp, ERR_MEMORY, "");
             return sp->ErrorCode;
        }
        mssblx->NodeOutflow = (double *) calloc(sp->Nobjects[NODE], sizeof(double));
        if ( mssblx->NodeOutflow == NULL )
        {
             report_writeErrorMsg(sp, ERR_MEMORY, "");
             return sp->ErrorCode;
        }
        for (j = 0; j < sp->Nobjects[NODE]; j++)
            mssblx->NodeInflow[j] = sp->Node[j].newVolume;
    }
    return sp->ErrorCode;
}

//=============================================================================

void massbal_close(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory used by mass balance system.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;
    TMassbalExport *mssblx = &sp->MassbalExport;

    FREE(mssbl->LoadingTotals);
    FREE(mssbl->QualTotals);
    FREE(mssbl->StepQualTotals);
    FREE(mssblx->NodeInflow);
    FREE(mssblx->NodeOutflow);
}

//=============================================================================

void massbal_report(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reports mass balance results.
//
{
    int    j;
    double gwArea = 0.0;

    TMassbalShared *mssbl = &sp->MassbalShared;
    TMassbalExport *mssblx = &sp->MassbalExport;

    if ( sp->Nobjects[SUBCATCH] > 0 )
    {
        if ( massbal_getRunoffError(sp) > MAX_RUNOFF_BALANCE_ERR ||
             sp->RptFlags.continuity == TRUE
           ) report_writeRunoffError(sp, &mssbl->RunoffTotals, mssblx->TotalArea);

        if ( sp->Nobjects[POLLUT] > 0 && !sp->IgnoreQuality )
        {
            if ( massbal_getLoadingError(sp) > MAX_RUNOFF_BALANCE_ERR ||
                 sp->RptFlags.continuity == TRUE
               ) report_writeLoadingError(sp, mssbl->LoadingTotals);
        }
    }

    if ( sp->Nobjects[AQUIFER] > 0  && !sp->IgnoreGwater )
    {
        if ( massbal_getGwaterError(sp) > MAX_RUNOFF_BALANCE_ERR ||
             sp->RptFlags.continuity == TRUE )
        {
            for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
            {
                if ( sp->Subcatch[j].groundwater ) gwArea += sp->Subcatch[j].area;
            }
            if ( gwArea > 0.0 ) report_writeGwaterError(sp, &mssbl->GwaterTotals, gwArea);
       }
    }

    if ( sp->Nobjects[NODE] > 0 && !sp->IgnoreRouting )
    {
        if ( massbal_getFlowError(sp) > MAX_FLOW_BALANCE_ERR ||
             sp->RptFlags.continuity == TRUE
           ) report_writeFlowError(sp, &mssbl->FlowTotals);
    
        if ( sp->Nobjects[POLLUT] > 0 && !sp->IgnoreQuality )
        {
            if ( massbal_getQualError(sp) > MAX_FLOW_BALANCE_ERR ||
                 sp->RptFlags.continuity == TRUE
               ) report_writeQualError(sp, mssbl->QualTotals);
        }
    }
}

//=============================================================================

double massbal_getBuildup(SWMM_Project *sp, int p)
//
//  Input:   p = pollutant index
//  Output:  returns total pollutant buildup (lbs or kg)
//  Purpose: computes current total buildup of a pollutant over study area.
//
{
    int    i, j;
    double load = 0.0;

    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        for (i = 0; i < sp->Nobjects[LANDUSE]; i++)
        {
            load += sp->Subcatch[j].landFactor[i].buildup[p];
        }
        load += sp->Subcatch[j].pondedQual[p] * sp->Pollut[p].mcf;
    }
    return load;
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void massbal_updateRunoffTotals(SWMM_Project *sp, int flowType, double v)
//
//  Input:   flowType = type of flow
//           v = flow volume (ft3)
//  Output:  none
//  Purpose: updates runoff totals after current time step.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    switch(flowType)
    {
    case RUNOFF_RAINFALL: mssbl->RunoffTotals.rainfall += v; break;
    case RUNOFF_EVAP:     mssbl->RunoffTotals.evap     += v; break;
    case RUNOFF_INFIL:    mssbl->RunoffTotals.infil    += v; break;
    case RUNOFF_RUNOFF:   mssbl->RunoffTotals.runoff   += v; break;
    case RUNOFF_DRAINS:   mssbl->RunoffTotals.drains   += v; break;
    case RUNOFF_RUNON:    mssbl->RunoffTotals.runon    += v; break;
    }
}

//=============================================================================

void massbal_updateGwaterTotals(SWMM_Project *sp, double vInfil,
        double vUpperEvap, double vLowerEvap, double vLowerPerc, double vGwater)
//
//  Input:   vInfil = volume depth of infiltrated water (ft)
//           vUpperEvap = volume depth of upper evaporation (ft)
//           vLowerEvap = volume depth of lower evaporation (ft)
//           vLowerPerc = volume depth of percolation to deep GW (ft)
//           vGwater = volume depth of groundwater outflow (ft)
//  Output:  none
//  Purpose: updates groundwater totals after current time step.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    mssbl->GwaterTotals.infil     += vInfil;
    mssbl->GwaterTotals.upperEvap += vUpperEvap;
    mssbl->GwaterTotals.lowerEvap += vLowerEvap;
    mssbl->GwaterTotals.lowerPerc += vLowerPerc;
    mssbl->GwaterTotals.gwater    += vGwater;
}

//=============================================================================

void massbal_initTimeStepTotals(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes routing totals for current time step.
//
{
    int j;

    TMassbalShared *mssbl = &sp->MassbalShared;

    mssbl->OldStepFlowTotals = mssbl->StepFlowTotals;
    mssbl->StepFlowTotals.dwInflow  = 0.0;
    mssbl->StepFlowTotals.wwInflow  = 0.0;
    mssbl->StepFlowTotals.gwInflow  = 0.0;
    mssbl->StepFlowTotals.iiInflow  = 0.0;
    mssbl->StepFlowTotals.exInflow  = 0.0;
    mssbl->StepFlowTotals.flooding  = 0.0;
    mssbl->StepFlowTotals.outflow   = 0.0;
    mssbl->StepFlowTotals.evapLoss  = 0.0;
    mssbl->StepFlowTotals.seepLoss  = 0.0;
    mssbl->StepFlowTotals.reacted   = 0.0;
    for (j=0; j<sp->Nobjects[POLLUT]; j++)
    {
        mssbl->StepQualTotals[j].dwInflow  = 0.0;
        mssbl->StepQualTotals[j].wwInflow  = 0.0;
        mssbl->StepQualTotals[j].gwInflow  = 0.0;
        mssbl->StepQualTotals[j].iiInflow  = 0.0;
        mssbl->StepQualTotals[j].exInflow  = 0.0;
        mssbl->StepQualTotals[j].flooding  = 0.0;
        mssbl->StepQualTotals[j].outflow   = 0.0;
        mssbl->StepQualTotals[j].reacted   = 0.0;
        mssbl->StepQualTotals[j].seepLoss  = 0.0;                                     //(5.1.008)
        mssbl->StepQualTotals[j].initStorage = 0.0;                                   //(5.1.010)
        mssbl->StepQualTotals[j].finalStorage = 0.0;                                  //(5.1.010)
    }
}

//=============================================================================

void massbal_addInflowFlow(SWMM_Project *sp, int type, double q)
//
//  Input:   type = type of inflow
//           q    = inflow rate (cfs)
//  Output:  none
//  Purpose: adds flow inflow to routing totals for current time step.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    switch (type)
    {
      case DRY_WEATHER_INFLOW: mssbl->StepFlowTotals.dwInflow += q; break;
      case WET_WEATHER_INFLOW: mssbl->StepFlowTotals.wwInflow += q; break;
      case GROUNDWATER_INFLOW: mssbl->StepFlowTotals.gwInflow += q; break;
      case RDII_INFLOW:        mssbl->StepFlowTotals.iiInflow += q; break;
      case EXTERNAL_INFLOW:    mssbl->StepFlowTotals.exInflow += q; break;
    }
}

//=============================================================================

void massbal_updateLoadingTotals(SWMM_Project *sp, int type, int p, double w)
//
//  Input:   type = type of inflow
//           p    = pollutant index
//           w    = mass loading
//  Output:  none
//  Purpose: adds inflow mass loading to loading totals for current time step.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    switch (type)
    {
      case BUILDUP_LOAD:     mssbl->LoadingTotals[p].buildup    += w; break;
      case DEPOSITION_LOAD:  mssbl->LoadingTotals[p].deposition += w; break;
      case SWEEPING_LOAD:    mssbl->LoadingTotals[p].sweeping   += w; break;
      case INFIL_LOAD:       mssbl->LoadingTotals[p].infil      += w; break;
      case BMP_REMOVAL_LOAD: mssbl->LoadingTotals[p].bmpRemoval += w; break;
      case RUNOFF_LOAD:      mssbl->LoadingTotals[p].runoff     += w; break;
      case FINAL_LOAD:       mssbl->LoadingTotals[p].finalLoad  += w; break;
    }
}

//=============================================================================

void massbal_addInflowQual(SWMM_Project *sp, int type, int p, double w)
//
//  Input:   type = type of inflow
//           p    = pollutant index
//           w    = mass flow rate (mass/sec)
//  Output:  none
//  Purpose: adds quality inflow to routing totals for current time step.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    if ( p < 0 || p >= sp->Nobjects[POLLUT] ) return;
    switch (type)
    {
      case DRY_WEATHER_INFLOW: mssbl->StepQualTotals[p].dwInflow += w; break;
      case WET_WEATHER_INFLOW: mssbl->StepQualTotals[p].wwInflow += w; break;
      case GROUNDWATER_INFLOW: mssbl->StepQualTotals[p].gwInflow += w; break;
      case EXTERNAL_INFLOW:    mssbl->StepQualTotals[p].exInflow += w; break;
      case RDII_INFLOW:        mssbl->StepQualTotals[p].iiInflow += w; break;
   }
}

//=============================================================================

////  This function was modified for release 5.1.007.  ////                    //(5.1.007)

void massbal_addOutflowFlow(SWMM_Project *sp, double q, int isFlooded)
//
//  Input:   q = outflow flow rate (cfs)
//           isFlooded = TRUE if outflow represents internal flooding
//  Output:  none
//  Purpose: adds flow outflow over current time step to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    if ( isFlooded ) mssbl->StepFlowTotals.flooding += q;
    else             mssbl->StepFlowTotals.outflow += q;
}

//=============================================================================

void massbal_addOutflowQual(SWMM_Project *sp, int p, double w, int isFlooded)
//
//  Input:   p = pollutant index
//           w = mass outflow rate (mass/sec)
//           isFlooded = TRUE if outflow represents internal flooding
//  Output:  none
//  Purpose: adds pollutant outflow over current time step to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    if ( p < 0 || p >= sp->Nobjects[POLLUT] ) return;
    if ( w >= 0.0 )
    {
        if ( isFlooded ) mssbl->StepQualTotals[p].flooding += w;
        else             mssbl->StepQualTotals[p].outflow += w;
    }
    else mssbl->StepQualTotals[p].exInflow -= w;
}

//=============================================================================

void massbal_addReactedMass(SWMM_Project *sp, int p, double w)
//
//  Input:   p = pollutant index
//           w = rate of mass reacted (mass/sec)
//  Output:  none
//  Purpose: adds mass reacted during current time step to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    if ( p < 0 || p >= sp->Nobjects[POLLUT] ) return;
    mssbl->StepQualTotals[p].reacted += w;
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

void massbal_addSeepageLoss(SWMM_Project *sp, int p, double w)
//
//  Input:   p = pollutant index
//           w = mass seepage rate (mass/sec)
//  Output:  none
//  Purpose: adds mass lost to seepage during current time step to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    if ( p < 0 || p >= sp->Nobjects[POLLUT] ) return;
    mssbl->StepQualTotals[p].seepLoss += w;
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

void massbal_addToFinalStorage(SWMM_Project *sp, int p, double w)
//
//  Input:   p = pollutant index
//           w = pollutant mass
//  Output:  none
//  Purpose: adds mass remaining on dry surface to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    if ( p < 0 || p >= sp->Nobjects[POLLUT] ) return;
    mssbl->StepQualTotals[p].finalStorage += w;
}

//=============================================================================

void massbal_addNodeLosses(SWMM_Project *sp, double evapLoss, double seepLoss)
//
//  Input:   evapLoss = evaporation loss from all nodes (ft3/sec)
//           seepLoss = seepage loss from all nodes (ft3/sec)
//  Output:  none
//  Purpose: adds node losses over current time step to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    mssbl->StepFlowTotals.evapLoss += evapLoss;
    mssbl->StepFlowTotals.seepLoss += seepLoss;
}

//=============================================================================

void massbal_addLinkLosses(SWMM_Project *sp, double evapLoss, double seepLoss)
//
//  Input:   evapLoss = evaporation loss from all links (ft3/sec)
//           infilLoss = infiltration loss from all links (ft3/sec)
//  Output:  none
//  Purpose: adds link losses over current time step to routing totals.
//
{
    TMassbalShared *mssbl = &sp->MassbalShared;

    mssbl->StepFlowTotals.evapLoss += evapLoss;
    mssbl->StepFlowTotals.seepLoss += seepLoss;
}

//=============================================================================

void massbal_updateRoutingTotals(SWMM_Project *sp, double tStep)
//
//  Input:   tStep = time step (sec)
//  Output:  none
//  Purpose: updates overall routing totals with totals from current time step.
//
{
    int j;

    TMassbalShared *mssbl = &sp->MassbalShared;
    TMassbalExport *mssblx = &sp->MassbalExport;

    mssbl->FlowTotals.dwInflow += mssbl->StepFlowTotals.dwInflow * tStep;
    mssbl->FlowTotals.wwInflow += mssbl->StepFlowTotals.wwInflow * tStep;
    mssbl->FlowTotals.gwInflow += mssbl->StepFlowTotals.gwInflow * tStep;
    mssbl->FlowTotals.iiInflow += mssbl->StepFlowTotals.iiInflow * tStep;
    mssbl->FlowTotals.exInflow += mssbl->StepFlowTotals.exInflow * tStep;
    mssbl->FlowTotals.flooding += mssbl->StepFlowTotals.flooding * tStep;
    mssbl->FlowTotals.outflow  += mssbl->StepFlowTotals.outflow * tStep;
    mssbl->FlowTotals.evapLoss += mssbl->StepFlowTotals.evapLoss * tStep;
    mssbl->FlowTotals.seepLoss += mssbl->StepFlowTotals.seepLoss * tStep;

    for (j = 0; j < sp->Nobjects[POLLUT]; j++)
    {
        mssbl->QualTotals[j].dwInflow += mssbl->StepQualTotals[j].dwInflow * tStep;
        mssbl->QualTotals[j].wwInflow += mssbl->StepQualTotals[j].wwInflow * tStep;
        mssbl->QualTotals[j].gwInflow += mssbl->StepQualTotals[j].gwInflow * tStep;
        mssbl->QualTotals[j].iiInflow += mssbl->StepQualTotals[j].iiInflow * tStep;
        mssbl->QualTotals[j].exInflow += mssbl->StepQualTotals[j].exInflow * tStep;
        mssbl->QualTotals[j].flooding += mssbl->StepQualTotals[j].flooding * tStep;
        mssbl->QualTotals[j].outflow  += mssbl->StepQualTotals[j].outflow * tStep;
        mssbl->QualTotals[j].reacted  += mssbl->StepQualTotals[j].reacted * tStep;
        mssbl->QualTotals[j].seepLoss += mssbl->StepQualTotals[j].seepLoss * tStep;          //(5.1.008)
        mssbl->QualTotals[j].finalStorage += mssbl->StepQualTotals[j].finalStorage;          //(5.1.010)
    }

    for ( j = 0; j < sp->Nobjects[NODE]; j++)
    {
        mssblx->NodeInflow[j] += sp->Node[j].inflow * tStep;
        if ( sp->Node[j].type == OUTFALL || 
            (sp->Node[j].degree == 0 && sp->Node[j].type != STORAGE) )                 //(5.1.012)
        {
            mssblx->NodeOutflow[j] += sp->Node[j].inflow * tStep;
        }
        else
        {
            mssblx->NodeOutflow[j] += sp->Node[j].outflow * tStep;
            if ( sp->Node[j].newVolume <= sp->Node[j].fullVolume ) 
                mssblx->NodeOutflow[j] += sp->Node[j].overflow * tStep;
        }
    }
}

//=============================================================================

double massbal_getStorage(SWMM_Project *sp, char isFinalStorage)
//
//  Input:   isFinalStorage = TRUE if at final time period
//  Output:  returns storage volume used (ft3)
//  Purpose: computes total system storage (nodes + links) filled
//
{
    int    j;
    double totalStorage = 0.0;
    double nodeStorage;

    TMassbalExport *mssblx = &sp->MassbalExport;

    // --- get volume in nodes
    for (j = 0; j < sp->Nobjects[NODE]; j++)
    {
        nodeStorage = sp->Node[j].newVolume;
        if ( isFinalStorage ) mssblx->NodeOutflow[j] += nodeStorage;
        totalStorage += nodeStorage;
    }

    // --- add contribution from minimum surface area (i.e., manhole diameter)
    //     to final storage under dynamic wave routing
    if ( isFinalStorage && sp->RouteModel == DW )
    {
        for (j = 0; j < sp->Nobjects[NODE]; j++)
        {
            if ( sp->Node[j].type != STORAGE &&
                 sp->Node[j].newDepth <= sp->Node[j].crownElev - sp->Node[j].invertElev )  //(5.1.007)
                totalStorage +=	sp->Node[j].newDepth * sp->MinSurfArea;
	}
    }

    // --- skip final link storage for Steady Flow routing 
    if ( isFinalStorage && sp->RouteModel == SF ) return totalStorage;

    // --- add on volume stored in links
    for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        totalStorage += sp->Link[j].newVolume;
    }
    return totalStorage;
}

//=============================================================================

void massbal_getSysFlows(SWMM_Project *sp, double f, double sysFlows[])
//
//  Input:   f = time weighting factor
//  Output:  sysFlows = array of total system flows
//  Purpose: retrieves time-weighted average of old and new system flows.
//
{
    double f1 = 1.0 - f;

    TMassbalShared *mssbl = &sp->MassbalShared;

    sysFlows[SYS_DWFLOW] = (f1 * mssbl->OldStepFlowTotals.dwInflow +
                             f * mssbl->StepFlowTotals.dwInflow) * UCF(sp, FLOW);
    sysFlows[SYS_GWFLOW] = (f1 * mssbl->OldStepFlowTotals.gwInflow +
                             f * mssbl->StepFlowTotals.gwInflow) * UCF(sp, FLOW);
    sysFlows[SYS_IIFLOW] = (f1 * mssbl->OldStepFlowTotals.iiInflow +
                             f * mssbl->StepFlowTotals.iiInflow) * UCF(sp, FLOW);
    sysFlows[SYS_EXFLOW] = (f1 * mssbl->OldStepFlowTotals.exInflow +
                             f * mssbl->StepFlowTotals.exInflow) * UCF(sp, FLOW);
    sysFlows[SYS_FLOODING] = (f1 * mssbl->OldStepFlowTotals.flooding +
                               f * mssbl->StepFlowTotals.flooding) * UCF(sp, FLOW);
    sysFlows[SYS_OUTFLOW] = (f1 * mssbl->OldStepFlowTotals.outflow +
                              f * mssbl->StepFlowTotals.outflow) * UCF(sp, FLOW);
    sysFlows[SYS_STORAGE] = (f1 * mssbl->OldStepFlowTotals.finalStorage +
                              f * mssbl->StepFlowTotals.finalStorage) * UCF(sp, VOLUME);
}

//=============================================================================

double massbal_getRunoffError(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: computes runoff mass balance error.
//
{
    int    j;
    double totalInflow;
    double totalOutflow;

    TMassbalShared *mssbl = &sp->MassbalShared;

    // --- find final storage on all subcatchments
    mssbl->RunoffTotals.finalStorage = 0.0;
    mssbl->RunoffTotals.finalSnowCover = 0.0;
    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        mssbl->RunoffTotals.finalStorage += subcatch_getStorage(sp, j);
        mssbl->RunoffTotals.finalSnowCover += snow_getSnowCover(sp, j);
    }

    // --- get snow removed from system
    mssbl->RunoffTotals.snowRemoved = sp->Snow.removed;

    // --- compute % difference between total inflow and outflow
    totalInflow  = mssbl->RunoffTotals.rainfall +
            mssbl->RunoffTotals.runon +                                        //(5.1.008)
            mssbl->RunoffTotals.initStorage +
            mssbl->RunoffTotals.initSnowCover;
    totalOutflow = mssbl->RunoffTotals.evap +
            mssbl->RunoffTotals.infil +
            mssbl->RunoffTotals.runoff +
            mssbl->RunoffTotals.drains +                                       //(5.1.008)
            mssbl->RunoffTotals.snowRemoved +
            mssbl->RunoffTotals.finalStorage +
            mssbl->RunoffTotals.finalSnowCover;
    mssbl->RunoffTotals.pctError = 0.0;
    if ( fabs(totalInflow - totalOutflow) < 1.0 )
    {
        mssbl->RunoffTotals.pctError = TINY;
    }
    else if ( totalInflow > 0.0 )
    {
        mssbl->RunoffTotals.pctError = 100.0 * (1.0 - totalOutflow / totalInflow);
    }
    else if ( totalOutflow > 0.0 )
    {
        mssbl->RunoffTotals.pctError = 100.0 * (totalInflow / totalOutflow - 1.0);
    }
    sp->RunoffError = mssbl->RunoffTotals.pctError;
    return mssbl->RunoffTotals.pctError;
}

//=============================================================================

double massbal_getLoadingError(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: computes runoff load mass balance error.
//
{
    int    j;
    double loadIn;
    double loadOut;
    double maxError = 0.0;

    TMassbalShared *mssbl = &sp->MassbalShared;

    for (j = 0; j < sp->Nobjects[POLLUT]; j++)
    {
        // --- get final pollutant loading remaining on land surface
        mssbl->LoadingTotals[j].finalLoad += massbal_getBuildup(sp, j);

        // --- compute total load added to study area
        loadIn = mssbl->LoadingTotals[j].initLoad +
                mssbl->LoadingTotals[j].buildup +
                mssbl->LoadingTotals[j].deposition;
    
        // --- compute total load removed from study area
        loadOut = mssbl->LoadingTotals[j].sweeping +
                mssbl->LoadingTotals[j].infil +
                mssbl->LoadingTotals[j].bmpRemoval +
                mssbl->LoadingTotals[j].runoff +
                mssbl->LoadingTotals[j].finalLoad;

        // --- compute mass balance error
        mssbl->LoadingTotals[j].pctError = 0.0;
        if ( fabs(loadIn - loadOut) < 0.001 )
        {
            mssbl->LoadingTotals[j].pctError = TINY;
        }
        else if ( loadIn > 0.0 )
        {
            mssbl->LoadingTotals[j].pctError = 100.0 * (1.0 - loadOut / loadIn);
        }
        else if ( loadOut > 0.0 )
        {
            mssbl->LoadingTotals[j].pctError = 100.0 * (loadIn / loadOut - 1.0);
        }
        maxError = MAX(maxError, mssbl->LoadingTotals[j].pctError);

        // --- report total counts as log10
        if ( sp->Pollut[j].units == COUNT )
        {
            mssbl->LoadingTotals[j].initLoad   = LOG10(mssbl->LoadingTotals[j].initLoad);
            mssbl->LoadingTotals[j].buildup    = LOG10(mssbl->LoadingTotals[j].buildup);
            mssbl->LoadingTotals[j].deposition = LOG10(mssbl->LoadingTotals[j].deposition);
            mssbl->LoadingTotals[j].sweeping   = LOG10(mssbl->LoadingTotals[j].sweeping);
            mssbl->LoadingTotals[j].infil      = LOG10(mssbl->LoadingTotals[j].infil);
            mssbl->LoadingTotals[j].bmpRemoval = LOG10(mssbl->LoadingTotals[j].bmpRemoval);
            mssbl->LoadingTotals[j].runoff     = LOG10(mssbl->LoadingTotals[j].runoff);
            mssbl->LoadingTotals[j].finalLoad  = LOG10(mssbl->LoadingTotals[j].finalLoad);
        }
    }
    return maxError;
}

//=============================================================================

double massbal_getGwaterError(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: computes groundwater mass balance error.
//
{
    int    j;
    double totalInflow;
    double totalOutflow;

    TMassbalShared *mssbl = &sp->MassbalShared;

    // --- find final storage in groundwater
    mssbl->GwaterTotals.finalStorage = 0.0;
    for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
    {
        mssbl->GwaterTotals.finalStorage += gwater_getVolume(sp, j) * sp->Subcatch[j].area;
    }

    // --- compute % difference between total inflow and outflow
    totalInflow  = mssbl->GwaterTotals.infil +
            mssbl->GwaterTotals.initStorage;
    totalOutflow = mssbl->GwaterTotals.upperEvap +
            mssbl->GwaterTotals.lowerEvap +
            mssbl->GwaterTotals.lowerPerc +
            mssbl->GwaterTotals.gwater +
            mssbl->GwaterTotals.finalStorage;
    mssbl->GwaterTotals.pctError = 0.0;
    if ( fabs(totalInflow - totalOutflow) < 1.0 )
    {
        mssbl->GwaterTotals.pctError = TINY;
    }
    else if ( totalInflow > 0.0 )
    {
        mssbl->GwaterTotals.pctError = 100.0 * (1.0 - totalOutflow / totalInflow);
    }
    else if ( totalOutflow > 0.0 )
    {
        mssbl->GwaterTotals.pctError = 100.0 * (totalInflow / totalOutflow - 1.0);
    }
    sp->GwaterError = mssbl->GwaterTotals.pctError;
    return mssbl->GwaterTotals.pctError;
}

//=============================================================================

////  The following function was re-written for release 5.1.008.  ////         //(5.1.008)

double massbal_getFlowError(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: computes flow routing mass balance error.
//
{
    double totalInflow;
    double totalOutflow;

    TMassbalShared *mssbl = &sp->MassbalShared;

    // --- get final volume of nodes and links
    mssbl->FlowTotals.finalStorage = massbal_getStorage(sp, TRUE);

    // --- add contributions to total inflow and outflow that are always positive
    totalInflow = mssbl->FlowTotals.initStorage + mssbl->FlowTotals.wwInflow +
            mssbl->FlowTotals.iiInflow;
    totalOutflow = mssbl->FlowTotals.finalStorage + mssbl->FlowTotals.flooding +
            mssbl->FlowTotals.evapLoss + mssbl->FlowTotals.seepLoss +
            mssbl->FlowTotals.reacted;

    // --- add on contributions that might be either positive or negative
    if ( mssbl->FlowTotals.dwInflow >= 0.0 )
        totalInflow += mssbl->FlowTotals.dwInflow;
    else
        totalOutflow -= mssbl->FlowTotals.dwInflow;

    if ( mssbl->FlowTotals.gwInflow >= 0.0 )
        totalInflow += mssbl->FlowTotals.gwInflow;
    else
        totalOutflow -= mssbl->FlowTotals.gwInflow;

    if ( mssbl->FlowTotals.exInflow >= 0.0 )
        totalInflow += mssbl->FlowTotals.exInflow;
    else
        totalOutflow -= mssbl->FlowTotals.exInflow;

    if ( mssbl->FlowTotals.outflow >= 0.0 )
        totalOutflow += mssbl->FlowTotals.outflow;
    else
        totalInflow -= mssbl->FlowTotals.outflow;

    // --- find percent difference between total inflow and outflow
    mssbl->FlowTotals.pctError = 0.0;
    if ( fabs(totalInflow - totalOutflow) < 1.0 )
    {
        mssbl->FlowTotals.pctError = TINY;
    }
    else if ( fabs(totalInflow) > 0.0 )
    {
        mssbl->FlowTotals.pctError = 100.0 * (1.0 - totalOutflow / totalInflow);
    }
    else if ( fabs(totalOutflow) > 0.0 )
    {
        mssbl->FlowTotals.pctError = 100.0 * (totalInflow / totalOutflow - 1.0);
    }
    sp->FlowError = mssbl->FlowTotals.pctError;
    return mssbl->FlowTotals.pctError;
}

//=============================================================================

double massbal_getQualError(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: computes water quality routing mass balance error.
//
{
    int    p;
    double maxQualError = 0.0;
    double totalInflow;
    double totalOutflow;
    double cf;

    TMassbalShared *mssbl = &sp->MassbalShared;

    // --- analyze each pollutant
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        // --- get final mass stored in nodes and links
        mssbl->QualTotals[p].finalStorage += massbal_getStoredMass(sp, p);                //(5.1.008)

        // --- compute % difference between total inflow and outflow
        totalInflow  = mssbl->QualTotals[p].dwInflow +
                mssbl->QualTotals[p].wwInflow +
                mssbl->QualTotals[p].gwInflow +
                mssbl->QualTotals[p].iiInflow +
                mssbl->QualTotals[p].exInflow +
                mssbl->QualTotals[p].initStorage;
        totalOutflow = mssbl->QualTotals[p].flooding +
                mssbl->QualTotals[p].outflow +
                mssbl->QualTotals[p].reacted +
                mssbl->QualTotals[p].seepLoss +                                //(5.1.008)
                mssbl->QualTotals[p].finalStorage;
        mssbl->QualTotals[p].pctError = 0.0;
        if ( fabs(totalInflow - totalOutflow) < 0.001 )
        {
            mssbl->QualTotals[p].pctError = TINY;
        }
        else if ( totalInflow > 0.0 )
        {
            mssbl->QualTotals[p].pctError = 100.0 * (1.0 - totalOutflow / totalInflow);
        }
        else if ( totalOutflow > 0.0 )
        {
            mssbl->QualTotals[p].pctError = 100.0 * (totalInflow / totalOutflow - 1.0);
        }

        // --- update max. error among all pollutants
        if ( fabs(mssbl->QualTotals[p].pctError) > fabs(maxQualError) )
        {
            maxQualError = mssbl->QualTotals[p].pctError;
        }

        // --- convert totals to reporting units (lbs, kg, or Log(Count))
        cf = LperFT3;
        if ( sp->Pollut[p].units == COUNT )
        {
            mssbl->QualTotals[p].dwInflow     = LOG10(cf * mssbl->QualTotals[p].dwInflow);
            mssbl->QualTotals[p].wwInflow     = LOG10(cf * mssbl->QualTotals[p].wwInflow);
            mssbl->QualTotals[p].gwInflow     = LOG10(cf * mssbl->QualTotals[p].gwInflow);
            mssbl->QualTotals[p].iiInflow     = LOG10(cf * mssbl->QualTotals[p].iiInflow);
            mssbl->QualTotals[p].exInflow     = LOG10(cf * mssbl->QualTotals[p].exInflow);
            mssbl->QualTotals[p].flooding     = LOG10(cf * mssbl->QualTotals[p].flooding);
            mssbl->QualTotals[p].outflow      = LOG10(cf * mssbl->QualTotals[p].outflow);
            mssbl->QualTotals[p].reacted      = LOG10(cf * mssbl->QualTotals[p].reacted);
            mssbl->QualTotals[p].seepLoss     = LOG10(cf * mssbl->QualTotals[p].seepLoss);   //(5.1.008)
            mssbl->QualTotals[p].initStorage  = LOG10(cf * mssbl->QualTotals[p].initStorage);
            mssbl->QualTotals[p].finalStorage = LOG10(cf * mssbl->QualTotals[p].finalStorage);
        }
        else
        {
            cf = cf * UCF(sp, MASS);
            if ( sp->Pollut[p].units == UG ) cf /= 1000.0;
            mssbl->QualTotals[p].dwInflow     *= cf;
            mssbl->QualTotals[p].wwInflow     *= cf;
            mssbl->QualTotals[p].gwInflow     *= cf;
            mssbl->QualTotals[p].iiInflow     *= cf;
            mssbl->QualTotals[p].exInflow     *= cf;
            mssbl->QualTotals[p].flooding     *= cf;
            mssbl->QualTotals[p].outflow      *= cf;
            mssbl->QualTotals[p].reacted      *= cf;
            mssbl->QualTotals[p].seepLoss     *= cf;
            mssbl->QualTotals[p].initStorage  *= cf;
            mssbl->QualTotals[p].finalStorage *= cf;
        }
    }
    sp->QualError = maxQualError;
    return maxQualError;
}
//=============================================================================

double massbal_getStepFlowError(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns fractional difference between total inflow and outflow.
//  Purpose: computes flow routing mass balance error at current time step.
//
{
    double totalInflow;
    double totalOutflow;

    TMassbalShared *mssbl = &sp->MassbalShared;

    // --- compute % difference between total inflow and outflow
    totalInflow  = mssbl->StepFlowTotals.dwInflow +
            mssbl->StepFlowTotals.wwInflow +
            mssbl->StepFlowTotals.gwInflow +
            mssbl->StepFlowTotals.iiInflow +
            mssbl->StepFlowTotals.exInflow;
    totalOutflow = mssbl->StepFlowTotals.flooding +
            mssbl->StepFlowTotals.outflow +
            mssbl->StepFlowTotals.evapLoss +
            mssbl->StepFlowTotals.seepLoss +
            mssbl->StepFlowTotals.reacted;
    if ( fabs(totalInflow) > 0.0 )                                             //(5.1.007)
        return 1.0 - totalOutflow / totalInflow;
    else if ( fabs(totalOutflow) > 0.0 )
        return totalInflow / totalOutflow - 1.0;                               //(5.1.007)
    else return 0.0;
}

//=============================================================================

double massbal_getStoredMass(SWMM_Project *sp, int p)
//
//  Input:   p = pollutant index
//  Output:  returns mass of pollutant.
//  Purpose: computes mass of pollutant stored in conveyance network.
//
{
    int j;
    double storedMass = 0.0;

    // --- get mass stored in nodes
    for (j = 0; j < sp->Nobjects[NODE]; j++)
        storedMass += sp->Node[j].newVolume * sp->Node[j].newQual[p];

    // --- get mass stored in links (except for Steady Flow routing)
    if ( sp->RouteModel != SF )                                                    //(5.1.011)
    {
        for (j = 0; j < sp->Nobjects[LINK]; j++)
            storedMass += sp->Link[j].newVolume * sp->Link[j].newQual[p];
    }
    return storedMass;
}

//=============================================================================

int massbal_getRoutingFlowTotal(SWMM_Project *sp, TRoutingTotals *RoutingTotal)
//
// Input:    element = element to return
// Return:   value
// Purpose:  Gets the routing total for toolkitAPI
//
{
	int errorcode = 0;

    TMassbalShared *mssbl = &sp->MassbalShared;

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

	else
	{
		memcpy(RoutingTotal, &mssbl->FlowTotals, sizeof(TRoutingTotals));
	}

	return errorcode;
}

int massbal_getRunoffTotal(SWMM_Project *sp, TRunoffTotals *runoffTot)
//
// Input:    element = element to return
// Return:   value
// Purpose:  Gets the runoff total for toolkitAPI
//
{
	int errorcode = 0;

    TMassbalShared *mssbl = &sp->MassbalShared;

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

	else
	{
		memcpy(runoffTot, &mssbl->RunoffTotals, sizeof(TRunoffTotals));
	}
	return errorcode;
}

double massbal_getTotalArea(SWMM_Project *sp)
//
// Return: Total Area for Runoff Surface
// Purpose: Used for Toolkit API Unit Conversion
{
    TMassbalExport *mssblx = &sp->MassbalExport;

	return mssblx->TotalArea;
}

int massbal_getNodeTotalInflow(SWMM_Project *sp, int index, double *value)
//
// Input:  NodeIndex
// Output: Volume
// Return: Error
// Purpose: Used for ToolkitAPI to pull total Node Inflow.
{
    int errorcode = 0;

    TMassbalExport *mssblx = &sp->MassbalExport;

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
    else
    {
		*value = mssblx->NodeInflow[index];
    }
    return errorcode;
}
