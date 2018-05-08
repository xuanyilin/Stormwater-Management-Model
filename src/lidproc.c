//-----------------------------------------------------------------------------
//   lidproc.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/12   (Build 5.1.001)
//             05/19/14   (Build 5.1.006)
//             09/15/14   (Build 5.1.007)
//             03/19/15   (Build 5.1.008)
//             04/30/15   (Build 5.1.009)
//             08/05/15   (Build 5.1.010)
//             08/01/16   (Build 5.1.011)
//             03/14/17   (Build 5.1.012)
//   Author:   L. Rossman (US EPA)
//
//   This module computes the hydrologic performance of an LID (Low Impact
//   Development) unit at a given point in time.
//
//   Build 5.1.007:
//   - Euler integration now applied to all LID types except Vegetative
//     Swale which continues to use successive approximation.
//   - LID layer flux routines were re-written to more accurately model
//     flooded conditions.
//
//   Build 5.1.008:
//   - MAX_STATE_VARS replaced with MAX_LAYERS.
//   - Optional soil layer added to Porous Pavement LID.
//   - Rooftop Disconnection added to types of LIDs.
//   - Separate accounting of drain flows added.
//   - Indicator for currently wet LIDs added.
//   - Detailed reporting procedure fixed.
//   - Possibile negative head on Bioretention Cell drain avoided.
//   - Bug in computing flow through Green Roof drainage mat fixed.
//
//   Build 5.1.009:
//   - Fixed typo in net flux rate for vegetative swale LID.
//
//   Build 5.1.010:
//   - New modified version of Green-Ampt used for surface layer infiltration.
//
//   Build 5.1.011:
//   - Re-named STOR_INFIL to STOR_EXFIL and StorageInfil to StorageExfil to
//     better reflect their meaning.
//   - Evaporation rates from sub-surface layers reduced by fraction of 
//     surface that is pervious (applies to block paver systems)
//   - Flux rate routines for LIDs with underdrains modified to produce more
//     physically meaningful results.
//   - Reporting of detailed results re-written.
//
//   Build 5.1.012:
//   - Modified upper limit for soil layer percolation.
//   - Modified upper limit on surface infiltration into rain gardens.
//   - Modified upper limit on drain flow for LIDs with storage layers.
//   - Used re-defined wasDry variable for LID reports to fix duplicate lines.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
#define STOPTOL  0.00328     // integration error tolerance in ft (= 1 mm)
#define MINFLOW  2.3e-8      // flow cutoff for dry conditions (= 0.001 in/hr)

//-----------------------------------------------------------------------------
//  Enumerations
//-----------------------------------------------------------------------------
enum LidLayerTypes {
    SURF,                    // surface layer
    SOIL,                    // soil layer
    STOR,                    // storage layer
    PAVE,                    // pavement layer
    DRAIN};                  // underdrain system

enum LidRptVars {
    SURF_INFLOW,             // inflow to surface layer
    TOTAL_EVAP,              // evaporation rate from all layers
    SURF_INFIL,              // infiltration into surface layer
    PAVE_PERC,               // percolation through pavement layer             //(5.1.008)
    SOIL_PERC,               // percolation through soil layer
    STOR_EXFIL,              // exfiltration out of storage layer              //(5.1.011)
    SURF_OUTFLOW,            // outflow from surface layer
    STOR_DRAIN,              // outflow from storage layer
    SURF_DEPTH,              // ponded depth on surface layer
    PAVE_DEPTH,              // water level in pavement layer                  //(5.1.011)
    SOIL_MOIST,              // moisture content of soil layer
    STOR_DEPTH,              // water level in storage layer
    MAX_RPT_VARS};

////  Added to release 5.1.008.  ////                                          //(5.1.008)
//-----------------------------------------------------------------------------
//  Imported variables 
//-----------------------------------------------------------------------------
//extern char HasWetLids;      // TRUE if any LIDs are wet (declared in runoff.c)

//-----------------------------------------------------------------------------
//  External Functions (declared in lid.h)
//-----------------------------------------------------------------------------
// lidproc_initWaterBalance  (called by lid_initState)
// lidproc_getOutflow        (called by evalLidUnit in lid.c)
// lidproc_saveResults       (called by evalLidUnit in lid.c)

//-----------------------------------------------------------------------------
// Local Functions
//-----------------------------------------------------------------------------
static void   barrelFluxRates(SWMM_Project *sp, double x[], double f[]);
static void   biocellFluxRates(SWMM_Project *sp, double x[], double f[]);
static void   greenRoofFluxRates(SWMM_Project *sp, double x[], double f[]);
static void   pavementFluxRates(SWMM_Project *sp, double x[], double f[]);
static void   trenchFluxRates(SWMM_Project *sp, double x[], double f[]);
static void   swaleFluxRates(SWMM_Project *sp, double x[], double f[]);
static void   roofFluxRates(SWMM_Project *sp, double x[], double f[]);                           //(5.1.008)

static double getSurfaceOutflowRate(SWMM_Project *sp, double depth);
static double getSurfaceOverflowRate(SWMM_Project *sp, double* surfaceDepth);
static double getPavementPermRate(SWMM_Project *sp);
static double getSoilPercRate(SWMM_Project *sp, double theta);                                   //(5.1.007)
static double getStorageExfilRate(SWMM_Project *sp);                                       //(5.1.011)
static double getStorageDrainRate(SWMM_Project *sp, double storageDepth, double soilTheta,       //(5.1.011)
              double paveDepth, double surfaceDepth);                          //(5.1.011)
static double getDrainMatOutflow(SWMM_Project *sp, double depth);
static void   getEvapRates(SWMM_Project *sp, double surfaceVol, double paveVol,                  //(5.1.008)
              double soilVol, double storageVol, double pervFrac);             //(5.1.011)

static void   updateWaterBalance(SWMM_Project *sp, TLidUnit *lidUnit,
        double inflow, double evap, double infil, double surfFlow,
        double drainFlow, double storage);

static int    modpuls_solve(SWMM_Project *sp, int n, double* x, double* xOld,
        double* xPrev, double* xMin, double* xMax, double* xTol, double* qOld,
        double* q, double dt, double omega,  //(5.1.007)
        void (*derivs)(SWMM_Project*, double*, double*));


//=============================================================================

void lidproc_initWaterBalance(TLidUnit *lidUnit, double initVol)
//
//  Purpose: initializes the water balance components of a LID unit.
//  Input:   lidUnit = a particular LID unit
//           initVol = initial water volume stored in the unit (ft)
//  Output:  none
//
{
    lidUnit->waterBalance.inflow = 0.0;
    lidUnit->waterBalance.evap = 0.0;
    lidUnit->waterBalance.infil = 0.0;
    lidUnit->waterBalance.surfFlow = 0.0;
    lidUnit->waterBalance.drainFlow = 0.0;
    lidUnit->waterBalance.initVol = initVol;
    lidUnit->waterBalance.finalVol = initVol;                                  //(5.1.008)
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double lidproc_getOutflow(SWMM_Project *sp, TLidUnit* lidUnit, TLidProc* lidProc,
        double inflow, double evap, double infil, double maxInfil, double tStep,
        double* lidEvap, double* lidInfil, double* lidDrain)
//
//  Purpose: computes runoff outflow from a single LID unit.
//  Input:   lidUnit  = ptr. to specific LID unit being analyzed
//           lidProc  = ptr. to generic LID process of the LID unit
//           inflow   = runoff rate captured by LID unit (ft/s)
//           evap     = potential evaporation rate (ft/s)
//           infil    = infiltration rate to native soil (ft/s)
//           maxInfil = max. infiltration rate to native soil (ft/s)
//           tStep    = time step (sec)
//  Output:  lidEvap  = evaporation rate for LID unit (ft/s)
//           lidInfil = infiltration rate for LID unit (ft/s)
//           lidDrain = drain flow for LID unit (ft/s)
//           returns surface runoff rate from the LID unit (ft/s)
//
{
    int    i;
    double x[MAX_LAYERS];        // layer moisture levels
    double xOld[MAX_LAYERS];     // work vector
    double xPrev[MAX_LAYERS];    // work vector
    double xMin[MAX_LAYERS];     // lower limit on moisture levels
    double xMax[MAX_LAYERS];     // upper limit on moisture levels
    double fOld[MAX_LAYERS];     // previously computed flux rates
    double f[MAX_LAYERS];        // newly computed flux rates

    // convergence tolerance on moisture levels (ft, moisture fraction , ft)
    double xTol[MAX_LAYERS] = {STOPTOL, STOPTOL, STOPTOL, STOPTOL};

    double omega = 0.0;          // integration time weighting

    //... define a pointer to function that computes flux rates through the LID
    void (*fluxRates) (SWMM_Project*, double *, double *) = NULL;

    TLidprocShared *ldprc = &sp->LidprocShared;

    //... save references to the LID process and LID unit
    ldprc->theLidProc = lidProc;
    ldprc->theLidUnit = lidUnit;

    //... save evap, max. infil. & time step to shared variables
    ldprc->EvapRate = evap;
    ldprc->MaxNativeInfil = maxInfil;
    ldprc->Tstep = tStep;

    //... store current moisture levels in vector x
    x[SURF] = ldprc->theLidUnit->surfaceDepth;
    x[SOIL] = ldprc->theLidUnit->soilMoisture;
    x[STOR] = ldprc->theLidUnit->storageDepth;
    x[PAVE] = ldprc->theLidUnit->paveDepth;

    //... initialize layer flux rates and moisture limits
    ldprc->SurfaceInflow  = inflow;
    ldprc->SurfaceInfil   = 0.0;
    ldprc->SurfaceEvap    = 0.0;
    ldprc->SurfaceOutflow = 0.0;
    ldprc->PaveEvap       = 0.0;
    ldprc->PavePerc       = 0.0;
    ldprc->SoilEvap       = 0.0;
    ldprc->SoilPerc       = 0.0;
    ldprc->StorageInflow  = 0.0;
    ldprc->StorageExfil   = 0.0;                                                      //(5.1.011)
    ldprc->StorageEvap    = 0.0;
    ldprc->StorageDrain   = 0.0;
    for (i = 0; i < MAX_LAYERS; i++)
    {
        f[i] = 0.0;
        fOld[i] = ldprc->theLidUnit->oldFluxRates[i];
        xMin[i] = 0.0;
        xMax[i] = BIG;
        ldprc->Xold[i] = x[i];
    }

    //... find Green-Ampt infiltration from surface layer
    if ( ldprc->theLidProc->lidType == POROUS_PAVEMENT )
        ldprc->SurfaceInfil = 0.0;
    else if ( ldprc->theLidUnit->soilInfil.Ks > 0.0 )
    {
        ldprc->SurfaceInfil =
            grnampt_getInfil(sp, &ldprc->theLidUnit->soilInfil, ldprc->Tstep,
                    ldprc->SurfaceInflow, ldprc->theLidUnit->surfaceDepth,
                             MOD_GREEN_AMPT);                                  //(5.1.010)
    }
    else ldprc->SurfaceInfil = infil;

    //... set moisture limits for soil & storage layers
    if ( ldprc->theLidProc->soil.thickness > 0.0 )
    {
        xMin[SOIL] = ldprc->theLidProc->soil.wiltPoint;
        xMax[SOIL] = ldprc->theLidProc->soil.porosity;
    }
    if ( ldprc->theLidProc->pavement.thickness > 0.0 )
    {
        xMax[PAVE] = ldprc->theLidProc->pavement.thickness;                           //(5.1.011)
    }
    if ( ldprc->theLidProc->storage.thickness > 0.0 )
    {
        xMax[STOR] = ldprc->theLidProc->storage.thickness;
    }
    if ( ldprc->theLidProc->lidType == GREEN_ROOF )
    {
        xMax[STOR] = ldprc->theLidProc->drainMat.thickness;
    }

    //... determine which flux rate function to use
    switch (ldprc->theLidProc->lidType)
    {
    case BIO_CELL:
    case RAIN_GARDEN:     fluxRates = &biocellFluxRates;   break;
    case GREEN_ROOF:      fluxRates = &greenRoofFluxRates; break;
    case INFIL_TRENCH:    fluxRates = &trenchFluxRates;    break;
    case POROUS_PAVEMENT: fluxRates = &pavementFluxRates;  break;
    case RAIN_BARREL:     fluxRates = &barrelFluxRates;    break;
    case ROOF_DISCON:     fluxRates = &roofFluxRates;      break;
    case VEG_SWALE:       fluxRates = &swaleFluxRates;
                          omega = 0.5;
                          break;
    default:              return 0.0;
    }

    //... update moisture levels and flux rates over the time step
    i = modpuls_solve(sp, MAX_LAYERS, x, xOld, xPrev, xMin, xMax, xTol,
                     fOld, f, tStep, omega, fluxRates);

/** For debugging only ********************************************
    if  (i == 0)
    {
        fprintf(Frpt.file,
        "\n  WARNING 09: integration failed to converge at %s %s",
            theDate, theTime);
        fprintf(Frpt.file,
        "\n              for LID %s placed in subcatchment %s.",
            theLidProc->ID, theSubcatch->ID);
    }
*******************************************************************/

    //... add any surface overflow to surface outflow
    if ( ldprc->theLidProc->surface.canOverflow || ldprc->theLidUnit->fullWidth == 0.0 )
    {
        ldprc->SurfaceOutflow += getSurfaceOverflowRate(sp, &x[SURF]);
    }

    //... save updated results
    ldprc->theLidUnit->surfaceDepth = x[SURF];
    ldprc->theLidUnit->paveDepth    = x[PAVE];                                        //(5.1.011)
    ldprc->theLidUnit->soilMoisture = x[SOIL];
    ldprc->theLidUnit->storageDepth = x[STOR];
    for (i = 0; i < MAX_LAYERS; i++) ldprc->theLidUnit->oldFluxRates[i] = f[i];

    //... assign values to LID unit evaporation, infiltration & drain flow
    *lidEvap = ldprc->SurfaceEvap + ldprc->PaveEvap + ldprc->SoilEvap + ldprc->StorageEvap;
    *lidInfil = ldprc->StorageExfil;
    *lidDrain = ldprc->StorageDrain;

    //... return surface outflow (per unit area) from unit
    return ldprc->SurfaceOutflow;
}

//=============================================================================

////  This function was re-written for release 5.1.011.  ////                  //(5.1.011)

void lidproc_saveResults(SWMM_Project *sp, TLidUnit* lidUnit, double ucfRainfall,
        double ucfRainDepth)
//
//  Purpose: updates the mass balance for an LID unit and saves
//           current flux rates to the LID report file.
//  Input:   lidUnit = ptr. to LID unit
//           ucfRainfall = units conversion factor for rainfall rate
//           ucfDepth = units conversion factor for rainfall depth
//  Output:  none
//
{
    double ucf;                        // units conversion factor
    double totalEvap;                  // total evaporation rate (ft/s)
    double totalVolume;                // total volume stored in LID (ft)
    double rptVars[MAX_RPT_VARS];      // array of reporting variables
    int    isDry = FALSE;              // true if current state of LID is dry
    char   timeStamp[24];              // date/time stamp
    double elapsedHrs;                 // elapsed hours

    TLidprocShared *ldprc = &sp->LidprocShared;
    TRunoffShared *rnff = &sp->RunoffShared;

    //... find total evap. rate and stored volume
    totalEvap = ldprc->SurfaceEvap + ldprc->PaveEvap + ldprc->SoilEvap +
            ldprc->StorageEvap;
    totalVolume = ldprc->SurfaceVolume + ldprc->PaveVolume + ldprc->SoilVolume +
            ldprc->StorageVolume;

    //... update mass balance totals
    updateWaterBalance(sp, ldprc->theLidUnit, ldprc->SurfaceInflow, totalEvap,
            ldprc->StorageExfil, ldprc->SurfaceOutflow, ldprc->StorageDrain,
            totalVolume);

    //... check if dry-weather conditions hold
    if ( ldprc->SurfaceInflow  < MINFLOW &&
            ldprc->SurfaceOutflow < MINFLOW &&
            ldprc->StorageDrain   < MINFLOW &&
            ldprc->StorageExfil   < MINFLOW &&
		    totalEvap      < MINFLOW
       ) isDry = TRUE;

    //... update status of HasWetLids
    if ( !isDry ) rnff->HasWetLids = TRUE;

    //... write results to LID report file                                     //(5.1.012)
    if ( lidUnit->rptFile )                                                    //(5.1.012)
    {
        //... convert rate results to original units (in/hr or mm/hr)
        ucf = ucfRainfall;
        rptVars[SURF_INFLOW]  = ldprc->SurfaceInflow*ucf;
        rptVars[TOTAL_EVAP]   = totalEvap*ucf;
        rptVars[SURF_INFIL]   = ldprc->SurfaceInfil*ucf;
        rptVars[PAVE_PERC]    = ldprc->PavePerc*ucf;
        rptVars[SOIL_PERC]    = ldprc->SoilPerc*ucf;
        rptVars[STOR_EXFIL]   = ldprc->StorageExfil*ucf;
        rptVars[SURF_OUTFLOW] = ldprc->SurfaceOutflow*ucf;
        rptVars[STOR_DRAIN]   = ldprc->StorageDrain*ucf;

        //... convert storage results to original units (in or mm)
        ucf = ucfRainDepth;
        rptVars[SURF_DEPTH] = ldprc->theLidUnit->surfaceDepth*ucf;
        rptVars[PAVE_DEPTH] = ldprc->theLidUnit->paveDepth;                           //(5.1.011)
        rptVars[SOIL_MOIST] = ldprc->theLidUnit->soilMoisture;
        rptVars[STOR_DEPTH] = ldprc->theLidUnit->storageDepth*ucf;

        //... if the current LID state is wet but the previous state was dry
        //    for more than one period then write the saved previous results   //(5.1.012)
        //    to the report file thus marking the end of a dry period          //(5.10012)
        if ( !isDry && ldprc->theLidUnit->rptFile->wasDry > 1)                        //(5.1.012)
        {
            fprintf(ldprc->theLidUnit->rptFile->file, "%s",
                    ldprc->theLidUnit->rptFile->results);
        }

        //... write the current results to a string which is saved between
        //    reporting periods
        elapsedHrs = sp->NewRunoffTime / 1000.0 / 3600.0;
        datetime_getTimeStamp(sp, M_D_Y, getDateTime(sp, sp->NewRunoffTime),
                24, timeStamp);
        sprintf(ldprc->theLidUnit->rptFile->results,
             "\n%20s\t %8.3f\t %8.3f\t %8.4f\t %8.3f\t %8.3f\t %8.3f\t %8.3f\t"
             "%8.3f\t %8.3f\t %8.3f\t %8.3f\t %8.3f\t %8.3f",
             timeStamp, elapsedHrs, rptVars[0], rptVars[1], rptVars[2],
             rptVars[3], rptVars[4], rptVars[5], rptVars[6], rptVars[7],
             rptVars[8], rptVars[9], rptVars[10], rptVars[11]);

        //... if the current LID state is dry
        if ( isDry )
        {
            //... if the previous state was wet then write the current
            //    results to file marking the start of a dry period
            if ( ldprc->theLidUnit->rptFile->wasDry == 0 )                            //(5.1.012)
            {
                fprintf(ldprc->theLidUnit->rptFile->file, "%s",
                        ldprc->theLidUnit->rptFile->results);
            }

            //... increment the number of successive dry periods               //(5.1.012)
            ldprc->theLidUnit->rptFile->wasDry++;                                     //(5.1.012)
        }

        //... if the current LID state is wet
        else
        {
            //... write the current results to the report file
			fprintf(ldprc->theLidUnit->rptFile->file, "%s",
			        ldprc->theLidUnit->rptFile->results);

            //... re-set the number of successive dry periods to 0             //(5.1.012)
			ldprc->theLidUnit->rptFile->wasDry = 0;                                   //(5.1.012)
        }
    }
}

//=============================================================================

////  New function for release 5.1.008.  ////                                  //(5.1.008)

void roofFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates for roof disconnection.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double surfaceDepth = x[SURF];

    TLidprocShared *ldprc = &sp->LidprocShared;

    getEvapRates(sp, surfaceDepth, 0.0, 0.0, 0.0, 1.0);                            //(5.1.011)
    ldprc->SurfaceVolume = surfaceDepth;
    ldprc->SurfaceInfil = 0.0;
    if ( ldprc->theLidProc->surface.alpha > 0.0 )
        ldprc->SurfaceOutflow = getSurfaceOutflowRate(sp, surfaceDepth);
    else getSurfaceOverflowRate(sp, &surfaceDepth);
    ldprc->StorageDrain = MIN(ldprc->theLidProc->drain.coeff/UCF(sp, RAINFALL),
            ldprc->SurfaceOutflow);
    ldprc->SurfaceOutflow -= ldprc->StorageDrain;
    f[SURF] = (ldprc->SurfaceInflow - ldprc->SurfaceEvap - ldprc->StorageDrain -
            ldprc->SurfaceOutflow);
}

//=============================================================================

////  This function was re-written for release 5.1.011.  ////                  //(5.1.011)

void greenRoofFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates from the layers of a green roof.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // Moisture level variables
    double surfaceDepth;
    double soilTheta;
    double storageDepth;

    // Intermediate variables
    double availVolume;
    double maxRate;

    TLidprocShared *ldprc = &sp->LidprocShared;

    // Green roof properties
    double soilThickness    = ldprc->theLidProc->soil.thickness;
    double storageThickness = ldprc->theLidProc->storage.thickness;
    double soilPorosity     = ldprc->theLidProc->soil.porosity;
    double storageVoidFrac  = ldprc->theLidProc->storage.voidFrac;
    double soilFieldCap     = ldprc->theLidProc->soil.fieldCap;
    double soilWiltPoint    = ldprc->theLidProc->soil.wiltPoint;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    soilTheta    = x[SOIL];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    ldprc->SurfaceVolume = surfaceDepth * ldprc->theLidProc->surface.voidFrac;
    ldprc->SoilVolume = soilTheta * soilThickness;
    ldprc->StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = ldprc->SoilVolume - soilWiltPoint * soilThickness;
    getEvapRates(sp, ldprc->SurfaceVolume, 0.0, availVolume, ldprc->StorageVolume, 1.0);
    if ( soilTheta >= soilPorosity ) ldprc->StorageEvap = 0.0;

    //... soil layer perc rate
    ldprc->SoilPerc = getSoilPercRate(sp, soilTheta);

    //... limit perc rate by available water
    availVolume = (soilTheta - soilFieldCap) * soilThickness;
    maxRate = MAX(availVolume, 0.0) / ldprc->Tstep - ldprc->SoilEvap;                        //(5.1.012)
    ldprc->SoilPerc = MIN(ldprc->SoilPerc, maxRate);
    ldprc->SoilPerc = MAX(ldprc->SoilPerc, 0.0);

    //... storage (drain mat) outflow rate
    ldprc->StorageExfil = 0.0;
    ldprc->StorageDrain = getDrainMatOutflow(sp, storageDepth);

    //... unit is full
    if ( soilTheta >= soilPorosity && storageDepth >= storageThickness )
    {
        //... outflow from both layers equals limiting rate
        maxRate = MIN(ldprc->SoilPerc, ldprc->StorageDrain);
        ldprc->SoilPerc = maxRate;
        ldprc->StorageDrain = maxRate;

        //... adjust inflow rate to soil layer
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
    }

    //... unit not full
    else
    {
        //... limit drainmat outflow by available storage volume
        maxRate = storageDepth * storageVoidFrac / ldprc->Tstep - ldprc->StorageEvap;        //(5.1.012)
        if ( storageDepth >= storageThickness ) maxRate += ldprc->SoilPerc;           //(5.1.012)
        maxRate = MAX(maxRate, 0.0);
        ldprc->StorageDrain = MIN(ldprc->StorageDrain, maxRate);

        //... limit soil perc inflow by unused storage volume
        maxRate = (storageThickness - storageDepth) * storageVoidFrac / ldprc->Tstep +
                ldprc->StorageDrain + ldprc->StorageEvap;
        ldprc->SoilPerc = MIN(ldprc->SoilPerc, maxRate);
                
        //... adjust surface infil. so soil porosity not exceeded
        maxRate = (soilPorosity - soilTheta) * soilThickness / ldprc->Tstep +
                ldprc->SoilPerc + ldprc->SoilEvap;
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
    }

    // ... find surface outflow rate
    ldprc->SurfaceOutflow = getSurfaceOutflowRate(sp, surfaceDepth);

    // ... compute overall layer flux rates
    f[SURF] = (ldprc->SurfaceInflow - ldprc->SurfaceEvap - ldprc->SurfaceInfil -
            ldprc->SurfaceOutflow) / ldprc->theLidProc->surface.voidFrac;
    f[SOIL] = (ldprc->SurfaceInfil - ldprc->SoilEvap - ldprc->SoilPerc) /
            ldprc->theLidProc->soil.thickness;
    f[STOR] = (ldprc->SoilPerc - ldprc->StorageEvap - ldprc->StorageDrain) /
            ldprc->theLidProc->storage.voidFrac;
}

//=============================================================================

////  This function was re-written for release 5.1.011.  ////                  //(5.1.011)

void biocellFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates from the layers of a bio-retention cell LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // Moisture level variables
    double surfaceDepth;
    double soilTheta;
    double storageDepth;

    // Intermediate variables
    double availVolume;
    double maxRate;

    TLidprocShared *ldprc = &sp->LidprocShared;

    // LID layer properties
    double soilThickness    = ldprc->theLidProc->soil.thickness;
    double soilPorosity     = ldprc->theLidProc->soil.porosity;
    double soilFieldCap     = ldprc->theLidProc->soil.fieldCap;
    double soilWiltPoint    = ldprc->theLidProc->soil.wiltPoint;
    double storageThickness = ldprc->theLidProc->storage.thickness;
    double storageVoidFrac  = ldprc->theLidProc->storage.voidFrac;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    soilTheta    = x[SOIL];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    ldprc->SurfaceVolume = surfaceDepth * ldprc->theLidProc->surface.voidFrac;
    ldprc->SoilVolume    = soilTheta * soilThickness;
    ldprc->StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = ldprc->SoilVolume - soilWiltPoint * soilThickness;
    getEvapRates(sp, ldprc->SurfaceVolume, 0.0, availVolume, ldprc->StorageVolume, 1.0);
    if ( soilTheta >= soilPorosity ) ldprc->StorageEvap = 0.0;

    //... soil layer perc rate
    ldprc->SoilPerc = getSoilPercRate(sp, soilTheta);

    //... limit perc rate by available water
    availVolume =  (soilTheta - soilFieldCap) * soilThickness;
    maxRate = MAX(availVolume, 0.0) / ldprc->Tstep - ldprc->SoilEvap;                        //(5.1.012)
    ldprc->SoilPerc = MIN(ldprc->SoilPerc, maxRate);
    ldprc->SoilPerc = MAX(ldprc->SoilPerc, 0.0);

    //... exfiltration rate out of storage layer
    ldprc->StorageExfil = getStorageExfilRate(sp);

    //... underdrain flow rate
    ldprc->StorageDrain = 0.0;
    if ( ldprc->theLidProc->drain.coeff > 0.0 )
    {
        ldprc->StorageDrain = getStorageDrainRate(sp, storageDepth, soilTheta,
                0.0, surfaceDepth);
    }

    //... special case of no storage layer present
    if ( storageThickness == 0.0 )
    {
        ldprc->StorageEvap = 0.0;
        maxRate = MIN(ldprc->SoilPerc, ldprc->StorageExfil);
        ldprc->SoilPerc = maxRate;
        ldprc->StorageExfil = maxRate;

////  Following code segment added to release 5.1.012  ////                    //(5.1.012)
        //... limit surface infil. by unused soil volume
        maxRate = (soilPorosity - soilTheta) * soilThickness / ldprc->Tstep +
                ldprc->SoilPerc + ldprc->SoilEvap;
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
//////////////////////////////////////////////////////////

	}

    //... storage & soil layers are full
    else if ( soilTheta >= soilPorosity && storageDepth >= storageThickness )
    {
        //... limiting rate is smaller of soil perc and storage outflow
        maxRate = ldprc->StorageExfil + ldprc->StorageDrain;
        if ( ldprc->SoilPerc < maxRate )
        {
            maxRate = ldprc->SoilPerc;
            if ( maxRate > ldprc->StorageExfil )
                ldprc->StorageDrain = maxRate - ldprc->StorageExfil;
            else
            {
                ldprc->StorageExfil = maxRate;
                ldprc->StorageDrain = 0.0;
            }
        }
        else ldprc->SoilPerc = maxRate;

        //... apply limiting rate to surface infil.
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
    }

    //... either layer not full
    else if ( storageThickness > 0.0 )
    {
        //... limit storage exfiltration by available storage volume
        maxRate = ldprc->SoilPerc - ldprc->StorageEvap +
                storageDepth*storageVoidFrac/ldprc->Tstep;
        ldprc->StorageExfil = MIN(ldprc->StorageExfil, maxRate);
        ldprc->StorageExfil = MAX(ldprc->StorageExfil, 0.0);

        //... limit underdrain flow by volume above drain offset
        if ( ldprc->StorageDrain > 0.0 )
        {
            maxRate = -ldprc->StorageExfil - ldprc->StorageEvap;                              //(5.1.012)
            if ( storageDepth >= storageThickness) maxRate += ldprc->SoilPerc;         //(5.1.012)
            if ( ldprc->theLidProc->drain.offset <= storageDepth )
            {
                maxRate += (storageDepth - ldprc->theLidProc->drain.offset) *
                           storageVoidFrac/ldprc->Tstep;
            }
            maxRate = MAX(maxRate, 0.0);
            ldprc->StorageDrain = MIN(ldprc->StorageDrain, maxRate);
        }

        //... limit soil perc by unused storage volume
        maxRate = ldprc->StorageExfil + ldprc->StorageDrain + ldprc->StorageEvap +
                  (storageThickness - storageDepth) *
                  storageVoidFrac/ldprc->Tstep;
        ldprc->SoilPerc = MIN(ldprc->SoilPerc, maxRate);

        //... limit surface infil. by unused soil volume
        maxRate = (soilPorosity - soilTheta) * soilThickness / ldprc->Tstep +
                ldprc->SoilPerc + ldprc->SoilEvap;
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
    }

    //... find surface layer outflow rate
    ldprc->SurfaceOutflow = getSurfaceOutflowRate(sp, surfaceDepth);

    //... compute overall layer flux rates
    f[SURF] = (ldprc->SurfaceInflow - ldprc->SurfaceEvap - ldprc->SurfaceInfil -
            ldprc->SurfaceOutflow) / ldprc->theLidProc->surface.voidFrac;
    f[SOIL] = (ldprc->SurfaceInfil - ldprc->SoilEvap - ldprc->SoilPerc) /
            ldprc->theLidProc->soil.thickness;
    if ( storageThickness == 0.0 ) f[STOR] = 0.0;
    else f[STOR] = (ldprc->SoilPerc - ldprc->StorageEvap - ldprc->StorageExfil -
            ldprc->StorageDrain) / ldprc->theLidProc->storage.voidFrac;
}

//=============================================================================

////  This function was re-written for release 5.1.011.  ////                  //(5.1.011)

void trenchFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates from the layers of an infiltration trench LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // Moisture level variables
    double surfaceDepth;
    double storageDepth;

    // Intermediate variables
    double availVolume;
    double maxRate;

    TLidprocShared *ldprc = &sp->LidprocShared;

    // Storage layer properties
    double storageThickness = ldprc->theLidProc->storage.thickness;
    double storageVoidFrac = ldprc->theLidProc->storage.voidFrac;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    ldprc->SurfaceVolume = surfaceDepth * ldprc->theLidProc->surface.voidFrac;
    ldprc->SoilVolume = 0.0;
    ldprc->StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = (storageThickness - storageDepth) * storageVoidFrac;
    getEvapRates(sp, ldprc->SurfaceVolume, 0.0, 0.0, ldprc->StorageVolume, 1.0);

    //... no storage evap if surface ponded
    if ( surfaceDepth > 0.0 ) ldprc->StorageEvap = 0.0;

    //... nominal storage inflow
    ldprc->StorageInflow = ldprc->SurfaceInflow + ldprc->SurfaceVolume / ldprc->Tstep;

    //... exfiltration rate out of storage layer
    ldprc->StorageExfil = getStorageExfilRate(sp);

    //... underdrain flow rate
    ldprc->StorageDrain = 0.0;
    if ( ldprc->theLidProc->drain.coeff > 0.0 )
    {
        ldprc->StorageDrain = getStorageDrainRate(sp, storageDepth, 0.0, 0.0, surfaceDepth);
    }

    //... limit storage exfiltration by available storage volume
    maxRate = ldprc->StorageInflow - ldprc->StorageEvap + storageDepth*storageVoidFrac/ldprc->Tstep;
    ldprc->StorageExfil = MIN(ldprc->StorageExfil, maxRate);
    ldprc->StorageExfil = MAX(ldprc->StorageExfil, 0.0);

    //... limit underdrain flow by volume above drain offset
    if ( ldprc->StorageDrain > 0.0 )
    {
        maxRate = -ldprc->StorageExfil - ldprc->StorageEvap;                                 //(5.1.012)
        if (storageDepth >= storageThickness ) maxRate += ldprc->StorageInflow;       //(5.1.012)
        if ( ldprc->theLidProc->drain.offset <= storageDepth )
        {
            maxRate += (storageDepth - ldprc->theLidProc->drain.offset) *
                       storageVoidFrac/ldprc->Tstep;
        }
        maxRate = MAX(maxRate, 0.0);
        ldprc->StorageDrain = MIN(ldprc->StorageDrain, maxRate);
    }

    //... limit storage inflow to not exceed storage layer capacity
    maxRate = (storageThickness - storageDepth)*storageVoidFrac/ldprc->Tstep +
            ldprc->StorageExfil + ldprc->StorageEvap + ldprc->StorageDrain;
    ldprc->StorageInflow = MIN(ldprc->StorageInflow, maxRate);

    //... equate surface infil to storage inflow
    ldprc->SurfaceInfil = ldprc->StorageInflow;

    //... find surface outflow rate
    ldprc->SurfaceOutflow = getSurfaceOutflowRate(sp, surfaceDepth);

    // ... find net fluxes for each layer
    f[SURF] = ldprc->SurfaceInflow - ldprc->SurfaceEvap - ldprc->StorageInflow -
            ldprc->SurfaceOutflow / ldprc->theLidProc->surface.voidFrac;
    f[STOR] = (ldprc->StorageInflow - ldprc->StorageEvap - ldprc->StorageExfil -
            ldprc->StorageDrain) / ldprc->theLidProc->storage.voidFrac;
    f[SOIL] = 0.0;
}

//=============================================================================

////  This function was re-written for release 5.1.011.  ////                  //(5.1.011)

void pavementFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates for the layers of a porous pavement LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    //... Moisture level variables
    double surfaceDepth;
    double paveDepth;
    double soilTheta;
    double storageDepth;

    TLidprocShared *ldprc = &sp->LidprocShared;

    //... Intermediate variables
    double pervFrac = (1.0 - ldprc->theLidProc->pavement.impervFrac);
    double storageInflow;    // inflow rate to storage layer (ft/s)
    double availVolume;
    double maxRate;

    //... LID layer properties
    double paveVoidFrac     = ldprc->theLidProc->pavement.voidFrac * pervFrac;
    double paveThickness    = ldprc->theLidProc->pavement.thickness;
    double soilThickness    = ldprc->theLidProc->soil.thickness;
    double soilPorosity     = ldprc->theLidProc->soil.porosity;
    double soilFieldCap     = ldprc->theLidProc->soil.fieldCap;
    double soilWiltPoint    = ldprc->theLidProc->soil.wiltPoint;
    double storageThickness = ldprc->theLidProc->storage.thickness;
    double storageVoidFrac  = ldprc->theLidProc->storage.voidFrac;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    paveDepth    = x[PAVE];
    soilTheta    = x[SOIL];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    ldprc->SurfaceVolume = surfaceDepth * ldprc->theLidProc->surface.voidFrac;
    ldprc->PaveVolume = paveDepth * paveVoidFrac;
    ldprc->SoilVolume = soilTheta * soilThickness;
    ldprc->StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = ldprc->SoilVolume - soilWiltPoint * soilThickness;
    getEvapRates(sp, ldprc->SurfaceVolume, ldprc->PaveVolume, availVolume,
            ldprc->StorageVolume, pervFrac);

    //... no storage evap if soil or pavement layer saturated
    if ( paveDepth >= paveThickness ||
       ( soilThickness > 0.0 && soilTheta >= soilPorosity )
       ) ldprc->StorageEvap = 0.0;

    //... find nominal rate of surface infiltration into pavement layer
    ldprc->SurfaceInfil = ldprc->SurfaceInflow + (ldprc->SurfaceVolume / ldprc->Tstep);

    //... find perc rate out of pavement layer
    ldprc->PavePerc = getPavementPermRate(sp);

    //... limit pavement perc by available water
    maxRate = ldprc->PaveVolume/ldprc->Tstep + ldprc->SurfaceInfil - ldprc->PaveEvap;
    maxRate = MAX(maxRate, 0.0);
    ldprc->PavePerc = MIN(ldprc->PavePerc, maxRate);

    //... find soil layer perc rate
    if ( soilThickness > 0.0 )
    {
        ldprc->SoilPerc = getSoilPercRate(sp, soilTheta);
        availVolume = (soilTheta - soilFieldCap) * soilThickness;
        maxRate = MAX(availVolume, 0.0) / ldprc->Tstep - ldprc->SoilEvap;                    //(5.1.012)
        ldprc->SoilPerc = MIN(ldprc->SoilPerc, maxRate);
        ldprc->SoilPerc = MAX(ldprc->SoilPerc, 0.0);
    }
    else ldprc->SoilPerc = ldprc->PavePerc;

    //... exfiltration rate out of storage layer
    ldprc->StorageExfil = getStorageExfilRate(sp);

    //... underdrain flow rate
    ldprc->StorageDrain = 0.0;
    if ( ldprc->theLidProc->drain.coeff > 0.0 )
    {
        ldprc->StorageDrain = getStorageDrainRate(sp, storageDepth, soilTheta,
                paveDepth, surfaceDepth);
    }

    //... check for adjacent saturated layers

    //... no soil layer, pavement & storage layers are full
    if ( soilThickness == 0.0 &&
         storageDepth >= storageThickness &&
         paveDepth >= paveThickness )
    {
        //... pavement outflow can't exceed storage outflow
        maxRate = ldprc->StorageEvap + ldprc->StorageDrain + ldprc->StorageExfil;
        if ( ldprc->PavePerc > maxRate ) ldprc->PavePerc = maxRate;

        //... storage outflow can't exceed pavement outflow
        else
        {
            //... use up available exfiltration capacity first
            ldprc->StorageExfil = MIN(ldprc->StorageExfil, ldprc->PavePerc);
            ldprc->StorageDrain = ldprc->PavePerc - ldprc->StorageExfil;
        }

        //... set soil perc to pavement perc
        ldprc->SoilPerc = ldprc->PavePerc;

        //... limit surface infil. by pavement perc
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, ldprc->PavePerc);
    }

    //... pavement, soil & storage layers are full
    else if ( soilThickness > 0 &&
              storageDepth >= storageThickness &&
              soilTheta >= soilPorosity &&
              paveDepth >= paveThickness )
    {
        //... find which layer has limiting flux rate
        maxRate = ldprc->StorageExfil + ldprc->StorageDrain;
        if ( ldprc->SoilPerc < maxRate) maxRate = ldprc->SoilPerc;
        else maxRate = MIN(maxRate, ldprc->PavePerc);

        //... use up available storage exfiltration capacity first
        if ( maxRate > ldprc->StorageExfil )
            ldprc->StorageDrain = maxRate - ldprc->StorageExfil;
        else
        {
            ldprc->StorageExfil = maxRate;
            ldprc->StorageDrain = 0.0;
        }
        ldprc->SoilPerc = maxRate;
        ldprc->PavePerc = maxRate;

        //... limit surface infil. by pavement perc
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, ldprc->PavePerc);
    }

    //... storage & soil layers are full
    else if ( soilThickness > 0.0 &&
              storageDepth >= storageThickness &&
              soilTheta >= soilPorosity )
    {
        //... soil perc can't exceed storage outflow
        maxRate = ldprc->StorageDrain + ldprc->StorageExfil;
        if ( ldprc->SoilPerc > maxRate ) ldprc->SoilPerc = maxRate;

        //... storage outflow can't exceed soil perc
        else
        {
            //... use up available exfiltration capacity first
            ldprc->StorageExfil = MIN(ldprc->StorageExfil, ldprc->SoilPerc);
            ldprc->StorageDrain = ldprc->SoilPerc - ldprc->StorageExfil;
        }

        //... limit surface infil. by available pavement volume
        availVolume = (paveThickness - paveDepth) * paveVoidFrac;
        maxRate = availVolume / ldprc->Tstep + ldprc->PavePerc + ldprc->PaveEvap;
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
    }

    //... soil and pavement layers are full
    else if ( soilThickness > 0.0 &&
              paveDepth >= paveThickness &&
              soilTheta >= soilPorosity )
    {
        ldprc->PavePerc = MIN(ldprc->PavePerc, ldprc->SoilPerc);
        ldprc->SoilPerc = ldprc->PavePerc;
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, ldprc->PavePerc);
    }

    //... no adjoining layers are full
    else
    {
        //... limit storage exfiltration by available storage volume
        //    (if no soil layer, SoilPerc is same as PavePerc)
        maxRate = ldprc->SoilPerc - ldprc->StorageEvap + ldprc->StorageVolume /
                ldprc->Tstep;
        maxRate = MAX(0.0, maxRate);
        ldprc->StorageExfil = MIN(ldprc->StorageExfil, maxRate);

        //... limit underdrain flow by volume above drain offset
        if ( ldprc->StorageDrain > 0.0 )
        {
            maxRate = -ldprc->StorageExfil - ldprc->StorageEvap;                             //(5.1.012)
            if (storageDepth >= storageThickness ) maxRate += ldprc->SoilPerc;        //(5.1.012)
            if ( ldprc->theLidProc->drain.offset <= storageDepth )
            {
                maxRate += (storageDepth - ldprc->theLidProc->drain.offset) *
                           storageVoidFrac/ldprc->Tstep;
            }
            maxRate = MAX(maxRate, 0.0);
            ldprc->StorageDrain = MIN(ldprc->StorageDrain, maxRate);
        }

        //... limit soil & pavement outflow by unused storage volume
        availVolume = (storageThickness - storageDepth) * storageVoidFrac;
        maxRate = availVolume/ldprc->Tstep + ldprc->StorageEvap +
                ldprc->StorageDrain + ldprc->StorageExfil;
        maxRate = MAX(maxRate, 0.0);
        if ( soilThickness > 0.0 )
        {
            ldprc->SoilPerc = MIN(ldprc->SoilPerc, maxRate);
            maxRate = (soilPorosity - soilTheta) * soilThickness / ldprc->Tstep +
                    ldprc->SoilPerc;
        }
        ldprc->PavePerc = MIN(ldprc->PavePerc, maxRate);

        //... limit surface infil. by available pavement volume
        availVolume = (paveThickness - paveDepth) * paveVoidFrac;
        maxRate = availVolume / ldprc->Tstep + ldprc->PavePerc + ldprc->PaveEvap;
        ldprc->SurfaceInfil = MIN(ldprc->SurfaceInfil, maxRate);
    }

    //... surface outflow
    ldprc->SurfaceOutflow = getSurfaceOutflowRate(sp, surfaceDepth);

    //... compute overall layer flux rates
    f[SURF] = ldprc->SurfaceInflow - ldprc->SurfaceEvap - ldprc->SurfaceInfil -
            ldprc->SurfaceOutflow;
    f[PAVE] = (ldprc->SurfaceInfil - ldprc->PaveEvap - ldprc->PavePerc) /
            paveVoidFrac;
    if ( ldprc->theLidProc->soil.thickness > 0.0)
    {
        f[SOIL] = (ldprc->PavePerc - ldprc->SoilEvap - ldprc->SoilPerc) /
                soilThickness;
        storageInflow = ldprc->SoilPerc;
    }
    else
    {
        f[SOIL] = 0.0;
        storageInflow = ldprc->PavePerc;
        ldprc->SoilPerc = 0.0;
    }
    f[STOR] = (storageInflow - ldprc->StorageEvap - ldprc->StorageExfil -
            ldprc->StorageDrain) / storageVoidFrac;
}

//=============================================================================

void swaleFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates from a vegetative swale LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double depth;            // depth of surface water in swale (ft)
    double topWidth;         // top width of full swale (ft)
    double botWidth;         // bottom width of swale (ft)
    double length;           // length of swale (ft)
    double surfInflow;       // inflow rate to swale (cfs)
    double surfWidth;        // top width at current water depth (ft)
    double surfArea;         // surface area of current water depth (ft2)
    double flowArea;         // x-section flow area (ft2)
    double lidArea;          // surface area of full swale (ft2)
    double hydRadius;        // hydraulic radius for current depth (ft)
    double slope;            // slope of swale side wall (run/rise)
    double volume;           // swale volume at current water depth (ft3)
    double dVdT;             // change in volume w.r.t. time (cfs)
    double dStore;           // depression storage depth (ft)
    double xDepth;           // depth above depression storage (ft)

    TLidprocShared *ldprc = &sp->LidprocShared;

    //... retrieve state variable from work vector
    depth = x[SURF];
    depth = MIN(depth, ldprc->theLidProc->surface.thickness);

    //... depression storage depth
    dStore = 0.0;

    //... get swale's bottom width
    //    (0.5 ft minimum to avoid numerical problems)
    slope = ldprc->theLidProc->surface.sideSlope;
    topWidth = ldprc->theLidUnit->fullWidth;
    topWidth = MAX(topWidth, 0.5);
    botWidth = topWidth - 2.0 * slope * ldprc->theLidProc->surface.thickness;
    if ( botWidth < 0.5 )
    {
        botWidth = 0.5;
        slope = 0.5 * (topWidth - 0.5) / ldprc->theLidProc->surface.thickness;
    }

    //... swale's length
    lidArea = ldprc->theLidUnit->area;
    length = lidArea / topWidth;

    //... top width, surface area and flow area of current ponded depth
    surfWidth = botWidth + 2.0 * slope * depth;
    surfArea = length * surfWidth;
    flowArea = (depth * (botWidth + slope * depth)) *
            ldprc->theLidProc->surface.voidFrac;

    //... wet volume and effective depth
    volume = length * flowArea;

    //... surface inflow into swale (cfs)
    surfInflow = ldprc->SurfaceInflow * lidArea;

    //... ET rate in cfs
    ldprc->SurfaceEvap = ldprc->EvapRate * surfArea;
    ldprc->SurfaceEvap = MIN(ldprc->SurfaceEvap, volume/ldprc->Tstep);

    //... infiltration rate to native soil in cfs
    ldprc->StorageExfil = ldprc->SurfaceInfil * surfArea;                                    //(5.1.011)

    //... no surface outflow if depth below depression storage
    xDepth = depth - dStore;
    if ( xDepth <= ZERO ) ldprc->SurfaceOutflow = 0.0;

    //... otherwise compute a surface outflow
    else
    {
        //... modify flow area to remove depression storage,
        flowArea -= (dStore * (botWidth + slope * dStore)) *
                ldprc->theLidProc->surface.voidFrac;
        if ( flowArea < ZERO ) ldprc->SurfaceOutflow = 0.0;
        else
        {
            //... compute hydraulic radius
            botWidth = botWidth + 2.0 * dStore * slope;
            hydRadius = botWidth + 2.0 * xDepth * sqrt(1.0 + slope*slope);
            hydRadius = flowArea / hydRadius;

            //... use Manning Eqn. to find outflow rate in cfs
            ldprc->SurfaceOutflow = ldprc->theLidProc->surface.alpha * flowArea *
                             pow(hydRadius, 2./3.);
        }
    }

    //... net flux rate (dV/dt) in cfs
    dVdT = surfInflow - ldprc->SurfaceEvap - ldprc->StorageExfil -
            ldprc->SurfaceOutflow;           //(5.1.011)

    //... when full, any net positive inflow becomes spillage
    if ( depth == ldprc->theLidProc->surface.thickness && dVdT > 0.0 )
    {
        ldprc->SurfaceOutflow += dVdT;
        dVdT = 0.0;
    }

    //... convert flux rates to ft/s
    ldprc->SurfaceEvap /= lidArea;
    ldprc->StorageExfil /= lidArea;                                                   //(5.1.011)
    ldprc->SurfaceOutflow /= lidArea;
    f[SURF] = dVdT / surfArea;
    f[SOIL] = 0.0;
    f[STOR] = 0.0;

    //... assign values to layer volumes
    ldprc->SurfaceVolume = volume / lidArea;
    ldprc->SoilVolume = 0.0;
    ldprc->StorageVolume = 0.0;
}

//=============================================================================

////  This function was re-written for release 5.1.007.  ////                  //(5.1.007)

void barrelFluxRates(SWMM_Project *sp, double x[], double f[])
//
//  Purpose: computes flux rates for a rain barrel LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double storageDepth = x[STOR];
	double head;
    double maxValue;

    TLidprocShared *ldprc = &sp->LidprocShared;

    //... assign values to layer volumes
    ldprc->SurfaceVolume = 0.0;
    ldprc->SoilVolume = 0.0;
    ldprc->StorageVolume = storageDepth;

    //... initialize flows
    ldprc->SurfaceInfil = 0.0;
    ldprc->SurfaceOutflow = 0.0;
    ldprc->StorageDrain = 0.0;

    //... compute outflow if time since last rain exceeds drain delay
    //    (dryTime is updated in lid.evalLidUnit at each time step)
    if ( ldprc->theLidProc->drain.delay == 0.0 ||
            ldprc->theLidUnit->dryTime >= ldprc->theLidProc->drain.delay )
	{
	    head = storageDepth - ldprc->theLidProc->drain.offset;
		if ( head > 0.0 )
	    {
		    ldprc->StorageDrain = getStorageDrainRate(sp, storageDepth, 0.0, 0.0, 0.0);
		    maxValue = (head/ldprc->Tstep);
		    ldprc->StorageDrain = MIN(ldprc->StorageDrain, maxValue);
		}
	}

    //... limit inflow to available storage
    ldprc->StorageInflow = ldprc->SurfaceInflow;
    maxValue = (ldprc->theLidProc->storage.thickness - storageDepth) /
            ldprc->Tstep + ldprc->StorageDrain;
    ldprc->StorageInflow = MIN(ldprc->StorageInflow, maxValue);
    ldprc->SurfaceInfil = ldprc->StorageInflow;

    //... assign values to layer flux rates
    f[SURF] = ldprc->SurfaceInflow - ldprc->StorageInflow;
    f[STOR] = ldprc->StorageInflow - ldprc->StorageDrain;
    f[SOIL] = 0.0;
}

//=============================================================================

double getSurfaceOutflowRate(SWMM_Project *sp, double depth)
//
//  Purpose: computes outflow rate from a LID's surface layer.
//  Input:   depth = depth of ponded water on surface layer (ft)
//  Output:  returns outflow from surface layer (ft/s)
//
//  Note: this function should not be applied to swales or rain barrels.
//
{
    double delta;
    double outflow;

    TLidprocShared *ldprc = &sp->LidprocShared;

    //... no outflow if ponded depth below storage depth
    delta = depth - ldprc->theLidProc->surface.thickness;
    if ( delta < 0.0 ) return 0.0;

    //... compute outflow from overland flow Manning equation
    outflow = ldprc->theLidProc->surface.alpha * pow(delta, 5.0/3.0) *
            ldprc->theLidUnit->fullWidth / ldprc->theLidUnit->area;
    outflow = MIN(outflow, delta / ldprc->Tstep);
    return outflow;
}

//=============================================================================

double getPavementPermRate(SWMM_Project *sp)
//
//  Purpose: computes reduced permeability of a pavement layer due to
//           clogging.
//  Input:   none
//  Output:  returns the reduced permeability of the pavement layer (ft/s).
//
{
    double permRate;
    double permReduction;

    TLidprocShared *ldprc = &sp->LidprocShared;

    permReduction = ldprc->theLidProc->pavement.clogFactor;
    if ( permReduction > 0.0 )
    {
        permReduction = ldprc->theLidUnit->waterBalance.inflow / permReduction;
        permReduction = MIN(permReduction, 1.0);
    }
    permRate = ldprc->theLidProc->pavement.kSat * (1.0 - permReduction);
    return permRate;
}

//=============================================================================

////  This function was modified for release 5.1.011.  ////                    //(5.1.011)

double getSoilPercRate(SWMM_Project *sp, double theta)
//
//  Purpose: computes percolation rate of water through a LID's soil layer.
//  Input:   theta = moisture content (fraction)
//  Output:  returns percolation rate within soil layer (ft/s)
//
{
    double delta;            // moisture deficit

    TLidprocShared *ldprc = &sp->LidprocShared;

    // ... no percolation if soil moisture <= field capacity
    if ( theta <= ldprc->theLidProc->soil.fieldCap ) return 0.0;

    // ... perc rate = unsaturated hydraulic conductivity
    delta = ldprc->theLidProc->soil.porosity - theta;
    return ldprc->theLidProc->soil.kSat * exp(-delta * ldprc->theLidProc->soil.kSlope);

}

//=============================================================================

double getStorageExfilRate(SWMM_Project *sp)                                                   //(5.1.011)
//
//  Purpose: computes exfiltration rate from storage zone into                 //(5.1.011)
//           native soil beneath a LID.
//  Input:   depth = depth of water storage zone (ft)
//  Output:  returns infiltration rate (ft/s)
//
{
    double infil = 0.0;
    double clogFactor = 0.0;

    TLidprocShared *ldprc = &sp->LidprocShared;

    if ( ldprc->theLidProc->storage.kSat == 0.0 ) return 0.0;
    if ( ldprc->MaxNativeInfil == 0.0 ) return 0.0;

    //... reduction due to clogging
    clogFactor = ldprc->theLidProc->storage.clogFactor;
    if ( clogFactor > 0.0 )
    {
        clogFactor = ldprc->theLidUnit->waterBalance.inflow / clogFactor;
        clogFactor = MIN(clogFactor, 1.0);
    }

    //... infiltration rate = storage Ksat reduced by any clogging
    infil = ldprc->theLidProc->storage.kSat * (1.0 - clogFactor);

    //... limit infiltration rate by any groundwater-imposed limit
    return MIN(infil, ldprc->MaxNativeInfil);
}

//=============================================================================

////  This function was modified for release 5.1.011.  ////                    //(5.1.011)

double  getStorageDrainRate(SWMM_Project *sp, double storageDepth,
        double soilTheta, double paveDepth, double surfaceDepth)
//
//  Purpose: computes underdrain flow rate in a LID's storage layer.
//  Input:   storageDepth = depth of water in storage layer (ft)
//           soilTheta    = moisture content of soil layer
//           paveDepth    = effective depth of water in pavement layer (ft)
//           surfaceDepth = depth of ponded water on surface layer (ft)
//  Output:  returns flow in underdrain (ft/s)
//
//  Note:    drain eqn. is evaluated in user's units.
//  Note:    head on drain is water depth in storage layer plus the
//           layers above it (soil, pavement, and surface in that order)
//           minus the drain outlet offset.
{
    double head = storageDepth;
    double outflow = 0.0;

    TLidprocShared *ldprc = &sp->LidprocShared;

    double paveThickness    = ldprc->theLidProc->pavement.thickness;
    double soilThickness    = ldprc->theLidProc->soil.thickness;
    double soilPorosity     = ldprc->theLidProc->soil.porosity;
    double soilFieldCap     = ldprc->theLidProc->soil.fieldCap;
    double storageThickness = ldprc->theLidProc->storage.thickness;

    // --- storage layer is full
    if ( storageDepth >= storageThickness )
    {
        // --- a soil layer exists
        if ( soilThickness > 0.0 )
        {
            // --- increase head by fraction of soil layer saturated
            if ( soilTheta > soilFieldCap )
            {
                head += (soilTheta - soilFieldCap) /
                        (soilPorosity - soilFieldCap) * soilThickness;

                // --- soil layer is saturated, increase head by water
                //     depth in layer above it
                if ( soilTheta >= soilPorosity )
                {
                    if ( paveThickness > 0.0 ) head += paveDepth;
                    else head += surfaceDepth;
                }
            }
        }

        // --- no soil layer so increase head by water level in pavement
        //     layer and possibly surface layer
        if ( paveThickness > 0.0 )
        {
            head += paveDepth;
            if ( paveDepth >= paveThickness ) head += surfaceDepth;
        }
    }

    // --- make head relative to drain offset
    head -= ldprc->theLidProc->drain.offset;

    // ... compute drain outflow from underdrain flow equation in user units
    //     (head in inches or mm, flow rate in in/hr or mm/hr)
    if ( head > ZERO )
    {
        head *= UCF(sp, RAINDEPTH);
        outflow = ldprc->theLidProc->drain.coeff *
                  pow(head, ldprc->theLidProc->drain.expon);
        outflow /= UCF(sp, RAINFALL);
    }
    return outflow;
}

//=============================================================================

////  This function was modified for release 5.1.007.  ////                    //(5.1.007)

double getDrainMatOutflow(SWMM_Project *sp, double depth)
//
//  Purpose: computes flow rate through a green roof's drainage mat.
//  Input:   depth = depth of water in drainage mat (ft)
//  Output:  returns flow in drainage mat (ft/s)
//
{
    TLidprocShared *ldprc = &sp->LidprocShared;

    //... default is to pass all inflow
    double result = ldprc->SoilPerc;

    //... otherwise use Manning eqn. if its parameters were supplied
    if ( ldprc->theLidProc->drainMat.alpha > 0.0 )
    {
        result = ldprc->theLidProc->drainMat.alpha * pow(depth, 5.0/3.0) *
                ldprc->theLidUnit->fullWidth / ldprc->theLidUnit->area *
                ldprc->theLidProc->drainMat.voidFrac;                                //(5.1.008)
    }
    return result;
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void getEvapRates(SWMM_Project *sp, double surfaceVol, double paveVol,
        double soilVol, double storageVol, double pervFrac)                                        //(5.1.011)
//
//  Purpose: computes surface, pavement, soil, and storage evaporation rates.
//  Input:   surfaceVol = volume/area of ponded water on surface layer (ft)
//           paveVol    = volume/area of water in pavement pores (ft)
//           soilVol    = volume/area of water in soil (or pavement) pores (ft)
//           storageVol = volume/area of water in storage layer (ft)
//           pervFrac   = fraction of surface layer that is pervious           //(5.1.011)
//  Output:  none
//
{
    double availEvap;

    TLidprocShared *ldprc = &sp->LidprocShared;

    //... surface evaporation flux
    availEvap = ldprc->EvapRate;
    ldprc->SurfaceEvap = MIN(availEvap, surfaceVol/ldprc->Tstep);
    ldprc->SurfaceEvap = MAX(0.0, ldprc->SurfaceEvap);
    availEvap = MAX(0.0, (availEvap - ldprc->SurfaceEvap));
    availEvap *= pervFrac;                                                     //(5.1.011)

    //... no subsurface evap if water is infiltrating
    if ( ldprc->SurfaceInfil > 0.0 )
    {
        ldprc->PaveEvap = 0.0;
        ldprc->SoilEvap = 0.0;
        ldprc->StorageEvap = 0.0;
    }
    else
    {
        //... pavement evaporation flux
        ldprc->PaveEvap = MIN(availEvap, paveVol / ldprc->Tstep);
        availEvap = MAX(0.0, (availEvap - ldprc->PaveEvap));

        //... soil evaporation flux
        ldprc->SoilEvap = MIN(availEvap, soilVol / ldprc->Tstep);
        availEvap = MAX(0.0, (availEvap - ldprc->SoilEvap));

        //... storage evaporation flux
        ldprc->StorageEvap = MIN(availEvap, storageVol / ldprc->Tstep);
    }
}

//=============================================================================

double getSurfaceOverflowRate(SWMM_Project *sp, double* surfaceDepth)
//
//  Purpose: finds surface overflow rate from a LID unit.
//  Input:   surfaceDepth = depth of water stored in surface layer (ft)
//  Output:  returns the overflow rate (ft/s)
//
{
    TLidprocShared *ldprc = &sp->LidprocShared;

    double delta = *surfaceDepth - ldprc->theLidProc->surface.thickness;
    if (  delta <= 0.0 ) return 0.0;
    *surfaceDepth = ldprc->theLidProc->surface.thickness;
    return delta * ldprc->theLidProc->surface.voidFrac / ldprc->Tstep;
}

//=============================================================================

void updateWaterBalance(SWMM_Project *sp, TLidUnit *lidUnit, double inflow,
        double evap, double infil, double surfFlow, double drainFlow,
        double storage)
//
//  Purpose: updates components of the water mass balance for a LID unit
//           over the current time step.
//  Input:   lidUnit   = a particular LID unit
//           inflow    = runon + rainfall to the LID unit (ft/s)
//           evap      = evaporation rate from the unit (ft/s)
//           infil     = infiltration out the bottom of the unit (ft/s)
//           surfFlow  = surface runoff from the unit (ft/s)
//           drainFlow = underdrain flow from the unit
//           storage   = volume of water stored in the unit (ft)
//  Output:  none
//
{
    TLidprocShared *ldprc = &sp->LidprocShared;

    lidUnit->waterBalance.inflow += inflow * ldprc->Tstep;
    lidUnit->waterBalance.evap += evap * ldprc->Tstep;
    lidUnit->waterBalance.infil += infil * ldprc->Tstep;
    lidUnit->waterBalance.surfFlow += surfFlow * ldprc->Tstep;
    lidUnit->waterBalance.drainFlow += drainFlow * ldprc->Tstep;
    lidUnit->waterBalance.finalVol = storage;
}

//=============================================================================

int modpuls_solve(SWMM_Project *sp, int n, double* x, double* xOld, double* xPrev,
                  double* xMin, double* xMax, double* xTol,
                  double* qOld, double* q, double dt, double omega,            //(5.1.007)
                  void (*derivs)(SWMM_Project*, double*, double*))
//
//  Purpose: solves system of equations dx/dt = q(x) for x at end of time step
//           dt using a modified Puls method.
//  Input:   n = number of state variables
//           x = vector of state variables
//           xOld = state variable values at start of time step
//           xPrev = state variable values from previous iteration
//           xMin = lower limits on state variables
//           xMax = upper limits on state variables
//           xTol = convergence tolerances on state variables
//           qOld = flux rates at start of time step
//           q = flux rates at end of time step
//           dt = time step (sec)
//           omega = time weighting parameter (use 0 for Euler method          //(5.1.007)
//                   or 0.5 for modified Puls method)                          //(5.1.007)
//           derivs = pointer to function that computes flux rates q as a
//                    function of state variables x
//  Output:  returns number of steps required for convergence (or 0 if
//           process doesn't converge)
//
{
    int i;
    int canStop;
    int steps = 1;
    int maxSteps = 20;

    //... initialize state variable values
    for (i=0; i<n; i++)
    {
        xOld[i] = x[i];
        xPrev[i] = x[i];
    }

    //... repeat until convergence achieved
    while (steps < maxSteps)
    {
        //... compute flux rates for current state levels
        canStop = 1;
        derivs(sp, x, q);

        //... update state levels based on current flux rates
        for (i=0; i<n; i++)
        {
            x[i] = xOld[i] + (omega*qOld[i] + (1.0 - omega)*q[i]) * dt;
            x[i] = MIN(x[i], xMax[i]);
            x[i] = MAX(x[i], xMin[i]);

            if ( omega > 0.0 &&                                                //(5.1.007)
                 fabs(x[i] - xPrev[i]) > xTol[i] ) canStop = 0;
            xPrev[i] = x[i];
        }

        //... return if process converges
        if (canStop) return steps;
        steps++;
    }

    //... no convergence so return 0
    return 0;
}
