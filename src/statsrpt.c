//-----------------------------------------------------------------------------
//   statsrpt.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14 (Build 5.1.001)
//             09/15/14 (Build 5.1.007)
//             03/19/15 (Build 5.1.008)
//             04/30/15 (Build 5.1.009)
//             08/01/16 (Build 5.1.011)
//   Author:   L. Rossman
//
//   Report writing functions for summary statistics.
//
//   Build 5.1.008:
//   - New Groundwater Summary table added.
//   - Reported Max. Depth added to Node Depth Summary table.
//
//   Build 5.1.009:
//   - Units on column heading in Node Inflow Summary table fixed.
//
//   Build 5.1.011:
//   - Redundant units conversion on max. reported node depth removed.
//   - Node Surcharge table only produced for dynamic wave routing.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "headers.h"
#include "lid.h"

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
extern TSubcatchStats* SubcatchStats;          // defined in STATS.C
extern TNodeStats*     NodeStats;
extern TLinkStats*     LinkStats;
extern TStorageStats*  StorageStats;
extern TOutfallStats*  OutfallStats;
extern TPumpStats*     PumpStats;
extern double          MaxOutfallFlow;
extern double          MaxRunoffFlow;
extern double*         NodeInflow;             // defined in MASSBAL.C
extern double*         NodeOutflow;            // defined in massbal.c

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
void    writeSubcatchRunoff(SWMM_Project *sp);
void    writeGroundwater(SWMM_Project *sp);                                                //(5.1.008)
void    writeSubcatchLoads(SWMM_Project *sp);
void    writeNodeDepths(SWMM_Project *sp);
void    writeNodeFlows(SWMM_Project *sp);
void    writeNodeSurcharge(SWMM_Project *sp);
void    writeNodeFlooding(SWMM_Project *sp);
void    writeStorageVolumes(SWMM_Project *sp);
void    writeOutfallLoads(SWMM_Project *sp);
void    writeLinkFlows(SWMM_Project *sp);
void    writeFlowClass(SWMM_Project *sp);
void    writeLinkSurcharge(SWMM_Project *sp);
void    writePumpFlows(SWMM_Project *sp);
void    writeLinkLoads(SWMM_Project *sp);

#define WRITE(x) (report_writeLine(sp,(x)))

static char   FlowFmt[6];
static double Vcf;

//=============================================================================

void statsrpt_writeReport(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reports simulation summary statistics.
//
{
    // --- set number of decimal places for reporting flow values
    if ( sp->FlowUnits == MGD || sp->FlowUnits == CMS ) strcpy(FlowFmt, "%9.3f");
    else strcpy(FlowFmt, "%9.2f");

    // --- volume conversion factor from ft3 to Mgal or Mliters
    if (sp->UnitSystem == US) Vcf = 7.48 / 1.0e6;
    else                  Vcf = 28.317 / 1.0e6;

    // --- report summary results for subcatchment runoff 
    if ( sp->Nobjects[SUBCATCH] > 0 )
    {
        if ( !sp->IgnoreRainfall ||
             (sp->Nobjects[SNOWMELT] > 0 && !sp->IgnoreSnowmelt) ||
             (sp->Nobjects[AQUIFER] > 0  && !sp->IgnoreGwater) )
        {
            writeSubcatchRunoff(sp);
            lid_writeWaterBalance(sp);
            if ( !sp->IgnoreGwater ) writeGroundwater(sp);                           //(5.1.008)
            if ( sp->Nobjects[POLLUT] > 0 && !sp->IgnoreQuality) writeSubcatchLoads(sp);
        }
    }

    // --- report summary results for flow routing
    if ( sp->Nobjects[LINK] > 0 && !sp->IgnoreRouting )
    {
        writeNodeDepths(sp);
        writeNodeFlows(sp);
        if ( sp->RouteModel == DW ) writeNodeSurcharge(sp);                          //(5.1.011)
        writeNodeFlooding(sp);
        writeStorageVolumes(sp);
        writeOutfallLoads(sp);
        writeLinkFlows(sp);
        writeFlowClass(sp);
        writeLinkSurcharge(sp);
        writePumpFlows(sp);
        if ( sp->Nobjects[POLLUT] > 0 && !sp->IgnoreQuality) writeLinkLoads(sp);
    }
}

//=============================================================================

void writeSubcatchRunoff(SWMM_Project *sp)
{
    int    j;
    double a, x, r;

    if ( sp->Nobjects[SUBCATCH] == 0 ) return;
    WRITE("");
    WRITE("***************************");
    WRITE("Subcatchment Runoff Summary");
    WRITE("***************************");
    WRITE("");
    fprintf(sp->Frpt.file,

"\n  --------------------------------------------------------------------------------------------------------"
"\n                            Total      Total      Total      Total      Total       Total     Peak  Runoff"
"\n                           Precip      Runon       Evap      Infil     Runoff      Runoff   Runoff   Coeff");
    if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file,
"\n  Subcatchment                 in         in         in         in         in    %8s      %3s",
        VolUnitsWords[sp->UnitSystem], FlowUnitWords[sp->FlowUnits]);
    else fprintf(sp->Frpt.file,
"\n  Subcatchment                 mm         mm         mm         mm         mm    %8s      %3s",
        VolUnitsWords[sp->UnitSystem], FlowUnitWords[sp->FlowUnits]);
    fprintf(sp->Frpt.file,
"\n  --------------------------------------------------------------------------------------------------------");

    for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
    {
        a = sp->Subcatch[j].area;
        if ( a == 0.0 ) continue;
        fprintf(sp->Frpt.file, "\n  %-20s", sp->Subcatch[j].ID);
        x = SubcatchStats[j].precip * UCF(sp, RAINDEPTH);
        fprintf(sp->Frpt.file, " %10.2f", x/a);
        x = SubcatchStats[j].runon * UCF(sp, RAINDEPTH); 
        fprintf(sp->Frpt.file, " %10.2f", x/a);
        x = SubcatchStats[j].evap * UCF(sp, RAINDEPTH);
        fprintf(sp->Frpt.file, " %10.2f", x/a);
        x = SubcatchStats[j].infil * UCF(sp, RAINDEPTH); 
        fprintf(sp->Frpt.file, " %10.2f", x/a);
        x = SubcatchStats[j].runoff * UCF(sp, RAINDEPTH);
        fprintf(sp->Frpt.file, " %10.2f", x/a);
        x = SubcatchStats[j].runoff * Vcf;
	fprintf(sp->Frpt.file, "%12.2f", x);
        x = SubcatchStats[j].maxFlow * UCF(sp, FLOW);
        fprintf(sp->Frpt.file, " %8.2f", x);
        r = SubcatchStats[j].precip + SubcatchStats[j].runon;
        if ( r > 0.0 ) r = SubcatchStats[j].runoff / r;
        fprintf(sp->Frpt.file, "%8.3f", r);
    }
    WRITE("");
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

void    writeGroundwater(SWMM_Project *sp)
{
    int i, j;
    int count = 0;
    double totalSeconds = sp->NewRunoffTime / 1000.;
    double x[9];

    if ( sp->Nobjects[SUBCATCH] == 0 ) return;
    for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
    {
        if ( sp->Subcatch[j].groundwater != NULL ) count++;
    }
    if ( count == 0 ) return;

    WRITE("");
    WRITE("*******************");
    WRITE("Groundwater Summary");
    WRITE("*******************");
    WRITE("");
    fprintf(sp->Frpt.file,

"\n  -----------------------------------------------------------------------------------------------------"
"\n                                            Total    Total  Maximum  Average  Average    Final    Final"
"\n                          Total    Total    Lower  Lateral  Lateral    Upper    Water    Upper    Water"
"\n                          Infil     Evap  Seepage  Outflow  Outflow   Moist.    Table   Moist.    Table");
    if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file,
"\n  Subcatchment               in       in       in       in      %3s                ft                ft",
        FlowUnitWords[sp->FlowUnits]);
    else fprintf(sp->Frpt.file,
"\n  Subcatchment               mm       mm       mm       mm      %3s                 m                 m",
        FlowUnitWords[sp->FlowUnits]);
    fprintf(sp->Frpt.file,
"\n  -----------------------------------------------------------------------------------------------------");

    for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
    {
        if ( sp->Subcatch[j].area == 0.0 || sp->Subcatch[j].groundwater == NULL ) continue;
        fprintf(sp->Frpt.file, "\n  %-20s", sp->Subcatch[j].ID);
        x[0] = sp->Subcatch[j].groundwater->stats.infil * UCF(sp, RAINDEPTH);
        x[1] = sp->Subcatch[j].groundwater->stats.evap * UCF(sp, RAINDEPTH);
        x[2] = sp->Subcatch[j].groundwater->stats.deepFlow * UCF(sp, RAINDEPTH);
        x[3] = sp->Subcatch[j].groundwater->stats.latFlow * UCF(sp, RAINDEPTH);
        x[4] = sp->Subcatch[j].groundwater->stats.maxFlow * UCF(sp, FLOW) * sp->Subcatch[j].area;
        x[5] = sp->Subcatch[j].groundwater->stats.avgUpperMoist / totalSeconds;
        x[6] = sp->Subcatch[j].groundwater->stats.avgWaterTable * UCF(sp, LENGTH) /
               totalSeconds;
        x[7] = sp->Subcatch[j].groundwater->stats.finalUpperMoist;
        x[8] = sp->Subcatch[j].groundwater->stats.finalWaterTable * UCF(sp, LENGTH);
        for (i = 0; i < 9; i++) fprintf(sp->Frpt.file, " %8.2f", x[i]);
    }
    WRITE("");
}

//=============================================================================

void writeSubcatchLoads(SWMM_Project *sp)
{
    int i, j, p;
    double x;
    double* totals; 
    char  units[15];
    char  subcatchLine[] = "--------------------";
    char  pollutLine[]   = "--------------";

    // --- create an array to hold total loads for each pollutant
    totals = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
    if ( totals )
    {
        // --- print the table headings 
        WRITE("");
        WRITE("****************************");
        WRITE("Subcatchment Washoff Summary");
        WRITE("****************************");
        WRITE("");
        fprintf(sp->Frpt.file, "\n  %s", subcatchLine);
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%s", pollutLine);
        fprintf(sp->Frpt.file, "\n                      ");
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%14s", Pollut[p].ID);
        fprintf(sp->Frpt.file, "\n  Subcatchment        ");
        for (p = 0; p < sp->Nobjects[POLLUT]; p++)
        {
            i = sp->UnitSystem;
            if ( Pollut[p].units == COUNT ) i = 2;
            strcpy(units, LoadUnitsWords[i]);
            fprintf(sp->Frpt.file, "%14s", units);
            totals[p] = 0.0;
        }
        fprintf(sp->Frpt.file, "\n  %s", subcatchLine);
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%s", pollutLine);

        // --- print the pollutant loadings from each subcatchment
        for ( j = 0; j < sp->Nobjects[SUBCATCH]; j++ )
        {
            fprintf(sp->Frpt.file, "\n  %-20s", sp->Subcatch[j].ID);
            for (p = 0; p < sp->Nobjects[POLLUT]; p++)
            {
                x = sp->Subcatch[j].totalLoad[p];
                totals[p] += x;
                if ( Pollut[p].units == COUNT ) x = LOG10(x);
				fprintf(sp->Frpt.file, "%14.3f", x);
            }
        }

        // --- print the total loading of each pollutant
        fprintf(sp->Frpt.file, "\n  %s", subcatchLine);
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%s", pollutLine);
        fprintf(sp->Frpt.file, "\n  System              ");
        for (p = 0; p < sp->Nobjects[POLLUT]; p++)
        {
            x = totals[p];
            if ( Pollut[p].units == COUNT ) x = LOG10(x);
			fprintf(sp->Frpt.file, "%14.3f", x);
        }
        free(totals);
        WRITE("");
    }
}

//=============================================================================

////  Function modified for release 5.1.008.  ////                             //(5.1.008)

void writeNodeDepths(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes simulation statistics for nodes to report file.
//
{
    int j, days, hrs, mins;

    if ( sp->Nobjects[LINK] == 0 ) return;

    WRITE("");
    WRITE("******************");
    WRITE("Node Depth Summary");
    WRITE("******************");
    WRITE("");

    fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------------------"
"\n                                 Average  Maximum  Maximum  Time of Max    Reported"
"\n                                   Depth    Depth      HGL   Occurrence   Max Depth");
    if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file,
"\n  Node                 Type         Feet     Feet     Feet  days hr:min        Feet");
    else fprintf(sp->Frpt.file,
"\n  Node                 Type       Meters   Meters   Meters  days hr:min      Meters");
    fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------------------");

    for ( j = 0; j < sp->Nobjects[NODE]; j++ )
    {
        fprintf(sp->Frpt.file, "\n  %-20s", sp->Node[j].ID);
        fprintf(sp->Frpt.file, " %-9s ", NodeTypeWords[sp->Node[j].type]);
        getElapsedTime(sp, NodeStats[j].maxDepthDate, &days, &hrs, &mins);
        fprintf(sp->Frpt.file, "%7.2f  %7.2f  %7.2f  %4d  %02d:%02d  %10.2f",
            NodeStats[j].avgDepth / sp->StepCount * UCF(sp, LENGTH),
            NodeStats[j].maxDepth * UCF(sp, LENGTH),
            (NodeStats[j].maxDepth + sp->Node[j].invertElev) * UCF(sp, LENGTH),
            days, hrs, mins, NodeStats[j].maxRptDepth);                        //(5.1.011)
    }
    WRITE("");
}

//=============================================================================

void writeNodeFlows(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes flow statistics for nodes to report file.
//
{
    int j;
    int days1, hrs1, mins1;

    WRITE("");
    WRITE("*******************");
    WRITE("Node Inflow Summary");
    WRITE("*******************");
    WRITE("");

    fprintf(sp->Frpt.file,
"\n  -------------------------------------------------------------------------------------------------"
"\n                                  Maximum  Maximum                  Lateral       Total        Flow"
"\n                                  Lateral    Total  Time of Max      Inflow      Inflow     Balance"
"\n                                   Inflow   Inflow   Occurrence      Volume      Volume       Error"
"\n  Node                 Type           %3s      %3s  days hr:min    %8s    %8s     Percent",
        FlowUnitWords[sp->FlowUnits], FlowUnitWords[sp->FlowUnits], VolUnitsWords[sp->UnitSystem],
        VolUnitsWords[sp->UnitSystem]);                                            //(5.1.009)
    fprintf(sp->Frpt.file,
"\n  -------------------------------------------------------------------------------------------------");

    for ( j = 0; j < sp->Nobjects[NODE]; j++ )
    {
        fprintf(sp->Frpt.file, "\n  %-20s", sp->Node[j].ID);
        fprintf(sp->Frpt.file, " %-9s", NodeTypeWords[sp->Node[j].type]);
        getElapsedTime(sp, NodeStats[j].maxInflowDate, &days1, &hrs1, &mins1);
        fprintf(sp->Frpt.file, FlowFmt, NodeStats[j].maxLatFlow * UCF(sp, FLOW));
        fprintf(sp->Frpt.file, FlowFmt, NodeStats[j].maxInflow * UCF(sp, FLOW));
        fprintf(sp->Frpt.file, "  %4d  %02d:%02d", days1, hrs1, mins1);
		fprintf(sp->Frpt.file, "%12.3g", NodeStats[j].totLatFlow * Vcf);
		fprintf(sp->Frpt.file, "%12.3g", NodeInflow[j] * Vcf);
    	if ( fabs(NodeOutflow[j]) < 1.0 )
            fprintf(sp->Frpt.file, "%12.3f %s",
                (NodeInflow[j]-NodeOutflow[j])*Vcf*1.0e6,
                VolUnitsWords2[sp->UnitSystem]);
	    else
            fprintf(sp->Frpt.file, "%12.3f", (NodeInflow[j]-NodeOutflow[j]) /
                                          NodeOutflow[j]*100.); 
    }
    WRITE("");
}

//=============================================================================

void writeNodeSurcharge(SWMM_Project *sp)
{
    int    j, n = 0;
    double t, d1, d2;

    WRITE("");
    WRITE("**********************");
    WRITE("Node Surcharge Summary");
    WRITE("**********************");
    WRITE("");

    for ( j = 0; j < sp->Nobjects[NODE]; j++ )
    {
        if ( sp->Node[j].type == OUTFALL ) continue;
        if ( NodeStats[j].timeSurcharged == 0.0 ) continue;
        t = MAX(0.01, (NodeStats[j].timeSurcharged / 3600.0));
        if ( n == 0 )
        {
            WRITE("Surcharging occurs when water rises above the top of the highest conduit.");
            fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------"
"\n                                               Max. Height   Min. Depth"
"\n                                   Hours       Above Crown    Below Rim");
    if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file,
"\n  Node                 Type      Surcharged           Feet         Feet");
    else fprintf(sp->Frpt.file,
"\n  Node                 Type      Surcharged         Meters       Meters");
    fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------");
            n = 1;
        }
        fprintf(sp->Frpt.file, "\n  %-20s", sp->Node[j].ID);
        fprintf(sp->Frpt.file, " %-9s", NodeTypeWords[sp->Node[j].type]);
        d1 = NodeStats[j].maxDepth + sp->Node[j].invertElev - sp->Node[j].crownElev;
        if ( d1 < 0.0 ) d1 = 0.0;
        d2 = sp->Node[j].fullDepth - NodeStats[j].maxDepth;
        if ( d2 < 0.0 ) d2 = 0.0;
        fprintf(sp->Frpt.file, "  %9.2f      %9.3f    %9.3f",
                t, d1*UCF(sp, LENGTH), d2*UCF(sp, LENGTH));
    }
    if ( n == 0 ) WRITE("No nodes were surcharged.");
    WRITE("");
}

//=============================================================================

void writeNodeFlooding(SWMM_Project *sp)
{
    int    j, n = 0;
    int    days, hrs, mins;
    double t;

    WRITE("");
    WRITE("*********************");
    WRITE("Node Flooding Summary");
    WRITE("*********************");
    WRITE("");

    for ( j = 0; j < sp->Nobjects[NODE]; j++ )
    {
        if ( sp->Node[j].type == OUTFALL ) continue;
        if ( NodeStats[j].timeFlooded == 0.0 ) continue;
        t = MAX(0.01, (NodeStats[j].timeFlooded / 3600.0));

        if ( n == 0 )
        {
            WRITE("Flooding refers to all water that overflows a node, whether it ponds or not.");
            fprintf(sp->Frpt.file,
"\n  --------------------------------------------------------------------------"
"\n                                                             Total   Maximum"
"\n                                 Maximum   Time of Max       Flood    Ponded"
"\n                        Hours       Rate    Occurrence      Volume");
            if ( sp->RouteModel == DW ) fprintf(sp->Frpt.file, "     Depth");
            else                    fprintf(sp->Frpt.file, "    Volume");
            fprintf(sp->Frpt.file,
"\n  Node                 Flooded       %3s   days hr:min    %8s",
                FlowUnitWords[sp->FlowUnits], VolUnitsWords[sp->UnitSystem]);
            if ( sp->RouteModel == DW )      fprintf(sp->Frpt.file, "    %6s",
                                         PondingUnitsWords[sp->UnitSystem]);
            else if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file, "  1000 ft3");
            else                         fprintf(sp->Frpt.file, "   1000 m3");
            fprintf(sp->Frpt.file,
"\n  --------------------------------------------------------------------------");
            n = 1;
        }
        fprintf(sp->Frpt.file, "\n  %-20s", sp->Node[j].ID);
        fprintf(sp->Frpt.file, " %7.2f ", t);
        fprintf(sp->Frpt.file, FlowFmt, NodeStats[j].maxOverflow * UCF(sp, FLOW));
        getElapsedTime(sp, NodeStats[j].maxOverflowDate, &days, &hrs, &mins);
        fprintf(sp->Frpt.file, "   %4d  %02d:%02d", days, hrs, mins);
		fprintf(sp->Frpt.file, "%12.3f", NodeStats[j].volFlooded * Vcf);
        if ( sp->RouteModel == DW )
            fprintf(sp->Frpt.file, " %9.3f",
                (NodeStats[j].maxDepth - sp->Node[j].fullDepth) * UCF(sp, LENGTH));
        else
            fprintf(sp->Frpt.file, " %9.3f", NodeStats[j].maxPondedVol /
                                         1000.0 * UCF(sp, VOLUME));
    }

    if ( n == 0 ) WRITE("No nodes were flooded.");
    WRITE("");
}

//=============================================================================

void writeStorageVolumes(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes simulation statistics for storage units to report file.
//
{
    int    j, k, days, hrs, mins;
    double avgVol, maxVol, pctAvgVol, pctMaxVol;
    double addedVol, pctEvapLoss, pctSeepLoss;

    if ( sp->Nnodes[STORAGE] > 0 )
    {
        WRITE("");
        WRITE("**********************");
        WRITE("Storage Volume Summary");
        WRITE("**********************");
        WRITE("");

        fprintf(sp->Frpt.file,
"\n  --------------------------------------------------------------------------------------------------"
"\n                         Average     Avg  Evap Exfil       Maximum     Max    Time of Max    Maximum"  //(5.1.007)
"\n                          Volume    Pcnt  Pcnt  Pcnt        Volume    Pcnt     Occurrence    Outflow");
        if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file,
"\n  Storage Unit          1000 ft3    Full  Loss  Loss      1000 ft3    Full    days hr:min        ");
        else fprintf(sp->Frpt.file,
"\n  Storage Unit           1000 m3    Full  Loss  Loss       1000 m3    Full    days hr:min        ");
        fprintf(sp->Frpt.file, "%3s", FlowUnitWords[sp->FlowUnits]);
        fprintf(sp->Frpt.file,
"\n  --------------------------------------------------------------------------------------------------");

        for ( j = 0; j < sp->Nobjects[NODE]; j++ )
        {
            if ( sp->Node[j].type != STORAGE ) continue;
            k = sp->Node[j].subIndex;
            fprintf(sp->Frpt.file, "\n  %-20s", sp->Node[j].ID);
            avgVol = StorageStats[k].avgVol / sp->StepCount;
            maxVol = StorageStats[k].maxVol;
            pctMaxVol = 0.0;
            pctAvgVol = 0.0;
            if ( sp->Node[j].fullVolume > 0.0 )
            {
                pctAvgVol = avgVol / sp->Node[j].fullVolume * 100.0;
                pctMaxVol = maxVol / sp->Node[j].fullVolume * 100.0;
            }
            pctEvapLoss = 0.0;
            pctSeepLoss = 0.0;
            addedVol = NodeInflow[j] + StorageStats[k].initVol;
            if ( addedVol > 0.0 )
            {
                pctEvapLoss = StorageStats[k].evapLosses / addedVol * 100.0;
                pctSeepLoss = StorageStats[k].exfilLosses / addedVol * 100.0;  //(5.1.007)
            }

            fprintf(sp->Frpt.file, "%10.3f    %4.0f  %4.0f  %4.0f    %10.3f    %4.0f",
                avgVol*UCF(sp, VOLUME)/1000.0, pctAvgVol, pctEvapLoss, pctSeepLoss,
                maxVol*UCF(sp, VOLUME)/1000.0, pctMaxVol);

            getElapsedTime(sp, StorageStats[k].maxVolDate, &days, &hrs, &mins);
            fprintf(sp->Frpt.file, "    %4d  %02d:%02d  ", days, hrs, mins);
            fprintf(sp->Frpt.file, FlowFmt, StorageStats[k].maxFlow*UCF(sp, FLOW));
        }
        WRITE("");
    }
}

//=============================================================================

void writeOutfallLoads(SWMM_Project *sp)
//
//  Input:   node
//  Output:  none
//  Purpose: writes simulation statistics for outfall nodess to report file.
//
{
    char    units[15];
    int     i, j, k, p;
    double  x;
    double  outfallCount, flowCount;
    double  flowSum, freqSum, volSum;
    double* totals;

    if ( sp->Nnodes[OUTFALL] > 0 )
    {
        // --- initial totals
        totals = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        for (p=0; p<sp->Nobjects[POLLUT]; p++) totals[p] = 0.0;
        flowSum = 0.0;
        freqSum = 0.0;
		volSum  = 0.0;

        // --- print table title
        WRITE("");
        WRITE("***********************");
        WRITE("Outfall Loading Summary");
        WRITE("***********************");
        WRITE("");

        // --- print table column headers
        fprintf(sp->Frpt.file,
 "\n  -----------------------------------------------------------"); 
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "--------------");
        fprintf(sp->Frpt.file,
 "\n                         Flow       Avg       Max       Total");
        for (p=0; p<sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file,"         Total");
        fprintf(sp->Frpt.file,
 "\n                         Freq      Flow      Flow      Volume");
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%14s", Pollut[p].ID);
        fprintf(sp->Frpt.file,
 "\n  Outfall Node           Pcnt       %3s       %3s    %8s",
            FlowUnitWords[sp->FlowUnits], FlowUnitWords[sp->FlowUnits],
			VolUnitsWords[sp->UnitSystem]);
        for (p = 0; p < sp->Nobjects[POLLUT]; p++)
        {
            i = sp->UnitSystem;
            if ( Pollut[p].units == COUNT ) i = 2;
            strcpy(units, LoadUnitsWords[i]);
            fprintf(sp->Frpt.file, "%14s", units);
        }
        fprintf(sp->Frpt.file,
 "\n  -----------------------------------------------------------");
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "--------------");

        // --- identify each outfall node
        for (j=0; j<sp->Nobjects[NODE]; j++)
        {
            if ( sp->Node[j].type != OUTFALL ) continue;
            k = sp->Node[j].subIndex;
            flowCount = OutfallStats[k].totalPeriods;

            // --- print node ID, flow freq., avg. flow, max. flow & flow vol.
            fprintf(sp->Frpt.file, "\n  %-20s", sp->Node[j].ID);
            x = 100.*flowCount/(double)sp->StepCount;
            fprintf(sp->Frpt.file, "%7.2f", x);
            freqSum += x;
            if ( flowCount > 0 )
                x = OutfallStats[k].avgFlow*UCF(sp, FLOW)/flowCount;
            else
                x = 0.0;
            flowSum += x;

            fprintf(sp->Frpt.file, " ");
            fprintf(sp->Frpt.file, FlowFmt, x);
            fprintf(sp->Frpt.file, " ");
            fprintf(sp->Frpt.file, FlowFmt, OutfallStats[k].maxFlow*UCF(sp, FLOW));
            fprintf(sp->Frpt.file, "%12.3f", NodeInflow[j] * Vcf);
            volSum += NodeInflow[j];

            // --- print load of each pollutant for outfall
            for (p=0; p<sp->Nobjects[POLLUT]; p++)
            {
                x = OutfallStats[k].totalLoad[p] * LperFT3 * Pollut[p].mcf;
                totals[p] += x;
                if ( Pollut[p].units == COUNT ) x = LOG10(x);
		fprintf(sp->Frpt.file, "%14.3f", x);
            }
        }

        // --- print total outfall loads
        outfallCount = sp->Nnodes[OUTFALL];
        fprintf(sp->Frpt.file,
 "\n  -----------------------------------------------------------"); 
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "--------------");

        fprintf(sp->Frpt.file, "\n  System              %7.2f ",
            freqSum/outfallCount);
        fprintf(sp->Frpt.file, FlowFmt, flowSum);
        fprintf(sp->Frpt.file, " ");
        fprintf(sp->Frpt.file, FlowFmt, MaxOutfallFlow*UCF(sp, FLOW));
		fprintf(sp->Frpt.file, "%12.3f", volSum * Vcf);

        for (p = 0; p < sp->Nobjects[POLLUT]; p++)
        {
            x = totals[p];
            if ( Pollut[p].units == COUNT ) x = LOG10(x);
            fprintf(sp->Frpt.file, "%14.3f", x);
        }
        WRITE("");
        free(totals);
    } 
}

//=============================================================================

void writeLinkFlows(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes simulation statistics for links to report file.
//
{
    int    j, k, days, hrs, mins;
    double v, fullDepth;

    if ( sp->Nobjects[LINK] == 0 ) return;
    WRITE("");
    WRITE("********************");
    WRITE("Link Flow Summary");
    WRITE("********************");
    WRITE("");

    fprintf(sp->Frpt.file,
"\n  -----------------------------------------------------------------------------"
"\n                                 Maximum  Time of Max   Maximum    Max/    Max/"
"\n                                  |Flow|   Occurrence   |Veloc|    Full    Full");
    if ( sp->UnitSystem == US ) fprintf(sp->Frpt.file,
"\n  Link                 Type          %3s  days hr:min    ft/sec    Flow   Depth",
        FlowUnitWords[sp->FlowUnits]);
    else fprintf(sp->Frpt.file,
"\n  Link                 Type          %3s  days hr:min     m/sec    Flow   Depth",
        FlowUnitWords[sp->FlowUnits]);
    fprintf(sp->Frpt.file,
"\n  -----------------------------------------------------------------------------");

    for ( j = 0; j < sp->Nobjects[LINK]; j++ )
    {
        // --- print link ID
        k = Link[j].subIndex;
        fprintf(sp->Frpt.file, "\n  %-20s", Link[j].ID);

        // --- print link type
        if ( Link[j].xsect.type == DUMMY ) fprintf(sp->Frpt.file, " DUMMY   ");
        else if ( Link[j].xsect.type == IRREGULAR ) fprintf(sp->Frpt.file, " CHANNEL ");
        else fprintf(sp->Frpt.file, " %-7s ", LinkTypeWords[Link[j].type]);

        // --- print max. flow & time of occurrence
        getElapsedTime(sp, LinkStats[j].maxFlowDate, &days, &hrs, &mins);
        fprintf(sp->Frpt.file, FlowFmt, LinkStats[j].maxFlow*UCF(sp, FLOW));
        fprintf(sp->Frpt.file, "  %4d  %02d:%02d", days, hrs, mins);

        // --- print max flow / flow capacity for pumps
        if ( Link[j].type == PUMP && Link[j].qFull > 0.0)
        {
            fprintf(sp->Frpt.file, "          ");
            fprintf(sp->Frpt.file, "  %6.2f",
                LinkStats[j].maxFlow / Link[j].qFull);
            continue;
        }

        // --- stop printing for dummy conduits
        if ( Link[j].xsect.type == DUMMY ) continue;

        // --- stop printing for outlet links (since they don't have xsections)
        if ( Link[j].type == OUTLET ) continue;

        // --- print max velocity & max/full flow for conduits
        if ( Link[j].type == CONDUIT )
        {
            v = LinkStats[j].maxVeloc*UCF(sp, LENGTH);
            if ( v > 50.0 ) fprintf(sp->Frpt.file, "    >50.00");
            else fprintf(sp->Frpt.file, "   %7.2f", v);
            fprintf(sp->Frpt.file, "  %6.2f", LinkStats[j].maxFlow / Link[j].qFull /
                                          (double)Conduit[k].barrels);
        }
        else fprintf(sp->Frpt.file, "                  ");

        // --- print max/full depth
        fullDepth = Link[j].xsect.yFull;
        if ( Link[j].type == ORIFICE &&
             Orifice[k].type == BOTTOM_ORIFICE ) fullDepth = 0.0;
        if ( fullDepth > 0.0 )
        {
            fprintf(sp->Frpt.file, "  %6.2f", LinkStats[j].maxDepth / fullDepth);
        }
        else fprintf(sp->Frpt.file, "        ");
    }
    WRITE("");
}

//=============================================================================

void writeFlowClass(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes flow classification fro each conduit to report file.
//
{
    int   i, j, k;

    if ( sp->RouteModel != DW ) return;
    WRITE("");
    WRITE("***************************");
    WRITE("Flow Classification Summary");
    WRITE("***************************");
    WRITE("");
    fprintf(sp->Frpt.file,
"\n  -------------------------------------------------------------------------------------"
"\n                      Adjusted    ---------- Fraction of Time in Flow Class ---------- "
"\n                       /Actual         Up    Down  Sub   Sup   Up    Down  Norm  Inlet "
"\n  Conduit               Length    Dry  Dry   Dry   Crit  Crit  Crit  Crit  Ltd   Ctrl  "
"\n  -------------------------------------------------------------------------------------");
    for ( j = 0; j < sp->Nobjects[LINK]; j++ )
    {
        if ( Link[j].type != CONDUIT ) continue;
        if ( Link[j].xsect.type == DUMMY ) continue;
        k = Link[j].subIndex;
        fprintf(sp->Frpt.file, "\n  %-20s", Link[j].ID);
        fprintf(sp->Frpt.file, "  %6.2f ", Conduit[k].modLength / Conduit[k].length);
        for ( i=0; i<MAX_FLOW_CLASSES; i++ )
        {
            fprintf(sp->Frpt.file, "  %4.2f",
                LinkStats[j].timeInFlowClass[i] /= sp->StepCount);
        }
        fprintf(sp->Frpt.file, "  %4.2f", LinkStats[j].timeNormalFlow /
                                      (sp->NewRoutingTime/1000.0));
        fprintf(sp->Frpt.file, "  %4.2f", LinkStats[j].timeInletControl /
                                      (sp->NewRoutingTime/1000.0)); 
    }
    WRITE("");
}

//=============================================================================

void writeLinkSurcharge(SWMM_Project *sp)
{
    int    i, j, n = 0;
    double t[5];

    WRITE("");
    WRITE("*************************");
    WRITE("Conduit Surcharge Summary");
    WRITE("*************************");
    WRITE("");
    for ( j = 0; j < sp->Nobjects[LINK]; j++ )
    {
        if ( Link[j].type != CONDUIT ||
			 Link[j].xsect.type == DUMMY ) continue; 
        t[0] = LinkStats[j].timeSurcharged / 3600.0;
        t[1] = LinkStats[j].timeFullUpstream / 3600.0;
        t[2] = LinkStats[j].timeFullDnstream / 3600.0;
        t[3] = LinkStats[j].timeFullFlow / 3600.0;
        if ( t[0] + t[1] + t[2] + t[3] == 0.0 ) continue;
        t[4] = LinkStats[j].timeCapacityLimited / 3600.0;
        for (i=0; i<5; i++) t[i] = MAX(0.01, t[i]);
        if (n == 0)
        {
            fprintf(sp->Frpt.file,
"\n  ----------------------------------------------------------------------------"
"\n                                                           Hours        Hours "
"\n                         --------- Hours Full --------   Above Full   Capacity"
"\n  Conduit                Both Ends  Upstream  Dnstream   Normal Flow   Limited"
"\n  ----------------------------------------------------------------------------");
            n = 1;
        }
        fprintf(sp->Frpt.file, "\n  %-20s", Link[j].ID);
        fprintf(sp->Frpt.file, "    %8.2f  %8.2f  %8.2f  %8.2f     %8.2f",
                t[0], t[1], t[2], t[3], t[4]);
    }
    if ( n == 0 ) WRITE("No conduits were surcharged.");
    WRITE("");
}

//=============================================================================

void writePumpFlows(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes simulation statistics for pumps to report file.
//
{
    int    j, k;
    double avgFlow, pctUtilized, pctOffCurve1, pctOffCurve2, totalSeconds;

    if ( sp->Nlinks[PUMP] == 0 ) return;

    WRITE("");
    WRITE("***************");
    WRITE("Pumping Summary");
    WRITE("***************");
    WRITE("");

    fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------------------------------------------"
"\n                                                  Min       Avg       Max     Total     Power    %% Time Off"
"\n                        Percent   Number of      Flow      Flow      Flow    Volume     Usage    Pump Curve"
"\n  Pump                 Utilized   Start-Ups       %3s       %3s       %3s  %8s     Kw-hr    Low   High"
"\n  ---------------------------------------------------------------------------------------------------------",
        FlowUnitWords[sp->FlowUnits], FlowUnitWords[sp->FlowUnits],
        FlowUnitWords[sp->FlowUnits], VolUnitsWords[sp->UnitSystem]);
    for ( j = 0; j < sp->Nobjects[LINK]; j++ )
    {
        if ( Link[j].type != PUMP ) continue;
        k = Link[j].subIndex;
        fprintf(sp->Frpt.file, "\n  %-20s", Link[j].ID);
        totalSeconds = sp->NewRoutingTime / 1000.0;
        pctUtilized = PumpStats[k].utilized / totalSeconds * 100.0;
        avgFlow = PumpStats[k].avgFlow;
        if ( PumpStats[k].totalPeriods > 0 )
            avgFlow /=  PumpStats[k].totalPeriods;
        fprintf(sp->Frpt.file, " %8.2f  %10d %9.2f %9.2f %9.2f %9.3f %9.2f",
            pctUtilized, PumpStats[k].startUps, PumpStats[k].minFlow*UCF(sp, FLOW),
            avgFlow*UCF(sp, FLOW), PumpStats[k].maxFlow*UCF(sp, FLOW), 
            PumpStats[k].volume*Vcf, PumpStats[k].energy);
        pctOffCurve1 = PumpStats[k].offCurveLow;
        pctOffCurve2 = PumpStats[k].offCurveHigh;
        if ( PumpStats[k].utilized > 0.0 )
        {
            pctOffCurve1 = pctOffCurve1 / PumpStats[k].utilized * 100.0;
            pctOffCurve2 = pctOffCurve2 / PumpStats[k].utilized * 100.0;
        }
        fprintf(sp->Frpt.file, " %6.1f %6.1f", pctOffCurve1, pctOffCurve2);
    }
    WRITE("");
}

//=============================================================================

void writeLinkLoads(SWMM_Project *sp)
{
    int i, j, p;
    double x;
    char  units[15];
    char  linkLine[] = "--------------------";
    char  pollutLine[]   = "--------------";

    // --- print the table headings 
    WRITE("");
    WRITE("***************************");
    WRITE("Link Pollutant Load Summary");
    WRITE("***************************");
    WRITE("");
    fprintf(sp->Frpt.file, "\n  %s", linkLine);
    for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%s", pollutLine);
    fprintf(sp->Frpt.file, "\n                      ");
    for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%14s", Pollut[p].ID);
    fprintf(sp->Frpt.file, "\n  Link                ");
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        i = sp->UnitSystem;
        if ( Pollut[p].units == COUNT ) i = 2;
        strcpy(units, LoadUnitsWords[i]);
        fprintf(sp->Frpt.file, "%14s", units);
    }
    fprintf(sp->Frpt.file, "\n  %s", linkLine);
    for (p = 0; p < sp->Nobjects[POLLUT]; p++) fprintf(sp->Frpt.file, "%s", pollutLine);

    // --- print the pollutant loadings carried by each link
    for ( j = 0; j < sp->Nobjects[LINK]; j++ )
    {
        fprintf(sp->Frpt.file, "\n  %-20s", Link[j].ID);
        for (p = 0; p < sp->Nobjects[POLLUT]; p++)
        {
            x = Link[j].totalLoad[p] * LperFT3 * Pollut[p].mcf;
            if ( Pollut[p].units == COUNT ) x = LOG10(x);
            if ( x < 10000. ) fprintf(sp->Frpt.file, "%14.3f", x);
            else fprintf(sp->Frpt.file, "%14.3e", x);
        }
    }
    WRITE("");
}
