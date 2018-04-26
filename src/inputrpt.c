//-----------------------------------------------------------------------------
//   inputrpt.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14 (Build 5.1.001)
//   Author:   L. Rossman
//
//   Report writing functions for input data summary.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <time.h>
#include "headers.h"
#include "lid.h"

#define WRITE(x) (report_writeLine(sp,(x)))

//=============================================================================

void inputrpt_writeInput(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: writes summary of input data to report file.
//
{
    int m;
    int i, k;
    int lidCount = 0;

    if ( sp->ErrorCode ) return;

    WRITE("");
    WRITE("*************");
    WRITE("Element Count");
    WRITE("*************");
    fprintf(sp->Frpt.file, "\n  Number of rain gages ...... %d", sp->Nobjects[GAGE]);
    fprintf(sp->Frpt.file, "\n  Number of subcatchments ... %d", sp->Nobjects[SUBCATCH]);
    fprintf(sp->Frpt.file, "\n  Number of nodes ........... %d", sp->Nobjects[NODE]);
    fprintf(sp->Frpt.file, "\n  Number of links ........... %d", sp->Nobjects[LINK]);
    fprintf(sp->Frpt.file, "\n  Number of pollutants ...... %d", sp->Nobjects[POLLUT]);
    fprintf(sp->Frpt.file, "\n  Number of land uses ....... %d", sp->Nobjects[LANDUSE]);

    if ( sp->Nobjects[POLLUT] > 0 )
    {
        WRITE("");
        WRITE("");
        WRITE("*****************");
        WRITE("Pollutant Summary");
        WRITE("*****************");
        fprintf(sp->Frpt.file,
    "\n                               Ppt.      GW         Kdecay");
        fprintf(sp->Frpt.file,
    "\n  Name                 Units   Concen.   Concen.    1/days    CoPollutant");
        fprintf(sp->Frpt.file,
    "\n  -----------------------------------------------------------------------");
        for (i = 0; i < sp->Nobjects[POLLUT]; i++)
        {
            fprintf(sp->Frpt.file, "\n  %-20s %5s%10.2f%10.2f%10.2f", Pollut[i].ID,
                QualUnitsWords[Pollut[i].units], Pollut[i].pptConcen,
                Pollut[i].gwConcen, Pollut[i].kDecay*SECperDAY);
            if ( Pollut[i].coPollut >= 0 )
                fprintf(sp->Frpt.file, "    %-s  (%.2f)",
                    Pollut[Pollut[i].coPollut].ID, Pollut[i].coFraction);
        }
    }

    if ( sp->Nobjects[LANDUSE] > 0 )
    {
        WRITE("");
        WRITE("");
        WRITE("***************");
        WRITE("Landuse Summary");
        WRITE("***************");
        fprintf(sp->Frpt.file,
    "\n                         Sweeping   Maximum      Last");
        fprintf(sp->Frpt.file,
    "\n  Name                   Interval   Removal     Swept");
        fprintf(sp->Frpt.file,
    "\n  ---------------------------------------------------");
        for (i=0; i<sp->Nobjects[LANDUSE]; i++)
        {
            fprintf(sp->Frpt.file, "\n  %-20s %10.2f%10.2f%10.2f", Landuse[i].ID,
                Landuse[i].sweepInterval, Landuse[i].sweepRemoval,
                Landuse[i].sweepDays0);
        }
    }

    if ( sp->Nobjects[GAGE] > 0 )
    {
        WRITE("");
        WRITE("");
        WRITE("****************");
        WRITE("Raingage Summary");
        WRITE("****************");
    fprintf(sp->Frpt.file,
"\n                                                      Data       Recording");
    fprintf(sp->Frpt.file,
"\n  Name                 Data Source                    Type       Interval ");
    fprintf(sp->Frpt.file,
"\n  ------------------------------------------------------------------------");
        for (i = 0; i < sp->Nobjects[GAGE]; i++)
        {
            if ( sp->Gage[i].tSeries >= 0 )
            {
                fprintf(sp->Frpt.file, "\n  %-20s %-30s ",
                    sp->Gage[i].ID, Tseries[sp->Gage[i].tSeries].ID);
                fprintf(sp->Frpt.file, "%-10s %3d min.",
                    RainTypeWords[sp->Gage[i].rainType],
                    (sp->Gage[i].rainInterval)/60);
            }
            else fprintf(sp->Frpt.file, "\n  %-20s %-30s",
                sp->Gage[i].ID, sp->Gage[i].fname);
        }
    }

    if ( sp->Nobjects[SUBCATCH] > 0 )
    {
        WRITE("");
        WRITE("");
        WRITE("********************");
        WRITE("Subcatchment Summary");
        WRITE("********************");
        fprintf(sp->Frpt.file,
"\n  Name                       Area     Width   %%Imperv    %%Slope Rain Gage            Outlet              ");
        fprintf(sp->Frpt.file,
"\n  -----------------------------------------------------------------------------------------------------------");
        for (i = 0; i < sp->Nobjects[SUBCATCH]; i++)
        {
            fprintf(sp->Frpt.file,"\n  %-20s %10.2f%10.2f%10.2f%10.4f %-20s ",
                sp->Subcatch[i].ID, sp->Subcatch[i].area*UCF(sp, LANDAREA),
                sp->Subcatch[i].width*UCF(sp, LENGTH),  sp->Subcatch[i].fracImperv*100.0,
                sp->Subcatch[i].slope*100.0, sp->Gage[sp->Subcatch[i].gage].ID);
            if ( sp->Subcatch[i].outNode >= 0 )
            {
                fprintf(sp->Frpt.file, "%-20s", sp->Node[sp->Subcatch[i].outNode].ID);
            }
            else if ( sp->Subcatch[i].outSubcatch >= 0 )
            {
                fprintf(sp->Frpt.file, "%-20s", sp->Subcatch[sp->Subcatch[i].outSubcatch].ID);
            }
            if ( sp->Subcatch[i].lidArea ) lidCount++;
        }
    }
    if ( lidCount > 0 ) lid_writeSummary(sp);

    if ( sp->Nobjects[NODE] > 0 )
    {
        WRITE("");
        WRITE("");
        WRITE("************");
        WRITE("Node Summary");
        WRITE("************");
        fprintf(sp->Frpt.file,
"\n                                           Invert      Max.    Ponded    External");
        fprintf(sp->Frpt.file,
"\n  Name                 Type                 Elev.     Depth      Area    Inflow  ");
        fprintf(sp->Frpt.file,
"\n  -------------------------------------------------------------------------------");
        for (i = 0; i < sp->Nobjects[NODE]; i++)
        {
            fprintf(sp->Frpt.file, "\n  %-20s %-16s%10.2f%10.2f%10.1f", sp->Node[i].ID,
                NodeTypeWords[sp->Node[i].type-JUNCTION],
                sp->Node[i].invertElev*UCF(sp, LENGTH),
                sp->Node[i].fullDepth*UCF(sp, LENGTH),
                sp->Node[i].pondedArea*UCF(sp, LENGTH)*UCF(sp, LENGTH));
            if ( sp->Node[i].extInflow || sp->Node[i].dwfInflow || sp->Node[i].rdiiInflow )
            {
                fprintf(sp->Frpt.file, "    Yes");
            }
        }
    }

    if ( sp->Nobjects[LINK] > 0 )
    {
        WRITE("");
        WRITE("");
        WRITE("************");
        WRITE("Link Summary");
        WRITE("************");
        fprintf(sp->Frpt.file,
"\n  Name             From Node        To Node          Type            Length    %%Slope Roughness");
        fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------------------------------");
        for (i = 0; i < sp->Nobjects[LINK]; i++)
        {
            // --- list end nodes in their original orientation
            if ( Link[i].direction == 1 )
                fprintf(sp->Frpt.file, "\n  %-16s %-16s %-16s ",
                    Link[i].ID, sp->Node[Link[i].node1].ID, sp->Node[Link[i].node2].ID);
            else
                fprintf(sp->Frpt.file, "\n  %-16s %-16s %-16s ",
                    Link[i].ID, sp->Node[Link[i].node2].ID, sp->Node[Link[i].node1].ID);

            // --- list link type
            if ( Link[i].type == PUMP )
            {
                k = Link[i].subIndex;
                fprintf(sp->Frpt.file, "%-5s PUMP  ",
                    PumpTypeWords[Pump[k].type]);
            }
            else fprintf(sp->Frpt.file, "%-12s",
                LinkTypeWords[Link[i].type-CONDUIT]);

            // --- list length, slope and roughness for conduit links
            if (Link[i].type == CONDUIT)
            {
                k = Link[i].subIndex;
                fprintf(sp->Frpt.file, "%10.1f%10.4f%10.4f",
                    Conduit[k].length*UCF(sp, LENGTH),
                    Conduit[k].slope*100.0*Link[i].direction,
                    Conduit[k].roughness);
            }
        }

        WRITE("");
        WRITE("");
        WRITE("*********************");
        WRITE("Cross Section Summary");
        WRITE("*********************");
        fprintf(sp->Frpt.file,
"\n                                        Full     Full     Hyd.     Max.   No. of     Full");
        fprintf(sp->Frpt.file,
"\n  Conduit          Shape               Depth     Area     Rad.    Width  Barrels     Flow");
        fprintf(sp->Frpt.file,
"\n  ---------------------------------------------------------------------------------------");
        for (i = 0; i < sp->Nobjects[LINK]; i++)
        {
            if (Link[i].type == CONDUIT)
            {
                k = Link[i].subIndex;
                fprintf(sp->Frpt.file, "\n  %-16s ", Link[i].ID);
                if ( Link[i].xsect.type == CUSTOM )
                    fprintf(sp->Frpt.file, "%-16s ", Curve[Link[i].xsect.transect].ID);
                else if ( Link[i].xsect.type == IRREGULAR )
                    fprintf(sp->Frpt.file, "%-16s ",
                    Transect[Link[i].xsect.transect].ID);
                else fprintf(sp->Frpt.file, "%-16s ",
                    XsectTypeWords[Link[i].xsect.type]);
                fprintf(sp->Frpt.file, "%8.2f %8.2f %8.2f %8.2f      %3d %8.2f",
                    Link[i].xsect.yFull*UCF(sp, LENGTH),
                    Link[i].xsect.aFull*UCF(sp, LENGTH)*UCF(sp, LENGTH),
                    Link[i].xsect.rFull*UCF(sp, LENGTH),
                    Link[i].xsect.wMax*UCF(sp, LENGTH),
                    Conduit[k].barrels,
                    Link[i].qFull*UCF(sp, FLOW));
            }
        }
    }

    if (sp->Nobjects[SHAPE] > 0)
    {
        WRITE("");
        WRITE("");
        WRITE("*************");
        WRITE("Shape Summary");
        WRITE("*************");
        for (i = 0; i < sp->Nobjects[SHAPE]; i++)
        {
            k = Shape[i].curve;
            fprintf(sp->Frpt.file, "\n\n  Shape %s", Curve[k].ID);
            fprintf(sp->Frpt.file, "\n  Area:  ");
            for ( m = 1; m < N_SHAPE_TBL; m++)
            {
                 if ( m % 5 == 1 ) fprintf(sp->Frpt.file,"\n          ");
                 fprintf(sp->Frpt.file, "%10.4f ", Shape[i].areaTbl[m]);
            }
            fprintf(sp->Frpt.file, "\n  Hrad:  ");
            for ( m = 1; m < N_SHAPE_TBL; m++)
            {
                 if ( m % 5 == 1 ) fprintf(sp->Frpt.file,"\n          ");
                 fprintf(sp->Frpt.file, "%10.4f ", Shape[i].hradTbl[m]);
            }
            fprintf(sp->Frpt.file, "\n  Width: ");
            for ( m = 1; m < N_SHAPE_TBL; m++)
            {
                 if ( m % 5 == 1 ) fprintf(sp->Frpt.file,"\n          ");
                 fprintf(sp->Frpt.file, "%10.4f ", Shape[i].widthTbl[m]);
            }
        }
    }
    WRITE("");

    if (sp->Nobjects[TRANSECT] > 0)
    {
        WRITE("");
        WRITE("");
        WRITE("****************");
        WRITE("Transect Summary");
        WRITE("****************");
        for (i = 0; i < sp->Nobjects[TRANSECT]; i++)
        {
            fprintf(sp->Frpt.file, "\n\n  Transect %s", Transect[i].ID);
            fprintf(sp->Frpt.file, "\n  Area:  ");
            for ( m = 1; m < N_TRANSECT_TBL; m++)
            {
                 if ( m % 5 == 1 ) fprintf(sp->Frpt.file,"\n          ");
                 fprintf(sp->Frpt.file, "%10.4f ", Transect[i].areaTbl[m]);
            }
            fprintf(sp->Frpt.file, "\n  Hrad:  ");
            for ( m = 1; m < N_TRANSECT_TBL; m++)
            {
                 if ( m % 5 == 1 ) fprintf(sp->Frpt.file,"\n          ");
                 fprintf(sp->Frpt.file, "%10.4f ", Transect[i].hradTbl[m]);
            }
            fprintf(sp->Frpt.file, "\n  Width: ");
            for ( m = 1; m < N_TRANSECT_TBL; m++)
            {
                 if ( m % 5 == 1 ) fprintf(sp->Frpt.file,"\n          ");
                 fprintf(sp->Frpt.file, "%10.4f ", Transect[i].widthTbl[m]);
            }
        }
    }
    WRITE("");
}
