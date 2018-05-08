//-----------------------------------------------------------------------------
//   link.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             09/15/14   (Build 5.1.007)
//             03/19/15   (Build 5.1.008)
//             08/05/15   (Build 5.1.010)
//             08/01/16   (Build 5.1.011)
//             03/14/17   (Build 5.1.012)
//   Author:   L. Rossman (EPA)
//             M. Tryby (EPA)
//
//   Conveyance system link functions
//
//   Build 5.1.007:
//   - Optional surcharging of weirs introduced.
//
//   Build 5.1.008:
//   - Bug in finding flow through surcharged weir fixed.
//   - Bug in finding if conduit is upstrm/dnstrm full fixed.
//   - Monthly conductivity adjustment applied to conduit seepage.
//   - Conduit seepage limited by conduit's flow rate.
//
//   Build 5.1.010:
//   - Support added for new ROADWAY_WEIR object.
//   - Time of last setting change initialized for links.
//
//   Build 5.1.011:
//   - Crest elevation of regulator links raised to downstream invert.
//   - Fixed converting roadWidth weir parameter to internal units.
//   - Weir shape parameter deprecated.
//   - Extra geometric parameters ignored for non-conduit open rectangular
//     cross sections.
//
//   Build 5.1.012:
//   - Conduit seepage rate now based on flow width, not wetted perimeter.
//   - Formula for side flow weir corrected.
//   - Crest length contraction adjustments corrected.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const double MIN_DELTA_Z = 0.001; // minimum elevation change for conduit
                                         // slopes (ft)

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  link_readParams        (called by parseLine in input.c)
//  link_readXsectParams   (called by parseLine in input.c)
//  link_readLossParams    (called by parseLine in input.c)
//  link_validate          (called by project_validate in project.c)
//  link_initState         (called by initObjects in swmm5.c)
//  link_setOldHydState    (called by routing_execute in routing.c)
//  link_setOldQualState   (called by routing_execute in routing.c)
//  link_setTargetSetting  (called by routing_execute in routing.c)
//  link_setSetting        (called by routing_execute in routing.c)
//  link_getResults        (called by output_saveLinkResults)
//  link_getLength         (called in dwflow.c, kinwave.c & flowrout.c)
//  link_getFroude         (called in dwflow.c)
//  link_getInflow         (called in flowrout.c & dynwave.c)
//  link_setOutfallDepth   (called in flowrout.c & dynwave.c)
//  link_getYcrit          (called by link_setOutfallDepth & in dwflow.c)
//  link_getYnorm          (called by conduit_initState, link_setOutfallDepth & in dwflow.c)
//  link_getVelocity       (called by link_getResults & stats_updateLinkStats)
//  link_getPower          (called by stats_updateLinkStats in stats.c)
//  link_getLossRate       (called in dwflow.c, kinwave.c & flowrout.c)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void   link_setParams(SWMM_Project *sp, int j, int type, int n1, int n2,
        int k, double x[]);
static void   link_convertOffsets(SWMM_Project *sp, int j);
static double link_getOffsetHeight(SWMM_Project *sp, int j, double offset,
        double elev);

static int    conduit_readParams(SWMM_Project *sp, int j, int k, char* tok[],
        int ntoks);
static void   conduit_validate(SWMM_Project *sp, int j, int k);
static void   conduit_initState(SWMM_Project *sp, int j, int k);
static void   conduit_reverse(SWMM_Project *sp, int j, int k);
static double conduit_getLength(SWMM_Project *sp, int j);
static double conduit_getLengthFactor(SWMM_Project *sp, int j, int k,
        double roughness);
static double conduit_getSlope(SWMM_Project *sp, int j);
static double conduit_getInflow(SWMM_Project *sp, int j);
static double conduit_getLossRate(SWMM_Project *sp, int j, double q, double tstep);              //(5.1.008)

static int    pump_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks);
static void   pump_validate(SWMM_Project *sp, int j, int k);
static void   pump_initState(SWMM_Project *sp, int j, int k);
static double pump_getInflow(SWMM_Project *sp, int j);

static int    orifice_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks);
static void   orifice_validate(SWMM_Project *sp, int j, int k);
static void   orifice_setSetting(SWMM_Project *sp, int j, double tstep);
static double orifice_getWeirCoeff(SWMM_Project *sp, int j, int k, double h);
static double orifice_getInflow(SWMM_Project *sp, int j);
static double orifice_getFlow(SWMM_Project *sp, int j, int k, double head,
        double f, int hasFlapGate);

static int    weir_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks);
static void   weir_validate(SWMM_Project *sp, int j, int k);
static void   weir_setSetting(SWMM_Project *sp, int j);                                          //(5.1.007)
static double weir_getInflow(SWMM_Project *sp, int j);
static double weir_getOpenArea(SWMM_Project *sp, int j, double y);
static void   weir_getFlow(SWMM_Project *sp, int j, int k, double head, double dir,
              int hasFlapGate, double* q1, double* q2);
static double weir_getOrificeFlow(SWMM_Project *sp, int j, double head, double y,
        double cOrif); //(5.1.007)
static double weir_getdqdh(SWMM_Project *sp, int k, double dir, double h,
        double q1, double q2);

static int    outlet_readParams(SWMM_Project *sp, int j, int k, char* tok[],
        int ntoks);
static double outlet_getFlow(SWMM_Project *sp, int k, double head);
static double outlet_getInflow(SWMM_Project *sp, int j);


//=============================================================================

int link_readParams(SWMM_Project *sp, int j, int type, int k, char* tok[], int ntoks)
//
//  Input:   j     = link index
//           type  = link type code
//           k     = link type index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads parameters for a specific type of link from a
//           tokenized line of input data.
//
{
    switch ( type )
    {
      case CONDUIT: return conduit_readParams(sp, j, k, tok, ntoks);
      case PUMP:    return pump_readParams(sp, j, k, tok, ntoks);
      case ORIFICE: return orifice_readParams(sp, j, k, tok, ntoks);
      case WEIR:    return weir_readParams(sp, j, k, tok, ntoks);
      case OUTLET:  return outlet_readParams(sp, j, k, tok, ntoks);
      default: return 0;
    }
}

//=============================================================================

int link_readXsectParams(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads a link's cross section parameters from a tokenized
//           line of input data.
//
{
    int    i, j, k;
    double x[4];

    // --- get index of link
    if ( ntoks < 6 ) return error_setInpError(sp, ERR_ITEMS, "");
    j = project_findObject(sp, LINK, tok[0]);
    if ( j < 0 ) return error_setInpError(sp, ERR_NAME, tok[0]);

    // --- get code of xsection shape
    k = findmatch(tok[1], XsectTypeWords);
    if ( k < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[1]);

    // --- assign default number of barrels to conduit
    if ( sp->Link[j].type == CONDUIT ) sp->Conduit[sp->Link[j].subIndex].barrels = 1;

    // --- assume link is not a culvert
    sp->Link[j].xsect.culvertCode = 0;

    // --- for irregular shape, find index of transect object
    if ( k == IRREGULAR )
    {
        i = project_findObject(sp, TRANSECT, tok[2]);
        if ( i < 0 ) return error_setInpError(sp, ERR_NAME, tok[2]);
        sp->Link[j].xsect.type = k;
        sp->Link[j].xsect.transect = i;
    }
    else
    {
        // --- parse max. depth & shape curve for a custom shape
        if ( k == CUSTOM )
        {
            if ( !getDouble(tok[2], &x[0]) || x[0] <= 0.0 )
               return error_setInpError(sp, ERR_NUMBER, tok[2]);
            i = project_findObject(sp, CURVE, tok[3]);
            if ( i < 0 ) return error_setInpError(sp, ERR_NAME, tok[3]);
            sp->Link[j].xsect.type = k;
            sp->Link[j].xsect.transect = i;
            sp->Link[j].xsect.yFull = x[0] / UCF(sp, LENGTH);
        }

        // --- parse and save geometric parameters
        else for (i = 2; i <= 5; i++)
        {
            if ( !getDouble(tok[i], &x[i-2]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i]);
        }

////  Following code segment added to release 5.1.011.  ////                   //(5.1.011)
        // --- ignore extra parameters for non-conduit open rectangular shapes 
        if ( sp->Link[j].type != CONDUIT && k == RECT_OPEN )
        {
            x[2] = 0.0;
            x[3] = 0.0;
        }
////
        if ( !xsect_setParams(sp, &sp->Link[j].xsect, k, x, UCF(sp, LENGTH)) )
        {
            return error_setInpError(sp, ERR_NUMBER, "");
        }

        // --- parse number of barrels if present
        if ( sp->Link[j].type == CONDUIT && ntoks >= 7 )
        {
            i = atoi(tok[6]);
            if ( i <= 0 ) return error_setInpError(sp, ERR_NUMBER, tok[6]);
            else sp->Conduit[sp->Link[j].subIndex].barrels = (char)i;
        }

        // --- parse culvert code if present
        if ( sp->Link[j].type == CONDUIT && ntoks >= 8 )
        {
            i = atoi(tok[7]);
            if ( i < 0 ) return error_setInpError(sp, ERR_NUMBER, tok[7]);
            else sp->Link[j].xsect.culvertCode = i;
        }

    }
    return 0;
}

//=============================================================================

int link_readLossParams(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads local loss parameters for a link from a tokenized
//           line of input data.
//
//  Format:  LinkID  cInlet  cOutlet  cAvg  FlapGate(YES/NO)  SeepRate         //(5.1.007)
//
{
    int    i, j, k;
    double x[3];
    double seepRate = 0.0;

    if ( ntoks < 4 ) return error_setInpError(sp, ERR_ITEMS, "");
    j = project_findObject(sp, LINK, tok[0]);
    if ( j < 0 ) return error_setInpError(sp, ERR_NAME, tok[0]);
    for (i=1; i<=3; i++)
    {
        if ( ! getDouble(tok[i], &x[i-1]) || x[i-1] < 0.0 )
        return error_setInpError(sp, ERR_NUMBER, tok[i]);
    }
    k = 0;
    if ( ntoks >= 5 )
    {
        k = findmatch(tok[4], NoYesWords);
        if ( k < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[4]);
    }
    if ( ntoks >= 6 )
    {
        if ( ! getDouble(tok[5], &seepRate) )
        return error_setInpError(sp, ERR_NUMBER, tok[5]);
    }
    sp->Link[j].cLossInlet   = x[0];
    sp->Link[j].cLossOutlet  = x[1];
    sp->Link[j].cLossAvg     = x[2];
    sp->Link[j].hasFlapGate  = k;
    sp->Link[j].seepRate     = seepRate / UCF(sp, RAINFALL);
    return 0;
}

//=============================================================================

void  link_setParams(SWMM_Project *sp, int j, int type, int n1, int n2, int k,
        double x[])
//
//  Input:   j   = link index
//           type = link type code
//           n1   = index of upstream node
//           n2   = index of downstream node
//           k    = index of link's sub-type
//           x    = array of parameter values
//  Output:  none
//  Purpose: sets parameters for a link.
//
{
    sp->Link[j].node1       = n1;
    sp->Link[j].node2       = n2;
    sp->Link[j].type        = type;
    sp->Link[j].subIndex    = k;
    sp->Link[j].offset1     = 0.0;
    sp->Link[j].offset2     = 0.0;
    sp->Link[j].q0          = 0.0;
    sp->Link[j].qFull       = 0.0;
    sp->Link[j].setting     = 1.0;
    sp->Link[j].targetSetting = 1.0;
    sp->Link[j].hasFlapGate = 0;
    sp->Link[j].qLimit      = 0.0;         // 0 means that no limit is defined
    sp->Link[j].direction   = 1;

    switch (type)
    {
      case CONDUIT:
        sp->Conduit[k].length    = x[0] / UCF(sp, LENGTH);
        sp->Conduit[k].modLength = sp->Conduit[k].length;
        sp->Conduit[k].roughness = x[1];
        sp->Link[j].offset1      = x[2] / UCF(sp, LENGTH);
        sp->Link[j].offset2      = x[3] / UCF(sp, LENGTH);
        sp->Link[j].q0           = x[4] / UCF(sp, FLOW);
        sp->Link[j].qLimit       = x[5] / UCF(sp, FLOW);
        break;

      case PUMP:
        sp->Pump[k].pumpCurve    = (int)x[0];
        sp->Link[j].hasFlapGate  = FALSE;
        sp->Pump[k].initSetting  = x[1];
        sp->Pump[k].yOn          = x[2] / UCF(sp, LENGTH);
        sp->Pump[k].yOff         = x[3] / UCF(sp, LENGTH);
        sp->Pump[k].xMin         = 0.0;
        sp->Pump[k].xMax         = 0.0;
        break;

      case ORIFICE:
        sp->Orifice[k].type      = (int)x[0];
        sp->Link[j].offset1      = x[1] / UCF(sp, LENGTH);
        sp->Link[j].offset2      = sp->Link[j].offset1;
        sp->Orifice[k].cDisch    = x[2];
        sp->Link[j].hasFlapGate  = (x[3] > 0.0) ? 1 : 0;
        sp->Orifice[k].orate     = x[4] * 3600.0;
        break;

      case WEIR:
        sp->Weir[k].type         = (int)x[0];
        sp->Link[j].offset1      = x[1] / UCF(sp, LENGTH);
        sp->Link[j].offset2      = sp->Link[j].offset1;
        sp->Weir[k].cDisch1      = x[2];
        sp->Link[j].hasFlapGate  = (x[3] > 0.0) ? 1 : 0;
        sp->Weir[k].endCon       = x[4];
        sp->Weir[k].cDisch2      = x[5];
        sp->Weir[k].canSurcharge = (int)x[6];                                      //(5.1.007)
        sp->Weir[k].roadWidth    = x[7] / UCF(sp, LENGTH);                             //(5.1.011)
        sp->Weir[k].roadSurface  = (int)x[8];                                      //(5.1.010)
//      sp->Weir[k].shape        = -(int)x[9];  //DELETED//                        //(5.1.011)
        break;

      case OUTLET:
        sp->Link[j].offset1      = x[0] / UCF(sp, LENGTH);
        sp->Link[j].offset2      = sp->Link[j].offset1;
        sp->Outlet[k].qCoeff     = x[1];
        sp->Outlet[k].qExpon     = x[2];
        sp->Outlet[k].qCurve     = (int)x[3];
        sp->Link[j].hasFlapGate  = (x[4] > 0.0) ? 1 : 0;
        sp->Outlet[k].curveType  = (int)x[5];

        xsect_setParams(sp, &sp->Link[j].xsect, DUMMY, NULL, 0.0);
        break;

    }
}

//=============================================================================

void  link_validate(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: validates a link's properties.
//
{
    int   n;

    if ( sp->LinkOffsets == ELEV_OFFSET ) link_convertOffsets(sp, j);
    switch ( sp->Link[j].type )
    {
      case CONDUIT: conduit_validate(sp, j, sp->Link[j].subIndex); break;
      case PUMP:    pump_validate(sp, j, sp->Link[j].subIndex);    break;
      case ORIFICE: orifice_validate(sp, j, sp->Link[j].subIndex); break;
      case WEIR:    weir_validate(sp, j, sp->Link[j].subIndex);    break;
    }

    // --- check if crest of regulator opening < invert of downstream node
    switch ( sp->Link[j].type )
    {
      case ORIFICE:
      case WEIR:
      case OUTLET:
          if ( sp->Node[sp->Link[j].node1].invertElev + sp->Link[j].offset1 <
               sp->Node[sp->Link[j].node2].invertElev )
          {
              sp->Link[j].offset1 = sp->Node[sp->Link[j].node2].invertElev -               //(5.1.011)
                                sp->Node[sp->Link[j].node1].invertElev;                //(5.1.011)
              report_writeWarningMsg(sp, WARN10, sp->Link[j].ID);
          }
    }    

    // --- force max. depth of end nodes to be >= link crown height
    //     at non-storage nodes

    // --- skip pumps and bottom orifices
    if ( sp->Link[j].type == PUMP ||
         (sp->Link[j].type == ORIFICE &&
          sp->Orifice[sp->Link[j].subIndex].type == BOTTOM_ORIFICE) ) return;

    // --- extend upstream node's full depth to link's crown elevation
    n = sp->Link[j].node1;
    if ( sp->Node[n].type != STORAGE )
    {
        sp->Node[n].fullDepth = MAX(sp->Node[n].fullDepth,
                            sp->Link[j].offset1 + sp->Link[j].xsect.yFull);
    }

    // --- do same for downstream node only for conduit links
    n = sp->Link[j].node2;
    if ( sp->Node[n].type != STORAGE && sp->Link[j].type == CONDUIT )
    {
        sp->Node[n].fullDepth = MAX(sp->Node[n].fullDepth,
                            sp->Link[j].offset2 + sp->Link[j].xsect.yFull);
    }
}

//=============================================================================

void link_convertOffsets(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: converts offset elevations to offset heights for a link.
//
{
    double elev;

    elev = sp->Node[sp->Link[j].node1].invertElev;
    sp->Link[j].offset1 = link_getOffsetHeight(sp, j, sp->Link[j].offset1, elev);
    if ( sp->Link[j].type == CONDUIT )
    {
        elev = sp->Node[sp->Link[j].node2].invertElev;
        sp->Link[j].offset2 = link_getOffsetHeight(sp, j, sp->Link[j].offset2, elev);
    }
    else sp->Link[j].offset2 = sp->Link[j].offset1;
}

//=============================================================================

double link_getOffsetHeight(SWMM_Project *sp, int j, double offset, double elev)
//
//  Input:   j = link index
//           offset = link elevation offset (ft)
//           elev   = node invert elevation (ft)
//  Output:  returns offset distance above node invert (ft)
//  Purpose: finds offset height for one end of a link.
//
{
    if ( offset <= MISSING || sp->Link[j].type == PUMP) return 0.0;
    offset -= elev;
    if ( offset >= 0.0 ) return offset;
    if ( offset >= -MIN_DELTA_Z ) return 0.0;
    report_writeWarningMsg(sp, WARN03, sp->Link[j].ID);
    return 0.0;
}

//=============================================================================

void link_initState(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: initializes a link's state variables at start of simulation.
//
{
    int   p;

    // --- initialize hydraulic state
    sp->Link[j].oldFlow   = sp->Link[j].q0;
    sp->Link[j].newFlow   = sp->Link[j].q0;
    sp->Link[j].oldDepth  = 0.0;
    sp->Link[j].newDepth  = 0.0;
    sp->Link[j].oldVolume = 0.0;
    sp->Link[j].newVolume = 0.0;
    sp->Link[j].setting   = 1.0;
    sp->Link[j].targetSetting = 1.0;
    sp->Link[j].timeLastSet = sp->StartDate;                                           //(5.1.010)
    sp->Link[j].inletControl  = FALSE;
    sp->Link[j].normalFlow    = FALSE;
    if ( sp->Link[j].type == CONDUIT ) conduit_initState(sp, j, sp->Link[j].subIndex);
    if ( sp->Link[j].type == PUMP    ) pump_initState(sp, j, sp->Link[j].subIndex);

    // --- initialize water quality state
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        sp->Link[j].oldQual[p] = 0.0;
        sp->Link[j].newQual[p] = 0.0;
		sp->Link[j].totalLoad[p] = 0.0;
    }
}

//=============================================================================

double  link_getInflow(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns link flow rate (cfs)
//  Purpose: finds total flow entering a link during current time step.
//
{
    if ( sp->Link[j].setting == 0 ) return 0.0;
    switch ( sp->Link[j].type )
    {
      case CONDUIT: return conduit_getInflow(sp, j);
      case PUMP:    return pump_getInflow(sp, j);
      case ORIFICE: return orifice_getInflow(sp, j);
      case WEIR:    return weir_getInflow(sp, j);
      case OUTLET:  return outlet_getInflow(sp, j);
      default:      return node_getOutflow(sp, sp->Link[j].node1, j);
    }
}

//=============================================================================

////  This function has been modified for release 5.1.008.  ////               //(5.1.008)

void link_setOldHydState(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: replaces link's old hydraulic state values with current ones.
//
{
    int k;

    sp->Link[j].oldDepth  = sp->Link[j].newDepth;
    sp->Link[j].oldFlow   = sp->Link[j].newFlow;
    sp->Link[j].oldVolume = sp->Link[j].newVolume;

    if ( sp->Link[j].type == CONDUIT )
    {
        k = sp->Link[j].subIndex;
        sp->Conduit[k].q1Old = sp->Conduit[k].q1;
        sp->Conduit[k].q2Old = sp->Conduit[k].q2;
    }
}

//=============================================================================

void link_setOldQualState(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: replaces link's old water quality state values with current ones.
//
{
    int p;
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        sp->Link[j].oldQual[p] = sp->Link[j].newQual[p];
        sp->Link[j].newQual[p] = 0.0;
    }
}

//=============================================================================

void link_setTargetSetting(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: updates a link's target setting.
//
{
    int k, n1;
    if ( sp->Link[j].type == PUMP )
    {
        k = sp->Link[j].subIndex;
        n1 = sp->Link[j].node1;
        sp->Link[j].targetSetting = sp->Link[j].setting;
        if ( sp->Pump[k].yOff > 0.0 &&
             sp->Link[j].setting > 0.0 &&
             sp->Node[n1].newDepth < sp->Pump[k].yOff ) sp->Link[j].targetSetting = 0.0;
        if ( sp->Pump[k].yOn > 0.0 &&
             sp->Link[j].setting == 0.0 &&
             sp->Node[n1].newDepth > sp->Pump[k].yOn )  sp->Link[j].targetSetting = 1.0;
    }
}

//=============================================================================

void link_setSetting(SWMM_Project *sp, int j, double tstep)
//
//  Input:   j = link index
//           tstep = time step over which setting is adjusted
//  Output:  none
//  Purpose: updates a link's setting as a result of a control action.
//
{
    if ( sp->Link[j].type == ORIFICE ) orifice_setSetting(sp, j, tstep);
    else if ( sp->Link[j].type == WEIR ) weir_setSetting(sp, j);                       //(5.1.007)
    else sp->Link[j].setting = sp->Link[j].targetSetting;
}

//=============================================================================

int link_setFlapGate(SWMM_Project *sp, int j, int n1, int n2, double q)
//
//  Input:   j = link index
//           n1 = index of node on upstream end of link
//           n2 = index of node on downstream end of link
//           q = signed flow value (value and units don't matter)
//  Output:  returns TRUE if there is reverse flow through a flap gate
//           associated with the link.
//  Purpose: based on the sign of the flow, determines if a flap gate
//           associated with the link should close or not.
//
{
    int    n = -1;

    // --- check for reverse flow through link's flap gate
    if ( sp->Link[j].hasFlapGate )
    {
        if ( q * (double)sp->Link[j].direction < 0.0 ) return TRUE;
    }

    // --- check for Outfall with flap gate node on inflow end of link
    if ( q < 0.0 ) n = n2;
    if ( q > 0.0 ) n = n1;
    if ( n >= 0 &&
         sp->Node[n].type == OUTFALL &&
         sp->Outfall[sp->Node[n].subIndex].hasFlapGate ) return TRUE;
    return FALSE;
}

//=============================================================================

void link_getResults(SWMM_Project *sp, int j, double f, float x[])
//
//  Input:   j = link index
//           f = time weighting factor
//  Output:  x = array of weighted results
//  Purpose: retrieves time-weighted average of old and new results for a link.
//
{
    int    p;                     // pollutant index
    double y,                     // depth
           q,                     // flow
           u,                     // velocity
           v,                     // volume
           c;                     // capacity, setting or concentration
    double f1 = 1.0 - f;

    y = f1*sp->Link[j].oldDepth + f*sp->Link[j].newDepth;
    q = f1*sp->Link[j].oldFlow + f*sp->Link[j].newFlow;
    v = f1*sp->Link[j].oldVolume + f*sp->Link[j].newVolume;
    u = link_getVelocity(sp, j, q, y);
    c = 0.0;
    if (sp->Link[j].type == CONDUIT)
    {
        if (sp->Link[j].xsect.type != DUMMY)
            c = xsect_getAofY(sp, &sp->Link[j].xsect, y) / sp->Link[j].xsect.aFull;
    }
    else c = sp->Link[j].setting;

    // --- override time weighting for pump flow between on/off states
    if (sp->Link[j].type == PUMP && sp->Link[j].oldFlow*sp->Link[j].newFlow == 0.0)
    {
        if ( f >= f1 ) q = sp->Link[j].newFlow;
        else           q = sp->Link[j].oldFlow;
    }

    y *= UCF(sp, LENGTH);
    v *= UCF(sp, VOLUME);
    q *= UCF(sp, FLOW) * (double)sp->Link[j].direction;
    u *= UCF(sp, LENGTH) * (double)sp->Link[j].direction;
    x[LINK_DEPTH]    = (float)y;
    x[LINK_FLOW]     = (float)q;
    x[LINK_VELOCITY] = (float)u;
    x[LINK_VOLUME]   = (float)v;
    x[LINK_CAPACITY] = (float)c;

    if ( !sp->IgnoreQuality ) for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        c = f1*sp->Link[j].oldQual[p] + f*sp->Link[j].newQual[p];
        x[LINK_QUAL+p] = (float)c;
    }
}

//=============================================================================

void link_setOutfallDepth(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: sets depth at outfall node connected to link j.
//
{
    int     k;                         // conduit index
    int     n;                         // outfall node index
    double  z;                         // invert offset height (ft)
    double  q;                         // flow rate (cfs)
    double  yCrit = 0.0;               // critical flow depth (ft)
    double  yNorm = 0.0;               // normal flow depth (ft)

    // --- find which end node of link is an outfall
    if ( sp->Node[sp->Link[j].node2].type == OUTFALL )
    {
        n = sp->Link[j].node2;
        z = sp->Link[j].offset2;
    }
    else if ( sp->Node[sp->Link[j].node1].type == OUTFALL )
    {
        n = sp->Link[j].node1;
        z = sp->Link[j].offset1;
    }
    else return;

    // --- find both normal & critical depth for current flow
    if ( sp->Link[j].type == CONDUIT )
    {
        k = sp->Link[j].subIndex;
        q = fabs(sp->Link[j].newFlow / sp->Conduit[k].barrels);
        yNorm = link_getYnorm(sp, j, q);
        yCrit = link_getYcrit(sp, j, q);
    }

    // --- set new depth at node
    node_setOutletDepth(sp, n, yNorm, yCrit, z);
}

//=============================================================================

double link_getYcrit(SWMM_Project *sp, int j, double q)
//
//  Input:   j = link index
//           q = link flow rate (cfs)
//  Output:  returns critical depth (ft)
//  Purpose: computes critical depth for given flow rate.
//
{
    return xsect_getYcrit(sp, &sp->Link[j].xsect, q);
}

//=============================================================================

double  link_getYnorm(SWMM_Project *sp, int j, double q)
//
//  Input:   j = link index
//           q = link flow rate (cfs)
//  Output:  returns normal depth (ft)
//  Purpose: computes normal depth for given flow rate.
//
{
    int    k;
    double s, a, y;

    if ( sp->Link[j].type != CONDUIT ) return 0.0;
    if ( sp->Link[j].xsect.type == DUMMY ) return 0.0;
    q = fabs(q);
    k = sp->Link[j].subIndex;
    if ( q > sp->Conduit[k].qMax ) q = sp->Conduit[k].qMax;
    if ( q <= 0.0 ) return 0.0;
    s = q / sp->Conduit[k].beta;
    a = xsect_getAofS(sp, &sp->Link[j].xsect, s);
    y = xsect_getYofA(sp, &sp->Link[j].xsect, a);
    return y;
}

//=============================================================================

double link_getLength(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns length (ft)
//  Purpose: finds true length of a link.
//
{
    if ( sp->Link[j].type == CONDUIT ) return conduit_getLength(sp, j);
    return 0.0;
}

//=============================================================================

double link_getVelocity(SWMM_Project *sp, int j, double flow, double depth)
//
//  Input:   j     = link index
//           flow  = link flow rate (cfs)
//           depth = link flow depth (ft)
//  Output:  returns flow velocity (fps)
//  Purpose: finds flow velocity given flow and depth.
//
{
    double area;
    double veloc = 0.0;
    int    k;

    if ( depth <= 0.01 ) return 0.0;
    if ( sp->Link[j].type == CONDUIT )
    {
        k = sp->Link[j].subIndex;
        flow /= sp->Conduit[k].barrels;
        area = xsect_getAofY(sp, &sp->Link[j].xsect, depth);
        if (area > FUDGE ) veloc = flow / area;
    }
    return veloc;
}

//=============================================================================

double link_getFroude(SWMM_Project *sp, int j, double v, double y)
//
//  Input:   j = link index
//           v = flow velocity (fps)
//           y = flow depth (ft)
//  Output:  returns Froude Number
//  Purpose: computes Froude Number for given velocity and flow depth
//
{
    TXsect*  xsect = &sp->Link[j].xsect;

    // --- return 0 if link is not a conduit
    if ( sp->Link[j].type != CONDUIT ) return 0.0;

    // --- return 0 if link empty or closed conduit is full
    if ( y <= FUDGE ) return 0.0;
    if ( !xsect_isOpen(xsect->type) &&
         xsect->yFull - y <= FUDGE ) return 0.0;

    // --- compute hydraulic depth
    y = xsect_getAofY(sp, xsect, y) / xsect_getWofY(sp, xsect, y);

    // --- compute Froude No.
    return fabs(v) / sqrt(GRAVITY * y);
}

//=============================================================================

double link_getPower(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns power consumed by link in kwatts
//  Purpose: computes power consumed by head loss (or head gain) of
//           water flowing through a link
//
{
    int    n1 = sp->Link[j].node1;
    int    n2 = sp->Link[j].node2;
    double dh = (sp->Node[n1].invertElev + sp->Node[n1].newDepth) -
                (sp->Node[n2].invertElev + sp->Node[n2].newDepth);
    double q =  fabs(sp->Link[j].newFlow);
    return fabs(dh) * q / 8.814 * KWperHP;
}

//=============================================================================

double link_getLossRate(SWMM_Project *sp, int j, double q, double tStep)                         //(5.1.008)
//
//  Input:   j = link index
//           q = flow rate (ft3/sec)                                           //(5.1.008)
//           tstep = time step (sec)
//  Output:  returns uniform loss rate in link (ft3/sec)
//  Purpose: computes rate at which flow volume is lost in a link due to       //(5.1.008)
//           evaporation and seepage.
//
{
    if ( sp->Link[j].type == CONDUIT ) return conduit_getLossRate(sp, j, q, tStep);    //(5.1.008)
    else return 0.0;
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

char  link_getFullState(double a1, double a2, double aFull)
//
//  Input:   a1 = upstream link area (ft2)
//           a2 = downstream link area (ft2)
//           aFull = area of full conduit
//  Output:  returns fullness state of a link
//  Purpose: determines if a link is upstream, downstream or completely full.
//  
{
    if ( a1 >= aFull )
    {
        if ( a2 >= aFull ) return ALL_FULL;
        else return UP_FULL;
    }
    if ( a2 >= aFull ) return DN_FULL;
    return 0;
}

//=============================================================================
//                    C O N D U I T   M E T H O D S
//=============================================================================

int  conduit_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = link index
//           k = conduit index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads conduit parameters from a tokenzed line of input.
//
{
    int    n1, n2;
    double x[6];
    char*  id;

    // --- check for valid ID and end node IDs
    if ( ntoks < 7 ) return error_setInpError(sp, ERR_ITEMS, "");
    id = project_findID(sp, LINK, tok[0]);                // link ID
    if ( id == NULL ) return error_setInpError(sp, ERR_NAME, tok[0]);
    n1 = project_findObject(sp, NODE, tok[1]);            // upstrm. node
    if ( n1 < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
    n2 = project_findObject(sp, NODE, tok[2]);            // dwnstrm. node
    if ( n2 < 0 ) return error_setInpError(sp, ERR_NAME, tok[2]);

    // --- parse length & Mannings N
    if ( !getDouble(tok[3], &x[0]) )
        return error_setInpError(sp, ERR_NUMBER, tok[3]);
    if ( !getDouble(tok[4], &x[1]) )
        return error_setInpError(sp, ERR_NUMBER, tok[4]);

    // --- parse offsets
    if ( sp->LinkOffsets == ELEV_OFFSET && *tok[5] == '*' ) x[2] = MISSING;
    else if ( !getDouble(tok[5], &x[2]) )
        return error_setInpError(sp, ERR_NUMBER, tok[5]);
    if ( sp->LinkOffsets == ELEV_OFFSET && *tok[6] == '*' ) x[3] = MISSING;
    else if ( !getDouble(tok[6], &x[3]) )
        return error_setInpError(sp, ERR_NUMBER, tok[6]);

   // --- parse optional parameters
    x[4] = 0.0;                                       // init. flow
    if ( ntoks >= 8 )
    {
        if ( !getDouble(tok[7], &x[4]) )
        return error_setInpError(sp, ERR_NUMBER, tok[7]);
    }
    x[5] = 0.0;
    if ( ntoks >= 9 )
    {
        if ( !getDouble(tok[8], &x[5]) )
        return error_setInpError(sp, ERR_NUMBER, tok[8]);
    }

    // --- add parameters to data base
    sp->Link[j].ID = id;
    link_setParams(sp, j, CONDUIT, n1, n2, k, x);
    return 0;
}

//=============================================================================

void  conduit_validate(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = conduit index
//  Output:  none
//  Purpose: validates a conduit's properties.
//
{
    double aa;
    double lengthFactor, roughness, slope;

    // --- a storage node cannot have a dummy outflow link
    if ( sp->Link[j].xsect.type == DUMMY && sp->RouteModel == DW )                     //(5.1.007)
    {
        if ( sp->Node[sp->Link[j].node1].type == STORAGE )
        {
            report_writeErrorMsg(sp, ERR_DUMMY_LINK, sp->Node[sp->Link[j].node1].ID);
            return;
        }
    }

    // --- if custom xsection, then set its parameters
    if ( sp->Link[j].xsect.type == CUSTOM )
        xsect_setCustomXsectParams(sp, &sp->Link[j].xsect);

    // --- if irreg. xsection, assign transect roughness to conduit
    if ( sp->Link[j].xsect.type == IRREGULAR )
    {
        xsect_setIrregXsectParams(sp, &sp->Link[j].xsect);
        sp->Conduit[k].roughness = sp->Transect[sp->Link[j].xsect.transect].roughness;
    }

    // --- if force main xsection, adjust units on D-W roughness height
    if ( sp->Link[j].xsect.type == FORCE_MAIN )
    {
        if ( sp->ForceMainEqn == D_W ) sp->Link[j].xsect.rBot /= UCF(sp, RAINDEPTH);
        if ( sp->Link[j].xsect.rBot <= 0.0 )
            report_writeErrorMsg(sp, ERR_XSECT, sp->Link[j].ID);
    }

    // --- check for valid length & roughness
    if ( sp->Conduit[k].length <= 0.0 )
        report_writeErrorMsg(sp, ERR_LENGTH, sp->Link[j].ID);
    if ( sp->Conduit[k].roughness <= 0.0 )
        report_writeErrorMsg(sp, ERR_ROUGHNESS, sp->Link[j].ID);
    if ( sp->Conduit[k].barrels <= 0 )
        report_writeErrorMsg(sp, ERR_BARRELS, sp->Link[j].ID);

    // --- check for valid xsection
    if ( sp->Link[j].xsect.type != DUMMY )
    {
        if ( sp->Link[j].xsect.type < 0 )
            report_writeErrorMsg(sp, ERR_NO_XSECT, sp->Link[j].ID);
        else if ( sp->Link[j].xsect.aFull <= 0.0 )
            report_writeErrorMsg(sp, ERR_XSECT, sp->Link[j].ID);
    }
    if ( sp->ErrorCode ) return;

    // --- check for negative offsets
    if ( sp->Link[j].offset1 < 0.0 )
    {
        report_writeWarningMsg(sp, WARN03, sp->Link[j].ID);
        sp->Link[j].offset1 = 0.0;
    }
	if ( sp->Link[j].offset2 < 0.0 )
    {
        report_writeWarningMsg(sp, WARN03, sp->Link[j].ID);
        sp->Link[j].offset2 = 0.0;
    }

    // --- adjust conduit offsets for partly filled circular xsection
    if ( sp->Link[j].xsect.type == FILLED_CIRCULAR )
    {
        sp->Link[j].offset1 += sp->Link[j].xsect.yBot;
        sp->Link[j].offset2 += sp->Link[j].xsect.yBot;
    }

    // --- compute conduit slope
    slope = conduit_getSlope(sp, j);
    sp->Conduit[k].slope = slope;

    // --- reverse orientation of conduit if using dynamic wave routing
    //     and slope is negative
    if ( sp->RouteModel == DW &&
         slope < 0.0 &&
         sp->Link[j].xsect.type != DUMMY )
    {
        conduit_reverse(sp, j, k);
    }

    // --- get equivalent Manning roughness for Force Mains
    //     for use when pipe is partly full
    roughness = sp->Conduit[k].roughness;
    if ( sp->RouteModel == DW && sp->Link[j].xsect.type == FORCE_MAIN )
    {
        roughness = forcemain_getEquivN(sp, j, k);
    }

    // --- adjust roughness for meandering natural channels
    if ( sp->Link[j].xsect.type == IRREGULAR )
    {
        lengthFactor = sp->Transect[sp->Link[j].xsect.transect].lengthFactor;
        roughness *= sqrt(lengthFactor);
    }

    // --- lengthen conduit if lengthening option is in effect
    lengthFactor = 1.0;
    if ( sp->RouteModel == DW &&
         sp->LengtheningStep > 0.0 &&
         sp->Link[j].xsect.type != DUMMY )
    {
        lengthFactor = conduit_getLengthFactor(sp, j, k, roughness);
    }

    if ( lengthFactor != 1.0 )
    {
        sp->Conduit[k].modLength = lengthFactor * conduit_getLength(sp, j);
        slope /= lengthFactor;
        roughness = roughness / sqrt(lengthFactor);
    }

    // --- compute roughness factor used when computing friction
    //     slope term in Dynamic Wave flow routing

    // --- special case for non-Manning Force Mains
    //     (roughness factor for full flow is saved in xsect.sBot)
    if ( sp->RouteModel == DW && sp->Link[j].xsect.type == FORCE_MAIN )
    {
        sp->Link[j].xsect.sBot =
            forcemain_getRoughFactor(sp, j, lengthFactor);
    }
    sp->Conduit[k].roughFactor = GRAVITY * SQR(roughness/PHI);

    // --- compute full flow through cross section
    if ( sp->Link[j].xsect.type == DUMMY ) sp->Conduit[k].beta = 0.0;
    else sp->Conduit[k].beta = PHI * sqrt(fabs(slope)) / roughness;
    sp->Link[j].qFull = sp->Link[j].xsect.sFull * sp->Conduit[k].beta;
    sp->Conduit[k].qMax = sp->Link[j].xsect.sMax * sp->Conduit[k].beta;

    // --- see if flow is supercritical most of time
    //     by comparing normal & critical velocities.
    //     (factor of 0.3 is for circular pipe 95% full)
    // NOTE: this factor was used in the past for a modified version of
    //       Kinematic Wave routing but is now deprecated.
    aa = sp->Conduit[k].beta / sqrt(32.2) *
         pow(sp->Link[j].xsect.yFull, 0.1666667) * 0.3;
    if ( aa >= 1.0 ) sp->Conduit[k].superCritical = TRUE;
    else             sp->Conduit[k].superCritical = FALSE;

    // --- set value of hasLosses flag
    if ( sp->Link[j].cLossInlet  == 0.0 &&
         sp->Link[j].cLossOutlet == 0.0 &&
         sp->Link[j].cLossAvg    == 0.0
       ) sp->Conduit[k].hasLosses = FALSE;
    else sp->Conduit[k].hasLosses = TRUE;
}

//=============================================================================

void conduit_reverse(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = conduit index
//  Output:  none
//  Purpose: reverses direction of a conduit
//
{
    int    i;
    double z;
    double cLoss;

    // --- reverse end nodes
    i = sp->Link[j].node1;
    sp->Link[j].node1 = sp->Link[j].node2;
    sp->Link[j].node2 = i;

    // --- reverse node offsets
    z = sp->Link[j].offset1;
    sp->Link[j].offset1 = sp->Link[j].offset2;
    sp->Link[j].offset2 = z;

    // --- reverse loss coeffs.
    cLoss = sp->Link[j].cLossInlet;
    sp->Link[j].cLossInlet = sp->Link[j].cLossOutlet;
    sp->Link[j].cLossOutlet = cLoss;

    // --- reverse direction & slope
    sp->Conduit[k].slope = -sp->Conduit[k].slope;
    sp->Link[j].direction *= (signed char)-1;

    // --- reverse initial flow value
    sp->Link[j].q0 = -sp->Link[j].q0;
}

//=============================================================================

double conduit_getLength(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns conduit's length (ft)
//  Purpose: finds true length of a conduit.
//
//  Note: for irregular natural channels, user inputs length of main
//        channel (for FEMA purposes) but program should use length
//        associated with entire flood plain. Transect.lengthFactor
//        is the ratio of these two lengths.
//
{
    int k = sp->Link[j].subIndex;
    int t;
    if ( sp->Link[j].xsect.type != IRREGULAR ) return sp->Conduit[k].length;
    t = sp->Link[j].xsect.transect;
    if ( t < 0 || t >= sp->Nobjects[TRANSECT] ) return sp->Conduit[k].length;
    return sp->Conduit[k].length / sp->Transect[t].lengthFactor;
}

//=============================================================================

double conduit_getLengthFactor(SWMM_Project *sp, int j, int k, double roughness)
//
//  Input:   j = link index
//           k = conduit index
//           roughness = conduit Manning's n
//  Output:  returns factor by which a conduit should be lengthened
//  Purpose: computes amount of conduit lengthing to improve numerical stability.
//
//  The following form of the Courant criterion is used:
//      L = t * v * (1 + Fr) / Fr
//  where L = conduit length, t = time step, v = velocity, & Fr = Froude No.
//  After substituting Fr = v / sqrt(gy), where y = flow depth, we get:
//    L = t * ( sqrt(gy) + v )
//
{
    double ratio;
    double yFull;
    double vFull;
    double tStep;

    // --- evaluate flow depth and velocity at full normal flow condition
    yFull = sp->Link[j].xsect.yFull;
    if ( xsect_isOpen(sp->Link[j].xsect.type) )
    {
        yFull = sp->Link[j].xsect.aFull / xsect_getWofY(sp, &sp->Link[j].xsect, yFull);
    }
    vFull = PHI / roughness * sp->Link[j].xsect.sFull *
            sqrt(fabs(sp->Conduit[k].slope)) / sp->Link[j].xsect.aFull;

    // --- determine ratio of Courant length to actual length
    if ( sp->LengtheningStep == 0.0 ) tStep = sp->RouteStep;
    else                          tStep = MIN(sp->RouteStep, sp->LengtheningStep);
    ratio = (sqrt(GRAVITY*yFull) + vFull) * tStep / conduit_getLength(sp, j);

    // --- return max. of 1.0 and ratio
    if ( ratio > 1.0 ) return ratio;
    else return 1.0;
}

//=============================================================================

double conduit_getSlope(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns conduit slope
//  Purpose: computes conduit slope.
//
{
    double elev1, elev2, delta, slope;
    double length = conduit_getLength(sp, j);

    // --- check that elevation drop > minimum allowable drop
    elev1 = sp->Link[j].offset1 + sp->Node[sp->Link[j].node1].invertElev;
    elev2 = sp->Link[j].offset2 + sp->Node[sp->Link[j].node2].invertElev;
    delta = fabs(elev1 - elev2);
    if ( delta < MIN_DELTA_Z )
    {
        report_writeWarningMsg(sp, WARN04, sp->Link[j].ID);
        delta = MIN_DELTA_Z;
    }

    // --- elevation drop cannot exceed conduit length
    if ( delta >= length )
    {
        report_writeWarningMsg(sp, WARN08, sp->Link[j].ID);
        slope = delta / length;
    }

    // --- slope = elev. drop / horizontal distance
    else slope = delta / sqrt(SQR(length) - SQR(delta));

    // -- check that slope exceeds minimum allowable slope
    if ( sp->MinSlope > 0.0 && slope < sp->MinSlope )
    {
        report_writeWarningMsg(sp, WARN05, sp->Link[j].ID);
        slope = sp->MinSlope;
        // keep min. slope positive for SF or KW routing
        if (sp->RouteModel == SF || sp->RouteModel == KW) return slope;
    }

    // --- change sign for adverse slope
    if ( elev1 < elev2 ) slope = -slope;
    return slope;
}

//=============================================================================

void  conduit_initState(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = conduit index
//  Output:  none
//  Purpose: sets initial conduit depth to normal depth of initial flow
//
{
    sp->Link[j].newDepth = link_getYnorm(sp, j, sp->Link[j].q0 / sp->Conduit[k].barrels);
    sp->Link[j].oldDepth = sp->Link[j].newDepth;
}

//=============================================================================

double conduit_getInflow(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns flow in link (cfs)
//  Purpose: finds inflow to conduit from upstream node.
//
{
    double qIn = node_getOutflow(sp, sp->Link[j].node1, j);
    if ( sp->Link[j].qLimit > 0.0 ) qIn = MIN(qIn, sp->Link[j].qLimit);
    return qIn;
}

//=============================================================================

double conduit_getLossRate(SWMM_Project *sp, int j, double q, double tStep)                      //(5.1.008)
//
//  Input:   j = link index
//           tStep = time step (sec)
//  Output:  returns rate of evaporation & seepage losses (ft3/sec)            //(5.1.008)
//  Purpose: computes volumetric rate of water evaporation & seepage           //(5.1.008)
//           from a conduit (per barrel).                                      //(5.1.008)
//
{
    TXsect *xsect;
    double depth = 0.5 * (sp->Link[j].oldDepth + sp->Link[j].newDepth);
    double length;
    double topWidth;
    //double wettedPerimeter;    //DEPRECATED                                  //(5.1.012)
    double maxLossRate;
    double evapLossRate = 0.0,
           seepLossRate = 0.0,
           totalLossRate = 0.0;

    if ( depth > FUDGE )
    {
        xsect = &sp->Link[j].xsect;
        length = conduit_getLength(sp, j);

        // --- find evaporation rate for open conduits
        if ( xsect_isOpen(xsect->type) && sp->Evap.rate > 0.0 )
        {
            topWidth = xsect_getWofY(sp, xsect, depth);
            evapLossRate = topWidth * length * sp->Evap.rate;
        }

        // --- compute seepage loss rate
        if ( sp->Link[j].seepRate > 0.0 )
        {
            // limit depth to depth at max width
            if ( depth >= xsect->ywMax ) depth = xsect->ywMax;
			
//// The following section was deprecated for release 5.1.012. ////            //(5.1.012)			
            // get wetted perimeter
//          wettedPerimeter = 0.0;
//          if ( depth > 0.0 )
//          {
//              wettedPerimeter = xsect_getAofY(xsect, depth) /
//                                xsect_getRofY(xsect, depth);
//          }
/////////////////////////////////////////////////////////////////

            // compute seepage loss rate across length of conduit
            seepLossRate = sp->Link[j].seepRate * xsect_getWofY(sp, xsect, depth) *    //(5.1.012)
                           length;
            seepLossRate *= sp->Adjust.hydconFactor;                               //(5.1.008)
        }

        // --- compute total loss rate
        totalLossRate = evapLossRate + seepLossRate;

        // --- total loss rate cannot exceed flow volume or flow rate          //(5.1.008)
        if ( totalLossRate > 0.0 )
        {
            maxLossRate = 0.5 * (sp->Link[j].oldVolume + sp->Link[j].newVolume) / tStep;
            maxLossRate = MIN(maxLossRate, fabs(q));                           //(5.1.008)
            if ( totalLossRate > maxLossRate )
            {
                evapLossRate = evapLossRate * maxLossRate / totalLossRate;
                seepLossRate = seepLossRate * maxLossRate / totalLossRate;
                totalLossRate = maxLossRate;
            }
        }
    }
    sp->Conduit[sp->Link[j].subIndex].evapLossRate = evapLossRate;
    sp->Conduit[sp->Link[j].subIndex].seepLossRate = seepLossRate;
    return totalLossRate;
}


//=============================================================================
//                        P U M P   M E T H O D S
//=============================================================================

int  pump_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = link index
//           k = pump index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads pump parameters from a tokenized line of input.
//
{
    int    m;
    int    n1, n2;
    double x[4];
    char*  id;

    // --- check for valid ID and end node IDs
    if ( ntoks < 3 ) return error_setInpError(sp, ERR_ITEMS, "");
    id = project_findID(sp, LINK, tok[0]);
    if ( id == NULL ) return error_setInpError(sp, ERR_NAME, tok[0]);
    n1 = project_findObject(sp, NODE, tok[1]);
    if ( n1 < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
    n2 = project_findObject(sp, NODE, tok[2]);
    if ( n2 < 0 ) return error_setInpError(sp, ERR_NAME, tok[2]);

    // --- parse curve name
    x[0] = -1.;
    if ( ntoks >= 4 )
    {
        if ( !strcomp(tok[3],"*") )
        {
            m = project_findObject(sp, CURVE, tok[3]);
            if ( m < 0 ) return error_setInpError(sp, ERR_NAME, tok[3]);
            x[0] = m;
        }
    }

    // --- parse init. status if present
    x[1] = 1.0;
    if ( ntoks >= 5 )
    {
        m = findmatch(tok[4], OffOnWords);
        if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[4]);
        x[1] = m;
    }

    // --- parse startup/shutoff depths if present
    x[2] = 0.0;
    if ( ntoks >= 6 )
    {
        if ( !getDouble(tok[5], &x[2]) || x[2] < 0.0)
        return error_setInpError(sp, ERR_NUMBER, tok[5]);
    }
    x[3] = 0.0;
    if ( ntoks >= 7 )
    {
        if ( !getDouble(tok[6], &x[3]) || x[3] < 0.0 )
        return error_setInpError(sp, ERR_NUMBER, tok[6]);
    }

    // --- add parameters to pump object
    sp->Link[j].ID = id;
    link_setParams(sp, j, PUMP, n1, n2, k, x);
    return 0;
}

//=============================================================================

void  pump_validate(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = pump index
//  Output:  none
//  Purpose: validates a pump's properties
//
{
    int    m, n1;
    double x, y;

    sp->Link[j].xsect.yFull = 0.0;

    // --- check for valid curve type
    m = sp->Pump[k].pumpCurve;
    if ( m < 0 )
    {
        sp->Pump[k].type = IDEAL_PUMP;
    }
    else
    {
        if ( sp->Curve[m].curveType < PUMP1_CURVE ||
             sp->Curve[m].curveType > PUMP4_CURVE )
            report_writeErrorMsg(sp, ERR_NO_CURVE, sp->Link[j].ID);

        // --- store pump curve type with pump's parameters
        else
        {
            sp->Pump[k].type = sp->Curve[m].curveType - PUMP1_CURVE;
            if ( table_getFirstEntry(sp, &sp->Curve[m], &x, &y) )
            {
                sp->Link[j].qFull = y;
                sp->Pump[k].xMin = x;
                sp->Pump[k].xMax = x;
                while ( table_getNextEntry(sp, &sp->Curve[m], &x, &y) )
                {
                    sp->Link[j].qFull = MAX(y, sp->Link[j].qFull);
                    sp->Pump[k].xMax = x;
                }
            }
            sp->Link[j].qFull /= UCF(sp, FLOW);
       }
    }

    // --- check that shutoff depth < startup depth
    if ( sp->Pump[k].yOn > 0.0 && sp->Pump[k].yOn <= sp->Pump[k].yOff )
        report_writeErrorMsg(sp, ERR_PUMP_LIMITS, sp->Link[j].ID);

    // --- assign wet well volume to inlet node of Type 1 pump
    if ( sp->Pump[k].type == TYPE1_PUMP )
    {
        n1 = sp->Link[j].node1;
        if ( sp->Node[n1].type != STORAGE )
            sp->Node[n1].fullVolume = MAX(sp->Node[n1].fullVolume,
                                      sp->Pump[k].xMax / UCF(sp, VOLUME));
    }

}

//=============================================================================

void  pump_initState(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = pump index
//  Output:  none
//  Purpose: initializes pump conditions at start of a simulation
//
{
    sp->Link[j].setting = sp->Pump[k].initSetting;
    sp->Link[j].targetSetting = sp->Pump[k].initSetting;
}

//=============================================================================

double pump_getInflow(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns pump flow (cfs)
//  Purpose: finds flow produced by a pump.
//
{
    int     k, m;
    int     n1, n2;
    double  vol, depth, head;
    double  qIn, qIn1, dh = 0.001;

    k = sp->Link[j].subIndex;
    m = sp->Pump[k].pumpCurve;
    n1 = sp->Link[j].node1;
    n2 = sp->Link[j].node2;

    // --- no flow if setting is closed
    sp->Link[j].flowClass = NO;
    sp->Link[j].setting = sp->Link[j].targetSetting;
    if ( sp->Link[j].setting == 0.0 ) return 0.0;

    // --- pump flow = node inflow for IDEAL_PUMP
    if ( sp->Pump[k].type == IDEAL_PUMP )
        qIn = sp->Node[n1].inflow + sp->Node[n1].overflow;

    // --- pumping rate depends on pump curve type
    else switch(sp->Curve[m].curveType)
    {
      case PUMP1_CURVE:
        vol = sp->Node[n1].newVolume * UCF(sp, VOLUME);
        qIn = table_intervalLookup(&sp->Curve[m], vol) / UCF(sp, FLOW);

        // --- check if off of pump curve
        if ( vol < sp->Pump[k].xMin || vol > sp->Pump[k].xMax )
            sp->Link[j].flowClass = YES;
        break;

      case PUMP2_CURVE:
        depth = sp->Node[n1].newDepth * UCF(sp, LENGTH);
        qIn = table_intervalLookup(&sp->Curve[m], depth) / UCF(sp, FLOW);

        // --- check if off of pump curve
        if ( depth < sp->Pump[k].xMin || depth > sp->Pump[k].xMax )
            sp->Link[j].flowClass = YES;
        break;

      case PUMP3_CURVE:
        head = ( (sp->Node[n2].newDepth + sp->Node[n2].invertElev) -
                 (sp->Node[n1].newDepth + sp->Node[n1].invertElev) );

		head = MAX(head, 0.0);

        qIn = table_lookup(&sp->Curve[m], head*UCF(sp, LENGTH)) / UCF(sp, FLOW);

        // --- compute dQ/dh (slope of pump curve) and
        //     reverse sign since flow decreases with increasing head
    	sp->Link[j].dqdh = -table_getSlope(&sp->Curve[m], head*UCF(sp, LENGTH)) *
                       UCF(sp, LENGTH) / UCF(sp, FLOW);

        // --- check if off of pump curve
        head *= UCF(sp, LENGTH);
        if ( head < sp->Pump[k].xMin || head > sp->Pump[k].xMax )
            sp->Link[j].flowClass = YES;
        break;

      case PUMP4_CURVE:
        depth = sp->Node[n1].newDepth;
        qIn = table_lookup(&sp->Curve[m], depth*UCF(sp, LENGTH)) / UCF(sp, FLOW);

        // --- compute dQ/dh (slope of pump curve)
        qIn1 = table_lookup(&sp->Curve[m], (depth+dh)*UCF(sp, LENGTH)) / UCF(sp, FLOW);
        sp->Link[j].dqdh = (qIn1 - qIn) / dh;

        // --- check if off of pump curve
        depth *= UCF(sp, LENGTH);
        if ( depth < sp->Pump[k].xMin ) sp->Link[j].flowClass = DN_DRY;
        if ( depth > sp->Pump[k].xMax ) sp->Link[j].flowClass = UP_DRY;
        break;

      default: qIn = 0.0;
    }

    // --- do not allow reverse flow through pump
    if ( qIn < 0.0 )  qIn = 0.0;
    return qIn * sp->Link[j].setting;
}


//=============================================================================
//                    O R I F I C E   M E T H O D S
//=============================================================================

int  orifice_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = link index
//           k = orifice index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads orifice parameters from a tokenized line of input.
//
{
    int    m;
    int    n1, n2;
    double x[5];
    char*  id;

    // --- check for valid ID and end node IDs
    if ( ntoks < 6 ) return error_setInpError(sp, ERR_ITEMS, "");
    id = project_findID(sp, LINK, tok[0]);
    if ( id == NULL ) return error_setInpError(sp, ERR_NAME, tok[0]);
    n1 = project_findObject(sp, NODE, tok[1]);
    if ( n1 < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
    n2 = project_findObject(sp, NODE, tok[2]);
    if ( n2 < 0 ) return error_setInpError(sp, ERR_NAME, tok[2]);

    // --- parse orifice parameters
    m = findmatch(tok[3], OrificeTypeWords);
    if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[3]);
    x[0] = m;                                              // type
    if ( sp->LinkOffsets == ELEV_OFFSET && *tok[4] == '*' ) x[1] = MISSING;
    else if ( ! getDouble(tok[4], &x[1]) )                 // crest height
        return error_setInpError(sp, ERR_NUMBER, tok[4]);
    if ( ! getDouble(tok[5], &x[2]) || x[2] < 0.0 )        // cDisch
        return error_setInpError(sp, ERR_NUMBER, tok[5]);
    x[3] = 0.0;
    if ( ntoks >= 7 )
    {
        m = findmatch(tok[6], NoYesWords);
        if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[6]);
        x[3] = m;                                          // flap gate
    }
    x[4] = 0.0;
    if ( ntoks >= 8 )
    {
        if ( ! getDouble(tok[7], &x[4]) || x[4] < 0.0 )    // orate
            return error_setInpError(sp, ERR_NUMBER, tok[7]);
    }

    // --- add parameters to orifice object
    sp->Link[j].ID = id;
    link_setParams(sp, j, ORIFICE, n1, n2, k, x);
    return 0;
}

//=============================================================================

void  orifice_validate(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = orifice index
//  Output:  none
//  Purpose: validates an orifice's properties
//
{
    int    err = 0;

    // --- check for valid xsection
    if ( sp->Link[j].xsect.type != RECT_CLOSED
    &&   sp->Link[j].xsect.type != CIRCULAR ) err = ERR_REGULATOR_SHAPE;
    if ( err > 0 )
    {
        report_writeErrorMsg(sp, err, sp->Link[j].ID);
        return;
    }

    // --- check for negative offset
    if ( sp->Link[j].offset1 < 0.0 ) sp->Link[j].offset1 = 0.0;

    // --- compute partial flow adjustment
    orifice_setSetting(sp, j, 0.0);

    // --- compute an equivalent length
    sp->Orifice[k].length = 2.0 * sp->RouteStep * sqrt(GRAVITY * sp->Link[j].xsect.yFull);
    sp->Orifice[k].length = MAX(200.0, sp->Orifice[k].length);
    sp->Orifice[k].surfArea = 0.0;
}

//=============================================================================

void  orifice_setSetting(SWMM_Project *sp, int j, double tstep)
//
//  Input:   j = link index
//           tstep = time step over which setting is adjusted (sec)
//  Output:  none
//  Purpose: updates an orifice's setting as a result of a control action.
//
{
    int    k = sp->Link[j].subIndex;
    double delta, step;
    double h, f;

    // --- case where adjustment rate is instantaneous
    if ( sp->Orifice[k].orate == 0.0 || tstep == 0.0)
        sp->Link[j].setting = sp->Link[j].targetSetting;

    // --- case where orifice setting depends on time step
    else
    {
        delta = sp->Link[j].targetSetting - sp->Link[j].setting;
        step = tstep / sp->Orifice[k].orate;
        if ( step + 0.001 >= fabs(delta) )
            sp->Link[j].setting = sp->Link[j].targetSetting;
        else sp->Link[j].setting += SGN(delta) * step;
    }

    // --- find effective orifice discharge coeff.
    h = sp->Link[j].setting * sp->Link[j].xsect.yFull;
    f = xsect_getAofY(sp, &sp->Link[j].xsect, h) * sqrt(2.0 * GRAVITY);
    sp->Orifice[k].cOrif = sp->Orifice[k].cDisch * f;

    // --- find equiv. discharge coeff. for when weir flow occurs
    sp->Orifice[k].cWeir = orifice_getWeirCoeff(sp, j, k, h) * f;
}

//=============================================================================

double orifice_getWeirCoeff(SWMM_Project *sp, int j, int k, double h)
//
//  Input:   j = link index
//           k = orifice index
//           h = height of orifice opening (ft)
//  Output:  returns a discharge coefficient (ft^1/2)
//  Purpose: computes the discharge coefficient for an orifice
//           at the critical depth where weir flow begins.
//
{
    double w, aOverL;

    // --- this is for bottom orifices
    if ( sp->Orifice[k].type == BOTTOM_ORIFICE )
    {
        // --- find critical height above opening where orifice flow
        //     turns into weir flow. It equals (Co/Cw)*(Area/Length)
        //     where Co is the orifice coeff., Cw is the weir coeff/sqrt(2g),
        //     Area is the area of the opening, and Length = circumference
        //     of the opening. For a basic sharp crested weir, Cw = 0.414.
        if (sp->Link[j].xsect.type == CIRCULAR) aOverL = h / 4.0;
        else
        {
            w = sp->Link[j].xsect.wMax;
            aOverL = (h*w) / (2.0*(h+w));
        }
        h = sp->Orifice[k].cDisch / 0.414 * aOverL;
        sp->Orifice[k].hCrit = h;
    }

    // --- this is for side orifices
    else
    {
        // --- critical height is simply height of opening
        sp->Orifice[k].hCrit = h;

        // --- head on orifice is distance to center line
        h = h / 2.0;
    }

    // --- return a coefficient for the critical depth
    return sp->Orifice[k].cDisch * sqrt(h);
}

//=============================================================================

double orifice_getInflow(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns orifice flow rate (cfs)
//  Purpose: finds the flow through an orifice.
//
{
    int    k, n1, n2;
    double head, h1, h2, y1, dir;
    double f;
    double hcrest = 0.0;
    double hcrown = 0.0;
    double hmidpt;
    double q, ratio;

    // --- get indexes of end nodes and link's orifice
    n1 = sp->Link[j].node1;
    n2 = sp->Link[j].node2;
    k  = sp->Link[j].subIndex;

    // --- find heads at upstream & downstream nodes
    if ( sp->RouteModel == DW )
    {
        h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
        h2 = sp->Node[n2].newDepth + sp->Node[n2].invertElev;
    }
    else
    {
        h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
        h2 = sp->Node[n1].invertElev;
    }
    dir = (h1 >= h2) ? +1.0 : -1.0;

    // --- exchange h1 and h2 for reverse flow
    y1 = sp->Node[n1].newDepth;
    if ( dir < 0.0 )
    {
        head = h1;
        h1 = h2;
        h2 = head;
        y1 = sp->Node[n2].newDepth;
    }

    // --- orifice is a bottom orifice (oriented in horizontal plane)
    if ( sp->Orifice[k].type == BOTTOM_ORIFICE )
    {
        // --- compute crest elevation
        hcrest = sp->Node[n1].invertElev + sp->Link[j].offset1;

        // --- compute head on orifice
        if (h1 < hcrest) head = 0.0;
        else if (h2 > hcrest) head = h1 - h2;
        else head = h1 - hcrest;

        // --- find fraction of critical height for which weir flow occurs
        f = head / sp->Orifice[k].hCrit;
        f = MIN(f, 1.0);
    }

    // --- otherwise orifice is a side orifice (oriented in vertical plane)
    else
    {
        // --- compute elevations of orifice crest and crown
        hcrest = sp->Node[n1].invertElev + sp->Link[j].offset1;
        hcrown = hcrest + sp->Link[j].xsect.yFull * sp->Link[j].setting;
        hmidpt = (hcrest + hcrown) / 2.0;

        // --- compute degree of inlet submergence
        if ( h1 < hcrown && hcrown > hcrest )
            f = (h1 - hcrest) / (hcrown - hcrest);
        else f = 1.0;

        // --- compute head on orifice
        if ( f < 1.0 )          head = h1 - hcrest;
        else if ( h2 < hmidpt ) head = h1 - hmidpt;
        else                    head = h1 - h2;
    }

    // --- return if head is negligible or flap gate closed
    if ( head <= FUDGE || y1 <= FUDGE ||
         link_setFlapGate(sp, j, n1, n2, dir) )
    {
        sp->Link[j].newDepth = 0.0;
        sp->Link[j].flowClass = DRY;
        sp->Orifice[k].surfArea = FUDGE * sp->Orifice[k].length;
        sp->Link[j].dqdh = 0.0;
        return 0.0;
    }

    // --- determine flow class
    sp->Link[j].flowClass = SUBCRITICAL;
    if ( hcrest > h2 )
    {
        if ( dir == 1.0 ) sp->Link[j].flowClass = DN_CRITICAL;
        else              sp->Link[j].flowClass = UP_CRITICAL;
    }

    // --- compute flow depth and surface area
    y1 = sp->Link[j].xsect.yFull * sp->Link[j].setting;
    if ( sp->Orifice[k].type == SIDE_ORIFICE )
    {
        sp->Link[j].newDepth = y1 * f;
        sp->Orifice[k].surfArea =
            xsect_getWofY(sp, &sp->Link[j].xsect, sp->Link[j].newDepth) *
            sp->Orifice[k].length;
    }
    else
    {
        sp->Link[j].newDepth = y1;
        sp->Orifice[k].surfArea = xsect_getAofY(sp, &sp->Link[j].xsect, y1);
    }

    // --- find flow through the orifice
    q = dir * orifice_getFlow(sp, j, k, head, f, sp->Link[j].hasFlapGate);

    // --- apply Villemonte eqn. to correct for submergence
    if ( f < 1.0 && h2 > hcrest )
    {
        ratio = (h2 - hcrest) / (h1 - hcrest);
        q *= pow( (1.0 - pow(ratio, 1.5)), 0.385);
    }
    return q;
}

//=============================================================================

double orifice_getFlow(SWMM_Project *sp, int j, int k,  double head, double f,
        int hasFlapGate)
//
//  Input:   j = link index
//           k = orifice index
//           head = head across orifice
//           f = fraction of critical depth filled
//           hasFlapGate = flap gate indicator
//  Output:  returns flow through an orifice
//  Purpose: computes flow through an orifice as a function of head.
//
{
    double area, q;
    double veloc, hLoss;

    // --- case where orifice is closed
    if ( head == 0.0 || f <= 0.0  )
    {
        sp->Link[j].dqdh = 0.0;
        return 0.0;
    }

    // --- case where inlet depth is below critical depth;
    //     orifice behaves as a weir
    else if ( f < 1.0 )
    {
        q = sp->Orifice[k].cWeir * pow(f, 1.5);
        sp->Link[j].dqdh = 1.5 * q / (f * sp->Orifice[k].hCrit);
    }

    // --- case where normal orifice flow applies
    else
    {
        q = sp->Orifice[k].cOrif * sqrt(head);
        sp->Link[j].dqdh = q / (2.0 * head);
    }

    // --- apply ARMCO adjustment for headloss from flap gate
    if ( hasFlapGate )
    {
        // --- compute velocity for current orifice flow
        area = xsect_getAofY(sp, &sp->Link[j].xsect,
                sp->Link[j].setting * sp->Link[j].xsect.yFull);
        veloc = q / area;

        // --- compute head loss from gate
        hLoss = (4.0 / GRAVITY) * veloc * veloc *
                 exp(-1.15 * veloc / sqrt(head) );

        // --- update head (for orifice flow)
        //     or critical depth fraction (for weir flow)
        if ( f < 1.0 )
        {
            f = f - hLoss/sp->Orifice[k].hCrit;
            if ( f < 0.0 ) f = 0.0;
        }
        else
        {
            head = head - hLoss;
            if ( head < 0.0 ) head = 0.0;
        }

        // --- make recursive call to this function, with hasFlapGate
        //     set to false, to find flow values at adjusted head value
        q = orifice_getFlow(sp, j, k, head, f, FALSE);
    }
    return q;
}

//=============================================================================
//                           W E I R   M E T H O D S
//=============================================================================

int   weir_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = link index
//           k = weir index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads weir parameters from a tokenized line of input.
//
{
    int    m;
    int    n1, n2;
    double x[9];                                                               //(5.1.010)
    char*  id;

    // --- check for valid ID and end node IDs
    if ( ntoks < 6 ) return error_setInpError(sp, ERR_ITEMS, "");
    id = project_findID(sp, LINK, tok[0]);
    if ( id == NULL ) return error_setInpError(sp, ERR_NAME, tok[0]);
    n1 = project_findObject(sp, NODE, tok[1]);
    if ( n1 < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
    n2 = project_findObject(sp, NODE, tok[2]);
    if ( n2 < 0 ) return error_setInpError(sp, ERR_NAME, tok[2]);

    // --- parse weir parameters
    m = findmatch(tok[3], WeirTypeWords);
    if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[3]);
    x[0] = m;                                              // type
    if ( sp->LinkOffsets == ELEV_OFFSET && *tok[4] == '*' ) x[1] = MISSING;
    else if ( ! getDouble(tok[4], &x[1]) )                 // height
        return error_setInpError(sp, ERR_NUMBER, tok[4]);
    if ( ! getDouble(tok[5], &x[2]) || x[2] < 0.0 )        // cDisch1
        return error_setInpError(sp, ERR_NUMBER, tok[5]);
    x[3] = 0.0;
    x[4] = 0.0;
    x[5] = 0.0;
    x[6] = 1.0;                                                                //(5.1.007)
    x[7] = 0.0;                                                                //(5.1.010)
    x[8] = 0.0;                                                                //(5.1.010)
    if ( ntoks >= 7 && *tok[6] != '*' )                                        //(5.1.007)
    {
        m = findmatch(tok[6], NoYesWords);
        if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[6]);
        x[3] = m;                                          // flap gate
    }
    if ( ntoks >= 8 && *tok[7] != '*' )                                        //(5.1.007)
    {
        if ( ! getDouble(tok[7], &x[4]) || x[4] < 0.0 )     // endCon
            return error_setInpError(sp, ERR_NUMBER, tok[7]);
    }
    if ( ntoks >= 9 && *tok[8] != '*' )                                        //(5.1.007)
    {
        if ( ! getDouble(tok[8], &x[5]) || x[5] < 0.0 )     // cDisch2
            return error_setInpError(sp, ERR_NUMBER, tok[8]);
    }

////  Following segment added for release 5.1.007.  ////                       //(5.1.007)
    if ( ntoks >= 10 && *tok[9] != '*' )
    {
        m = findmatch(tok[9], NoYesWords);
        if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[9]);
        x[6] = m;                                           // canSurcharge
    }
////

////  Following segment added for release 5.1.010.  ////                       //(5.1.010)
    if ( (m = (int)x[0]) == ROADWAY_WEIR )
    {
        if ( ntoks >= 11 )                                  // road width
        {
            if ( ! getDouble(tok[10], &x[7]) || x[7] < 0.0 ) 
                return error_setInpError(sp, ERR_NUMBER, tok[10]);
        }
        if ( ntoks >= 12 )                                  // road surface
        {
            if ( strcomp(tok[11], "PAVED") ) x[8] = 1.0;
            else if ( strcomp(tok[11], "GRAVEL") ) x[8] = 2.0;
        }
    }
////

    // --- add parameters to weir object
    sp->Link[j].ID = id;
    link_setParams(sp, j, WEIR, n1, n2, k, x);
    return 0;
}

//=============================================================================

void  weir_validate(SWMM_Project *sp, int j, int k)
//
//  Input:   j = link index
//           k = weir index
//  Output:  none
//  Purpose: validates a weir's properties
//
{
    int    err = 0;
    double q, q1, q2, head;                                                    //(5.1.008)
 
    // --- check for valid cross section
    switch ( sp->Weir[k].type)
    {
      case TRANSVERSE_WEIR:
      case SIDEFLOW_WEIR:
      case ROADWAY_WEIR:                                                       //(5.1.010)
        if ( sp->Link[j].xsect.type != RECT_OPEN ) err = ERR_REGULATOR_SHAPE;
        sp->Weir[k].slope = 0.0;
        break;

      case VNOTCH_WEIR:
        if ( sp->Link[j].xsect.type != TRIANGULAR ) err = ERR_REGULATOR_SHAPE;
        else
        {
            sp->Weir[k].slope = sp->Link[j].xsect.sBot;
        }
        break;

      case TRAPEZOIDAL_WEIR:
        if ( sp->Link[j].xsect.type != TRAPEZOIDAL ) err = ERR_REGULATOR_SHAPE;
        else
        {
            sp->Weir[k].slope = sp->Link[j].xsect.sBot;
        }
        break;
    }
    if ( err > 0 )
    {
        report_writeErrorMsg(sp, err, sp->Link[j].ID);
        return;
    }

    // --- check for negative offset
    if ( sp->Link[j].offset1 < 0.0 ) sp->Link[j].offset1 = 0.0;

    // --- compute an equivalent length
    sp->Weir[k].length = 2.0 * sp->RouteStep * sqrt(GRAVITY * sp->Link[j].xsect.yFull);
    sp->Weir[k].length = MAX(200.0, sp->Weir[k].length);
    sp->Weir[k].surfArea = 0.0;

////  Following code segment added to release 5.1.008.  ////                   //(5.1.008)

    // --- find flow through weir when water level equals weir height
    head = sp->Link[j].xsect.yFull;
    weir_getFlow(sp, j, k, head, 1.0, FALSE, &q1, &q2);
    q = q1 + q2;
////

    // --- compute equivalent orifice coeff. (for CFS flow units)
    head = head / 2.0;  // head seen by equivalent orifice
    sp->Weir[k].cSurcharge = q / sqrt(head); 
}

//=============================================================================

////  New function added to release 5.1.007.  ////                             //(5.1.007)

void weir_setSetting(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  none
//  Purpose: updates a weir's setting as a result of a control action.
//
{
    int    k = sp->Link[j].subIndex;
    double h, q, q1, q2;

    // --- adjust weir setting
    sp->Link[j].setting = sp->Link[j].targetSetting;
    if ( !sp->Weir[k].canSurcharge ) return;
    if ( sp->Weir[k].type == ROADWAY_WEIR ) return;                                //(5.1.010)

    // --- find orifice coeff. for surcharged flow
    if ( sp->Link[j].setting == 0.0 ) sp->Weir[k].cSurcharge = 0.0;
    else
    {
        // --- find flow through weir when water level equals weir height
        h = sp->Link[j].setting * sp->Link[j].xsect.yFull;
        weir_getFlow(sp, j, k, h, 1.0, FALSE, &q1, &q2);
        q = q1 + q2;

        // --- compute equivalent orifice coeff. (for CFS flow units)
        h = h / 2.0;  // head seen by equivalent orifice
        sp->Weir[k].cSurcharge = q / sqrt(h);
    }
}

//=============================================================================

double weir_getInflow(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  returns weir flow rate (cfs)
//  Purpose: finds the flow over a weir.
//
{
    int    n1;          // index of upstream node
    int    n2;          // index of downstream node
    int    k;           // index of weir
    double q1;          // flow through central part of weir (cfs)
    double q2;          // flow through end sections of weir (cfs)
    double head;        // head on weir (ft)
    double h1;          // upstrm nodal head (ft)
    double h2;          // downstrm nodal head (ft)
    double hcrest;      // head at weir crest (ft)
    double hcrown;      // head at weir crown (ft)
    double y;           // water depth in weir (ft)
    double dir;         // direction multiplier
    double ratio;
    double weirPower[] = {1.5,       // transverse weir
                          5./3.,     // side flow weir
                          2.5,       // v-notch weir
                          1.5};      // trapezoidal weir

    n1 = sp->Link[j].node1;
    n2 = sp->Link[j].node2;
    k  = sp->Link[j].subIndex;
    if ( sp->RouteModel == DW )
    {
        h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
        h2 = sp->Node[n2].newDepth + sp->Node[n2].invertElev;
    }
    else
    {
        h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
        h2 = sp->Node[n1].invertElev;
    }
    dir = (h1 > h2) ? +1.0 : -1.0;

    // --- exchange h1 and h2 for reverse flow
    if ( dir < 0.0 )
    {
        head = h1;
        h1 = h2;
        h2 = head;
    }

    // --- find head of weir's crest and crown
    hcrest = sp->Node[n1].invertElev + sp->Link[j].offset1;
    hcrown = hcrest + sp->Link[j].xsect.yFull;

////  Added to release 5.1.010.  ////                                          //(5.1.010)
    // --- treat a roadway weir as a special case
    if ( sp->Weir[k].type == ROADWAY_WEIR )
        return roadway_getInflow(sp, j, dir, hcrest, h1, h2);
////

    // --- adjust crest ht. for partially open weir
    hcrest += (1.0 - sp->Link[j].setting) * sp->Link[j].xsect.yFull;

    // --- compute head relative to weir crest
    head = h1 - hcrest;

    // --- return if head is negligible or flap gate closed
    sp->Link[j].dqdh = 0.0;
    if ( head <= FUDGE || hcrest >= hcrown ||
         link_setFlapGate(sp, j, n1, n2, dir) )
    {
        sp->Link[j].newDepth = 0.0;
        sp->Link[j].flowClass = DRY;
        return 0.0;
    }

    // --- determine flow class
    sp->Link[j].flowClass = SUBCRITICAL;
    if ( hcrest > h2 )
    {
        if ( dir == 1.0 ) sp->Link[j].flowClass = DN_CRITICAL;
        else              sp->Link[j].flowClass = UP_CRITICAL;
    }

    // --- compute new equivalent surface area
    y = sp->Link[j].xsect.yFull - (hcrown - MIN(h1, hcrown));
    sp->Weir[k].surfArea = xsect_getWofY(sp, &sp->Link[j].xsect, y) * sp->Weir[k].length;

////  New section added to release 5.1.007.  ////                              //(5.1.007)

    // --- head is above crown
    if ( h1 >= hcrown )
    {
        // --- use equivalent orifice if weir can surcharge
        if ( sp->Weir[k].canSurcharge )
        {
            y = (hcrest + hcrown) / 2.0;
            if ( h2 < y ) head = h1 - y;
            else          head = h1 - h2;
            y = hcrown - hcrest;
            q1 = weir_getOrificeFlow(sp, j, head, y, sp->Weir[k].cSurcharge);
            sp->Link[j].newDepth = y;
            return dir * q1;
        }

        // --- otherwise limit head to height of weir opening
        else head = hcrown - hcrest;
    }
////

    // --- use weir eqn. to find flows through central (q1)
    //     and end sections (q2) of weir
    weir_getFlow(sp, j, k, head, dir, sp->Link[j].hasFlapGate, &q1, &q2);

    // --- apply Villemonte eqn. to correct for submergence
    if ( h2 > hcrest )
    {
        ratio = (h2 - hcrest) / (h1 - hcrest);
        q1 *= pow( (1.0 - pow(ratio, weirPower[sp->Weir[k].type])), 0.385);
        if ( q2 > 0.0 )
            q2 *= pow( (1.0 - pow(ratio, weirPower[VNOTCH_WEIR])), 0.385);
    }

    // --- return total flow through weir
    sp->Link[j].newDepth = MIN((h1 - hcrest), sp->Link[j].xsect.yFull);
    return dir * (q1 + q2);
}

//=============================================================================

void weir_getFlow(SWMM_Project *sp, int j, int k,  double head, double dir, int hasFlapGate,
                  double* q1, double* q2)
//
//  Input:   j    = link index
//           k    = weir index
//           head = head across weir (ft)
//           dir  = flow direction indicator
//           hasFlapGate = flap gate indicator
//  Output:  q1 = flow through central portion of weir (cfs)
//           q2 = flow through end sections of weir (cfs)
//  Purpose: computes flow over weir given head.
//
{
    double length;
    double h;
    double y;
    double hLoss;
    double area;
    double veloc;
    int    wType;

    // --- q1 = flow through central portion of weir,
    //     q2 = flow through end sections of trapezoidal weir
    *q1 = 0.0;
    *q2 = 0.0;
    sp->Link[j].dqdh = 0.0;
    if ( head <= 0.0 ) return;

    // --- convert weir length & head to original units
    length = sp->Link[j].xsect.wMax * UCF(sp, LENGTH);
    h = head * UCF(sp, LENGTH);

////  Following code segment re-located.  ////                                 //(5.1.012)
    // --- reduce length when end contractions present
    //length -= 0.1 * sp->Weir[k].endCon * h;
    //length = MAX(length, 0.0);
/////////////////////////////////////////////

    // --- use appropriate formula for weir flow
    wType = sp->Weir[k].type;
    if ( wType == VNOTCH_WEIR &&
         sp->Link[j].setting < 1.0 ) wType = TRAPEZOIDAL_WEIR;
    switch (wType)
    {
      case TRANSVERSE_WEIR:

        // --- reduce length when end contractions present                     //(5.1.012)
        length -= 0.1 * sp->Weir[k].endCon * h;                                    //(5.1.012)
        length = MAX(length, 0.0);                                             //(5.1.012)
        *q1 = sp->Weir[k].cDisch1 * length * pow(h, 1.5);
        break;

      case SIDEFLOW_WEIR:

        // --- reduce length when end contractions present                     //(5.1.012)
        length -= 0.1 * sp->Weir[k].endCon * h;                                    //(5.1.012)
        length = MAX(length, 0.0);                                             //(5.1.012)

        // --- weir behaves as a transverse weir under reverse flow
        if ( dir < 0.0 )
            *q1 = sp->Weir[k].cDisch1 * length * pow(h, 1.5);
        else

////   Corrected formula  ////                                                 //(5.1.012)
// (see Metcalf & Eddy, Inc., Wastewater Engineering, McGraw-Hill, 1972 p. 164).
            *q1 = sp->Weir[k].cDisch1 * pow(length, 0.83) * pow(h, 1.67);

        break;

      case VNOTCH_WEIR:
        *q1 = sp->Weir[k].cDisch1 * sp->Weir[k].slope * pow(h, 2.5);
        break;

      case TRAPEZOIDAL_WEIR:
        y = (1.0 - sp->Link[j].setting) * sp->Link[j].xsect.yFull;
        length = xsect_getWofY(sp, &sp->Link[j].xsect, y) * UCF(sp, LENGTH);

////  End contractions don't apply to trapezoidal weirs ////                   //(5.1.012)
        //length -= 0.1 * sp->Weir[k].endCon * h;                                  //(5.1.012)
        //length = MAX(length, 0.0);                                           //(5.1.012)

        *q1 = sp->Weir[k].cDisch1 * length * pow(h, 1.5);
        *q2 = sp->Weir[k].cDisch2 * sp->Weir[k].slope * pow(h, 2.5);
    }

    // --- convert CMS flows to CFS
    if ( sp->UnitSystem == SI )
    {
        *q1 /= M3perFT3;
        *q2 /= M3perFT3;
    }

    // --- apply ARMCO adjustment for headloss from flap gate
    if ( hasFlapGate )
    {
        // --- compute flow area & velocity for current weir flow
        area = weir_getOpenArea(sp, j, head);
        if ( area > TINY )
        {
            veloc = (*q1 + *q2) / area;

            // --- compute headloss and subtract from original head
            hLoss = (4.0 / GRAVITY) * veloc * veloc *
                     exp(-1.15 * veloc / sqrt(head) );
            head = head - hLoss;
            if ( head < 0.0 ) head = 0.0;

            // --- make recursive call to this function, with hasFlapGate
            //     set to false, to find flow values at adjusted head value
            weir_getFlow(sp, j, k, head, dir, FALSE, q1, q2);
        }
    }
    sp->Link[j].dqdh = weir_getdqdh(sp, k, dir, head, *q1, *q2);
}

//=============================================================================

double weir_getOrificeFlow(SWMM_Project *sp, int j, double head, double y, double cOrif)
//
//  Input:   j = link index
//           head = head across weir (ft)
//           y = height of upstream water level above weir crest (ft)
//           cOrif = orifice flow coefficient
//  Output:  returns flow through weir
//  Purpose: finds flow through a surcharged weir using the orifice equation.
//
{
    double a, q, v, hloss;

    // --- evaluate the orifice flow equation
    q = cOrif * sqrt(head);

    // --- apply Armco adjustment if weir has a flap gate
    if ( sp->Link[j].hasFlapGate )
    {
        a = weir_getOpenArea(sp, j, y);
        if ( a > 0.0 )
        {
            v = q / a;
            hloss = (4.0 / GRAVITY) * v * v * exp(-1.15 * v / sqrt(y) );
            head -= hloss;
            head = MAX(head, 0.0);
            q = cOrif * sqrt(head);
        }
    }
    if ( head > 0.0 ) sp->Link[j].dqdh = q / (2.0 * head);
    else sp->Link[j].dqdh = 0.0;
    return q;
}

//=============================================================================

double weir_getOpenArea(SWMM_Project *sp, int j, double y)
//
//  Input:   j = link index
//           y = depth of water above weir crest (ft)
//  Output:  returns area between weir crest and y (ft2)
//  Purpose: finds flow area through a weir.
//
{
    double z, zy;

    // --- find offset of weir crest due to control setting
    z = (1.0 - sp->Link[j].setting) * sp->Link[j].xsect.yFull;

    // --- ht. of crest + ht of water above crest
    zy = z + y;
    zy = MIN(zy, sp->Link[j].xsect.yFull);

    // --- return difference between area of offset + water depth
    //     and area of just the offset
    return xsect_getAofY(sp, &sp->Link[j].xsect, zy) -
           xsect_getAofY(sp, &sp->Link[j].xsect, z);
}

//=============================================================================

double  weir_getdqdh(SWMM_Project *sp, int k, double dir, double h, double q1,
        double q2)
{
    double q1h;
    double q2h;

    if ( fabs(h) < FUDGE ) return 0.0;
    q1h = fabs(q1/h);
    q2h = fabs(q2/h);

    switch (sp->Weir[k].type)
    {
      case TRANSVERSE_WEIR: return 1.5 * q1h;

      case SIDEFLOW_WEIR:
        // --- weir behaves as a transverse weir under reverse flow
        if ( dir < 0.0 ) return 1.5 * q1h;
        else return 1.67 * q1h;                                                //(5.1.012)

      case VNOTCH_WEIR:
        if ( q2h == 0.0 ) return 2.5 * q1h;  // Fully open
        else return 1.5 * q1h + 2.5 * q2h;   // Partly open

      case TRAPEZOIDAL_WEIR: return 1.5 * q1h + 2.5 * q2h;
    }
    return 0.0;
}


//=============================================================================
//               O U T L E T    D E V I C E    M E T H O D S
//=============================================================================

int outlet_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = link index
//           k = outlet index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads outlet parameters from a tokenized  line of input.
//
{
    int    i, m, n;
    int    n1, n2;
    double x[6];
    char*  id;
    char*  s;

    // --- check for valid ID and end node IDs
    if ( ntoks < 6 ) return error_setInpError(sp, ERR_ITEMS, "");
    id = project_findID(sp, LINK, tok[0]);
    if ( id == NULL ) return error_setInpError(sp, ERR_NAME, tok[0]);
    n1 = project_findObject(sp, NODE, tok[1]);
    if ( n1 < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
    n2 = project_findObject(sp, NODE, tok[2]);
    if ( n2 < 0 ) return error_setInpError(sp, ERR_NAME, tok[2]);

    // --- get height above invert
    if ( sp->LinkOffsets == ELEV_OFFSET && *tok[3] == '*' ) x[0] = MISSING;
    else
    {
        if ( ! getDouble(tok[3], &x[0]) )
            return error_setInpError(sp, ERR_NUMBER, tok[3]);
	if ( sp->LinkOffsets == DEPTH_OFFSET && x[0] < 0.0 ) x[0] = 0.0;
    }

    // --- see if outlet flow relation is tabular or functional
    m = findmatch(tok[4], RelationWords);
    if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[4]);
    x[1] = 0.0;
    x[2] = 0.0;
    x[3] = -1.0;
    x[4] = 0.0;

    // --- see if rating curve is head or depth based
    x[5] = NODE_DEPTH;                                //default is depth-based
    s = strtok(tok[4], "/");                          //parse token for
    s = strtok(NULL, "/");                            //  qualifier term
    if ( strcomp(s, w_HEAD) ) x[5] = NODE_HEAD;       //check if its "HEAD"

    // --- get params. for functional outlet device
    if ( m == FUNCTIONAL )
    {
        if ( ntoks < 7 ) return error_setInpError(sp, ERR_ITEMS, "");
        if ( ! getDouble(tok[5], &x[1]) )
            return error_setInpError(sp, ERR_NUMBER, tok[5]);
        if ( ! getDouble(tok[6], &x[2]) )
            return error_setInpError(sp, ERR_NUMBER, tok[6]);
        n = 7;
    }

    // --- get name of outlet rating curve
    else
    {
        i = project_findObject(sp, CURVE, tok[5]);
        if ( i < 0 ) return error_setInpError(sp, ERR_NAME, tok[5]);
        x[3] = i;
        n = 6;
    }

    // --- check if flap gate specified
    if ( ntoks > n)
    {
        i = findmatch(tok[n], NoYesWords);
        if ( i < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[n]);
        x[4] = i;
    }

    // --- add parameters to outlet object
    sp->Link[j].ID = id;
    link_setParams(sp, j, OUTLET, n1, n2, k, x);
    return 0;
}

//=============================================================================

double outlet_getInflow(SWMM_Project *sp, int j)
//
//  Input:   j = link index
//  Output:  outlet flow rate (cfs)
//  Purpose: finds the flow through an outlet.
//
{
    int    k, n1, n2;
    double head, hcrest, h1, h2, y1, dir;

    // --- get indexes of end nodes
    n1 = sp->Link[j].node1;
    n2 = sp->Link[j].node2;
    k  = sp->Link[j].subIndex;

    // --- find heads at upstream & downstream nodes
    if ( sp->RouteModel == DW )
    {
        h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
        h2 = sp->Node[n2].newDepth + sp->Node[n2].invertElev;
    }
    else
    {
        h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
        h2 = sp->Node[n1].invertElev;
    }
    dir = (h1 >= h2) ? +1.0 : -1.0;

    // --- exchange h1 and h2 for reverse flow
    y1 = sp->Node[n1].newDepth;
    if ( dir < 0.0 )
    {
        y1 = h1;
        h1 = h2;
        h2 = y1;
        y1 = sp->Node[n2].newDepth;
    }

    // --- for a NODE_DEPTH rating curve the effective head across the
    //     outlet is the depth above the crest elev. while for a NODE_HEAD
    //     curve it is the difference between upstream & downstream heads
    hcrest = sp->Node[n1].invertElev + sp->Link[j].offset1;
    if ( sp->Outlet[k].curveType == NODE_HEAD && sp->RouteModel == DW )
        head = h1 - MAX(h2, hcrest);
    else head = h1 - hcrest;

    // --- no flow if either no effective head difference,
    //     no upstream water available, or closed flap gate
    if ( head <= FUDGE || y1 <= FUDGE ||
         link_setFlapGate(sp, j, n1, n2, dir) )
    {
        sp->Link[j].newDepth = 0.0;
        sp->Link[j].flowClass = DRY;
        return 0.0;
    }

    // --- otherwise use rating curve to compute flow
    sp->Link[j].newDepth = head;
    sp->Link[j].flowClass = SUBCRITICAL;
    return dir * sp->Link[j].setting * outlet_getFlow(sp, k, head);
}

//=============================================================================

double outlet_getFlow(SWMM_Project *sp, int k, double head)
//
//  Input:   k    = outlet index
//           head = head across outlet (ft)
//  Output:  returns outlet flow rate (cfs)
//  Purpose: computes flow rate through an outlet given head.
//
{
    int    m;
    double h;

    // --- convert head to original units
    h = head * UCF(sp, LENGTH);

    // --- look-up flow in rating curve table if provided
    m = sp->Outlet[k].qCurve;
    if ( m >= 0 ) return table_lookup(&sp->Curve[m], h) / UCF(sp, FLOW);

    // --- otherwise use function to find flow
    else return sp->Outlet[k].qCoeff * pow(h, sp->Outlet[k].qExpon) / UCF(sp, FLOW);
}
