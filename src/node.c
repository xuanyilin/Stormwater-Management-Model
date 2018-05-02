//-----------------------------------------------------------------------------
//   node.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             09/15/14   (Build 5.1.007)
//             04/02/15   (Build 5.1.008)
//             08/05/15   (Build 5.1.010)
//   Author:   L. Rossman
//
//   Conveyance system node functions.
//
//   Build 5.1.007:
//   - Ponded area property for storage nodes deprecated.
//   - Support for Green-Ampt seepage from bottom and sides of storage node added.
//   - Storage node evap. & seepage losses now computed at start of each routing
//     time step.
//
//   Build 5.1.008:
//   - Support added for sending outfall discharge to a subctchment.
//
//   Build 5.1.010:
//   - Storage losses now based on node's new volume instead of old volume.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "findroot.h"

//-----------------------------------------------------------------------------                  
//  Local Declarations
//-----------------------------------------------------------------------------
typedef struct
{
    int     k;                  // storage unit index
    double  v;                  // storage unit volume (ft3)
} TStorageVol;

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  node_readParams        (called from readNode in input.c)
//  node_validate          (called from project_validate) 
//  node_initState         (called from project_init)
//  node_setOldHydState    (called from routing_execute)
//  node_setOldQualState   (called from routing_execute)
//  node_initInflow        (called from routing_execute)
//  node_setOutletDepth    (called from routing_execute)
//  node_getLosses         (called from routing_execute)                       //(5.1.008)
//  node_getSystemOutflow  (called from removeOutflows in routing.c)
//  node_getResults        (called from output_saveNodeResults)
//  node_getPondedArea     (called from initNodeStates in dynwave.c)
//  node_getOutflow        (called from link_getInflow & conduit_getInflow)
//  node_getMaxOutflow     (called from flowrout.c and dynwave.c)
//  node_getSurfArea
//  node_getDepth
//  node_getVolume
//  node_getPondedDepth    removed                                             //(5.1.008)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void   node_setParams(SWMM_Project *sp, int j, int type, int k, double x[]);
static int    junc_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks);

static int    outfall_readParams(SWMM_Project *sp, int j, int k, char* tok[],
        int ntoks);
static void   outfall_setOutletDepth(SWMM_Project *sp, int j, double yNorm,
        double yCrit, double z);

static int    storage_readParams(SWMM_Project *sp, int j, int k, char* tok[],
        int ntoks);
static double storage_getDepth(SWMM_Project *sp, int j, double v);
static double storage_getVolume(SWMM_Project *sp, int j, double d);
static double storage_getSurfArea(SWMM_Project *sp, int j, double d);
static void   storage_getVolDiff(SWMM_Project *sp, double y, double* f,
        double* df, void* p);
static double storage_getOutflow(SWMM_Project *sp, int j, int i);
static double storage_getLosses(SWMM_Project *sp, int j, double tStep);

static int    divider_readParams(SWMM_Project *sp, int j, int k, char* tok[],
        int ntoks);
static void   divider_validate(SWMM_Project *sp, int j);
static double divider_getOutflow(SWMM_Project *sp, int j, int link);


//=============================================================================

int node_readParams(SWMM_Project *sp, int j, int type, int k, char* tok[], int ntoks)
//
//  Input:   j = node index
//           type = node type code
//           k = index of node type
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads node properties from a tokenized line of input.
//
{
    switch ( type )
    {
      case JUNCTION: return junc_readParams(sp, j, k, tok, ntoks);
      case OUTFALL:  return outfall_readParams(sp, j, k, tok, ntoks);
      case STORAGE:  return storage_readParams(sp, j, k, tok, ntoks);
      case DIVIDER:  return divider_readParams(sp, j, k, tok, ntoks);
      default:       return 0;
    }
}

//=============================================================================

void  node_setParams(SWMM_Project *sp, int j, int type, int k, double x[])
//
//  Input:   j = node index
//           type = node type code
//           k = index of node type
//           x[] = array of property values
//  Output:  none
//  Purpose: assigns property values to a node.
//
{
    sp->Node[j].type       = type;
    sp->Node[j].subIndex   = k;
    sp->Node[j].invertElev = x[0] / UCF(sp, LENGTH);
    sp->Node[j].crownElev  = sp->Node[j].invertElev;
    sp->Node[j].initDepth  = 0.0;
    sp->Node[j].newVolume  = 0.0;
    sp->Node[j].fullVolume = 0.0;
    sp->Node[j].fullDepth  = 0.0;
    sp->Node[j].surDepth   = 0.0;
    sp->Node[j].pondedArea = 0.0;
    sp->Node[j].degree     = 0;
    switch (type)
    {
      case JUNCTION:
        sp->Node[j].fullDepth = x[1] / UCF(sp, LENGTH);
        sp->Node[j].initDepth = x[2] / UCF(sp, LENGTH);
        sp->Node[j].surDepth  = x[3] / UCF(sp, LENGTH);
        sp->Node[j].pondedArea = x[4] / (UCF(sp, LENGTH)*UCF(sp, LENGTH));
        break;

      case OUTFALL:
        sp->Outfall[k].type        = (int)x[1];
        sp->Outfall[k].outfallStage= x[2] / UCF(sp, LENGTH);
        sp->Outfall[k].tideCurve   = (int)x[3];
        sp->Outfall[k].stageSeries = (int)x[4];
        sp->Outfall[k].hasFlapGate = (char)x[5];

////  Following code segment added to release 5.1.008.  ////                   //(5.1.008)

        sp->Outfall[k].routeTo     = (int)x[6];
        sp->Outfall[k].wRouted     = NULL;
        if ( sp->Outfall[k].routeTo >= 0 )
        {
            sp->Outfall[k].wRouted =
                (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        }
////
        break;

      case STORAGE:
        sp->Node[j].fullDepth  = x[1] / UCF(sp, LENGTH);
        sp->Node[j].initDepth  = x[2] / UCF(sp, LENGTH);
        sp->Storage[k].aCoeff  = x[3];
        sp->Storage[k].aExpon  = x[4];
        sp->Storage[k].aConst  = x[5];
        sp->Storage[k].aCurve  = (int)x[6];
        // x[7] (ponded depth) is deprecated.                                  //(5.1.007)
        sp->Storage[k].fEvap   = x[8];
        break;

      case DIVIDER:
        sp->Divider[k].link      = (int)x[1];
        sp->Divider[k].type      = (int)x[2];
        sp->Divider[k].flowCurve = (int)x[3];
        sp->Divider[k].qMin      = x[4] / UCF(sp, FLOW);
        sp->Divider[k].dhMax     = x[5];
        sp->Divider[k].cWeir     = x[6];
        sp->Node[j].fullDepth    = x[7] / UCF(sp, LENGTH);
        sp->Node[j].initDepth    = x[8] / UCF(sp, LENGTH);
        sp->Node[j].surDepth     = x[9] / UCF(sp, LENGTH);
        sp->Node[j].pondedArea   = x[10] / (UCF(sp, LENGTH)*UCF(sp, LENGTH));
        break;
    }
}

//=============================================================================

void  node_validate(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: validates a node's properties.
//
{
    TDwfInflow* inflow;

    // --- see if full depth was increased to accommodate conduit crown
    if ( sp->Node[j].fullDepth > sp->Node[j].oldDepth && sp->Node[j].oldDepth > 0.0 )
    {
        report_writeWarningMsg(sp, WARN02, sp->Node[j].ID);
    }

    // --- check that initial depth does not exceed max. depth
    if ( sp->Node[j].initDepth > sp->Node[j].fullDepth + sp->Node[j].surDepth )
        report_writeErrorMsg(sp, ERR_NODE_DEPTH, sp->Node[j].ID);

    if ( sp->Node[j].type == DIVIDER ) divider_validate(sp, j);

    // --- initialize dry weather inflows
    inflow = sp->Node[j].dwfInflow;
    while (inflow)
    {
        inflow_initDwfInflow(sp, inflow);
        inflow = inflow->next;
    }
}

//=============================================================================

void node_initState(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: initializes a node's state variables at start of simulation.
//
{
    int p, k;                                                                  //(5.1.007)

    // --- initialize depth
    sp->Node[j].oldDepth = sp->Node[j].initDepth;
    sp->Node[j].newDepth = sp->Node[j].oldDepth;
    sp->Node[j].crownElev = sp->Node[j].invertElev;

    sp->Node[j].fullVolume = node_getVolume(sp, j, sp->Node[j].fullDepth);
    sp->Node[j].oldVolume = node_getVolume(sp, j, sp->Node[j].oldDepth);
    sp->Node[j].newVolume = sp->Node[j].oldVolume;

    // --- initialize water quality state
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        sp->Node[j].oldQual[p]  = 0.0;
        sp->Node[j].newQual[p]  = 0.0;
    }

    // --- initialize any inflow
    sp->Node[j].oldLatFlow = 0.0;
    sp->Node[j].newLatFlow = 0.0;
    sp->Node[j].losses = 0.0;                                                      //(5.1.007)


////  Following code section added to release 5.1.007.  ////                   //(5.1.007)

    // --- initialize storage nodes
    if ( sp->Node[j].type == STORAGE )
    {
        // --- set hydraulic residence time to 0
        k = sp->Node[j].subIndex;
        sp->Storage[k].hrt = 0.0;

        // --- initialize exfiltration properties
        if ( sp->Storage[k].exfil ) exfil_initState(sp, k);
    }
////

////  Following code section added to release 5.1.008.  ////                   //(5.1.008)

    // --- initialize flow stream routed from outfall onto a subcatchment
    if ( sp->Node[j].type == OUTFALL )
    {
        k = sp->Node[j].subIndex;
        if ( sp->Outfall[k].routeTo >= 0 )
        {
            sp->Outfall[k].vRouted = 0.0;
            for (p = 0; p < sp->Nobjects[POLLUT]; p++) sp->Outfall[k].wRouted[p] = 0.0;
        }
    }
////
}

//=============================================================================

void node_setOldHydState(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: replaces a node's old hydraulic state values with new ones.
//
{
    sp->Node[j].oldDepth    = sp->Node[j].newDepth;
    sp->Node[j].oldLatFlow  = sp->Node[j].newLatFlow;
    sp->Node[j].oldVolume   = sp->Node[j].newVolume;
}

//=============================================================================

void node_setOldQualState(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: replaces a node's old water quality state values with new ones.
//
{
    int p;
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        sp->Node[j].oldQual[p] = sp->Node[j].newQual[p];
        sp->Node[j].newQual[p] = 0.0;
    }
}

//=============================================================================

void node_initInflow(SWMM_Project *sp, int j, double tStep)
//
//  Input:   j = node index
//           tStep = time step (sec)
//  Output:  none
//  Purpose: initializes a node's inflow at start of next time step.
//
{
    // --- initialize inflow & outflow
    sp->Node[j].oldFlowInflow = sp->Node[j].inflow;
    sp->Node[j].oldNetInflow  = sp->Node[j].inflow - sp->Node[j].outflow;
    sp->Node[j].inflow = sp->Node[j].newLatFlow;
    sp->Node[j].outflow = sp->Node[j].losses;                                          //(5.1.007)

    // --- set overflow to any excess stored volume
    if ( sp->Node[j].newVolume > sp->Node[j].fullVolume )
        sp->Node[j].overflow = (sp->Node[j].newVolume - sp->Node[j].fullVolume) / tStep;
    else sp->Node[j].overflow = 0.0;
}

//=============================================================================

double node_getDepth(SWMM_Project *sp, int j, double v)
//
//  Input:   j = node index
//           v = volume (ft3)
//  Output:  returns depth of water at a node (ft)
//  Purpose: computes a node's water depth from its volume.
//
{
    switch ( sp->Node[j].type )
    {
      case STORAGE: return storage_getDepth(sp, j, v);
      default:      return 0.0;
    }
}

//=============================================================================

double node_getVolume(SWMM_Project *sp, int j, double d)
//
//  Input:   j = node index
//           d = water depth (ft)
//  Output:  returns volume of water at a node (ft3)
//  Purpose: computes volume stored at a node from its water depth.
//
{
    switch ( sp->Node[j].type )
    {
      case STORAGE: return storage_getVolume(sp, j, d);

      default:
        if ( sp->Node[j].fullDepth > 0.0 )
            return sp->Node[j].fullVolume * (d / sp->Node[j].fullDepth);
        else return 0.0;
    }
}

//=============================================================================

double  node_getSurfArea(SWMM_Project *sp, int j, double d)
//
//  Input:   j = node index
//           d = water depth (ft)
//  Output:  returns surface area of water at a node (ft2)
//  Purpose: computes surface area of water stored at a node from water depth.
//
{
    switch (sp->Node[j].type)
    {
      case STORAGE: return storage_getSurfArea(sp, j, d);
      default:      return 0.0;        
    }
}

//=============================================================================

double node_getOutflow(SWMM_Project *sp, int j, int k)
//
//  Input:   j = node index
//           k = link index
//  Output:  returns flow rate (cfs)
//  Purpose: computes outflow from node available for inflow into a link.
//
{
    switch ( sp->Node[j].type )
    {
      case DIVIDER: return divider_getOutflow(sp, j, k);
      case STORAGE: return storage_getOutflow(sp, j, k);
      default:      return sp->Node[j].inflow + sp->Node[j].overflow;
    }
}

//=============================================================================

double node_getMaxOutflow(SWMM_Project *sp, int j, double q, double tStep)
//
//  Input:   j = node index
//           q = original outflow rate (cfs)
//           tStep = time step (sec)
//  Output:  returns modified flow rate (cfs)
//  Purpose: limits outflow rate from a node with storage volume.
//
{
    double qMax;
    if ( sp->Node[j].fullVolume > 0.0 )
    {
        qMax = sp->Node[j].inflow + sp->Node[j].oldVolume / tStep;
        if ( q > qMax ) q = qMax;
    }
    return MAX(0.0, q);
}

//=============================================================================

double node_getSystemOutflow(SWMM_Project *sp, int j, int *isFlooded)
//
//  Input:   j = node index
//           isFlooded = TRUE if node becomes flooded
//  Output:  returns flow rate lost from system (cfs)
//  Purpose: computes flow rate at outfalls and flooded nodes.
//
{
    double outflow = 0.0;;

    // --- assume there is no flooding
    *isFlooded = FALSE;

    // --- if node is an outfall
    if ( sp->Node[j].type == OUTFALL )
    {
        // --- node receives inflow from outfall conduit
        if ( sp->Node[j].outflow == 0.0 ) outflow = sp->Node[j].inflow;

        // --- node sends flow into outfall conduit
        //     (therefore it has a negative outflow)
        else

////  Following code segment modified for release 5.1.007.  ////               //(5.1.007)
        {
            if ( sp->Node[j].inflow == 0.0 )
            {
                outflow = -sp->Node[j].outflow;
                sp->Node[j].inflow = fabs(outflow);
            }
        }

        // --- set overflow and volume to 0
        sp->Node[j].overflow = 0.0;
        sp->Node[j].newVolume = 0.0;
    }

    // --- node is a terminal node under Steady or Kin. Wave routing
    else if ( sp->RouteModel != DW &&
              sp->Node[j].degree == 0 &&
              sp->Node[j].type != STORAGE
            )
    {
        if ( sp->Node[j].outflow == 0.0 ) outflow = sp->Node[j].inflow;
        sp->Node[j].overflow = 0.0;
        sp->Node[j].newVolume = 0.0;
    }

    // --- otherwise node is an interior node and any
    //     overflow is considered as system outflow and flooding
    else 
    {
        if ( sp->Node[j].newVolume <= sp->Node[j].fullVolume)
            outflow = sp->Node[j].overflow;
        if ( outflow > 0.0 ) *isFlooded = TRUE;
    }
    return outflow;
}

//=============================================================================

void node_getResults(SWMM_Project *sp, int j, double f, float x[])
//
//  Input:   j = node index
//           f = weighting factor
//           x[] = array of nodal reporting variables
//  Output:  none
//  Purpose: computes weighted average of old and new results at a node.
//
{
    int    p;
    double z;
    double f1 = 1.0 - f;

    z = (f1 * sp->Node[j].oldDepth + f * sp->Node[j].newDepth) * UCF(sp, LENGTH);
    x[NODE_DEPTH] = (float)z;
    z = sp->Node[j].invertElev * UCF(sp, LENGTH);
    x[NODE_HEAD] = x[NODE_DEPTH] + (float)z;
    z = (f1*sp->Node[j].oldVolume + f*sp->Node[j].newVolume) * UCF(sp, VOLUME);
    x[NODE_VOLUME]  = (float)z;
    z = (f1*sp->Node[j].oldLatFlow + f*sp->Node[j].newLatFlow) * UCF(sp, FLOW); 
    x[NODE_LATFLOW] = (float)z;
    z = (f1*sp->Node[j].oldFlowInflow + f*sp->Node[j].inflow) * UCF(sp, FLOW);
    x[NODE_INFLOW] = (float)z;
    z = sp->Node[j].overflow * UCF(sp, FLOW);
    x[NODE_OVERFLOW] = (float)z;

    if ( !sp->IgnoreQuality ) for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        z = f1*sp->Node[j].oldQual[p] + f*sp->Node[j].newQual[p];
        x[NODE_QUAL+p] = (float)z;
    }
}

//=============================================================================

void   node_setOutletDepth(SWMM_Project *sp, int j, double yNorm, double yCrit, double z)
//
//  Input:   j = node index
//           yNorm = normal flow depth (ft)
//           yCrit = critical flow depth (ft)
//           z = offset of connecting outfall link from node invert (ft)
//  Output:  none
//  Purpose: sets water depth at a node that serves as an outlet point.
//
{
    switch (sp->Node[j].type)
    {
      // --- do nothing if outlet is a storage unit
      case STORAGE:
        return;

      // --- if outlet is a designated outfall then use outfall's specs
      case OUTFALL:
        outfall_setOutletDepth(sp, j, yNorm, yCrit, z);
        break;

      // --- for all other nodes, use min. of critical & normal depths
      default:
        if ( z > 0.0 ) sp->Node[j].newDepth = 0.0;
        else sp->Node[j].newDepth = MIN(yNorm, yCrit);
    }
}

//=============================================================================

double node_getPondedArea(SWMM_Project *sp, int j, double d)
//
//  Input:   j = node index
//           d = water depth (ft)
//  Output:  returns surface area of water at a node (ft2)
//  Purpose: computes surface area of water at a node based on depth.
//
{
    double a;

    // --- use regular getSurfArea function if node not flooded
    if ( d <= sp->Node[j].fullDepth || sp->Node[j].pondedArea == 0.0 )
    {
        return node_getSurfArea(sp, j, d);
    }

    // --- compute ponded depth
    d = d - sp->Node[j].fullDepth;

    // --- use ponded area for flooded node
    a = sp->Node[j].pondedArea;
    if ( a <= 0.0 ) a = node_getSurfArea(sp, j, sp->Node[j].fullDepth);
    return a;
}

//=============================================================================

double node_getLosses(SWMM_Project *sp, int j, double tStep)
//
//  Input:   j = node index
//           tStep = time step (sec)
//  Output:  returns water loss rate at node (ft3)
//  Purpose: computes the rates of evaporation and infiltration over a given
//           time step for a node.
//
{
    if ( sp->Node[j].type == STORAGE ) return storage_getLosses(sp, j, tStep);
    else return 0.0;
}

//=============================================================================
//                   J U N C T I O N   M E T H O D S
//=============================================================================

int junc_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = node index
//           k = junction index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error message
//  Purpose: reads a junction's properties from a tokenized line of input.
//
//  Format of input line is:
//     nodeID  elev  maxDepth  initDepth  surDepth  aPond 
{
    int    i;
    double x[6];
    char*  id;

    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");
    id = project_findID(NODE, tok[0]);
    if ( id == NULL ) return error_setInpError(ERR_NAME, tok[0]);

    // --- parse invert elev., max. depth, init. depth, surcharged depth,
    //     & ponded area values
    for ( i = 1; i <= 5; i++ )
    {
        x[i-1] = 0.0;
        if ( i < ntoks )
        {
            if ( ! getDouble(tok[i], &x[i-1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
    }

    // --- check for non-negative values (except for invert elev.)
    for ( i = 1; i <= 4; i++ )
    {
        if ( x[i] < 0.0 ) return error_setInpError(ERR_NUMBER, tok[i+1]);
    }

    // --- add parameters to junction object
    sp->Node[j].ID = id;
    node_setParams(sp, j, JUNCTION, k, x);
    return 0;
}


//=============================================================================
//                   S T O R A G E   M E T H O D S
//=============================================================================

int storage_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = node index
//           k = storage unit index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error message
//  Purpose: reads a storage unit's properties from a tokenized line of input.
//
//  Format of input line is:
//     nodeID  elev  maxDepth  initDepth  FUNCTIONAL a1 a2 a0 aPond fEvap (infil)
//     nodeID  elev  maxDepth  initDepth  TABULAR    curveID  aPond fEvap (infil)
//
{
    int    i, m, n;
    double x[9];
    char*  id;

    // --- get ID name
    if ( ntoks < 6 ) return error_setInpError(ERR_ITEMS, "");
    id = project_findID(NODE, tok[0]);
    if ( id == NULL ) return error_setInpError(ERR_NAME, tok[0]);

    // --- get invert elev, max. depth, & init. depth
    for ( i = 1; i <= 3; i++ )
    {
        if ( ! getDouble(tok[i], &x[i-1]) )
            return error_setInpError(ERR_NUMBER, tok[i]);
    }

    // --- get surf. area relation type
    m = findmatch(tok[4], RelationWords);
    if ( m < 0 ) return error_setInpError(ERR_KEYWORD, tok[4]);
    x[3] = 0.0;                        // a1 
    x[4] = 0.0;                        // a2
    x[5] = 0.0;                        // a0
    x[6] = -1.0;                       // curveID
    x[7] = 0.0;                        // aPond
    x[8] = 0.0;                        // fEvap

    // --- get surf. area function coeffs.
    if ( m == FUNCTIONAL )
    {
        for (i=5; i<=7; i++)
        {
            if ( i < ntoks )
            {
                if ( ! getDouble(tok[i], &x[i-2]) )
                    return error_setInpError(ERR_NUMBER, tok[i]);
            }
        }
        n = 8;
    }

    // --- get surf. area curve name
    else
    {
        m = project_findObject(CURVE, tok[5]);
        if ( m < 0 ) return error_setInpError(ERR_NAME, tok[5]);
        x[6] = m;
        n = 6;
    }

    // --- ignore next token if present (deprecated ponded area property)      //(5.1.007) 
    if ( ntoks > n)
    {
        if ( ! getDouble(tok[n], &x[7]) )
            return error_setInpError(ERR_NUMBER, tok[n]);
        n++;
    }

    // --- get evaporation fraction if present
    if ( ntoks > n )
    {
        if ( ! getDouble(tok[n], &x[8]) )
            return error_setInpError(ERR_NUMBER, tok[n]);
        n++;
    }

    // --- add parameters to storage unit object
    sp->Node[j].ID = id;
    node_setParams(sp, j, STORAGE, k, x);

    // --- read exfiltration parameters if present
    if ( ntoks > n ) return exfil_readStorageParams(sp, k, tok, ntoks, n);         //(5.1.007)
    return 0;
}

//=============================================================================

double storage_getDepth(SWMM_Project *sp, int j, double v)
//
//  Input:   j = node index
//           v = volume (ft3)
//  Output:  returns depth of water at a storage node (ft)
//  Purpose: computes a storage node's water depth from its volume.
//
{
    int    k = sp->Node[j].subIndex;
    int    i = sp->Storage[k].aCurve;
    double d, e;
	TStorageVol storageVol;

    // --- return max depth if a max. volume has been computed
    //     and volume is > max. volume
    if ( sp->Node[j].fullVolume > 0.0
    &&   v >= sp->Node[j].fullVolume ) return sp->Node[j].fullDepth;
    if ( v == 0.0 ) return 0.0;

    // --- use tabular area v. depth curve
    if ( i >= 0 )
        return table_getInverseArea(&sp->Curve[i], v*UCF(sp, VOLUME)) / UCF(sp, LENGTH);

    // --- use functional area v. depth relation
    else
    {
        v *= UCF(sp, VOLUME);
        if ( sp->Storage[k].aExpon == 0.0 )
        {
            d = v / (sp->Storage[k].aConst + sp->Storage[k].aCoeff);
        }
        else if ( sp->Storage[k].aConst == 0.0 )
        {
            e = 1.0 / (sp->Storage[k].aExpon + 1.0);
            d = pow(v / (sp->Storage[k].aCoeff * e), e);
        }
        else
        {
            storageVol.k = k;
            storageVol.v = v;
            d = v / (sp->Storage[k].aConst + sp->Storage[k].aCoeff);
            findroot_Newton(sp, 0.0, sp->Node[j].fullDepth*UCF(sp, LENGTH), &d,
                            0.001, storage_getVolDiff, &storageVol);            
        }
        d /= UCF(sp, LENGTH);
        if ( d > sp->Node[j].fullDepth ) d = sp->Node[j].fullDepth;
        return d;
    }
}

//=============================================================================

void  storage_getVolDiff(SWMM_Project *sp, double y, double* f, double* df,
        void* p)
//
//  Input:   y = depth of water (ft)
//  Output:  f = volume of water (ft3)
//           df = dVolume/dDepth (ft2)
//  Purpose: computes volume and its derivative with respect to depth
//           at storage node Kstar using the node's area versus depth function.
//
{
    int    k;
    double e, v;
    TStorageVol* storageVol;
		
    // ... cast void pointer p to a TStorageVol object
    storageVol = (TStorageVol *)p;
    k = storageVol->k;

    // ... find storage volume at depth y
    e = sp->Storage[k].aExpon + 1.0;
    v = sp->Storage[k].aConst * y + sp->Storage[k].aCoeff / e * pow(y, e);

    // ... compute difference between this volume and target volume
    //     as well as its derivative w.r.t. y
    *f = v - storageVol->v;
    *df = sp->Storage[k].aConst + sp->Storage[k].aCoeff * pow(y, e-1.0);
}

//=============================================================================

double storage_getVolume(SWMM_Project *sp, int j, double d)
//
//  Input:   j = node index
//           d = depth (ft)
//  Output:  returns volume of stored water (ft3)
//  Purpose: computes a storage node's water volume from its depth.
//
{
    int    k = sp->Node[j].subIndex;
    int    i = sp->Storage[k].aCurve;
    double v;

    // --- return full volume if depth >= max. depth
    if ( d == 0.0 ) return 0.0;
    if ( d >= sp->Node[j].fullDepth
    &&   sp->Node[j].fullVolume > 0.0 ) return sp->Node[j].fullVolume;

    // --- use table integration if area v. depth table exists
    if ( i >= 0 )
        return table_getArea(&sp->Curve[i], d*UCF(sp, LENGTH)) / UCF(sp, VOLUME);

    // --- otherwise use functional area v. depth relation
    else
    {
        d *= UCF(sp, LENGTH);
        v = sp->Storage[k].aConst * d;
        v += sp->Storage[k].aCoeff / (sp->Storage[k].aExpon+1.0) *
             pow(d, sp->Storage[k].aExpon+1.0);
        return v / UCF(sp, VOLUME);
    }
}

//=============================================================================

double storage_getSurfArea(SWMM_Project *sp, int j, double d)
//
//  Input:   j = node index
//           d = depth (ft)
//  Output:  returns surface area (ft2)
//  Purpose: computes a storage node's surface area from its water depth.
//
{
    double area;
    int k = sp->Node[j].subIndex;
    int i = sp->Storage[k].aCurve;
    if ( i >= 0 )
        area = table_lookupEx(&sp->Curve[i], d*UCF(sp, LENGTH));
    else
    {
        if ( sp->Storage[k].aCoeff <= 0.0 ) area = sp->Storage[k].aConst;
        else if ( sp->Storage[k].aExpon == 0.0 )
            area = sp->Storage[k].aConst + sp->Storage[k].aCoeff;
        else area = sp->Storage[k].aConst + sp->Storage[k].aCoeff *
                    pow(d*UCF(sp, LENGTH), sp->Storage[k].aExpon);
    }
    return area / UCF(sp, LENGTH) / UCF(sp, LENGTH);
}

//=============================================================================

double storage_getOutflow(SWMM_Project *sp, int j, int i)
//
//  Input:   j = node index
//           i = link index
//  Output:  returns flow from storage node into conduit link (cfs)
//  Purpose: finds outflow from a storage node into its connecting conduit link
//           ( non-conduit links have their own getInflow functions).
//
{
    int    k;
    double a, y;

    // --- link must be a conduit
    if ( sp->Link[i].type != CONDUIT ) return 0.0;

    // --- find depth of water in conduit
    y = sp->Node[j].newDepth - sp->Link[i].offset1;

    // --- return 0 if conduit empty or full flow if full
    if ( y <= 0.0 ) return 0.0;
    if ( y >= sp->Link[i].xsect.yFull ) return sp->Link[i].qFull;

    // --- if partially full, return normal flow
    k = sp->Link[i].subIndex;
    a = xsect_getAofY(sp, &sp->Link[i].xsect, y);
    return sp->Conduit[k].beta * xsect_getSofA(sp, &sp->Link[i].xsect, a);
}

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

double storage_getLosses(SWMM_Project *sp, int j, double tStep)
//
//  Input:   j = node index
//           tStep = time step (sec)
//  Output:  returns evaporation + infiltration rate (cfs)
//  Purpose: computes combined rate of water evaporated & infiltrated from
//           a storage node.
//
{
    int    k;
    double depth;
    double area;
    double evapRate = 0.0;
    double exfilRate = 0.0;
    double totalLoss = 0.0;
    double lossRatio;
    TExfil* exfil;

    // --- if node has some stored volume
    if ( sp->Node[j].newVolume > FUDGE )                                           //(5.1.010)
    {
        // --- get node's evap. rate (ft/s) &  exfiltration object
        k = sp->Node[j].subIndex;
        evapRate = sp->Evap.rate * sp->Storage[k].fEvap;
        exfil = sp->Storage[k].exfil;

        // --- if either of these apply
        if ( evapRate > 0.0 || exfil != NULL) 
        {
            // --- obtain storage depth & surface area 
            depth = sp->Node[j].newDepth;                                          //(5.1.010)
            area = storage_getSurfArea(sp, j, depth);

            // --- compute evap rate over this area (cfs)
            evapRate = area * evapRate;

            // --- find exfiltration rate (cfs) through bottom and side banks
            if ( exfil != NULL )
            {
                exfilRate = exfil_getLoss(sp, exfil, tStep, depth, area);
            }

            // --- total loss over time step cannot exceed stored volume
            totalLoss = (evapRate + exfilRate) * tStep;
            if ( totalLoss > sp->Node[j].newVolume )                               //(5.1.010)
            {
                lossRatio = sp->Node[j].newVolume / totalLoss;                     //(5.1.010)
                evapRate *= lossRatio;
                exfilRate *= lossRatio; 
            }
        }
    }

    // --- save evap & infil losses at the node
    sp->Storage[sp->Node[j].subIndex].evapLoss = evapRate * tStep;
    sp->Storage[sp->Node[j].subIndex].exfilLoss = exfilRate * tStep;
    return evapRate + exfilRate;
}


//=============================================================================
//                   D I V I D E R   M E T H O D S
//=============================================================================

int divider_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = node index
//           k = divider index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error message
//  Purpose: reads a flow divider's properties from a tokenized line of input.
//
//  Format of input line is:
//    nodeID  elev  divLink  TABULAR  curveID (optional params)
//    nodeID  elev  divLink  OVERFLOW (optional params)
//    nodeID  elev  divLink  CUTOFF  qCutoff (optional params)
//    nodeID  elev  divLink  WEIR    qMin  dhMax  cWeir (optional params)
//  where optional params are:
//    maxDepth  initDepth  surDepth  aPond    
//
{
    int    i, m, m1, m2, n;
    double x[11];
    char*  id;

    // --- get ID name
    if ( ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
    id = project_findID(NODE, tok[0]);
    if ( id == NULL ) return error_setInpError(ERR_NAME, tok[0]);

    // --- get invert elev.
    if ( ! getDouble(tok[1], &x[0]) ) return error_setInpError(ERR_NUMBER, tok[1]);

    // --- initialize parameter values
    for ( i=1; i<11; i++) x[i] = 0.0;

    // --- check if no diverted link supplied
    if ( strlen(tok[2]) == 0 || strcmp(tok[2], "*") == 0 ) x[1] = -1.0;

    // --- otherwise get index of diverted link
    else
    {
        m1 = project_findObject(LINK, tok[2]);
        if ( m1 < 0 ) return error_setInpError(ERR_NAME, tok[2]);
        x[1] = m1;
    }
    
    // --- get divider type
	n = 4;
    m1 = findmatch(tok[3], DividerTypeWords);
    if ( m1 < 0 ) return error_setInpError(ERR_KEYWORD, tok[3]);
    x[2] = m1;

    // --- get index of flow diversion curve for Tabular divider
    x[3] = -1;
    if ( m1 == TABULAR_DIVIDER )
    {
        if ( ntoks < 5 ) return error_setInpError(ERR_ITEMS, "");
        m2 = project_findObject(CURVE, tok[4]);
        if ( m2 < 0 ) return error_setInpError(ERR_NAME, tok[4]);
        x[3] = m2;
        n = 5;
    }

    // --- get cutoff flow for Cutoff divider
    if ( m1 == CUTOFF_DIVIDER )
    {
        if ( ntoks < 5 ) return error_setInpError(ERR_ITEMS, "");
        if ( ! getDouble(tok[4], &x[4]) )
            return error_setInpError(ERR_NUMBER, tok[4]);
        n = 5;
    }

    // --- get qmin, dhMax, & cWeir for Weir divider
    if ( m1 == WEIR_DIVIDER )
    {
        if ( ntoks < 7 ) return error_setInpError(ERR_ITEMS, "");
        for (i=4; i<7; i++)
             if ( ! getDouble(tok[i], &x[i]) )
                 return error_setInpError(ERR_NUMBER, tok[i]);
        n = 7;
    }

    // --- no parameters needed for Overflow divider
    if ( m1 == OVERFLOW_DIVIDER ) n = 4;

    // --- retrieve optional full depth, init. depth, surcharged depth
    //      & ponded area
    m = 7;
    for (i=n; i<ntoks && m<11; i++)
    {
        if ( ! getDouble(tok[i], &x[m]) )
        {
            return error_setInpError(ERR_NUMBER, tok[i]);
        }
        m++;
    }
 
    // --- add parameters to data base
    sp->Node[j].ID = id;
    node_setParams(sp, j, DIVIDER, k, x);
    return 0;
}

//=============================================================================

void  divider_validate(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: validates a flow divider's properties.
//
{
    int i, k;

    // --- check that diverted link is attached to divider
    k = sp->Node[j].subIndex;
    i = sp->Divider[k].link;
    if ( i < 0 || ( sp->Link[i].node1 != j && sp->Link[i].node2 != j) )
    {
        report_writeErrorMsg(sp, ERR_DIVIDER_LINK, sp->Node[j].ID);
    }

    // --- validate parameters supplied for weir-type divider
    if ( sp->Divider[k].type == WEIR_DIVIDER )
    {
        if ( sp->Divider[k].dhMax <= 0.0 || sp->Divider[k].cWeir <= 0.0 )
            report_writeErrorMsg(sp, ERR_WEIR_DIVIDER, sp->Node[j].ID);
        else
        {
            // --- find flow when weir is full
            sp->Divider[k].qMax = sp->Divider[k].cWeir * pow(sp->Divider[k].dhMax, 1.5)
                              / UCF(sp, FLOW);
            if ( sp->Divider[k].qMin > sp->Divider[k].qMax )
                report_writeErrorMsg(sp, ERR_WEIR_DIVIDER, sp->Node[j].ID);
        }
    }
}

//=============================================================================

double divider_getOutflow(SWMM_Project *sp, int j, int k)
//
//  Input:   j = node index
//           k = index of diversion link
//  Output:  returns diverted flow rate (cfs)
//  Purpose: computes flow sent through divider node into its diversion link.
//
//  NOTE: requires that links be previously sorted so that the non-diversion
//        link always gets evaluated before the diversion link
{
    int    i;                     // index of divider node
    int    m;                     // index of diverted flow table
    double qIn;                   // inflow to divider
    double qOut;                  // diverted outflow
    double f;                     // fraction of weir divider full

    qIn = sp->Node[j].inflow + sp->Node[j].overflow;
    i = sp->Node[j].subIndex;
    switch ( sp->Divider[i].type )
    {
      case CUTOFF_DIVIDER:
        if ( qIn <= sp->Divider[i].qMin ) qOut = 0.0;
        else qOut = qIn - sp->Divider[i].qMin;
        break;

      case OVERFLOW_DIVIDER:
        // --- outflow sent into non-diversion link is simply node's inflow
        if ( k != sp->Divider[i].link ) qOut = qIn;

        // --- diversion link receives any excess of node's inflow and
        //     outflow sent previously into non-diversion link
        else qOut = qIn - sp->Node[j].outflow;
        if ( qOut < FLOW_TOL ) qOut = 0.0;
        return qOut;

      case WEIR_DIVIDER:
        // --- no flow if inflow < qMin
        if ( qIn <= sp->Divider[i].qMin ) qOut = 0.0;

        // --- otherwise use weir eqn.
        else
        {
            // --- find fractional depth of flow over weir
            f = (qIn - sp->Divider[i].qMin) /
                (sp->Divider[i].qMax - sp->Divider[i].qMin);

            // --- if weir surcharged, use orifice eqn.
            if ( f > 1.0 ) qOut = sp->Divider[i].qMax * sqrt(f);
            
            // --- otherwise use weir eqn.
            else qOut = sp->Divider[i].cWeir *
                        pow(f*sp->Divider[i].dhMax, 1.5) / UCF(sp, FLOW);
        }
        break;

      case TABULAR_DIVIDER:
        m = sp->Divider[i].flowCurve;
        if ( m >= 0 )
            qOut = table_lookup(&sp->Curve[m], qIn * UCF(sp, FLOW)) / UCF(sp, FLOW);
        else qOut = 0.0;
        break;

      default: qOut = 0.0;
    }

    // --- make sure outflow doesn't exceed inflow
    if ( qOut > qIn ) qOut = qIn;

    // --- if link k not the diversion link, then re-define qOut as 
    //     the undiverted flow
    if ( k != sp->Divider[i].link )
    {
        qOut = qIn - qOut;
    }
    return qOut;
}


//=============================================================================
//                    O U T F A L L   M E T H O D S
//=============================================================================

int outfall_readParams(SWMM_Project *sp, int j, int k, char* tok[], int ntoks)
//
//  Input:   j = node index
//           k = outfall index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error message
//  Purpose: reads an outfall's properties from a tokenized line of input.
//
//  Format of input line is:
//    nodeID  elev  FIXED  outfallStage (flapGate) (routeTo)
//    nodeID  elev  TIDAL  curveID (flapGate) (routeTo)
//    nodeID  elev  TIMESERIES  tseriesID (flapGate) (routTo)
//    nodeID  elev  FREE (flapGate) (routeTo)
//    nodeID  elev  NORMAL (flapGate) (routeTo)
//
{
    int    i, m, n;
    double x[7];                                                               //(5.1.008)
    char*  id;

    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");
    id = project_findID(NODE, tok[0]);                      // node ID
    if ( id == NULL )
        return error_setInpError(ERR_NAME, tok[0]);
    if ( ! getDouble(tok[1], &x[0]) )                       // invert elev. 
        return error_setInpError(ERR_NUMBER, tok[1]);
    i = findmatch(tok[2], OutfallTypeWords);               // outfall type
    if ( i < 0 ) return error_setInpError(ERR_KEYWORD, tok[2]);
    x[1] = i;                                              // outfall type
    x[2] = 0.0;                                            // fixed stage
    x[3] = -1.;                                            // tidal curve
    x[4] = -1.;                                            // tide series
    x[5] = 0.;                                             // flap gate
    x[6] = -1.;                                            // route to subcatch//(5.1.008)

    n = 4;
    if ( i >= STAGED_OUTFALL )
    {
        if ( ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
        n = 5;
        switch ( i )
        {
        case STAGED_OUTFALL:                                // fixed stage
          if ( ! getDouble(tok[3], &x[2]) )
              return error_setInpError(ERR_NUMBER, tok[3]);
          break;
        case TIDAL_OUTFALL:                                // tidal curve
          m = project_findObject(CURVE, tok[3]);              
          if ( m < 0 ) return error_setInpError(ERR_NAME, tok[3]);
          x[3] = m;
          break;
        case TIMESERIES_OUTFALL:                           // stage time series
          m = project_findObject(TSERIES, tok[3]);            
          if ( m < 0 ) return error_setInpError(ERR_NAME, tok[3]);
          x[4] = m;
          sp->Tseries[m].refersTo = TIMESERIES_OUTFALL;
        }
    }
    if ( ntoks == n )
    {
        m = findmatch(tok[n-1], NoYesWords);               // flap gate
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, tok[n-1]);
        x[5] = m;
    }

////  Added for release 5.1.008.  ////                                         //(5.1.008)
    if ( ntoks == n+1)
    {
        m = project_findObject(SUBCATCH, tok[n]);
        if ( m < 0 ) return error_setInpError(ERR_NAME, tok[n]);
        x[6] = m;
    }
////

    sp->Node[j].ID = id;
    node_setParams(sp, j, OUTFALL, k, x);
    return 0;
}

//=============================================================================

void outfall_setOutletDepth(SWMM_Project *sp, int j, double yNorm, double yCrit, double z)
//
//  Input:   j = node index
//           yNorm = normal flow depth in outfall conduit (ft)
//           yCrit = critical flow depth in outfall conduit (ft)
//           z = height to outfall conduit invert (ft)
//  Output:  none
//  Purpose: sets water depth at an outfall node.
//
{
    double   x, y;                     // x,y values in table
    double   yNew;                     // new depth above invert elev. (ft)
    double   stage;                    // water elevation at outfall (ft)
    int      k;                        // table index
    int      i = sp->Node[j].subIndex;     // outfall index
    DateTime currentDate;              // current date/time in days

    switch ( sp->Outfall[i].type )
    {
      case FREE_OUTFALL:
        if ( z > 0.0 ) sp->Node[j].newDepth = 0.0;
        else sp->Node[j].newDepth = MIN(yNorm, yCrit);
        return;

      case NORMAL_OUTFALL:
        if ( z > 0.0 ) sp->Node[j].newDepth = 0.0;
        else sp->Node[j].newDepth = yNorm;
        return;

      case STAGED_OUTFALL:
        stage = sp->Outfall[i].outfallStage;
        break;

      case TIDAL_OUTFALL:
        k = sp->Outfall[i].tideCurve;
        table_getFirstEntry(sp, &sp->Curve[k], &x, &y);
        currentDate = sp->NewRoutingTime / MSECperDAY;
        x += ( currentDate - floor(currentDate) ) * 24.0;
        stage = table_lookup(&sp->Curve[k], x) / UCF(sp, LENGTH);
        break;

      case TIMESERIES_OUTFALL:
        k = sp->Outfall[i].stageSeries;
        currentDate = sp->StartDateTime + sp->NewRoutingTime / MSECperDAY;
        stage = table_tseriesLookup(sp, &sp->Tseries[k], currentDate, TRUE) /
                UCF(sp, LENGTH);
        break;
      default: stage = sp->Node[j].invertElev;
    }

    // --- now determine depth at node given outfall stage elev.
 
    // --- let critical flow depth be min. of critical & normal depth
    yCrit = MIN(yCrit, yNorm);

    // --- if elev. of critical depth is below outfall stage elev. then
    //     the outfall stage determines node depth
    if ( yCrit + z + sp->Node[j].invertElev < stage )
    {
        yNew = stage - sp->Node[j].invertElev;
    }

    // --- otherwise if the outfall conduit lies above the outfall invert
    else if ( z > 0.0 )
    {
        // --- if the outfall stage lies below the bottom of the outfall
        //     conduit then the result is distance from node invert to stage
        if ( stage < sp->Node[j].invertElev + z )
            yNew = MAX(0.0, (stage - sp->Node[j].invertElev));

        // --- otherwise stage lies between bottom of conduit and critical
        //     depth in conduit so result is elev. of critical depth
        else yNew = z + yCrit;
    }

    // --- and for case where there is no conduit offset and outfall stage
    //     lies below critical depth, then node depth = critical depth 
    else yNew = yCrit;
    sp->Node[j].newDepth = yNew;
}

//=============================================================================
