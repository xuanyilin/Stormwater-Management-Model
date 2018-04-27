//-----------------------------------------------------------------------------
//   dynwave.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (5.1.001)
//             03/28/14   (5.1.002)
//             09/15/14   (5.1.007)
//             03/19/15   (5.1.008)
//             08/01/16   (5.1.011)
//   Author:   L. Rossman (EPA)
//             M. Tryby (EPA)
//             R. Dickinson (CDM)
//
//   Dynamic wave flow routing functions.
//
//   This module solves the dynamic wave flow routing equations using
//   Picard Iterations (i.e., a method of successive approximations)
//   to solve the explicit form of the continuity and momentum equations
//   for conduits.
//
//   Build 5.1.002:
//   - Only non-ponded nodal surface area is saved for use in
//     surcharge algorithm.
//
//   Build 5.1.007:
//   - Node losses added to node outflow variable instead of treated
//     as a separate item when computing change in node flow volume.
//
//   Build 5.1.008:
//   - Module-specific constants moved here from project.c.
//   - Support added for user-specified minimum variable time step.
//   - Node crown elevations found here instead of in flowrout.c module.
//   - OpenMP use to parallelize findLinkFlows() & findNodeDepths().
//   - Bug in finding complete list of capacity limited links fixed.
//
//   Build 5.1.011:
//   - Added test for failed memory allocation.
//   - Fixed illegal array index bug for Ideal Pumps.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include "headers.h"
#include <stdlib.h>
#include <math.h>
#if defined(_OPENMP)
  #include <omp.h>                                                             //(5.1.008)
#endif

//-----------------------------------------------------------------------------
//     Constants 
//-----------------------------------------------------------------------------
static const double MINTIMESTEP =  0.001;   // min. time step (sec)            //(5.1.008)
static const double OMEGA       =  0.5;     // under-relaxation parameter

//  Constants moved here from project.c  //                                    //(5.1.008)
const double DEFAULT_SURFAREA  = 12.566; // Min. nodal surface area (~4 ft diam.)
const double DEFAULT_HEADTOL   = 0.005;  // Default head tolerance (ft)
const int    DEFAULT_MAXTRIALS = 8;      // Max. trials per time step


//-----------------------------------------------------------------------------
//  Data Structures
//-----------------------------------------------------------------------------
typedef struct 
{
    char    converged;                 // TRUE if iterations for a node done
    double  newSurfArea;               // current surface area (ft2)
    double  oldSurfArea;               // previous surface area (ft2)
    double  sumdqdh;                   // sum of dqdh from adjoining links
    double  dYdT;                      // change in depth w.r.t. time (ft/sec)
} TXnode;

//-----------------------------------------------------------------------------
//  Shared Variables
//-----------------------------------------------------------------------------
static double  VariableStep;           // size of variable time step (sec)
static TXnode* Xnode;                  // extended nodal information

static double  Omega;                  // actual under-relaxation parameter
static int     Steps;                  // number of Picard iterations

//-----------------------------------------------------------------------------
//  Function declarations
//-----------------------------------------------------------------------------
static void   initRoutingStep(SWMM_Project *sp);
static void   initNodeStates(SWMM_Project *sp);
static void   findBypassedLinks(SWMM_Project *sp);
static void   findLimitedLinks(SWMM_Project *sp);

static void   findLinkFlows(SWMM_Project *sp, double dt);
static int    isTrueConduit(SWMM_Project *sp, int link);
static void   findNonConduitFlow(SWMM_Project *sp, int link, double dt);
static void   findNonConduitSurfArea(SWMM_Project *sp, int link);
static double getModPumpFlow(SWMM_Project *sp, int link, double q, double dt);
static void   updateNodeFlows(SWMM_Project *sp, int link);

static int    findNodeDepths(SWMM_Project *sp, double dt);
static void   setNodeDepth(SWMM_Project *sp, int node, double dt);
static double getFloodedDepth(SWMM_Project *sp, int node, int canPond, double dV,
        double yNew, double yMax, double dt);

static double getVariableStep(SWMM_Project *sp, double maxStep);
static double getLinkStep(SWMM_Project *sp, double tMin, int *minLink);
static double getNodeStep(SWMM_Project *sp, double tMin, int *minNode);

//=============================================================================

////  This function was modified for release 5.1.008.  ////                    //(5.1.008)

void dynwave_init(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes dynamic wave routing method.
//
{
    int i, j;
    double z;

    VariableStep = 0.0;
    Xnode = (TXnode *) calloc(sp->Nobjects[NODE], sizeof(TXnode));

////  Added to release 5.1.011.  ////                                          //(5.1.011)
    if ( Xnode == NULL )
    {
        report_writeErrorMsg(sp, ERR_MEMORY,
            " Not enough memory for dynamic wave routing.");
        return;
    }
//////////////////////////////////////

    // --- initialize node surface areas & crown elev.
    for (i = 0; i < sp->Nobjects[NODE]; i++ )
    {
        Xnode[i].newSurfArea = 0.0;
        Xnode[i].oldSurfArea = 0.0;
        sp->Node[i].crownElev = sp->Node[i].invertElev;
    }

    // --- update node crown elev. & initialize links
    for (i = 0; i < sp->Nobjects[LINK]; i++)
    {
        j = sp->Link[i].node1;
        z = sp->Node[j].invertElev + sp->Link[i].offset1 + sp->Link[i].xsect.yFull;
        sp->Node[j].crownElev = MAX(sp->Node[j].crownElev, z);
        j = sp->Link[i].node2;
        z = sp->Node[j].invertElev + sp->Link[i].offset2 + sp->Link[i].xsect.yFull;
        sp->Node[j].crownElev = MAX(sp->Node[j].crownElev, z);
        sp->Link[i].flowClass = DRY;
        sp->Link[i].dqdh = 0.0;
    }
}

//=============================================================================

void  dynwave_close()
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory allocated for dynamic wave routing method.
//
{
    FREE(Xnode);
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

void dynwave_validate(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: adjusts dynamic wave routing options.
//
{
    if ( sp->MinRouteStep > sp->RouteStep ) sp->MinRouteStep = sp->RouteStep;
    if ( sp->MinRouteStep < MINTIMESTEP ) sp->MinRouteStep = MINTIMESTEP;
	if ( sp->MinSurfArea == 0.0 ) sp->MinSurfArea = DEFAULT_SURFAREA;
	else sp->MinSurfArea /= UCF(sp, LENGTH) * UCF(sp, LENGTH);
    if ( sp->HeadTol == 0.0 ) sp->HeadTol = DEFAULT_HEADTOL;
    else sp->HeadTol /= UCF(sp, LENGTH);
	if ( sp->MaxTrials == 0 ) sp->MaxTrials = DEFAULT_MAXTRIALS;
}

//=============================================================================

double dynwave_getRoutingStep(SWMM_Project *sp, double fixedStep)
//
//  Input:   fixedStep = user-supplied fixed time step (sec)
//  Output:  returns routing time step (sec)
//  Purpose: computes variable routing time step if applicable.
//
{
    // --- use user-supplied fixed step if variable step option turned off
    //     or if its smaller than the min. allowable variable time step
    if ( sp->CourantFactor == 0.0 ) return fixedStep;
    if ( fixedStep < MINTIMESTEP ) return fixedStep;

    // --- at start of simulation (when current variable step is zero)
    //     use the minimum allowable time step
    if ( VariableStep == 0.0 )
    {
        VariableStep = sp->MinRouteStep;                                           //(5.1.008)
    }

    // --- otherwise compute variable step based on current flow solution
    else VariableStep = getVariableStep(sp, fixedStep);

    // --- adjust step to be a multiple of a millisecond
    VariableStep = floor(1000.0 * VariableStep) / 1000.0;
    return VariableStep;
}

//=============================================================================

int dynwave_execute(SWMM_Project *sp, double tStep)
//
//  Input:   links = array of topo sorted links indexes
//           tStep = time step (sec)
//  Output:  returns number of iterations used
//  Purpose: routes flows through drainage network over current time step.
//
{
    int converged;

    // --- initialize
    if ( sp->ErrorCode ) return 0;
    Steps = 0;
    converged = FALSE;
    Omega = OMEGA;
    initRoutingStep(sp);

    // --- keep iterating until convergence 
    while ( Steps < sp->MaxTrials )
    {
        // --- execute a routing step & check for nodal convergence
        initNodeStates(sp);
        findLinkFlows(sp, tStep);
        converged = findNodeDepths(sp, tStep);
        Steps++;
        if ( Steps > 1 )
        {
            if ( converged ) break;

            // --- check if link calculations can be skipped in next step
            findBypassedLinks(sp);
        }
    }
    if ( !converged ) sp->NonConvergeCount++;

    //  --- identify any capacity-limited conduits
    findLimitedLinks(sp);
    return Steps;
}

//=============================================================================

void   initRoutingStep(SWMM_Project *sp)
{
    int i;
    for (i = 0; i < sp->Nobjects[NODE]; i++)
    {
        Xnode[i].converged = FALSE;
        Xnode[i].dYdT = 0.0;
    }
    for (i = 0; i < sp->Nobjects[LINK]; i++)
    {
        sp->Link[i].bypassed = FALSE;
        sp->Link[i].surfArea1 = 0.0;
        sp->Link[i].surfArea2 = 0.0;
    }

    // --- a2 preserves conduit area from solution at last time step
    for ( i = 0; i < sp->Nlinks[CONDUIT]; i++) sp->Conduit[i].a2 = sp->Conduit[i].a1;
}

//=============================================================================

void initNodeStates(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes node's surface area, inflow & outflow
//
{
    int i;

    for (i = 0; i < sp->Nobjects[NODE]; i++)
    {
        // --- initialize nodal surface area
        if ( sp->AllowPonding )
        {
            Xnode[i].newSurfArea = node_getPondedArea(sp, i, sp->Node[i].newDepth);
        }
        else
        {
            Xnode[i].newSurfArea = node_getSurfArea(sp, i, sp->Node[i].newDepth);
        }
        if ( Xnode[i].newSurfArea < sp->MinSurfArea )
        {
            Xnode[i].newSurfArea = sp->MinSurfArea;
        }

////  Following code section modified for release 5.1.007  ////                //(5.1.007)
        // --- initialize nodal inflow & outflow
        sp->Node[i].inflow = 0.0;
        sp->Node[i].outflow = sp->Node[i].losses;
        if ( sp->Node[i].newLatFlow >= 0.0 )
        {    
            sp->Node[i].inflow += sp->Node[i].newLatFlow;
        }
        else
        {    
            sp->Node[i].outflow -= sp->Node[i].newLatFlow;
        }
        Xnode[i].sumdqdh = 0.0;
    }
}

//=============================================================================

void   findBypassedLinks(SWMM_Project *sp)
{
    int i;
    for (i = 0; i < sp->Nobjects[LINK]; i++)
    {
        if ( Xnode[sp->Link[i].node1].converged &&
             Xnode[sp->Link[i].node2].converged )
             sp->Link[i].bypassed = TRUE;
        else sp->Link[i].bypassed = FALSE;
    }
}

//=============================================================================

void  findLimitedLinks(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: determines if a conduit link is capacity limited.
//
{
    int    j, n1, n2, k;
    double h1, h2;

    for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        // ---- check only non-dummy conduit links
        if ( !isTrueConduit(sp, j) ) continue;                                     //(5.1.008)

        // --- check that upstream end is full
        k = sp->Link[j].subIndex;
        sp->Conduit[k].capacityLimited = FALSE;
        if ( sp->Conduit[k].a1 >= sp->Link[j].xsect.aFull )
        {
            // --- check if HGL slope > conduit slope
            n1 = sp->Link[j].node1;
            n2 = sp->Link[j].node2;
            h1 = sp->Node[n1].newDepth + sp->Node[n1].invertElev;
            h2 = sp->Node[n2].newDepth + sp->Node[n2].invertElev;
            if ( (h1 - h2) > fabs(sp->Conduit[k].slope) * sp->Conduit[k].length )
                sp->Conduit[k].capacityLimited = TRUE;
        }
    }
}

//=============================================================================

void findLinkFlows(SWMM_Project *sp, double dt)
{
    int i;

    // --- find new flow in each non-dummy conduit
#pragma omp parallel num_threads(sp->NumThreads)                                   //(5.1.008)
{
    #pragma omp for                                                            //(5.1.008)
    for ( i = 0; i < sp->Nobjects[LINK]; i++)
    {
        if ( isTrueConduit(sp, i) && !sp->Link[i].bypassed )
            dwflow_findConduitFlow(sp, i, Steps, Omega, dt);
    }
}

    // --- update inflow/outflows for nodes attached to non-dummy conduits
    for ( i = 0; i < sp->Nobjects[LINK]; i++)
    {
        if ( isTrueConduit(sp, i) ) updateNodeFlows(sp, i);
    }

    // --- find new flows for all dummy conduits, pumps & regulators
    for ( i = 0; i < sp->Nobjects[LINK]; i++)
    {
        if ( !isTrueConduit(sp, i) )
        {	
            if ( !sp->Link[i].bypassed ) findNonConduitFlow(sp, i, dt);
            updateNodeFlows(sp, i);
        }
    }
}

//=============================================================================

int isTrueConduit(SWMM_Project *sp, int j)
{
    return ( sp->Link[j].type == CONDUIT && sp->Link[j].xsect.type != DUMMY );
}

//=============================================================================

void findNonConduitFlow(SWMM_Project *sp, int i, double dt)
//
//  Input:   i = link index
//           dt = time step (sec)
//  Output:  none
//  Purpose: finds new flow in a non-conduit-type link
//
{
    double qLast;                      // previous link flow (cfs)
    double qNew;                       // new link flow (cfs)

    // --- get link flow from last iteration
    qLast = sp->Link[i].newFlow;
    sp->Link[i].dqdh = 0.0;

    // --- get new inflow to link from its upstream node
    //     (link_getInflow returns 0 if flap gate closed or pump is offline)
    qNew = link_getInflow(sp, i);
    if ( sp->Link[i].type == PUMP ) qNew = getModPumpFlow(sp, i, qNew, dt);

    // --- find surface area at each end of link
    findNonConduitSurfArea(sp, i);

    // --- apply under-relaxation with flow from previous iteration;
    // --- do not allow flow to change direction without first being 0
    if ( Steps > 0 && sp->Link[i].type != PUMP ) 
    {
        qNew = (1.0 - Omega) * qLast + Omega * qNew;
        if ( qNew * qLast < 0.0 ) qNew = 0.001 * SGN(qNew);
    }
    sp->Link[i].newFlow = qNew;
}

//=============================================================================

double getModPumpFlow(SWMM_Project *sp, int i, double q, double dt)
//
//  Input:   i = link index
//           q = pump flow from pump curve (cfs)
//           dt = time step (sec)
//  Output:  returns modified pump flow rate (cfs)
//  Purpose: modifies pump curve pumping rate depending on amount of water
//           available at pump's inlet node.
//
{
    int    j = sp->Link[i].node1;          // pump's inlet node index
    int    k = sp->Link[i].subIndex;       // pump's index
    double newNetInflow;               // inflow - outflow rate (cfs)
    double netFlowVolume;              // inflow - outflow volume (ft3)
    double y;                          // node depth (ft)

    if ( q == 0.0 ) return q;

    // --- case where inlet node is a storage node: 
    //     prevent node volume from going negative
    if ( sp->Node[j].type == STORAGE ) return node_getMaxOutflow(sp, j, q, dt);

    // --- case where inlet is a non-storage node
    switch ( sp->Pump[k].type )
    {
      // --- for Type1 pump, a volume is computed for inlet node,
      //     so make sure it doesn't go negative
      case TYPE1_PUMP:
        return node_getMaxOutflow(sp, j, q, dt);

      // --- for other types of pumps, if pumping rate would make depth
      //     at upstream node negative, then set pumping rate = inflow
      case TYPE2_PUMP:
      case TYPE4_PUMP:
      case TYPE3_PUMP:
         newNetInflow = sp->Node[j].inflow - sp->Node[j].outflow - q;
         netFlowVolume = 0.5 * (sp->Node[j].oldNetInflow + newNetInflow ) * dt;
         y = sp->Node[j].oldDepth + netFlowVolume / Xnode[j].newSurfArea;
         if ( y <= 0.0 ) return sp->Node[j].inflow;
    }
    return q;
}

//=============================================================================

void  findNonConduitSurfArea(SWMM_Project *sp, int i)
//
//  Input:   i = link index
//  Output:  none
//  Purpose: finds the surface area contributed by a non-conduit
//           link to its upstream and downstream nodes.
//
{
    if ( sp->Link[i].type == ORIFICE )
    {
        sp->Link[i].surfArea1 = sp->Orifice[sp->Link[i].subIndex].surfArea / 2.;
    }

    // --- no surface area for weirs to maintain SWMM 4 compatibility
/*
    else if ( sp->Link[i].type == WEIR )
    {
        Xlink[i].surfArea1 = Weir[sp->Link[i].subIndex].surfArea / 2.;
    }
*/

    else sp->Link[i].surfArea1 = 0.0;
    sp->Link[i].surfArea2 = sp->Link[i].surfArea1;
    if ( sp->Link[i].flowClass == UP_CRITICAL ||
        sp->Node[sp->Link[i].node1].type == STORAGE ) sp->Link[i].surfArea1 = 0.0;
    if ( sp->Link[i].flowClass == DN_CRITICAL ||
        sp->Node[sp->Link[i].node2].type == STORAGE ) sp->Link[i].surfArea2 = 0.0;
}

//=============================================================================

void updateNodeFlows(SWMM_Project *sp, int i)
//
//  Input:   i = link index
//           q = link flow rate (cfs)
//  Output:  none
//  Purpose: updates cumulative inflow & outflow at link's end nodes.
//
{
    int    k;                                                                  //(5.1.011)
    int    barrels = 1;
    int    n1 = sp->Link[i].node1;
    int    n2 = sp->Link[i].node2;
    double q = sp->Link[i].newFlow;
    double uniformLossRate = 0.0;

    // --- compute any uniform seepage loss from a conduit
    if ( sp->Link[i].type == CONDUIT )
    {
        k = sp->Link[i].subIndex;
        uniformLossRate = sp->Conduit[k].evapLossRate + sp->Conduit[k].seepLossRate; 
        barrels = sp->Conduit[k].barrels;
    }

    // --- update total inflow & outflow at upstream/downstream nodes
    if ( q >= 0.0 )
    {
        sp->Node[n1].outflow += q + uniformLossRate;
        sp->Node[n2].inflow  += q;
    }
    else
    {
        sp->Node[n1].inflow   -= q;
        sp->Node[n2].outflow  -= q - uniformLossRate;
    }

    // --- add surf. area contributions to upstream/downstream nodes
    Xnode[sp->Link[i].node1].newSurfArea += sp->Link[i].surfArea1 * barrels;
    Xnode[sp->Link[i].node2].newSurfArea += sp->Link[i].surfArea2 * barrels;

    // --- update summed value of dqdh at each end node
    Xnode[sp->Link[i].node1].sumdqdh += sp->Link[i].dqdh;
    if ( sp->Link[i].type == PUMP )
    {
        k = sp->Link[i].subIndex;
        if ( sp->Pump[k].type != TYPE4_PUMP )                                      //(5.1.011)
        {
            Xnode[n2].sumdqdh += sp->Link[i].dqdh;
        }
    }
    else Xnode[n2].sumdqdh += sp->Link[i].dqdh;
}

//=============================================================================

int findNodeDepths(SWMM_Project *sp, double dt)
{
    int i;
    int converged;      // convergence flag
    double yOld;        // previous node depth (ft)

    // --- compute outfall depths based on flow in connecting link
    for ( i = 0; i < sp->Nobjects[LINK]; i++ ) link_setOutfallDepth(sp, i);

    // --- compute new depth for all non-outfall nodes and determine if
    //     depth change from previous iteration is below tolerance
    converged = TRUE;
#pragma omp parallel num_threads(sp->NumThreads)                                   //(5.1.008)
{
    #pragma omp for private(yOld)                                              //(5.1.008)
    for ( i = 0; i < sp->Nobjects[NODE]; i++ )
    {
        if ( sp->Node[i].type == OUTFALL ) continue;
        yOld = sp->Node[i].newDepth;
        setNodeDepth(sp, i, dt);
        Xnode[i].converged = TRUE;
        if ( fabs(yOld - sp->Node[i].newDepth) > sp->HeadTol )
        {
            converged = FALSE;
            Xnode[i].converged = FALSE;
        }
    }
}                                                                              //(5.1.008)
    return converged;
}

//=============================================================================

void setNodeDepth(SWMM_Project *sp, int i, double dt)
//
//  Input:   i  = node index
//           dt = time step (sec)
//  Output:  none
//  Purpose: sets depth at non-outfall node after current time step.
//
{
    int     canPond;                   // TRUE if node can pond overflows
    int     isPonded;                  // TRUE if node is currently ponded 
    double  dQ;                        // inflow minus outflow at node (cfs)
    double  dV;                        // change in node volume (ft3)
    double  dy;                        // change in node depth (ft)
    double  yMax;                      // max. depth at node (ft)
    double  yOld;                      // node depth at previous time step (ft)
    double  yLast;                     // previous node depth (ft)
    double  yNew;                      // new node depth (ft)
    double  yCrown;                    // depth to node crown (ft)
    double  surfArea;                  // node surface area (ft2)
    double  denom;                     // denominator term
    double  corr;                      // correction factor
    double  f;                         // relative surcharge depth

    // --- see if node can pond water above it
    canPond = (sp->AllowPonding && sp->Node[i].pondedArea > 0.0);
    isPonded = (canPond && sp->Node[i].newDepth > sp->Node[i].fullDepth);

    // --- initialize values
    yCrown = sp->Node[i].crownElev - sp->Node[i].invertElev;
    yOld = sp->Node[i].oldDepth;
    yLast = sp->Node[i].newDepth;
    sp->Node[i].overflow = 0.0;
    surfArea = Xnode[i].newSurfArea;

    // --- determine average net flow volume into node over the time step
    dQ = sp->Node[i].inflow - sp->Node[i].outflow;
    dV = 0.5 * (sp->Node[i].oldNetInflow + dQ) * dt;

    // --- if node not surcharged, base depth change on surface area        
    if ( yLast <= yCrown || sp->Node[i].type == STORAGE || isPonded )
    {
        dy = dV / surfArea;
        yNew = yOld + dy;

        // --- save non-ponded surface area for use in surcharge algorithm     //(5.1.002)
        if ( !isPonded ) Xnode[i].oldSurfArea = surfArea;                      //(5.1.002)

        // --- apply under-relaxation to new depth estimate
        if ( Steps > 0 )
        {
            yNew = (1.0 - Omega) * yLast + Omega * yNew;
        }

        // --- don't allow a ponded node to drop much below full depth
        if ( isPonded && yNew < sp->Node[i].fullDepth )
            yNew = sp->Node[i].fullDepth - FUDGE;
    }

    // --- if node surcharged, base depth change on dqdh
    //     NOTE: depth change is w.r.t depth from previous
    //     iteration; also, do not apply under-relaxation.
    else
    {
        // --- apply correction factor for upstream terminal nodes
        corr = 1.0;
        if ( sp->Node[i].degree < 0 ) corr = 0.6;

        // --- allow surface area from last non-surcharged condition
        //     to influence dqdh if depth close to crown depth
        denom = Xnode[i].sumdqdh;
        if ( yLast < 1.25 * yCrown )
        {
            f = (yLast - yCrown) / yCrown;
            denom += (Xnode[i].oldSurfArea/dt -
                      Xnode[i].sumdqdh) * exp(-15.0 * f);
        }

        // --- compute new estimate of node depth
        if ( denom == 0.0 ) dy = 0.0;
        else dy = corr * dQ / denom;
        yNew = yLast + dy;
        if ( yNew < yCrown ) yNew = yCrown - FUDGE;

        // --- don't allow a newly ponded node to rise much above full depth
        if ( canPond && yNew > sp->Node[i].fullDepth )
            yNew = sp->Node[i].fullDepth + FUDGE;
    }

    // --- depth cannot be negative
    if ( yNew < 0 ) yNew = 0.0;

    // --- determine max. non-flooded depth
    yMax = sp->Node[i].fullDepth;
    if ( canPond == FALSE ) yMax += sp->Node[i].surDepth;

    // --- find flooded depth & volume
    if ( yNew > yMax )
    {
        yNew = getFloodedDepth(sp, i, canPond, dV, yNew, yMax, dt);
    }
    else sp->Node[i].newVolume = node_getVolume(sp, i, yNew);

    // --- compute change in depth w.r.t. time
    Xnode[i].dYdT = fabs(yNew - yOld) / dt;

    // --- save new depth for node
    sp->Node[i].newDepth = yNew;
}

//=============================================================================

double getFloodedDepth(SWMM_Project *sp, int i, int canPond, double dV,
        double yNew, double yMax, double dt)
//
//  Input:   i  = node index
//           canPond = TRUE if water can pond over node
//           isPonded = TRUE if water is currently ponded
//           dV = change in volume over time step (ft3)
//           yNew = current depth at node (ft)
//           yMax = max. depth at node before ponding (ft)
//           dt = time step (sec)
//  Output:  returns depth at node when flooded (ft)
//  Purpose: computes depth, volume and overflow for a flooded node.
//
{
    if ( canPond == FALSE )
    {
        sp->Node[i].overflow = dV / dt;
        sp->Node[i].newVolume = sp->Node[i].fullVolume;
        yNew = yMax;
    }
    else
    {
        sp->Node[i].newVolume = MAX((sp->Node[i].oldVolume+dV), sp->Node[i].fullVolume);
        sp->Node[i].overflow = (sp->Node[i].newVolume - 
            MAX(sp->Node[i].oldVolume, sp->Node[i].fullVolume)) / dt;
    }
    if ( sp->Node[i].overflow < FUDGE ) sp->Node[i].overflow = 0.0;
    return yNew;

}

//=============================================================================

double getVariableStep(SWMM_Project *sp, double maxStep)
//
//  Input:   maxStep = user-supplied max. time step (sec)
//  Output:  returns time step (sec)
//  Purpose: finds time step that satisfies stability criterion but
//           is no greater than the user-supplied max. time step.
//
{
    int    minLink = -1;                // index of link w/ min. time step
    int    minNode = -1;                // index of node w/ min. time step
    double tMin;                        // allowable time step (sec)
    double tMinLink;                    // allowable time step for links (sec)
    double tMinNode;                    // allowable time step for nodes (sec)

    // --- find stable time step for links & then nodes
    tMin = maxStep;
    tMinLink = getLinkStep(sp, tMin, &minLink);
    tMinNode = getNodeStep(sp, tMinLink, &minNode);

    // --- use smaller of the link and node time step
    tMin = tMinLink;
    if ( tMinNode < tMin )
    {
        tMin = tMinNode ;
        minLink = -1;
    }

    // --- update count of times the minimum node or link was critical
    stats_updateCriticalTimeCount(minNode, minLink);

    // --- don't let time step go below an absolute minimum
    if ( tMin < sp->MinRouteStep ) tMin = sp->MinRouteStep;                            //(5.1.008)
    return tMin;
}

//=============================================================================

double getLinkStep(SWMM_Project *sp, double tMin, int *minLink)
//
//  Input:   tMin = critical time step found so far (sec)
//  Output:  minLink = index of link with critical time step;
//           returns critical time step (sec)
//  Purpose: finds critical time step for conduits based on Courant criterion.
//
{
    int    i;                           // link index
    int    k;                           // conduit index
    double q;                           // conduit flow (cfs)
    double t;                           // time step (sec)
    double tLink = tMin;                // critical link time step (sec)

    // --- examine each conduit link
    for ( i = 0; i < sp->Nobjects[LINK]; i++ )
    {
        if ( sp->Link[i].type == CONDUIT )
        {
            // --- skip conduits with negligible flow, area or Fr
            k = sp->Link[i].subIndex;
            q = fabs(sp->Link[i].newFlow) / sp->Conduit[k].barrels;
            if ( q <= 0.05 * sp->Link[i].qFull
            ||   sp->Conduit[k].a1 <= FUDGE
            ||   sp->Link[i].froude <= 0.01 
               ) continue;

            // --- compute time step to satisfy Courant condition
            t = sp->Link[i].newVolume / sp->Conduit[k].barrels / q;
            t = t * sp->Conduit[k].modLength / link_getLength(sp, i);
            t = t * sp->Link[i].froude / (1.0 + sp->Link[i].froude) * sp->CourantFactor;

            // --- update critical link time step
            if ( t < tLink )
            {
                tLink = t;
                *minLink = i;
            }
        }
    }
    return tLink;
}

//=============================================================================

double getNodeStep(SWMM_Project *sp, double tMin, int *minNode)
//
//  Input:   tMin = critical time step found so far (sec)
//  Output:  minNode = index of node with critical time step;
//           returns critical time step (sec)
//  Purpose: finds critical time step for nodes based on max. allowable
//           projected change in depth.
//
{
    int    i;                           // node index
    double maxDepth;                    // max. depth allowed at node (ft)
    double dYdT;                        // change in depth per unit time (ft/sec)
    double t1;                          // time needed to reach depth limit (sec)
    double tNode = tMin;                // critical node time step (sec)

    // --- find smallest time so that estimated change in nodal depth
    //     does not exceed safety factor * maxdepth
    for ( i = 0; i < sp->Nobjects[NODE]; i++ )
    {
        // --- see if node can be skipped
        if ( sp->Node[i].type == OUTFALL ) continue;
        if ( sp->Node[i].newDepth <= FUDGE) continue;
        if ( sp->Node[i].newDepth  + FUDGE >=
             sp->Node[i].crownElev - sp->Node[i].invertElev ) continue;

        // --- define max. allowable depth change using crown elevation
        maxDepth = (sp->Node[i].crownElev - sp->Node[i].invertElev) * 0.25;
        if ( maxDepth < FUDGE ) continue;
        dYdT = Xnode[i].dYdT;
        if (dYdT < FUDGE ) continue;

        // --- compute time to reach max. depth & compare with critical time
        t1 = maxDepth / dYdT;
        if ( t1 < tNode )
        {
            tNode = t1;
            *minNode = i;
        }
    }
    return tNode;
}
