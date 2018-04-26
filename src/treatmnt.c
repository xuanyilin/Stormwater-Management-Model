//-----------------------------------------------------------------------------
//   treatmnt.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             03/19/15   (Build 5.1.008)
//   Author:   L. Rossman
//
//   Pollutant treatment functions.
//
//   Build 5.1.008:
//   - A bug in evaluating recursive calls to treatment functions was fixed. 
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const int PVMAX = 5;            // number of process variables
enum   ProcessVarType {pvHRT,          // hydraulic residence time
                       pvDT,           // time step duration
                       pvFLOW,         // flow rate
                       pvDEPTH,        // water height above invert
                       pvAREA};        // storage surface area

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static int     ErrCode;                // treatment error code
static int     J;                      // index of node being analyzed
static double  Dt;                     // curent time step (sec)
static double  Q;                      // node inflow (cfs)
static double  V;                      // node volume (ft3)
static double* R;                      // array of pollut. removals
static double* Cin;                    // node inflow concentrations
//static TTreatment* Treatment; // defined locally in treatmnt_treat()         //(5.1.008)

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  treatmnt_open           (called from routing_open)
//  treatment_close         (called from routing_close)
//  treatmnt_readExpression (called from parseLine in input.c)
//  treatmnt_delete         (called from deleteObjects in project.c)
//  treatmnt_setInflow      (called from qualrout_execute)
//  treatmnt_treat          (called from findNodeQual in qualrout.c)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    createTreatment(SWMM_Project *sp, int node);
static double getRemoval(SWMM_Project *sp, int pollut);
static int    getVariableIndex(SWMM_Project *sp, char* s);
static double getVariableValue(SWMM_Project *sp, int varCode);


//=============================================================================

int  treatmnt_open(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: allocates memory for computing pollutant removals by treatment.
//
{
    R = NULL;
    Cin = NULL;
    if ( sp->Nobjects[POLLUT] > 0 )
    {
        R = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        Cin = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        if ( R == NULL || Cin == NULL)
        {
            report_writeErrorMsg(sp, ERR_MEMORY, "");
            return FALSE;
        }
    }
    return TRUE;
}

//=============================================================================

void treatmnt_close(void)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: frees memory used for computing pollutant removals by treatment.
//
{
    FREE(R);
    FREE(Cin);
}

//=============================================================================

int  treatmnt_readExpression(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads a treatment expression from a tokenized line of input.
//
{
    char  s[MAXLINE+1];
    char* expr;
    int   i, j, k, p;
    MathExpr* equation;                // ptr. to a math. expression

    // --- retrieve node & pollutant
    if ( ntoks < 3 ) return error_setInpError(ERR_ITEMS, "");
    j = project_findObject(NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(ERR_NAME, tok[0]);
    p = project_findObject(POLLUT, tok[1]);
    if ( p < 0 ) return error_setInpError(ERR_NAME, tok[1]);

    // --- concatenate remaining tokens into a single string
    strcpy(s, tok[2]);
    for ( i=3; i<ntoks; i++)
    {
        strcat(s, " ");
        strcat(s, tok[i]);
    }

    // --- check treatment type
    if      ( UCHAR(s[0]) == 'R' ) k = 0;
    else if ( UCHAR(s[0]) == 'C' ) k = 1;
    else return error_setInpError(ERR_KEYWORD, tok[2]);

    // --- start treatment expression after equals sign
    expr = strchr(s, '=');
    if ( expr == NULL ) return error_setInpError(ERR_KEYWORD, "");
    else expr++;

    // --- create treatment objects at node j if they don't already exist
    if ( sp->Node[j].treatment == NULL )
    {
        if ( !createTreatment(sp, j) ) return error_setInpError(ERR_MEMORY, "");
    }

    // --- create a parsed expression tree from the string expr
    //     (getVariableIndex is the function that converts a treatment
    //      variable's name into an index number) 
    equation = mathexpr_create(sp, expr, getVariableIndex);
    if ( equation == NULL )
        return error_setInpError(ERR_TREATMENT_EXPR, "");

    // --- save the treatment parameters in the node's treatment object
    sp->Node[j].treatment[p].treatType = k;
    sp->Node[j].treatment[p].equation = equation;
    return 0;
}

//=============================================================================

void treatmnt_delete(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: deletes the treatment objects for each pollutant at a node.
//
{
    int p;
    if ( sp->Node[j].treatment )
    {
        for (p=0; p<sp->Nobjects[POLLUT]; p++)
            mathexpr_delete(sp->Node[j].treatment[p].equation);
        free(sp->Node[j].treatment);
    }
    sp->Node[j].treatment = NULL;
}

//=============================================================================

void  treatmnt_setInflow(SWMM_Project *sp, double qIn, double wIn[])
//
//  Input:   j = node index
//           qIn = flow inflow rate (cfs)
//           wIn = pollutant mass inflow rate (mass/sec)
//  Output:  none
//  Purpose: computes and saves array of inflow concentrations to a node.
//
{
    int    p;
    if ( qIn > 0.0 )
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) Cin[p] = wIn[p]/qIn;
    else
        for (p = 0; p < sp->Nobjects[POLLUT]; p++) Cin[p] = 0.0;
}

//=============================================================================

void  treatmnt_treat(SWMM_Project *sp, int j, double q, double v, double tStep)
//
//  Input:   j     = node index
//           q     = inflow to node (cfs)
//           v     = volume of node (ft3)
//           tStep = routing time step (sec)
//  Output:  none
//  Purpose: updates pollutant concentrations at a node after treatment.
//
{
    int    p;                          // pollutant index
    double cOut;                       // concentration after treatment
    double massLost;                   // mass lost by treatment per time step
    TTreatment* treatment;             // pointer to treatment object          //(5.1.008)

    // --- set locally shared variables for node j
    if ( sp->Node[j].treatment == NULL ) return;
    ErrCode = 0;
    J  = j;                            // current node
    Dt = tStep;                        // current time step
    Q  = q;                            // current inflow rate
    V  = v;                            // current node volume

    // --- initialze each removal to indicate no value 
    for ( p = 0; p < sp->Nobjects[POLLUT]; p++) R[p] = -1.0;

    // --- determine removal of each pollutant
    for ( p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        // --- removal is zero if there is no treatment equation
        treatment = &sp->Node[j].treatment[p];                                     //(5.1.008)
        if ( treatment->equation == NULL ) R[p] = 0.0;                         //(5.1.008)

        // --- no removal for removal-type expression when there is no inflow 
	    else if ( treatment->treatType == REMOVAL && q <= ZERO ) R[p] = 0.0;   //(5.1.008)

        // --- otherwise evaluate the treatment expression to find R[p]
        else getRemoval(sp, p);
    }

    // --- check for error condition
    if ( ErrCode == ERR_CYCLIC_TREATMENT )
    {
         report_writeErrorMsg(sp, ERR_CYCLIC_TREATMENT, sp->Node[J].ID);
    }

    // --- update nodal concentrations and mass balances
    else for ( p = 0; p < sp->Nobjects[POLLUT]; p++ )
    {
        if ( R[p] == 0.0 ) continue;
        treatment = &sp->Node[j].treatment[p];                                     //(5.1.008)

        // --- removal-type treatment equations get applied to inflow stream

        if ( treatment->treatType == REMOVAL )                                 //(5.1.008)
        {
            // --- if no pollutant in inflow then cOut is current nodal concen.
            if ( Cin[p] == 0.0 ) cOut = sp->Node[j].newQual[p];

            // ---  otherwise apply removal to influent concen.
            else cOut = (1.0 - R[p]) * Cin[p];

            // --- cOut can't be greater than mixture concen. at node
            //     (i.e., in case node is a storage unit) 
            cOut = MIN(cOut, sp->Node[j].newQual[p]);
        }

        // --- concentration-type equations get applied to nodal concentration
        else
        {
            cOut = (1.0 - R[p]) * sp->Node[j].newQual[p];
        }

        // --- mass lost must account for any initial mass in storage 
        massLost = (Cin[p]*q*tStep + sp->Node[j].oldQual[p]*sp->Node[j].oldVolume - 
                   cOut*(q*tStep + sp->Node[j].oldVolume)) / tStep; 
        massLost = MAX(0.0, massLost); 

        // --- add mass loss to mass balance totals and revise nodal concentration
        massbal_addReactedMass(sp, p, massLost);
        sp->Node[j].newQual[p] = cOut;
    }
}

//=============================================================================

int  createTreatment(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: creates a treatment object for each pollutant at a node.
//
{
    int p;
    sp->Node[j].treatment = (TTreatment *) calloc(sp->Nobjects[POLLUT],
                                              sizeof(TTreatment));
    if ( sp->Node[j].treatment == NULL )
    {
        return FALSE;
    }
    for (p = 0; p < sp->Nobjects[POLLUT]; p++)
    {
        sp->Node[j].treatment[p].equation = NULL;
    }
    return TRUE;
}

//=============================================================================

int  getVariableIndex(SWMM_Project *sp, char* s)
//
//  Input:   s = name of a process variable or pollutant
//  Output:  returns index of process variable or pollutant
//  Purpose: finds position of process variable/pollutant in list of names.
//
{
    // --- check for a process variable first
    int k;
    int m = PVMAX;                     // PVMAX is number of process variables

    k = findmatch(s, ProcessVarWords);
    if ( k >= 0 ) return k;

    // --- then check for a pollutant concentration
    k = project_findObject(POLLUT, s);
    if ( k >= 0 ) return (k + m);

    // --- finally check for a pollutant removal
    if ( UCHAR(s[0]) == 'R' && s[1] == '_')
    {
        k = project_findObject(POLLUT, s+2);
        if ( k >= 0 ) return (sp->Nobjects[POLLUT] + k + m);
    }
    return -1;
}

//=============================================================================

double getVariableValue(SWMM_Project *sp, int varCode)
//
//  Input:   varCode = code number of process variable or pollutant
//  Output:  returns current value of variable
//  Purpose: finds current value of a process variable or pollutant concen.,
//           making reference to the node being evaluated which is stored in
//           shared variable J.
//
{
    int    p;
    double a1, a2, y;
    TTreatment* treatment;                                                     //(5.1.008)

    // --- variable is a process variable
    if ( varCode < PVMAX )
    {
        switch ( varCode )
        {
          case pvHRT:                                 // HRT in hours
            if ( sp->Node[J].type == STORAGE )
            {
                return Storage[sp->Node[J].subIndex].hrt / 3600.0;
            }
            else return 0.0;

          case pvDT:
            return Dt;                                // time step in seconds

          case pvFLOW:
            return Q * UCF(sp, FLOW);                     // flow in user's units

          case pvDEPTH:
            y = (sp->Node[J].oldDepth + sp->Node[J].newDepth) / 2.0;
            return y * UCF(sp, LENGTH);                   // depth in ft or m

          case pvAREA:
            a1 = node_getSurfArea(sp, J, sp->Node[J].oldDepth);
            a2 = node_getSurfArea(sp, J, sp->Node[J].newDepth);
            return (a1 + a2) / 2.0 * UCF(sp, LENGTH) * UCF(sp, LENGTH);
            
          default: return 0.0;
        }
    }

    // --- variable is a pollutant concentration
    else if ( varCode < PVMAX + sp->Nobjects[POLLUT] )
    {
        p = varCode - PVMAX;
        treatment = &sp->Node[J].treatment[p];                                     //(5.1.008)
        if ( treatment->treatType == REMOVAL ) return Cin[p];                  //(5.1.008)
        return sp->Node[J].newQual[p];
    }

    // --- variable is a pollutant removal
    else
    {
        p = varCode - PVMAX - sp->Nobjects[POLLUT];
        if ( p >= sp->Nobjects[POLLUT] ) return 0.0;
        return getRemoval(sp, p);
    }
}

//=============================================================================

double  getRemoval(SWMM_Project *sp, int p)
//
//  Input:   p = pollutant index
//  Output:  returns fractional removal of pollutant
//  Purpose: computes removal of a specific pollutant
//
{
    double c0 = sp->Node[J].newQual[p];    // initial node concentration
    double r;                          // removal value
    TTreatment* treatment;                                                     //(5.1.008)

    // --- case where removal already being computed for another pollutant
    if ( R[p] > 1.0 || ErrCode )
    {
        ErrCode = 1;
        return 0.0;
    }

    // --- case where removal already computed
    if ( R[p] >= 0.0 && R[p] <= 1.0 ) return R[p];

    // --- set R[p] to value > 1 to show that value is being sought
    //     (prevents infinite recursive calls in case two removals
    //     depend on each other)
    R[p] = 10.0;

    // --- case where current concen. is zero
    if ( c0 == 0.0 )
    {
        R[p] = 0.0;
        return 0.0;
    }

    // --- apply treatment eqn.
    treatment = &sp->Node[J].treatment[p];                                         //(5.1.008)
    r = mathexpr_eval(sp, treatment->equation, getVariableValue);              //(5.1.008)
    r = MAX(0.0, r);

    // --- case where treatment eqn. is for removal
    if ( treatment->treatType == REMOVAL )                                     //(5.1.008)
    {
        r = MIN(1.0, r);
        R[p] = r;
    }

    // --- case where treatment eqn. is for effluent concen.
    else
    {
        r = MIN(c0, r);
        R[p] = 1.0 - r/c0;
    }
    return R[p];
}

//=============================================================================
