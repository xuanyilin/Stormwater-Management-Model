//-----------------------------------------------------------------------------
//   controls.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/21/14 (Build 5.1.001)
//             03/19/15 (Build 5.1.008)
//             04/30/15 (Build 5.1.009)
//             08/05/15 (Build 5.1.010)
//             08/01/16 (Build 5.1.011)
//   Author:   L. Rossman
//
//   Rule-based controls functions.
//
//   Control rules have the format:
//     RULE name
//     IF <premise>
//     AND / OR <premise>
//     etc.
//     THEN <action>
//     AND  <action>
//     etc.
//     ELSE <action>
//     AND  <action>
//     etc.
//     PRIORITY <p>
//
//   <premise> consists of:
//      <variable> <relational operator> value / <variable>
//   where <variable> is <object type> <id name> <attribute>
//   E.g.: Node 123 Depth > 4.5
//         Node 456 Depth < Node 123 Depth
//
//   <action> consists of:
//      <variable> = setting
//   E.g.: Pump abc status = OFF
//         Weir xyz setting = 0.5
//
//  Build 5.1.008:
//  - Support added for r.h.s. variables in rule premises.
//  - Node volume added as a premise variable.
//
//  Build 5.1.009:
//  - Fixed problem with parsing a RHS premise variable.
//
//  Build 5.1.010:
//  - Support added for link TIMEOPEN & TIMECLOSED premises.
//
//  Build 5.1.011:
//  - Support added for DAYOFYEAR attribute.
//  - Modulated controls no longer included in reported control actions.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"
#include "controls.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
enum RuleState    {r_RULE, r_IF, r_AND, r_OR, r_THEN, r_ELSE, r_PRIORITY,
                   r_ERROR};
enum RuleObject   {r_NODE, r_LINK, r_CONDUIT, r_PUMP, r_ORIFICE, r_WEIR,
	               r_OUTLET, r_SIMULATION};
enum RuleAttrib   {r_DEPTH, r_HEAD, r_VOLUME, r_INFLOW, r_FLOW, r_STATUS,      //(5.1.008)
                   r_SETTING, r_TIMEOPEN, r_TIMECLOSED, r_TIME, r_DATE,        //(5.1.010)
                   r_CLOCKTIME, r_DAYOFYEAR, r_DAY, r_MONTH};                  //(5.1.011)
enum RuleRelation {EQ, NE, LT, LE, GT, GE};
enum RuleSetting  {r_CURVE, r_TIMESERIES, r_PID, r_NUMERIC};

static char* ObjectWords[] =
    {"NODE", "LINK", "CONDUIT", "PUMP", "ORIFICE", "WEIR", "OUTLET",
	 "SIMULATION", NULL};
static char* AttribWords[] =
    {"DEPTH", "HEAD", "VOLUME", "INFLOW", "FLOW", "STATUS", "SETTING",         //(5.1.008)
     "TIMEOPEN", "TIMECLOSED","TIME", "DATE", "CLOCKTIME", "DAYOFYEAR",        //(5.1.011)
     "DAY", "MONTH", NULL};                                                    //(5.1.011)
static char* RelOpWords[] = {"=", "<>", "<", "<=", ">", ">=", NULL};
static char* StatusWords[]  = {"OFF", "ON", NULL};
static char* ConduitWords[] = {"CLOSED", "OPEN", NULL};
static char* SettingTypeWords[] = {"CURVE", "TIMESERIES", "PID", NULL};


//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//     controls_create
//     controls_delete
//     controls_addRuleClause
//     controls_evaluate

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
int    addPremise(SWMM_Project *sp, int r, int type, char* Tok[], int nToks);
int    getPremiseVariable(SWMM_Project *sp, char* tok[], int* k, struct TVariable* v);
int    getPremiseValue(SWMM_Project *sp, char* token, int attrib, double* value);
int    addAction(SWMM_Project *sp, int r, char* Tok[], int nToks);

int    evaluatePremise(SWMM_Project *sp, struct TPremise* p, double tStep);
double getVariableValue(SWMM_Project *sp, struct TVariable v);
int    compareTimes(SWMM_Project *sp, double lhsValue, int relation, double rhsValue,
       double halfStep);
int    compareValues(SWMM_Project *sp, double lhsValue, int relation, double rhsValue);

void   updateActionList(SWMM_Project *sp, struct TAction* a);
int    executeActionList(SWMM_Project *sp, DateTime currentTime);
void   clearActionList(SWMM_Project *sp);
void   deleteActionList(SWMM_Project *sp);
void   deleteRules(SWMM_Project *sp);

int    findExactMatch(char *s, char *keyword[]);
int    setActionSetting(SWMM_Project *sp, char* tok[], int nToks, int* curve, int* tseries,
       int* attrib, double* value);
void   updateActionValue(SWMM_Project *sp, struct TAction* a, DateTime currentTime,
        double dt);
double getPIDSetting(SWMM_Project *sp, struct TAction* a, double dt);

//=============================================================================

int  controls_create(SWMM_Project *sp, int n)
//
//  Input:   n = total number of control rules
//  Output:  returns error code
//  Purpose: creates an array of control rules.
//
{
   int r;

   TControlsShared *cntrl = &sp->ControlsShared;

   cntrl->ActionList = NULL;
   cntrl->InputState = r_PRIORITY;
   cntrl->RuleCount = n;
   if ( n == 0 ) return 0;
   cntrl->Rules = (struct TRule *) calloc(cntrl->RuleCount, sizeof(struct TRule));
   if (cntrl->Rules == NULL) return ERR_MEMORY;
   for ( r=0; r < cntrl->RuleCount; r++ )
   {
       cntrl->Rules[r].ID = NULL;
       cntrl->Rules[r].firstPremise = NULL;
       cntrl->Rules[r].lastPremise = NULL;
       cntrl->Rules[r].thenActions = NULL;
       cntrl->Rules[r].elseActions = NULL;
       cntrl->Rules[r].priority = 0.0;
   }
   return 0;
}

//=============================================================================

void controls_delete(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: deletes all control rules.
//
{
    TControlsShared *cntrl = &sp->ControlsShared;

   if ( cntrl->RuleCount == 0 ) return;
   deleteActionList(sp);
   deleteRules(sp);
}

//=============================================================================

int  controls_addRuleClause(SWMM_Project *sp, int r, int keyword, char* tok[],
        int nToks)
//
//  Input:   r = rule index
//           keyword = the clause's keyword code (IF, THEN, etc.)
//           tok = an array of string tokens that comprises the clause
//           nToks = number of tokens
//  Output:  returns an error  code
//  Purpose: addd a new clause to a control rule.
//
{
    TControlsShared *cntrl = &sp->ControlsShared;

    switch (keyword)
    {
      case r_RULE:
        if ( cntrl->Rules[r].ID == NULL )
            cntrl->Rules[r].ID = project_findID(sp, CONTROL, tok[1]);
        cntrl->InputState = r_RULE;
        if ( nToks > 2 ) return ERR_RULE;
        return 0;

      case r_IF:
        if ( cntrl->InputState != r_RULE ) return ERR_RULE;
        cntrl->InputState = r_IF;
        return addPremise(sp, r, r_AND, tok, nToks);

      case r_AND:
        if ( cntrl->InputState == r_IF ) return addPremise(sp, r, r_AND, tok, nToks);
        else if ( cntrl->InputState == r_THEN || cntrl->InputState == r_ELSE )
            return addAction(sp, r, tok, nToks);
        else return ERR_RULE;

      case r_OR:
        if ( cntrl->InputState != r_IF ) return ERR_RULE;
        return addPremise(sp, r, r_OR, tok, nToks);

      case r_THEN:
        if ( cntrl->InputState != r_IF ) return ERR_RULE;
        cntrl->InputState = r_THEN;
        return addAction(sp, r, tok, nToks);

      case r_ELSE:
        if ( cntrl->InputState != r_THEN ) return ERR_RULE;
        cntrl->InputState = r_ELSE;
        return addAction(sp, r, tok, nToks);

      case r_PRIORITY:
        if ( cntrl->InputState != r_THEN && cntrl->InputState != r_ELSE ) return ERR_RULE;
        cntrl->InputState = r_PRIORITY;
        if ( !getDouble(tok[1], &cntrl->Rules[r].priority) ) return ERR_NUMBER;
        if ( nToks > 2 ) return ERR_RULE;
        return 0;
    }
    return 0;
}

//=============================================================================

int controls_evaluate(SWMM_Project *sp, DateTime currentTime,
        DateTime elapsedTime, double tStep)
//
//  Input:   currentTime = current simulation date/time
//           elapsedTime = decimal days since start of simulation
//           tStep = simulation time step (days)
//  Output:  returns number of new actions taken
//  Purpose: evaluates all control rules at current time of the simulation.
//
{
    int    r;                          // control rule index
    int    result;                     // TRUE if rule premises satisfied
    struct TPremise* p;                // pointer to rule premise clause
    struct TAction*  a;                // pointer to rule action clause

    TControlsShared *cntrl = &sp->ControlsShared;

    // --- save date and time to shared variables
    cntrl->CurrentDate = floor(currentTime);
    cntrl->CurrentTime = currentTime - floor(currentTime);
    sp->ElapsedTime = elapsedTime;

    // --- evaluate each rule
    if ( cntrl->RuleCount == 0 ) return 0;
    clearActionList(sp);
    for (r = 0; r < cntrl->RuleCount; r++)
    {
        // --- evaluate rule's premises
        result = TRUE;
        p = cntrl->Rules[r].firstPremise;
        while (p)
        {
            if ( p->type == r_OR )
            {
                if ( result == FALSE )
                    result = evaluatePremise(sp, p, tStep);
            }
            else
            {
                if ( result == FALSE ) break;
                result = evaluatePremise(sp, p, tStep);
            }
            p = p->next;
        }    

        // --- if premises true, add THEN clauses to action list
        //     else add ELSE clauses to action list
        if ( result == TRUE ) a = cntrl->Rules[r].thenActions;
        else                  a = cntrl->Rules[r].elseActions;
        while (a)
        {
            updateActionValue(sp, a, currentTime, tStep);
            updateActionList(sp, a);
            a = a->next;
        }
    }

    // --- execute actions on action list
    if ( cntrl->ActionList ) return executeActionList(sp, currentTime);
    else return 0;
}

//=============================================================================

//  This function was revised to add support for r.h.s. premise variables. //  //(5.1.008)

int  addPremise(SWMM_Project *sp, int r, int type, char* tok[], int nToks)
//
//  Input:   r = control rule index
//           type = type of premise (IF, AND, OR)
//           tok = array of string tokens containing premise statement
//           nToks = number of string tokens
//  Output:  returns an error code
//  Purpose: adds a new premise to a control rule.
//
{
    int    relation, n, err = 0;
    double value = MISSING;
    struct TPremise* p;
    struct TVariable v1;
    struct TVariable v2;

    TControlsShared *cntrl = &sp->ControlsShared;

    // --- check for minimum number of tokens
    if ( nToks < 5 ) return ERR_ITEMS;

    // --- get LHS variable
    n = 1;
    err = getPremiseVariable(sp, tok, &n, &v1);
    if ( err > 0 ) return err;

    // --- get relational operator
    n++;
    relation = findExactMatch(tok[n], RelOpWords);
    if ( relation < 0 ) return error_setInpError(ERR_KEYWORD, tok[n]);
    n++;

    // --- initialize RHS variable
    v2.attribute = -1;
    v2.link = -1;
    v2.node = -1;

    // --- check that more tokens remain
    if ( n >= nToks ) return error_setInpError(ERR_ITEMS, "");
        
    // --- see if a RHS variable is supplied
    if ( findmatch(tok[n], ObjectWords) >= 0 && n + 3 >= nToks )
    {
        err = getPremiseVariable(sp, tok, &n, &v2);
        if ( err > 0 ) return ERR_RULE;                                        //(5.1.009)
        if ( v1.attribute != v2.attribute)                                     //(5.1.009)
            report_writeWarningMsg(sp, WARN11, cntrl->Rules[r].ID);            //(5.1.009)
    }

    // --- otherwise get value to which LHS variable is compared to
    else
    {
        err = getPremiseValue(sp, tok[n], v1.attribute, &value);
        n++;
    }
    if ( err > 0 ) return err;

    // --- make sure another clause is not on same line
    if ( n < nToks && findmatch(tok[n], RuleKeyWords) >= 0 ) return ERR_RULE;

    // --- create the premise object
    p = (struct TPremise *) malloc(sizeof(struct TPremise));
    if ( !p ) return ERR_MEMORY;
    p->type      = type;
    p->lhsVar    = v1;
    p->rhsVar    = v2;
    p->relation  = relation;
    p->value     = value;
    p->next      = NULL;
    if ( cntrl->Rules[r].firstPremise == NULL )
    {
        cntrl->Rules[r].firstPremise = p;
    }
    else
    {
        cntrl->Rules[r].lastPremise->next = p;
    }
    cntrl->Rules[r].lastPremise = p;
    return 0;
}

//=============================================================================

int getPremiseVariable(SWMM_Project *sp, char* tok[], int* k, struct TVariable* v)
//
//  Input:   tok = array of string tokens containing premise statement
//           k = index of current token
//  Output:  returns an error code; updates k to new current token and
//           places identity of specified variable in v
//  Purpose: parses a variable (e.g., Node 123 Depth) specified in a
//           premise clause of a control rule.
//
{
    int    n = *k;
    int    node = -1;
    int    link = -1;
    int    obj, attrib;

    // --- get object type
    obj = findmatch(tok[n], ObjectWords);
    if ( obj < 0 ) return error_setInpError(ERR_KEYWORD, tok[n]);

    // --- get object index from its name
    n++;
    switch (obj)
    {
      case r_NODE:
        node = project_findObject(sp, NODE, tok[n]);
        if ( node < 0 ) return error_setInpError(ERR_NAME, tok[n]);
        break;

      case r_LINK:
      case r_CONDUIT:
      case r_PUMP:
      case r_ORIFICE:
      case r_WEIR:
      case r_OUTLET:
        link = project_findObject(sp, LINK, tok[n]);
        if ( link < 0 ) return error_setInpError(ERR_NAME, tok[n]);
        break;
      default: n--;
    }
    n++;

    // --- get attribute index from its name
    attrib = findmatch(tok[n], AttribWords);
    if ( attrib < 0 ) return error_setInpError(ERR_KEYWORD, tok[n]);

    // --- check that attribute belongs to object type
    if ( obj == r_NODE ) switch (attrib)
    {
      case r_DEPTH:
      case r_HEAD:
      case r_VOLUME:                                                           //(5.1.008)
      case r_INFLOW: break;
      default: return error_setInpError(ERR_KEYWORD, tok[n]);
    }

////  Added to release 5.1.010.  ////                                          //(5.1.010)
    // --- check for link TIMEOPEN & TIMECLOSED attributes
    else if ( link >= 0  &&
            ( (attrib == r_TIMEOPEN ||
               attrib == r_TIMECLOSED)
            ))
    {
 
    }
////

    else if ( obj == r_LINK || obj == r_CONDUIT ) switch (attrib)
    {
      case r_STATUS:
      case r_DEPTH:
      case r_FLOW: break;
      default: return error_setInpError(ERR_KEYWORD, tok[n]);
    }
    else if ( obj == r_PUMP ) switch (attrib)
    {
      case r_FLOW:
      case r_STATUS: break;
      default: return error_setInpError(ERR_KEYWORD, tok[n]);
    }
    else if ( obj == r_ORIFICE || obj == r_WEIR ||
              obj == r_OUTLET ) switch (attrib)
    {
      case r_SETTING: break;
      default: return error_setInpError(ERR_KEYWORD, tok[n]);
    }
    else switch (attrib)
    {
      case r_TIME:
      case r_DATE:
      case r_CLOCKTIME:
      case r_DAY:
      case r_MONTH:
      case r_DAYOFYEAR: break;                                                 //(5.1.011)
      default: return error_setInpError(ERR_KEYWORD, tok[n]);
    }

    // --- populate variable structure
    v->node      = node;
    v->link      = link;
    v->attribute = attrib;
    *k = n;
    return 0;
}

//=============================================================================

int getPremiseValue(SWMM_Project *sp, char* token, int attrib, double* value)
//
//  Input:   token = a string token
//           attrib = index of a node/link attribute
//  Output:  value = attribute value;
//           returns an error code;
//  Purpose: parses the numerical value of a particular node/link attribute
//           in the premise clause of a control rule.
//
{
    char   strDate[25];                                                        //(5.1.011)
    switch (attrib)
    {
      case r_STATUS:
        *value = findmatch(token, StatusWords);
		if ( *value < 0.0 ) *value = findmatch(token, ConduitWords);
        if ( *value < 0.0 ) return error_setInpError(ERR_KEYWORD, token);
        break;

      case r_TIME:
      case r_CLOCKTIME:
      case r_TIMEOPEN:                                                         //(5.1.010)
      case r_TIMECLOSED:                                                       //(5.1.010)
        if ( !datetime_strToTime(token, value) )
            return error_setInpError(ERR_DATETIME, token);
        break;

      case r_DATE:
        if ( !datetime_strToDate(sp, token, value) )
            return error_setInpError(ERR_DATETIME, token);
        break;

      case r_DAY:
        if ( !getDouble(token, value) ) 
            return error_setInpError(ERR_NUMBER, token);
        if ( *value < 1.0 || *value > 7.0 )
             return error_setInpError(ERR_DATETIME, token);
        break;

      case r_MONTH:
        if ( !getDouble(token, value) )
            return error_setInpError(ERR_NUMBER, token);
        if ( *value < 1.0 || *value > 12.0 )
             return error_setInpError(ERR_DATETIME, token);
        break;

////  This code block added to release 5.1.011.  ////                          //(5.1.011)
      case r_DAYOFYEAR:
        strncpy(strDate, token, 6);
        strcat(strDate, "/1947");
        if ( datetime_strToDate(sp, strDate, value) )
        {
            *value = datetime_dayOfYear(*value);
        }
        else if ( !getDouble(token, value) || *value < 1 || *value > 365 )
            return error_setInpError(ERR_DATETIME, token);
        break;
////////////////////////////////////////////////////
       
      default: if ( !getDouble(token, value) )
          return error_setInpError(ERR_NUMBER, token);
    }
    return 0;
}

//=============================================================================

int  addAction(SWMM_Project *sp, int r, char* tok[], int nToks)
//
//  Input:   r = control rule index
//           tok = array of string tokens containing action statement
//           nToks = number of string tokens
//  Output:  returns an error code
//  Purpose: adds a new action to a control rule.
//
{
    int    obj, link, attrib;
    int    curve = -1, tseries = -1;
    int    n;
    int    err;

    TControlsShared *cntrl = &sp->ControlsShared;

    double values[] = {1.0, 0.0, 0.0};

    struct TAction* a;

    // --- check for proper number of tokens
    if ( nToks < 6 ) return error_setInpError(ERR_ITEMS, "");

    // --- check for valid object type
    obj = findmatch(tok[1], ObjectWords);
    if ( obj != r_LINK && obj != r_CONDUIT && obj != r_PUMP && 
         obj != r_ORIFICE && obj != r_WEIR && obj != r_OUTLET )
        return error_setInpError(ERR_KEYWORD, tok[1]);

    // --- check that object name exists and is of correct type
    link = project_findObject(sp, LINK, tok[2]);
    if ( link < 0 ) return error_setInpError(ERR_NAME, tok[2]);
    switch (obj)
    {
      case r_CONDUIT:
	if ( sp->Link[link].type != CONDUIT )
	    return error_setInpError(ERR_NAME, tok[2]);
	break;
      case r_PUMP:
        if ( sp->Link[link].type != PUMP )
            return error_setInpError(ERR_NAME, tok[2]);
        break;
      case r_ORIFICE:
        if ( sp->Link[link].type != ORIFICE )
            return error_setInpError(ERR_NAME, tok[2]);
        break;
      case r_WEIR:
        if ( sp->Link[link].type != WEIR )
            return error_setInpError(ERR_NAME, tok[2]);
        break;
      case r_OUTLET:
        if ( sp->Link[link].type != OUTLET )
            return error_setInpError(ERR_NAME, tok[2]);
        break;
    }

    // --- check for valid attribute name
    attrib = findmatch(tok[3], AttribWords);
    if ( attrib < 0 ) return error_setInpError(ERR_KEYWORD, tok[3]);

    // --- get control action setting
    if ( obj == r_CONDUIT )
    {
        if ( attrib == r_STATUS )
        {
            values[0] = findmatch(tok[5], ConduitWords);
            if ( values[0] < 0.0 )
                return error_setInpError(ERR_KEYWORD, tok[5]);
        }
        else return error_setInpError(ERR_KEYWORD, tok[3]);
    }

    else if ( obj == r_PUMP )
    {
        if ( attrib == r_STATUS )
        {
            values[0] = findmatch(tok[5], StatusWords);
            if ( values[0] < 0.0 )
                return error_setInpError(ERR_KEYWORD, tok[5]);
        }
        else if ( attrib == r_SETTING )
        {
            err = setActionSetting(sp, tok, nToks, &curve, &tseries,
                                   &attrib, values);
            if ( err > 0 ) return err;
        }
        else return error_setInpError(ERR_KEYWORD, tok[3]);
    }

    else if ( obj == r_ORIFICE || obj == r_WEIR || obj == r_OUTLET )
    {
        if ( attrib == r_SETTING )
        {
           err = setActionSetting(sp, tok, nToks, &curve, &tseries,
                                  &attrib, values);
           if ( err > 0 ) return err;
           if (  attrib == r_SETTING
           && (values[0] < 0.0 || values[0] > 1.0) ) 
               return error_setInpError(ERR_NUMBER, tok[5]);
        }
        else return error_setInpError(ERR_KEYWORD, tok[3]);
    }
    else return error_setInpError(ERR_KEYWORD, tok[1]);

    // --- check if another clause is on same line
    n = 6;
    if ( curve >= 0 || tseries >= 0 ) n = 7;
    if ( attrib == r_PID ) n = 9;
    if ( n < nToks && findmatch(tok[n], RuleKeyWords) >= 0 ) return ERR_RULE;

    // --- create the action object
    a = (struct TAction *) malloc(sizeof(struct TAction));
    if ( !a ) return ERR_MEMORY;
    a->rule      = r;
    a->link      = link;
    a->attribute = attrib;
    a->curve     = curve;
    a->tseries   = tseries;
    a->value     = values[0];
    if ( attrib == r_PID )
    {
        a->kp = values[0];
        a->ki = values[1];
        a->kd = values[2];
        a->e1 = 0.0;
        a->e2 = 0.0;
    }
    if ( cntrl->InputState == r_THEN )
    {
        a->next = cntrl->Rules[r].thenActions;
        cntrl->Rules[r].thenActions = a;
    }
    else
    {
        a->next = cntrl->Rules[r].elseActions;
        cntrl->Rules[r].elseActions = a;
    }
    return 0;
}

//=============================================================================

int  setActionSetting(SWMM_Project *sp, char* tok[], int nToks, int* curve,
        int* tseries, int* attrib, double values[])
//
//  Input:   tok = array of string tokens containing action statement
//           nToks = number of string tokens
//  Output:  curve = index of controller curve
//           tseries = index of controller time series
//           attrib = r_PID if PID controller used
//           values = values of control settings
//           returns an error code
//  Purpose: identifies how control actions settings are determined.
//
{
    int k, m;

    // --- see if control action is determined by a Curve or Time Series
    if (nToks < 6) return error_setInpError(ERR_ITEMS, "");
    k = findmatch(tok[5], SettingTypeWords);
    if ( k >= 0 && nToks < 7 ) return error_setInpError(ERR_ITEMS, "");
    switch (k)
    {

    // --- control determined by a curve - find curve index
    case r_CURVE:
        m = project_findObject(sp, CURVE, tok[6]);
        if ( m < 0 ) return error_setInpError(ERR_NAME, tok[6]);
        *curve = m;
        break;

    // --- control determined by a time series - find time series index
    case r_TIMESERIES:
        m = project_findObject(sp, TSERIES, tok[6]);
        if ( m < 0 ) return error_setInpError(ERR_NAME, tok[6]);
        *tseries = m;
        sp->Tseries[m].refersTo = CONTROL;
        break;

    // --- control determined by PID controller 
    case r_PID:
        if (nToks < 9) return error_setInpError(ERR_ITEMS, "");
        for (m=6; m<=8; m++)
        {
            if ( !getDouble(tok[m], &values[m-6]) )
                return error_setInpError(ERR_NUMBER, tok[m]);
        }
        *attrib = r_PID;
        break;

    // --- direct numerical control is used
    default:
        if ( !getDouble(tok[5], &values[0]) )
            return error_setInpError(ERR_NUMBER, tok[5]);
    }
    return 0;
}

//=============================================================================

void  updateActionValue(SWMM_Project *sp, struct TAction* a, DateTime currentTime,
        double dt)
//
//  Input:   a = an action object
//           currentTime = current simulation date/time (days)
//           dt = time step (days)
//  Output:  none
//  Purpose: updates value of actions found from Curves or Time Series.
//
{
    TControlsShared *cntrl = &sp->ControlsShared;

    if ( a->curve >= 0 )
    {
        a->value = table_lookup(&sp->Curve[a->curve], cntrl->ControlValue);
    }
    else if ( a->tseries >= 0 )
    {
        a->value = table_tseriesLookup(sp, &sp->Tseries[a->tseries], currentTime, TRUE);
    }
    else if ( a->attribute == r_PID )
    {
        a->value = getPIDSetting(sp, a, dt);
    }
}

//=============================================================================

double getPIDSetting(SWMM_Project *sp, struct TAction* a, double dt)
//
//  Input:   a = an action object
//           dt = current time step (days)
//  Output:  returns a new link setting 
//  Purpose: computes a new setting for a link subject to a PID controller.
//
//  Note:    a->kp = gain coefficient,
//           a->ki = integral time (minutes)
//           a->k2 = derivative time (minutes)
//           a->e1 = error from previous time step
//           a->e2 = error from two time steps ago
{
    double e0, setting;
	double p, i, d, update;
	double tolerance = 0.0001;

	TControlsShared *cntrl = &sp->ControlsShared;

	// --- convert time step from days to minutes
	dt *= 1440.0;

    // --- determine relative error in achieving controller set point
    e0 = cntrl->SetPoint - cntrl->ControlValue;
    if ( fabs(e0) > TINY )
    {
        if ( cntrl->SetPoint != 0.0 ) e0 = e0/cntrl->SetPoint;
        else                   e0 = e0/cntrl->ControlValue;
    }

	// --- reset previous errors to 0 if controller gets stuck
	if (fabs(e0 - a->e1) < tolerance)
	{
		a->e2 = 0.0;
		a->e1 = 0.0;
	}

    // --- use the recursive form of the PID controller equation to
    //     determine the new setting for the controlled link
	p = (e0 - a->e1);
	if ( a->ki == 0.0 ) i = 0.0;
	else i = e0 * dt / a->ki;
	d = a->kd * (e0 - 2.0*a->e1 + a->e2) / dt;
	update = a->kp * (p + i + d);
	if ( fabs(update) < tolerance ) update = 0.0;
	setting = sp->Link[a->link].targetSetting + update;

	// --- update previous errors
    a->e2 = a->e1;
    a->e1 = e0;

    // --- check that new setting lies within feasible limits
    if ( setting < 0.0 ) setting = 0.0;
    if (sp->Link[a->link].type != PUMP && setting > 1.0 ) setting = 1.0;
    return setting;
}

//=============================================================================

void updateActionList(SWMM_Project *sp, struct TAction* a)
//
//  Input:   a = an action object
//  Output:  none
//  Purpose: adds a new action to the list of actions to be taken.
//
{
    struct TActionList* listItem;
    struct TAction* a1;

    TControlsShared *cntrl = &sp->ControlsShared;

    double priority = cntrl->Rules[a->rule].priority;

    // --- check if link referred to in action is already listed
    listItem = cntrl->ActionList;
    while ( listItem )
    {
        a1 = listItem->action;
        if ( !a1 ) break;
        if ( a1->link == a->link )
        {
            // --- replace old action if new action has higher priority
            if ( priority > cntrl->Rules[a1->rule].priority ) listItem->action = a;
            return;
        }
        listItem = listItem->next;
    }

    // --- action not listed so add it to ActionList
    if ( !listItem )
    {
        listItem = (struct TActionList *) malloc(sizeof(struct TActionList));
        listItem->next = cntrl->ActionList;
        cntrl->ActionList = listItem;
    }
    listItem->action = a;
}

//=============================================================================

int executeActionList(SWMM_Project *sp, DateTime currentTime)
//
//  Input:   currentTime = current date/time of the simulation
//  Output:  returns number of new actions taken
//  Purpose: executes all actions required by fired control rules.
//
{
    struct TActionList* listItem;
    struct TActionList* nextItem;
    struct TAction* a1;
    int count = 0;

    TControlsShared *cntrl = &sp->ControlsShared;

    listItem = cntrl->ActionList;
    while ( listItem )
    {
        a1 = listItem->action;
        if ( !a1 ) break;
        if ( a1->link >= 0 )
        {
            if ( sp->Link[a1->link].targetSetting != a1->value )
            {
                sp->Link[a1->link].targetSetting = a1->value;
                if ( sp->RptFlags.controls && a1->curve < 0                        //(5.1.011)
                     && a1->tseries < 0 && a1->attribute != r_PID )            //(5.1.011)
                    report_writeControlAction(sp, currentTime, sp->Link[a1->link].ID,
                                              a1->value, cntrl->Rules[a1->rule].ID);
                count++;
            }
        }
        nextItem = listItem->next;
        listItem = nextItem;
    }
    return count;
}

//=============================================================================

////  This function was re-written for release 5.1.011.  ////                  //(5.1.011)

int evaluatePremise(SWMM_Project *sp, struct TPremise* p, double tStep)
//
//  Input:   p = a control rule premise condition
//           tStep = current time step (days)
//  Output:  returns TRUE if the condition is true or FALSE otherwise
//  Purpose: evaluates the truth of a control rule premise condition.
//
{
    double lhsValue, rhsValue;
    int    result = FALSE;

    TControlsShared *cntrl = &sp->ControlsShared;

    lhsValue = getVariableValue(sp, p->lhsVar);
    if ( p->value == MISSING ) rhsValue = getVariableValue(sp, p->rhsVar);
    else                       rhsValue = p->value;
    if ( lhsValue == MISSING || rhsValue == MISSING ) return FALSE;
    switch (p->lhsVar.attribute)
    {
    case r_TIME:
    case r_CLOCKTIME:
        return compareTimes(sp, lhsValue, p->relation, rhsValue, tStep/2.0);
    case r_TIMEOPEN:
    case r_TIMECLOSED:
        result = compareTimes(sp, lhsValue, p->relation, rhsValue, tStep/2.0);
        cntrl->ControlValue = lhsValue * 24.0;  // convert time from days to hours
        return result;
    default:
        return compareValues(sp, lhsValue, p->relation, rhsValue);
    }
}

//=============================================================================

double getVariableValue(SWMM_Project *sp, struct TVariable v)
{
    int i = v.node;
    int j = v.link;

    TControlsShared *cntrl = &sp->ControlsShared;

    switch ( v.attribute )
    {
      case r_TIME:
        return sp->ElapsedTime;
        
      case r_DATE:
        return cntrl->CurrentDate;

      case r_CLOCKTIME:
        return cntrl->CurrentTime;

      case r_DAY:
        return datetime_dayOfWeek(cntrl->CurrentDate);

      case r_MONTH:
        return datetime_monthOfYear(cntrl->CurrentDate);

      case r_DAYOFYEAR:                                                        //(5.1.011)
        return datetime_dayOfYear(cntrl->CurrentDate);                                //(5.1.011)

      case r_STATUS:
        if ( j < 0 ||
            (sp->Link[j].type != CONDUIT && sp->Link[j].type != PUMP) ) return MISSING;
        else return sp->Link[j].setting;
        
      case r_SETTING:
        if ( j < 0 || (sp->Link[j].type != ORIFICE && sp->Link[j].type != WEIR) )
            return MISSING;
        else return sp->Link[j].setting;

      case r_FLOW:
        if ( j < 0 ) return MISSING;
        else return sp->Link[j].direction*sp->Link[j].newFlow*UCF(sp, FLOW);

      case r_DEPTH:
        if ( j >= 0 ) return sp->Link[j].newDepth*UCF(sp, LENGTH);
        else if ( i >= 0 )
            return sp->Node[i].newDepth*UCF(sp, LENGTH);
        else return MISSING;

      case r_HEAD:
        if ( i < 0 ) return MISSING;
        return (sp->Node[i].newDepth + sp->Node[i].invertElev) * UCF(sp, LENGTH);

      case r_VOLUME:                                                           //(5.1.008)
        if ( i < 0 ) return MISSING;
        return (sp->Node[i].newVolume * UCF(sp, VOLUME));

      case r_INFLOW:
        if ( i < 0 ) return MISSING;
        else return sp->Node[i].newLatFlow*UCF(sp, FLOW);

////  This section added to release 5.1.010.  ////                             //(5.1.010)
      case r_TIMEOPEN:
          if ( j < 0 ) return MISSING;
          if ( sp->Link[j].setting <= 0.0 ) return MISSING;
          return cntrl->CurrentDate + cntrl->CurrentTime - sp->Link[j].timeLastSet;

      case r_TIMECLOSED:
          if ( j < 0 ) return MISSING;
          if ( sp->Link[j].setting > 0.0 ) return MISSING;
          return cntrl->CurrentDate + cntrl->CurrentTime - sp->Link[j].timeLastSet;
////

      default: return MISSING;
    }
}

//=============================================================================

int compareTimes(SWMM_Project *sp, double lhsValue, int relation, double rhsValue, double halfStep)
//
//  Input:   lhsValue = date/time value on left hand side of relation
//           relation = relational operator code (see RuleRelation enumeration)
//           rhsValue = date/time value on right hand side of relation 
//           halfStep = 1/2 the current time step (days)
//  Output:  returns TRUE if time relation is satisfied
//  Purpose: evaluates the truth of a relation between two date/times.
//
{
    if ( relation == EQ )
    {
        if ( lhsValue >= rhsValue - halfStep
        &&   lhsValue < rhsValue + halfStep ) return TRUE;
        return FALSE;
    }
    else if ( relation == NE )
    {
        if ( lhsValue < rhsValue - halfStep
        ||   lhsValue >= rhsValue + halfStep ) return TRUE;
        return FALSE;
    }
    else return compareValues(sp, lhsValue, relation, rhsValue);
}

//=============================================================================

int compareValues(SWMM_Project *sp, double lhsValue, int relation, double rhsValue)
//  Input:   lhsValue = value on left hand side of relation
//           relation = relational operator code (see RuleRelation enumeration)
//           rhsValue = value on right hand side of relation 
//  Output:  returns TRUE if relation is satisfied
//  Purpose: evaluates the truth of a relation between two values.
{
    TControlsShared *cntrl = &sp->ControlsShared;

    cntrl->SetPoint = rhsValue;
    cntrl->ControlValue = lhsValue;
    switch (relation)
    {
      case EQ: if ( lhsValue == rhsValue ) return TRUE; break;
      case NE: if ( lhsValue != rhsValue ) return TRUE; break;
      case LT: if ( lhsValue <  rhsValue ) return TRUE; break;
      case LE: if ( lhsValue <= rhsValue ) return TRUE; break;
      case GT: if ( lhsValue >  rhsValue ) return TRUE; break;
      case GE: if ( lhsValue >= rhsValue ) return TRUE; break;
    }
    return FALSE;
}

//=============================================================================

void clearActionList(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: clears the list of actions to be executed.
//
{
    TControlsShared *cntrl = &sp->ControlsShared;

    struct TActionList* listItem;
    listItem = cntrl->ActionList;
    while ( listItem )
    {
        listItem->action = NULL;
        listItem = listItem->next;
    }
}

//=============================================================================

void  deleteActionList(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: frees the memory used to hold the list of actions to be executed.
//
{
    TControlsShared *cntrl = &sp->ControlsShared;

    struct TActionList* listItem;
    struct TActionList* nextItem;
    listItem = cntrl->ActionList;
    while ( listItem )
    {
        nextItem = listItem->next;
        free(listItem);
        listItem = nextItem;
    }
    cntrl->ActionList = NULL;
}

//=============================================================================

void  deleteRules(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: frees the memory used for all of the control rules.
//
{
   struct TPremise* p;
   struct TPremise* pnext;
   struct TAction*  a;
   struct TAction*  anext;
   int r;

   TControlsShared *cntrl = &sp->ControlsShared;

   for (r = 0; r < cntrl->RuleCount; r++)
   {
      p = cntrl->Rules[r].firstPremise;
      while ( p )
      {
         pnext = p->next;
         free(p);
         p = pnext;
      }
      a = cntrl->Rules[r].thenActions;
      while (a )
      {
         anext = a->next;
         free(a);
         a = anext;
      }
      a = cntrl->Rules[r].elseActions;
      while (a )
      {
         anext = a->next;
         free(a);
         a = anext;
      }
   }
   FREE(cntrl->Rules);
   cntrl->RuleCount = 0;
}

//=============================================================================

int  findExactMatch(char *s, char *keyword[])
//
//  Input:   s = character string
//           keyword = array of keyword strings
//  Output:  returns index of keyword which matches s or -1 if no match found  
//  Purpose: finds exact match between string and array of keyword strings.
//
{
   int i = 0;
   while (keyword[i] != NULL)
   {
      if ( strcomp(s, keyword[i]) ) return(i);
      i++;
   }
   return(-1);
}

//=============================================================================
