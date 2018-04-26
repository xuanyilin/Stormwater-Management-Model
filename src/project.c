//-----------------------------------------------------------------------------
//   project.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/19/14  (Build 5.1.000)
//             04/14/14  (Build 5.1.004)
//             09/15/14  (Build 5.1.007)
//             03/19/15  (Build 5.1.008)
//             04/30/15  (Build 5.1.009)
//             08/01/16  (Build 5.1.011)
//             03/14/17  (Build 5.1.012)
//   Author:   L. Rossman
//
//   Project management functions.
//
//   This module provides project-related services such as:
//   o opening a new project and reading its input data
//   o allocating and freeing memory for project objects
//   o setting default values for object properties and options
//   o initializing the internal state of all objects
//   o managing hash tables for identifying objects by ID name
//
//   Build 5.1.004:
//   - Ignore RDII option added.
//
//   Build 5.1.007:
//   - Default monthly adjustments for climate variables included.
//   - User-supplied GW flow equaitions initialized to NULL.
//   - Storage node exfiltration object initialized to NULL.
//   - Freeing of memory used for storage node exfiltration included.
//
//   Build 5.1.008:
//   - Constants used for dynamic wave routing moved to dynwave.c.
//   - Input processing of minimum time step & number of
//     parallel threads for dynamic wave routing added.
//   - Default values of hyd. conductivity adjustments added.
//   - Freeing of memory used for outfall pollutant load added.
//
//   Build 5.1.009:
//   - Fixed bug in computing total duration introduced in 5.1.008.
//
//   Build 5.1.011:
//   - Memory management of hydraulic event dates array added.
//
//   Build 5.1.012:
//   - Minimum conduit slope option initialized to 0 (none).
//   - NO/YES no longer accepted as options for NORMAL_FLOW_LIMITED.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>                                                              //(5.1.008)
#if defined(_OPENMP)
  #include <omp.h>                                                             //(5.1.008)
#else
  int omp_get_num_threads(void) { return 1;}
#endif
#include "headers.h"
#include "lid.h" 
#include "hash.h"
#include "mempool.h"

#include "swmm5.h"
//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
////  Constants for DYNWAVE flow routing moved to dynwave.c.  ////             //(5.1.008)

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static HTtable* Htable[MAX_OBJ_TYPES]; // Hash tables for object ID names
static char     MemPoolAllocated;      // TRUE if memory pool allocated 

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  project_open           (called from swmm_open in swmm5.c)
//  project_close          (called from swmm_close in swmm5.c)
//  project_readInput      (called from swmm_open in swmm5.c)
//  project_readOption     (called from readOption in input.c)
//  project_validate       (called from swmm_open in swmm5.c)
//  project_init           (called from swmm_start in swmm5.c)
//  project_addObject      (called from addObject in input.c)
//  project_createMatrix   (called from openFileForInput in iface.c)
//  project_freeMatrix     (called from iface_closeRoutingFiles)
//  project_findObject
//  project_findID

//-----------------------------------------------------------------------------
//  Function declarations
//-----------------------------------------------------------------------------
static void initPointers(SWMM_Project *sp);
static void setDefaults(SWMM_Project *sp);
static void openFiles(SWMM_Project *sp, char *f1, char *f2, char *f3);
static void createObjects(SWMM_Project *sp);
static void deleteObjects(SWMM_Project *sp);
static void createHashTables(SWMM_Project *sp);
static void deleteHashTables(void);


//=============================================================================

void project_open(SWMM_Project *sp, char *f1, char *f2, char *f3)
//
//  Input:   f1 = pointer to name of input file
//           f2 = pointer to name of report file
//           f3 = pointer to name of binary output file
//  Output:  none
//  Purpose: opens a new SWMM project.
//
{
    initPointers(sp);
    setDefaults(sp);
    openFiles(sp, f1, f2, f3);
}

//=============================================================================

void project_readInput(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: retrieves project data from input file.
//
{
    // --- create hash tables for fast retrieval of objects by ID names
    createHashTables(sp);

    // --- count number of objects in input file and create them
    input_countObjects(sp);
    createObjects(sp);

    // --- read project data from input file
    input_readData(sp);
    if ( sp->ErrorCode ) return;

    // --- establish starting & ending date/time
    sp->StartDateTime = sp->StartDate + sp->StartTime;
    sp->EndDateTime   = sp->EndDate + sp->EndTime;
    sp->ReportStart   = sp->ReportStartDate + sp->ReportStartTime;
    sp->ReportStart   = MAX(sp->ReportStart, sp->StartDateTime);

    // --- check for valid starting & ending date/times
    if ( sp->EndDateTime <= sp->StartDateTime )
    {
        report_writeErrorMsg(sp, ERR_START_DATE, "");
    }
    else if ( sp->EndDateTime <= sp->ReportStart )
    {
        report_writeErrorMsg(sp, ERR_REPORT_DATE, "");
    }
    else
    {
////  Following code segment was modified for release 5.1.009.  ////           //(5.1.009)
////
        // --- compute total duration of simulation in seconds
        sp->TotalDuration = floor((sp->EndDateTime - sp->StartDateTime) * SECperDAY);

        // --- reporting step must be <= total duration
        if ( (double)sp->ReportStep > sp->TotalDuration )
        {
            sp->ReportStep = (int)(sp->TotalDuration);
        }

        // --- reporting step can't be < routing step
        if ( (double)sp->ReportStep < sp->RouteStep )
        {
            report_writeErrorMsg(sp, ERR_REPORT_STEP, "");
        }

        // --- convert total duration to milliseconds
        sp->TotalDuration *= 1000.0;
    }
////
}

//=============================================================================

void project_validate(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: checks validity of project data.
//
{
    int i;
    int j;
    int err;

    // --- validate Curves and TimeSeries
    for ( i=0; i<sp->Nobjects[CURVE]; i++ )
    {
         err = table_validate(&Curve[i]);
         if ( err ) report_writeErrorMsg(sp, ERR_CURVE_SEQUENCE, Curve[i].ID);
    }
    for ( i=0; i<sp->Nobjects[TSERIES]; i++ )
    {
        err = table_validate(&Tseries[i]);
        if ( err ) report_writeTseriesErrorMsg(sp, err, &Tseries[i]);
    }

    // --- validate hydrology objects
    //     (NOTE: order is important !!!!)
    climate_validate(sp);
    lid_validate(sp);
    if ( sp->Nobjects[SNOWMELT] == 0 ) sp->IgnoreSnowmelt = TRUE;
    if ( sp->Nobjects[AQUIFER]  == 0 ) sp->IgnoreGwater   = TRUE;
    for ( i=0; i<sp->Nobjects[GAGE]; i++ )     gage_validate(sp, i);
    for ( i=0; i<sp->Nobjects[AQUIFER]; i++ )  gwater_validateAquifer(sp, i);
    for ( i=0; i<sp->Nobjects[SUBCATCH]; i++ ) subcatch_validate(sp, i);
    for ( i=0; i<sp->Nobjects[SNOWMELT]; i++ ) snow_validateSnowmelt(sp, i);

    // --- compute geometry tables for each shape curve
    j = 0;
    for ( i=0; i<sp->Nobjects[CURVE]; i++ )
    {
        if ( Curve[i].curveType == SHAPE_CURVE )
        {
            Curve[i].refersTo = j;
            Shape[j].curve = i;
            if ( !shape_validate(&Shape[j], &Curve[i]) )
                report_writeErrorMsg(sp, ERR_CURVE_SEQUENCE, Curve[i].ID);
            j++;
        }
    }

    // --- validate links before nodes, since the latter can
    //     result in adjustment of node depths
    for ( i=0; i<sp->Nobjects[NODE]; i++) sp->Node[i].oldDepth = sp->Node[i].fullDepth;
    for ( i=0; i<sp->Nobjects[LINK]; i++) link_validate(sp, i);
    for ( i=0; i<sp->Nobjects[NODE]; i++) node_validate(sp, i);

    // --- adjust time steps if necessary
    if ( sp->DryStep < sp->WetStep )
    {
        report_writeWarningMsg(sp, WARN06, "");
        sp->DryStep = sp->WetStep;
    }
    if ( sp->RouteStep > (double)sp->WetStep )
    {
        report_writeWarningMsg(sp, WARN07, "");
        sp->RouteStep = sp->WetStep;
    }

    // --- adjust individual reporting flags to match global reporting flag
    if ( sp->RptFlags.subcatchments == ALL )
        for (i=0; i<sp->Nobjects[SUBCATCH]; i++) sp->Subcatch[i].rptFlag = TRUE;
    if ( sp->RptFlags.nodes == ALL )
        for (i=0; i<sp->Nobjects[NODE]; i++) sp->Node[i].rptFlag = TRUE;
    if ( sp->RptFlags.links == ALL )
        for (i=0; i<sp->Nobjects[LINK]; i++) sp->Link[i].rptFlag = TRUE;

    // --- validate dynamic wave options
    if ( sp->RouteModel == DW ) dynwave_validate(sp);                                //(5.1.008)

#pragma omp parallel                                                           //(5.1.008)
{
    if ( sp->NumThreads == 0 ) sp->NumThreads = omp_get_num_threads();                 //(5.1.008)
    else sp->NumThreads = MIN(sp->NumThreads, omp_get_num_threads());                  //(5.1.008)
}
    if ( sp->Nobjects[LINK] < 4 * sp->NumThreads ) sp->NumThreads = 1;                     //(5.1.008)

}

//=============================================================================

void project_close(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: closes a SWMM project.
//
{
    deleteObjects(sp);
    deleteHashTables();
}

//=============================================================================

int  project_init(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: initializes the internal state of all objects.
// 
{
    int j;
    climate_initState(sp);
    lid_initState(sp);
    for (j=0; j<sp->Nobjects[TSERIES]; j++)  table_tseriesInit(&Tseries[j]);
    for (j=0; j<sp->Nobjects[GAGE]; j++)     gage_initState(sp, j);
    for (j=0; j<sp->Nobjects[SUBCATCH]; j++) subcatch_initState(sp, j);
    for (j=0; j<sp->Nobjects[NODE]; j++)     node_initState(sp, j);
    for (j=0; j<sp->Nobjects[LINK]; j++)     link_initState(sp, j);
    return sp->ErrorCode;
}

//=============================================================================

int   project_addObject(int type, char *id, int n)
//
//  Input:   type = object type
//           id   = object ID string
//           n    = object index
//  Output:  returns 0 if object already added, 1 if not, -1 if hashing fails
//  Purpose: adds an object ID to a hash table
//
{
    int  result;
    int  len;
    char *newID;

    // --- do nothing if object already placed in hash table
    if ( project_findObject(type, id) >= 0 ) return 0;

    // --- use memory from the hash tables' common memory pool to store
    //     a copy of the object's ID string
    len = strlen(id) + 1;
    newID = (char *) Alloc(len*sizeof(char));
    strcpy(newID, id);

    // --- insert object's ID into the hash table for that type of object
    result = HTinsert(Htable[type], newID, n);
    if ( result == 0 ) result = -1;
    return result;
}

//=============================================================================

int DLLEXPORT  project_findObject(int type, char *id)
//
//  Input:   type = object type
//           id   = object ID
//  Output:  returns index of object with given ID, or -1 if ID not found
//  Purpose: uses hash table to find index of an object with a given ID.
//
{
    return HTfind(Htable[type], id);
}




//=============================================================================

char  *project_findID(int type, char *id)
//
//  Input:   type = object type
//           id   = ID name being sought
//  Output:  returns pointer to location where object's ID string is stored
//  Purpose: uses hash table to find address of given string entry.
//
{
    return HTfindKey(Htable[type], id);
}

//=============================================================================

double ** project_createMatrix(int nrows, int ncols)
//
//  Input:   nrows = number of rows (0-based)
//           ncols = number of columns (0-based)
//  Output:  returns a pointer to a matrix
//  Purpose: allocates memory for a matrix of doubles.
//
{
    int i,j;
    double **a;

    // --- allocate pointers to rows
    a = (double **) malloc(nrows * sizeof(double *));
    if ( !a ) return NULL;
    
    // --- allocate rows and set pointers to them
    a[0] = (double *) malloc (nrows * ncols * sizeof(double));
    if ( !a[0] ) return NULL;
    for ( i = 1; i < nrows; i++ ) a[i] = a[i-1] + ncols;

    for ( i = 0; i < nrows; i++)
    {
        for ( j = 0; j < ncols; j++) a[i][j] = 0.0;
    }
    
    // --- return pointer to array of pointers to rows
    return a;
}

//=============================================================================

void project_freeMatrix(double **a)
//
//  Input:   a = matrix of floats
//  Output:  none
//  Purpose: frees memory allocated for a matrix of doubles.
//
{
    if ( a != NULL )
    {
        if ( a[0] != NULL ) free( a[0] );
        free( a );
    }
}

//=============================================================================

int project_readOption(SWMM_Project *sp, char* s1, char* s2)
//
//  Input:   s1 = option keyword
//           s2 = string representation of option's value
//  Output:  returns error code
//  Purpose: reads a project option from a pair of string tokens.
//
//  NOTE:    all project options have default values assigned in setDefaults().
//
{
    int      k, m, h, s;
    double   tStep;
    char     strDate[25];
    DateTime aTime;
    DateTime aDate;

    // --- determine which option is being read
    k = findmatch(s1, OptionWords);
    if ( k < 0 ) return error_setInpError(ERR_KEYWORD, s1);
    switch ( k )
    {
      // --- choice of flow units
      case FLOW_UNITS:
        m = findmatch(s2, FlowUnitWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        sp->FlowUnits = m;
        if ( sp->FlowUnits <= MGD ) sp->UnitSystem = US;
        else                    sp->UnitSystem = SI;
        break;

      // --- choice of infiltration modeling method
      case INFIL_MODEL:
        m = findmatch(s2, InfilModelWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        sp->InfilModel = m;
        break;

      // --- choice of flow routing method
      case ROUTE_MODEL:
        m = findmatch(s2, RouteModelWords);
        if ( m < 0 ) m = findmatch(s2, OldRouteModelWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        if ( m == NO_ROUTING ) sp->IgnoreRouting = TRUE;
        else sp->RouteModel = m;
        if ( sp->RouteModel == EKW ) sp->RouteModel = KW;
        break;

      // --- simulation start date
      case START_DATE:
        if ( !datetime_strToDate(s2, &sp->StartDate) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        break;

      // --- simulation start time of day
      case START_TIME:
        if ( !datetime_strToTime(s2, &sp->StartTime) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        break;

      // --- simulation ending date
      case END_DATE:
        if ( !datetime_strToDate(s2, &sp->EndDate) ) 
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        break;

      // --- simulation ending time of day
      case END_TIME:
        if ( !datetime_strToTime(s2, &sp->EndTime) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        break;

      // --- reporting start date
      case REPORT_START_DATE:
        if ( !datetime_strToDate(s2, &sp->ReportStartDate) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        break;

      // --- reporting start time of day
      case REPORT_START_TIME:
        if ( !datetime_strToTime(s2, &sp->ReportStartTime) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        break;

      // --- day of year when street sweeping begins or when it ends
      //     (year is arbitrarily set to 1947 so that the dayOfYear
      //      function can be applied)
      case SWEEP_START:
      case SWEEP_END:
        strcpy(strDate, s2);
        strcat(strDate, "/1947");
        if ( !datetime_strToDate(strDate, &aDate) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        m = datetime_dayOfYear(aDate);
        if ( k == SWEEP_START ) sp->ReportStep = m;
        else sp->SweepEnd = m;
        break;

      // --- number of antecedent dry days
      case START_DRY_DAYS:
        sp->StartDryDays = atof(s2);
        if ( sp->StartDryDays < 0.0 )
        {
            return error_setInpError(ERR_NUMBER, s2);
        }
        break;

      // --- runoff or reporting time steps
      //     (input is in hrs:min:sec format, time step saved as seconds)
      case WET_STEP:
      case DRY_STEP:
      case REPORT_STEP:
        if ( !datetime_strToTime(s2, &aTime) )
        {
            return error_setInpError(ERR_DATETIME, s2);
        }
        datetime_decodeTime(aTime, &h, &m, &s);
        h += 24*(int)aTime;
        s = s + 60*m + 3600*h;
        if ( s <= 0 ) return error_setInpError(ERR_NUMBER, s2);
        switch ( k )
        {
          case WET_STEP:     sp->WetStep = s;     break;
          case DRY_STEP:     sp->DryStep = s;     break;
          case REPORT_STEP:  sp->ReportStep = s;  break;
        }
        break;

      // --- type of damping applied to inertial terms of dynamic wave routing
      case INERT_DAMPING:
        m = findmatch(s2, InertDampingWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        else sp->InertDamping = m;
        break;

      // --- Yes/No options (NO = 0, YES = 1)
      case ALLOW_PONDING:
      case SLOPE_WEIGHTING:
      case SKIP_STEADY_STATE:
      case IGNORE_RAINFALL:
      case IGNORE_SNOWMELT:
      case IGNORE_GWATER:
      case IGNORE_ROUTING:
      case IGNORE_QUALITY:
      case IGNORE_RDII:                                                        //(5.1.004)
        m = findmatch(s2, NoYesWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        switch ( k )
        {
          case ALLOW_PONDING:     sp->AllowPonding    = m;  break;
          case SLOPE_WEIGHTING:   sp->SlopeWeighting  = m;  break;
          case SKIP_STEADY_STATE: sp->SkipSteadyState = m;  break;
          case IGNORE_RAINFALL:   sp->IgnoreRainfall  = m;  break;
          case IGNORE_SNOWMELT:   sp->IgnoreSnowmelt  = m;  break;
          case IGNORE_GWATER:     sp->IgnoreGwater    = m;  break;
          case IGNORE_ROUTING:    sp->IgnoreRouting   = m;  break;
          case IGNORE_QUALITY:    sp->IgnoreQuality   = m;  break;
          case IGNORE_RDII:       sp->IgnoreRDII      = m;  break;                 //(5.1.004)
        }
        break;

      case NORMAL_FLOW_LTD: 
        m = findmatch(s2, NormalFlowWords); 
        //if ( m < 0 ) m = findmatch(s2, NoYesWords);   DEPRECATED             //(5.1.012)
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        sp->NormalFlowLtd = m;
        break;

      case FORCE_MAIN_EQN:
        m = findmatch(s2, ForceMainEqnWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        sp->ForceMainEqn = m;
        break;

      case LINK_OFFSETS:
        m = findmatch(s2, LinkOffsetWords);
        if ( m < 0 ) return error_setInpError(ERR_KEYWORD, s2);
        sp->LinkOffsets = m;
        break;

// TODO: This option is a no op. It should be deprecated.
//      // --- compatibility option for selecting solution method for
//      //     dynamic wave flow routing (NOT CURRENTLY USED)
//      case COMPATIBILITY:
//        if      ( strcomp(s2, "3") ) Compatibility = SWMM3;
//        else if ( strcomp(s2, "4") ) Compatibility = SWMM4;
//        else if ( strcomp(s2, "5") ) Compatibility = SWMM5;
//        else return error_setInpError(ERR_KEYWORD, s2);
//        break;

      // --- routing or lengthening time step (in decimal seconds)
      //     (lengthening time step is used in Courant stability formula
      //     to artificially lengthen conduits for dynamic wave flow routing
      //     (a value of 0 means that no lengthening is used))
      case ROUTE_STEP:
      case LENGTHENING_STEP:
        if ( !getDouble(s2, &tStep) )
        {
            if ( !datetime_strToTime(s2, &aTime) )
            {
                return error_setInpError(ERR_NUMBER, s2);
            }
            else
            {
                datetime_decodeTime(aTime, &h, &m, &s);
                h += 24*(int)aTime;
                s = s + 60*m + 3600*h;
                tStep = s;
            }
        }
        if ( k == ROUTE_STEP )
        {
            if ( tStep <= 0.0 ) return error_setInpError(ERR_NUMBER, s2);
            sp->RouteStep = tStep;
        }
        else sp->LengtheningStep = MAX(0.0, tStep);
        break;

////  Following code section added to release 5.1.008.  ////                   //(5.1.008)

     // --- minimum variable time step for dynamic wave routing
      case MIN_ROUTE_STEP:
        if ( !getDouble(s2, &sp->MinRouteStep) || sp->MinRouteStep < 0.0 )
            return error_setInpError(ERR_NUMBER, s2);
        break;

      case NUM_THREADS:
        m = atoi(s2);
        if ( m < 0 ) return error_setInpError(ERR_NUMBER, s2);
        sp->NumThreads = m;
        break;
 ////

      // --- safety factor applied to variable time step estimates under
      //     dynamic wave flow routing (value of 0 indicates that variable
      //     time step option not used)
      case VARIABLE_STEP:
        if ( !getDouble(s2, &sp->CourantFactor) )
            return error_setInpError(ERR_NUMBER, s2);
        if ( sp->CourantFactor < 0.0 || sp->CourantFactor > 2.0 )
            return error_setInpError(ERR_NUMBER, s2);
        break;

      // --- minimum surface area (ft2 or sq. meters) associated with nodes
      //     under dynamic wave flow routing 
      case MIN_SURFAREA:
        sp->MinSurfArea = atof(s2);
        break;

      // --- minimum conduit slope (%)
      case MIN_SLOPE:
        if ( !getDouble(s2, &sp->MinSlope) )
            return error_setInpError(ERR_NUMBER, s2);
        if ( sp->MinSlope < 0.0 || sp->MinSlope >= 100 )
            return error_setInpError(ERR_NUMBER, s2);
        sp->MinSlope /= 100.0;
        break;

      // --- maximum trials / time step for dynamic wave routing
      case MAX_TRIALS:
        m = atoi(s2);
        if ( m < 0 ) return error_setInpError(ERR_NUMBER, s2);
        sp->MaxTrials = m;
        break;

      // --- head convergence tolerance for dynamic wave routing
      case HEAD_TOL:
        if ( !getDouble(s2, &sp->HeadTol) )
        {
            return error_setInpError(ERR_NUMBER, s2);
        }
        break;

      // --- steady state tolerance on system inflow - outflow
      case SYS_FLOW_TOL:
        if ( !getDouble(s2, &sp->SysFlowTol) )
        {
            return error_setInpError(ERR_NUMBER, s2);
        }
        sp->SysFlowTol /= 100.0;
        break;

      // --- steady state tolerance on nodal lateral inflow
      case LAT_FLOW_TOL:
        if ( !getDouble(s2, &sp->LatFlowTol) )
        {
            return error_setInpError(ERR_NUMBER, s2);
        }
        sp->LatFlowTol /= 100.0;
        break;

      case TEMPDIR: // sp->Temporary Directory
        sstrncpy(sp->TempDir, s2, MAXFNAME);
        break;

    }
    return 0;
}

//=============================================================================

void initPointers(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: assigns NULL to all dynamic arrays for a new project.
//
{
    sp->Gage     = NULL;
    sp->Subcatch = NULL;
    sp->Node     = NULL;
    sp->Outfall  = NULL;
    sp->Divider  = NULL;
    sp->Storage  = NULL;
    sp->Link     = NULL;
    Conduit  = NULL;
    Pump     = NULL;
    Orifice  = NULL;
    Weir     = NULL;
    Outlet   = NULL;
    Pollut   = NULL;
    Landuse  = NULL;
    Pattern  = NULL;
    Curve    = NULL;
    Tseries  = NULL;
    Transect = NULL;
    Shape    = NULL;
    sp->Aquifer    = NULL;
    sp->UnitHyd    = NULL;
    sp->Snowmelt   = NULL;
    Event      = NULL;                                                         //(5.1.011)
    MemPoolAllocated = FALSE;
}

//=============================================================================

void setDefaults(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: assigns default values to project variables.
//
{
   int i, j;

   // Project title & temp. file path
   for (i = 0; i < MAXTITLE; i++) strcpy(sp->Title[i], "");
   strcpy(sp->TempDir, "");

   // Interface files
   sp->Frain.mode      = SCRATCH_FILE;     // Use scratch rainfall file
   sp->Fclimate.mode   = NO_FILE;
   sp->Frunoff.mode    = NO_FILE;
   sp->Frdii.mode      = NO_FILE;
   sp->Fhotstart1.mode = NO_FILE;
   sp->Fhotstart2.mode = NO_FILE;
   sp->Finflows.mode   = NO_FILE;
   sp->Foutflows.mode  = NO_FILE;
   sp->Frain.file      = NULL;
   sp->Fclimate.file   = NULL;
   sp->Frunoff.file    = NULL;
   sp->Frdii.file      = NULL;
   sp->Fhotstart1.file = NULL;
   sp->Fhotstart2.file = NULL;
   sp->Finflows.file   = NULL;
   sp->Foutflows.file  = NULL;
   sp->Fout.file   = NULL;
   sp->Fout.mode   = NO_FILE;

   // Analysis options
   sp->UnitSystem      = US;               // US unit system
   sp->FlowUnits       = CFS;              // CFS flow units
   sp->InfilModel      = HORTON;           // Horton infiltration method
   sp->RouteModel      = KW;               // Kin. wave flow routing method
   sp->AllowPonding    = FALSE;            // No ponding at nodes
   sp->InertDamping    = SOME;             // Partial inertial damping
   sp->NormalFlowLtd   = BOTH;             // Default normal flow limitation
   sp->ForceMainEqn    = H_W;              // Hazen-Williams eqn. for force mains
   sp->LinkOffsets     = DEPTH_OFFSET;     // Use depth for link offsets
   sp->LengtheningStep = 0;                // No lengthening of conduits
   sp->CourantFactor   = 0.0;              // No variable time step 
   sp->MinSurfArea     = 0.0;              // Force use of default min. surface area
   sp->MinSlope        = 0.0;              // No user supplied minimum conduit slope //(5.1.012)
   sp->SkipSteadyState = FALSE;            // Do flow routing in steady state periods
   sp->IgnoreRainfall  = FALSE;            // Analyze rainfall/runoff
   sp->IgnoreRDII      = FALSE;            // Analyze RDII                         //(5.1.004)
   sp->IgnoreSnowmelt  = FALSE;            // Analyze snowmelt 
   sp->IgnoreGwater    = FALSE;            // Analyze groundwater 
   sp->IgnoreRouting   = FALSE;            // Analyze flow routing
   sp->IgnoreQuality   = FALSE;            // Analyze water quality
   sp->WetStep         = 300;              // Runoff wet time step (secs)
   sp->DryStep         = 3600;             // Runoff dry time step (secs)
   sp->RouteStep       = 300.0;            // Routing time step (secs)
   sp->MinRouteStep    = 0.5;              // Minimum variable time step (sec)     //(5.1.008)
   sp->ReportStep      = 900;              // Reporting time step (secs)
   sp->StartDryDays    = 0.0;              // Antecedent dry days
   sp->MaxTrials       = 0;                // Force use of default max. trials 
   sp->HeadTol         = 0.0;              // Force use of default head tolerance
   sp->SysFlowTol      = 0.05;             // System flow tolerance for steady state
   sp->LatFlowTol      = 0.05;             // Lateral flow tolerance for steady state
   sp->NumThreads      = 0;                // Number of parallel threads to use
   sp->NumEvents       = 0;                // Number of detailed routing events    //(5.1.011)

   // Deprecated options
   sp->SlopeWeighting  = TRUE;             // Use slope weighting
//   Compatibility   = SWMM4;            // Use SWMM 4 up/dn weighting method

   // Starting & ending date/time
   sp->StartDate   = datetime_encodeDate(2004, 1, 1);
   sp->StartTime       = datetime_encodeTime(0,0,0);
   sp->StartDateTime   = sp->StartDate + sp->StartTime;
   sp->EndDate         = sp->StartDate;
   sp->EndTime         = 0.0;
   sp->ReportStartDate = NO_DATE;
   sp->ReportStartTime = NO_DATE;
   sp->ReportStep      = 1;
   sp->SweepEnd        = 365;

   // Reporting options
   sp->RptFlags.input         = FALSE;
   sp->RptFlags.continuity    = TRUE;
   sp->RptFlags.flowStats     = TRUE;
   sp->RptFlags.controls      = FALSE;
   sp->RptFlags.subcatchments = FALSE;
   sp->RptFlags.nodes         = FALSE;
   sp->RptFlags.links         = FALSE;
   sp->RptFlags.nodeStats     = FALSE;

   // sp->Temperature data
   sp->Temp.dataSource  = NO_TEMP;
   sp->Temp.tSeries     = -1;
   sp->Temp.ta          = 70.0;
   sp->Temp.elev        = 0.0;
   sp->Temp.anglat      = 40.0;
   sp->Temp.dtlong      = 0.0;
   sp->Temp.tmax        = MISSING;

   // Wind speed data
   sp->Wind.type = MONTHLY_WIND;
   for ( i=0; i<12; i++ ) sp->Wind.aws[i] = 0.0;

   // Snowmelt parameters
   sp->Snow.snotmp      = 34.0;
   sp->Snow.tipm        = 0.5;
   sp->Snow.rnm         = 0.6;

   // Snow areal depletion curves for pervious and impervious surfaces
   for ( i=0; i<2; i++ )
   {
       for ( j=0; j<10; j++) sp->Snow.adc[i][j] = 1.0;
   }

   // Evaporation rates
   sp->Evap.type = CONSTANT_EVAP;
   for (i=0; i<12; i++)
   {
       sp->Evap.monthlyEvap[i] = 0.0;
       sp->Evap.panCoeff[i]    = 1.0;
   }
   sp->Evap.recoveryPattern = -1;
   sp->Evap.recoveryFactor  = 1.0;
   sp->Evap.tSeries = -1;
   sp->Evap.dryOnly = FALSE;

////  Following code segment added to release 5.1.007.  ////                   //(5.1.007)
////
   // Climate adjustments
   for (i = 0; i < 12; i++)
   {
       sp->Adjust.temp[i] = 0.0;   // additive adjustments
       sp->Adjust.evap[i] = 0.0;   // additive adjustments
       sp->Adjust.rain[i] = 1.0;   // multiplicative adjustments
       sp->Adjust.hydcon[i] = 1.0; // hyd. conductivity adjustments                //(5.1.008)
   }
   sp->Adjust.rainFactor = 1.0;
   sp->Adjust.hydconFactor = 1.0;                                                  //(5.1.008)
////
}

//=============================================================================

void openFiles(SWMM_Project *sp, char *f1, char *f2, char *f3)
//
//  Input:   f1 = name of input file
//           f2 = name of report file
//           f3 = name of binary output file
//  Output:  none
//  Purpose: opens a project's input and report files.
//
{
    // --- initialize file pointers to NULL
    sp->Finp.file = NULL;
    sp->Frpt.file = NULL;
    sp->Fout.file = NULL;

    // --- save file names
    sstrncpy(sp->Finp.name, f1, MAXFNAME);
    sstrncpy(sp->Frpt.name, f2, MAXFNAME);
    sstrncpy(sp->Fout.name, f3, MAXFNAME);

    // --- check that file names are not identical
    if (strcomp(f1, f2) || strcomp(f1, f3) || strcomp(f2, f3))
    {
        writecon(FMT11);
        sp->ErrorCode = ERR_FILE_NAME;
        return;
    }

    // --- open input and report files
    if ((sp->Finp.file = fopen(f1,"rt")) == NULL)
    {
        writecon(FMT12);
        writecon(f1);
        sp->ErrorCode = ERR_INP_FILE;
        return;
    }
    if ((sp->Frpt.file = fopen(f2,"wt")) == NULL)
    {
       writecon(FMT13);
       sp->ErrorCode = ERR_RPT_FILE;
       return;
    }
}

//=============================================================================

void createObjects(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: allocates memory for project's objects.
//
//  NOTE: number of each type of object has already been determined in
//        project_readInput().
//
{
    int j, k;

    // --- allocate memory for each category of object
    if ( sp->ErrorCode ) return;
    sp->Gage     = (TGage *)     calloc(sp->Nobjects[GAGE],     sizeof(TGage));
    sp->Subcatch = (TSubcatch *) calloc(sp->Nobjects[SUBCATCH], sizeof(TSubcatch));
    sp->Node     = (TNode *)     calloc(sp->Nobjects[NODE],     sizeof(TNode));
    sp->Outfall  = (TOutfall *)  calloc(sp->Nnodes[OUTFALL],    sizeof(TOutfall));
    sp->Divider  = (TDivider *)  calloc(sp->Nnodes[DIVIDER],    sizeof(TDivider));
    sp->Storage  = (TStorage *)  calloc(sp->Nnodes[STORAGE],    sizeof(TStorage));
    sp->Link     = (TLink *)     calloc(sp->Nobjects[LINK],     sizeof(TLink));
    Conduit  = (TConduit *)  calloc(sp->Nlinks[CONDUIT],    sizeof(TConduit));
    Pump     = (TPump *)     calloc(sp->Nlinks[PUMP],       sizeof(TPump));
    Orifice  = (TOrifice *)  calloc(sp->Nlinks[ORIFICE],    sizeof(TOrifice));
    Weir     = (TWeir *)     calloc(sp->Nlinks[WEIR],       sizeof(TWeir));
    Outlet   = (TOutlet *)   calloc(sp->Nlinks[OUTLET],     sizeof(TOutlet));
    Pollut   = (TPollut *)   calloc(sp->Nobjects[POLLUT],   sizeof(TPollut));
    Landuse  = (TLanduse *)  calloc(sp->Nobjects[LANDUSE],  sizeof(TLanduse));
    Pattern  = (TPattern *)  calloc(sp->Nobjects[TIMEPATTERN],  sizeof(TPattern));
    Curve    = (TTable *)    calloc(sp->Nobjects[CURVE],    sizeof(TTable));
    Tseries  = (TTable *)    calloc(sp->Nobjects[TSERIES],  sizeof(TTable));
    sp->Aquifer  = (TAquifer *)  calloc(sp->Nobjects[AQUIFER],  sizeof(TAquifer));
    sp->UnitHyd  = (TUnitHyd *)  calloc(sp->Nobjects[UNITHYD],  sizeof(TUnitHyd));
    sp->Snowmelt = (TSnowmelt *) calloc(sp->Nobjects[SNOWMELT], sizeof(TSnowmelt));
    Shape    = (TShape *)    calloc(sp->Nobjects[SHAPE],    sizeof(TShape));

////  Added to release 5.1.011.  ////                                          //(5.1.011)
    // --- create array of detailed routing event periods
    Event = (TEvent *) calloc(sp->NumEvents+1, sizeof(TEvent));
    Event[sp->NumEvents].start = BIG;
    Event[sp->NumEvents].end = BIG + 1.0;
////

    // --- create LID objects
    lid_create(sp, sp->Nobjects[LID], sp->Nobjects[SUBCATCH]);

    // --- create control rules
    sp->ErrorCode = controls_create(sp->Nobjects[CONTROL]);
    if ( sp->ErrorCode ) return;

    // --- create cross section transects
    sp->ErrorCode = transect_create(sp->Nobjects[TRANSECT]);
    if ( sp->ErrorCode ) return;

    // --- allocate memory for infiltration data
    infil_create(sp, sp->Nobjects[SUBCATCH], sp->InfilModel);

    // --- allocate memory for water quality state variables
    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        sp->Subcatch[j].initBuildup =
                              (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Subcatch[j].oldQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Subcatch[j].newQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Subcatch[j].pondedQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Subcatch[j].totalLoad  = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
    }
    for (j = 0; j < sp->Nobjects[NODE]; j++)
    {
        sp->Node[j].oldQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Node[j].newQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Node[j].extInflow = NULL;
        sp->Node[j].dwfInflow = NULL;
        sp->Node[j].rdiiInflow = NULL;
        sp->Node[j].treatment = NULL;
    }
    for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        sp->Link[j].oldQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Link[j].newQual = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        sp->Link[j].totalLoad = (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
    }

    // --- allocate memory for land use buildup/washoff functions
    for (j = 0; j < sp->Nobjects[LANDUSE]; j++)
    {
        Landuse[j].buildupFunc =
            (TBuildup *) calloc(sp->Nobjects[POLLUT], sizeof(TBuildup));
        Landuse[j].washoffFunc =
            (TWashoff *) calloc(sp->Nobjects[POLLUT], sizeof(TWashoff));
    }

    // --- allocate memory for subcatchment landuse factors
    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        sp->Subcatch[j].landFactor =
            (TLandFactor *) calloc(sp->Nobjects[LANDUSE], sizeof(TLandFactor));
        for (k = 0; k < sp->Nobjects[LANDUSE]; k++)
        {
            sp->Subcatch[j].landFactor[k].buildup =
                (double *) calloc(sp->Nobjects[POLLUT], sizeof(double));
        }
    }

    // --- initialize buildup & washoff functions
    for (j = 0; j < sp->Nobjects[LANDUSE]; j++)
    {
        for (k = 0; k < sp->Nobjects[POLLUT]; k++)
        {
            Landuse[j].buildupFunc[k].funcType = NO_BUILDUP;
            Landuse[j].buildupFunc[k].normalizer = PER_AREA;
            Landuse[j].washoffFunc[k].funcType = NO_WASHOFF;
        }
    }

    // --- initialize rain gage properties
    for (j = 0; j < sp->Nobjects[GAGE]; j++)
    {
        sp->Gage[j].tSeries = -1;
        strcpy(sp->Gage[j].fname, "");
    }

    // --- initialize subcatchment properties
    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        sp->Subcatch[j].outSubcatch = -1;
        sp->Subcatch[j].outNode     = -1;
        sp->Subcatch[j].infil       = -1;
        sp->Subcatch[j].groundwater = NULL;
        sp->Subcatch[j].gwLatFlowExpr = NULL;                                      //(5.1.007)
        sp->Subcatch[j].gwDeepFlowExpr = NULL;                                     //(5.1.007)
        sp->Subcatch[j].snowpack    = NULL;
        sp->Subcatch[j].lidArea     = 0.0;
        for (k = 0; k < sp->Nobjects[POLLUT]; k++)
        {
            sp->Subcatch[j].initBuildup[k] = 0.0;
        }
    }

    // --- initialize RDII unit hydrograph properties
    for ( j = 0; j < sp->Nobjects[UNITHYD]; j++ ) rdii_initUnitHyd(sp, j);

    // --- initialize snowmelt properties
    for ( j = 0; j < sp->Nobjects[SNOWMELT]; j++ ) snow_initSnowmelt(sp, j);

    // --- initialize storage node exfiltration                                //(5.1.007)
    for (j = 0; j < sp->Nnodes[STORAGE]; j++) sp->Storage[j].exfil = NULL;             //(5.1.007)

    // --- initialize link properties
    for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        sp->Link[j].xsect.type   = -1;
        sp->Link[j].cLossInlet   = 0.0;
        sp->Link[j].cLossOutlet  = 0.0;
        sp->Link[j].cLossAvg     = 0.0;
        sp->Link[j].hasFlapGate  = FALSE;
    }
    for (j = 0; j < sp->Nlinks[PUMP]; j++) Pump[j].pumpCurve  = -1;

    // --- initialize reporting flags
    for (j = 0; j < sp->Nobjects[SUBCATCH]; j++) sp->Subcatch[j].rptFlag = FALSE;
    for (j = 0; j < sp->Nobjects[NODE]; j++) sp->Node[j].rptFlag = FALSE;
    for (j = 0; j < sp->Nobjects[LINK]; j++) sp->Link[j].rptFlag = FALSE;

    //  --- initialize curves, time series, and time patterns
    for (j = 0; j < sp->Nobjects[CURVE]; j++)   table_init(&Curve[j]);
    for (j = 0; j < sp->Nobjects[TSERIES]; j++) table_init(&Tseries[j]);
    for (j = 0; j < sp->Nobjects[TIMEPATTERN]; j++) inflow_initDwfPattern(j);
}

//=============================================================================

void deleteObjects(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory allocated for a project's objects.
//
//  NOTE: care is taken to first free objects that are properties of another
//        object before the latter is freed (e.g., we must free a
//        subcatchment's land use factors before freeing the subcatchment).
//
{
    int j, k;

    // --- free memory for landuse factors & groundwater
    if ( sp->Subcatch ) for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        for (k = 0; k < sp->Nobjects[LANDUSE]; k++)
        {
            FREE(sp->Subcatch[j].landFactor[k].buildup);
        }
        FREE(sp->Subcatch[j].landFactor);
        FREE(sp->Subcatch[j].groundwater);
        gwater_deleteFlowExpression(sp, j);
        FREE(sp->Subcatch[j].snowpack);
    }

    // --- free memory for buildup/washoff functions
    if ( Landuse ) for (j = 0; j < sp->Nobjects[LANDUSE]; j++)
    {
        FREE(Landuse[j].buildupFunc);
        FREE(Landuse[j].washoffFunc)
    }

    // --- free memory for water quality state variables
    if ( sp->Subcatch ) for (j = 0; j < sp->Nobjects[SUBCATCH]; j++)
    {
        FREE(sp->Subcatch[j].initBuildup);
        FREE(sp->Subcatch[j].oldQual);
        FREE(sp->Subcatch[j].newQual);
        FREE(sp->Subcatch[j].pondedQual);
        FREE(sp->Subcatch[j].totalLoad);
    }
    if ( sp->Node ) for (j = 0; j < sp->Nobjects[NODE]; j++)
    {
        FREE(sp->Node[j].oldQual);
        FREE(sp->Node[j].newQual);
    }
    if ( sp->Link ) for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        FREE(sp->Link[j].oldQual);
        FREE(sp->Link[j].newQual);
        FREE(sp->Link[j].totalLoad);
    }

    // --- free memory used for rainfall infiltration
    infil_delete();

////  Added for release 5.1.007.  ////                                         //(5.1.007)
////
    // --- free memory used for storage exfiltration
    if ( sp->Node ) for (j = 0; j < sp->Nnodes[STORAGE]; j++)
    {
        if ( sp->Storage[j].exfil )
        {
            FREE(sp->Storage[j].exfil->btmExfil);
            FREE(sp->Storage[j].exfil->bankExfil);
            FREE(sp->Storage[j].exfil);
        }
    }
////

    // --- free memory used for outfall pollutants loads                       //(5.1.008)
    if ( sp->Node ) for (j = 0; j < sp->Nnodes[OUTFALL]; j++)                          //(5.1.008)
        FREE(sp->Outfall[j].wRouted);                                              //(5.1.008)

    // --- free memory used for nodal inflows & treatment functions
    if ( sp->Node ) for (j = 0; j < sp->Nobjects[NODE]; j++)
    {
        inflow_deleteExtInflows(sp, j);
        inflow_deleteDwfInflows(sp, j);
        rdii_deleteRdiiInflow(sp, j);
        treatmnt_delete(sp, j);
    }

    // --- delete table entries for curves and time series
    if ( Tseries ) for (j = 0; j < sp->Nobjects[TSERIES]; j++)
        table_deleteEntries(&Tseries[j]);
    if ( Curve ) for (j = 0; j < sp->Nobjects[CURVE]; j++)
        table_deleteEntries(&Curve[j]);

    // --- delete cross section transects
    transect_delete();

    // --- delete control rules
    controls_delete();

    // --- delete LIDs
    lid_delete();

    // --- now free each major category of object
    FREE(sp->Gage);
    FREE(sp->Subcatch);
    FREE(sp->Node);
    FREE(sp->Outfall);
    FREE(sp->Divider);
    FREE(sp->Storage);
    FREE(sp->Link);
    FREE(Conduit);
    FREE(Pump);
    FREE(Orifice);
    FREE(Weir);
    FREE(Outlet);
    FREE(Pollut);
    FREE(Landuse);
    FREE(Pattern);
    FREE(Curve);
    FREE(Tseries);
    FREE(sp->Aquifer);
    FREE(sp->UnitHyd);
    FREE(sp->Snowmelt);
    FREE(Shape);
    FREE(Event);                                                               //(5.1.011)
}

//=============================================================================

void createHashTables(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: allocates memory for object ID hash tables
//
{   int j;
    MemPoolAllocated = FALSE;
    for (j = 0; j < MAX_OBJ_TYPES ; j++)
    {
        Htable[j] = HTcreate();
        if ( Htable[j] == NULL ) report_writeErrorMsg(sp, ERR_MEMORY, "");
    }

    // --- initialize memory pool used to store object ID's
    if ( AllocInit() == NULL ) report_writeErrorMsg(sp, ERR_MEMORY, "");
    else MemPoolAllocated = TRUE;
}

//=============================================================================

void deleteHashTables()
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory allocated for object ID hash tables
//
{
    int j;
    for (j = 0; j < MAX_OBJ_TYPES; j++)
    {
        if ( Htable[j] != NULL ) HTfree(Htable[j]);
    }

    // --- free object ID memory pool
    if ( MemPoolAllocated ) AllocFreePool();
}

//=============================================================================












