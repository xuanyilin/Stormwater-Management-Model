//-----------------------------------------------------------------------------
//   rdii.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//             04/04/14   (Build 5.1.003)
//             04/14/14   (Build 5.1.004)
//             09/15/14   (Build 5.1.007)
//   Author:   L. Rossman (EPA)
//             R. Dickinson (CDM)
//
//   RDII processing functions.
//
//   Note: RDII means rainfall dependent infiltration/inflow,
//         UH means unit hydrograph.
//
//   Build 5.1.007:
//   - Ignore RDII option implemented.
//   - Rainfall climate adjustment implemented.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "headers.h"

//-----------------------------------------------------------------------------
// Definition of 4-byte integer, 4-byte real and 8-byte real types
//-----------------------------------------------------------------------------
#define INT4  int
#define REAL4 float
#define REAL8 double
#define FILE_STAMP "SWMM5-RDII"

//-----------------------------------------------------------------------------
// Constants
//-----------------------------------------------------------------------------
const double ZERO_RDII = 0.0001;       // Minimum non-zero RDII inflow (cfs)
const char   FileStamp[] = FILE_STAMP;

//-----------------------------------------------------------------------------
// Data Structures
//-----------------------------------------------------------------------------
enum FileTypes {BINARY, TEXT};         // File mode types


//-----------------------------------------------------------------------------
// Imported Variables
//-----------------------------------------------------------------------------
//#ifdef __cplusplus
//extern const double Qcf[];             // flow units conversion factors
                                       // (see swmm5.c)
//#else
extern const double Qcf[];                   // flow units conversion factors
                                       // (see swmm5.c)
//#endif

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  rdii_readRdiiInflow     (called from parseLine in input.c)
//  rdii_deleteRdiiInflow   (called from deleteObjects in project.c)
//  rdii_initUnitHyd        (called from createObjects in project.c)
//  rdii_readUnitHydParams  (called from parseLine in input.c)
//  rdii_openRdii           (called from rain_open)
//  rdii_closeRdii          (called from rain_close)
//  rdii_getNumRdiiFlows    (called from addRdiiInflows in routing.c)
//  rdii_getRdiiFlow        (called from addRdiiInflows in routing.c)

//-----------------------------------------------------------------------------
// Function Declarations
//-----------------------------------------------------------------------------
// --- functions used to create a RDII file
static int    readOldUHFormat(SWMM_Project *sp, int j, int m, char* tok[], int ntoks);
static void   setUnitHydParams(SWMM_Project *sp, int j, int i, int m, double x[]);
static void   createRdiiFile(SWMM_Project *sp);
static int    getNumRdiiNodes(SWMM_Project *sp);
static void   validateRdii(SWMM_Project *sp);

static void   openRdiiProcessor(SWMM_Project *sp);
static int    allocRdiiMemory(SWMM_Project *sp);
static int    getRainInterval(SWMM_Project *sp, int i);
static int    getMaxPeriods(SWMM_Project *sp, int i, int k);
static void   initGageData(SWMM_Project *sp);
static void   initUnitHydData(SWMM_Project *sp);
static int    openNewRdiiFile(SWMM_Project *sp);
static void   getRainfall(SWMM_Project *sp, DateTime currentDate);

static double applyIA(SWMM_Project *sp, int j, int k, DateTime aDate, double dt,
              double rainDepth);
static void   updateDryPeriod(SWMM_Project *sp, int j, int k, double rain,
        int gageInterval);
static void   getUnitHydRdii(SWMM_Project *sp, DateTime currentDate);
static double getUnitHydConvol(SWMM_Project *sp, int j, int k, int gageInterval);
static double getUnitHydOrd(SWMM_Project *sp, int j, int m, int k, double t);

static int    getNodeRdii(SWMM_Project *sp);
static void   saveRdiiFlows(SWMM_Project *sp, DateTime currentDate);
static void   closeRdiiProcessor(SWMM_Project *sp);
static void   freeRdiiMemory(SWMM_Project *sp);

// --- functions used to read an existing RDII file
static int   readRdiiFileHeader(SWMM_Project *sp);
static void  readRdiiFlows(SWMM_Project *sp);

static void  openRdiiTextFile(SWMM_Project *sp);
static int   readRdiiTextFileHeader(SWMM_Project *sp);
static void  readRdiiTextFlows(SWMM_Project *sp);

//=============================================================================
//                   Management of RDII-Related Data
//=============================================================================

int rdii_readRdiiInflow(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads properties of an RDII inflow from a line of input.
//
{
    int    j, k;
    double a;
    TRdiiInflow* inflow;

    // --- check for proper number of items
    if ( ntoks < 3 ) return error_setInpError(sp, ERR_ITEMS, "");

    // --- check that node receiving RDII exists
    j = project_findObject(sp, NODE, tok[0]);
    if ( j < 0 ) return error_setInpError(sp, ERR_NAME, tok[0]);

    // --- check that RDII unit hydrograph exists
    k = project_findObject(sp, UNITHYD, tok[1]);
    if ( k < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);

    // --- read in sewer area value
    if ( !getDouble(tok[2], &a) || a < 0.0 )
        return error_setInpError(sp, ERR_NUMBER, tok[2]);

    // --- create the RDII inflow object if it doesn't already exist
    inflow = sp->Node[j].rdiiInflow;
    if ( inflow == NULL )
    {
        inflow = (TRdiiInflow *) malloc(sizeof(TRdiiInflow));
        if ( !inflow ) return error_setInpError(sp, ERR_MEMORY, "");
    }

    // --- assign UH & area to inflow object
    inflow->unitHyd = k;
    inflow->area = a / UCF(sp, LANDAREA);

    // --- assign inflow object to node
    sp->Node[j].rdiiInflow = inflow;
    return 0;
}

//=============================================================================

void rdii_initUnitHyd(SWMM_Project *sp, int j)
//
//  Input:   j = UH group index
//  Output:  none
//  Purpose: initializes properties of a unit hydrograph group.
//
{
    int i;                             // individual UH index
    int m;                             // month index

    for ( m=0; m<12; m++)
    {
        for (i=0; i<3; i++)
        {
            sp->UnitHyd[j].iaMax[m][i]   = 0.0;
            sp->UnitHyd[j].iaRecov[m][i] = 0.0;
            sp->UnitHyd[j].iaInit[m][i]  = 0.0;
            sp->UnitHyd[j].r[m][i]       = 0.0;
            sp->UnitHyd[j].tPeak[m][i]   = 0;
            sp->UnitHyd[j].tBase[m][i]   = 0;
        }
    }
}

//=============================================================================

int rdii_readUnitHydParams(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads parameters of an RDII unit hydrograph from a line of input.
//
{
    int i, j, k, m, g;
    double x[6];

    // --- check that RDII UH object exists in database
    j = project_findObject(sp, UNITHYD, tok[0]);
    if ( j < 0 ) return error_setInpError(sp, ERR_NAME, tok[0]);

    // --- assign UH ID to name in hash table
    if ( sp->UnitHyd[j].ID == NULL )
        sp->UnitHyd[j].ID = project_findID(sp, UNITHYD, tok[0]);

    // --- line has 2 tokens; assign rain gage to UH object
    if ( ntoks == 2 )
    {
        g = project_findObject(sp, GAGE, tok[1]);
        if ( g < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
        sp->UnitHyd[j].rainGage = g;
        return 0;
    }
    else if ( ntoks < 6 ) return error_setInpError(sp, ERR_ITEMS, "");

    // --- find which month UH params apply to
    m = datetime_findMonth(tok[1]);
    if ( m == 0 )
    {
        if ( !match(tok[1], w_ALL) )
            return error_setInpError(sp, ERR_KEYWORD, tok[1]);
    }

    // --- find type of UH being specified
    k = findmatch(tok[2], UHTypeWords);

    // --- if no type match, try using older UH line format
    if ( k < 0 ) return readOldUHFormat(sp, j, m, tok, ntoks);

    // --- read the R-T-K parameters
    for ( i = 0; i < 3; i++ )
    {
        if ( ! getDouble(tok[i+3], &x[i]) )
            return error_setInpError(sp, ERR_NUMBER, tok[i+3]);
    }

    // --- read the IA parameters if present
    for (i = 3; i < 6; i++)
    {
        x[i] = 0.0;
        if ( ntoks > i+3 )
        {
            if ( ! getDouble(tok[i+3], &x[i]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i+2]);
        }
    }

    // --- save UH params
    setUnitHydParams(sp, j, k, m, x);
    return 0;
}

//=============================================================================

int readOldUHFormat(SWMM_Project *sp, int j, int m, char* tok[], int ntoks)
//
//  Input:   j = unit hydrograph index
//           m = month of year (0 = all months)
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads parameters of a set of RDII unit hydrographs from a line of
//           input.
//
{
    int    i, k;
    double p[9], x[6];

    // --- check for proper number of tokens
    if ( ntoks < 11 ) return error_setInpError(sp, ERR_ITEMS, "");

    // --- read 3 sets of r-t-k values
    for ( i = 0; i < 9; i++ )
    {
        if ( ! getDouble(tok[i+2], &p[i]) )
            return error_setInpError(sp, ERR_NUMBER, tok[i+2]);
    }

    // --- read initial abstraction parameters
    for (i = 0; i < 3; i++)
    {
        x[i+3] = 0.0;
        if ( ntoks > i+11 )
        {
            if ( ! getDouble(tok[i+11], &x[i+3]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i+11]);
        }
    }

    // --- save UH parameters
    for ( k = 0; k < 3; k++)
    {
        for ( i = 0; i < 3; i++)
        {
            x[i] = p[3*k + i];
            setUnitHydParams(sp, j, k, m, x);
        }
    }
    return 0;
}

//=============================================================================

void setUnitHydParams(SWMM_Project *sp, int j, int i, int m, double x[])
//
//  Input:   j = unit hydrograph index
//           i = type of UH response (short, medium or long term)
//           m = month of year (0 = all months)
//           x = array of UH parameters
//  Output:  none
//  Purpose: assigns parameters to a unit hydrograph for a specified month of year.
//
{
    int    m1, m2;                     // start/end month indexes
    double t,                          // UH time to peak (hrs)
           k,                          // UH k-value
           tBase;                      // UH base time (hrs)

    // --- find range of months that share same parameter values
    if ( m == 0 )
    {
        m1 = 0;
        m2 = 11;
    }
    else
    {
        m1 = m-1;
        m2 = m1;
    }

    // --- for each month in the range
    for (m=m1; m<=m2; m++)
    {
        // --- set UH response ratio, time to peak, & base time
        sp->UnitHyd[j].r[m][i] = x[0];
        t = x[1];
        k = x[2];
        tBase = t * (1.0 + k);                              // hours
        sp->UnitHyd[j].tPeak[m][i] = (long)(t * 3600.);         // seconds
        sp->UnitHyd[j].tBase[m][i] = (long)(tBase * 3600.);     // seconds

        // -- set initial abstraction parameters
        sp->UnitHyd[j].iaMax[m][i]   = x[3];
        sp->UnitHyd[j].iaRecov[m][i] = x[4];
        sp->UnitHyd[j].iaInit[m][i]  = x[5];
    }
}

//=============================================================================

void rdii_deleteRdiiInflow(SWMM_Project *sp, int j)
//
//  Input:   j = node index
//  Output:  none
//  Purpose: deletes the RDII inflow object for a node.
//
{
    if ( sp->Node[j].rdiiInflow )
    {
        free(sp->Node[j].rdiiInflow);
        sp->Node[j].rdiiInflow = NULL;
    }
}


//=============================================================================
//                 Reading Inflow Data From a RDII File
//=============================================================================

void rdii_openRdii(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: opens an exisiting RDII interface file or creates a new one.
//
{
    char  fStamp[] = FILE_STAMP;

    TRdiiShared *rd = &sp->RdiiShared;

    rd->RdiiNodeIndex = NULL;
    rd->RdiiNodeFlow = NULL;
    rd->NumRdiiNodes = 0;
    rd->RdiiStartDate = NO_DATE;

    // --- create the RDII file if existing file not being used
    if ( sp->IgnoreRDII ) return;                                                  //(5.1.004)
    if ( sp->Frdii.mode != USE_FILE ) createRdiiFile(sp);
    if ( sp->Frdii.mode == NO_FILE || sp->ErrorCode ) return;

    // --- try to open the RDII file in binary mode
    sp->Frdii.file = fopen(sp->Frdii.name, "rb");
    if ( sp->Frdii.file == NULL)
    {
        if ( sp->Frdii.mode == SCRATCH_FILE )
        {
            report_writeErrorMsg(sp, ERR_RDII_FILE_SCRATCH, "");
        }
        else
        {
            report_writeErrorMsg(sp, ERR_RDII_FILE_OPEN, sp->Frdii.name);
        }
        return;
    }

    // --- check for valid file stamp
    fread(fStamp, sizeof(char), strlen(FileStamp), sp->Frdii.file);
    if ( strcmp(fStamp, FileStamp) == 0 )
    {
        rd->RdiiFileType = BINARY;
        sp->ErrorCode = readRdiiFileHeader(sp);
    }

    // --- if stamp invalid try to open the file in text mode
    else
    {
        fclose(sp->Frdii.file);
        rd->RdiiFileType = TEXT;
        openRdiiTextFile(sp);
    }

    // --- catch any error
    if ( sp->ErrorCode )
    {
        report_writeErrorMsg(sp, sp->ErrorCode, sp->Frdii.name);
    }

    // --- read the first set of RDII flows form the file
    else readRdiiFlows(sp);
}

//=============================================================================

void openRdiiTextFile(SWMM_Project *sp)
{
    // --- try to open the RDII file in text mode
    sp->Frdii.file = fopen(sp->Frdii.name, "rt");
    if ( sp->Frdii.file == NULL)
    {
        if ( sp->Frdii.mode == SCRATCH_FILE )
        {
            report_writeErrorMsg(sp, ERR_RDII_FILE_SCRATCH, "");
        }
        else
        {
            report_writeErrorMsg(sp, ERR_RDII_FILE_OPEN, sp->Frdii.name);
        }
        return;
    }

    // --- read header records from file
    sp->ErrorCode = readRdiiTextFileHeader(sp);
    if ( sp->ErrorCode )
    {
        report_writeErrorMsg(sp, sp->ErrorCode, sp->Frdii.name);
    }
}

//=============================================================================

void rdii_closeRdii(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: closes the RDII interface file.
//
{
    TRdiiShared *rd = &sp->RdiiShared;

    if ( sp->Frdii.file ) fclose(sp->Frdii.file);
    if ( sp->Frdii.mode == SCRATCH_FILE ) remove(sp->Frdii.name);
    FREE(rd->RdiiNodeIndex);
    FREE(rd->RdiiNodeFlow);
}

//=============================================================================

int rdii_getNumRdiiFlows(SWMM_Project *sp, DateTime aDate)
//
//  Input:   aDate = current date/time
//  Output:  returns 0 if no RDII flow or number of nodes with RDII inflows
//  Purpose: finds number of RDII inflows at a specified date.
//
{
    TRdiiShared *rd = &sp->RdiiShared;

    // --- default result is 0 indicating no RDII inflow at specified date
    if ( rd->NumRdiiNodes == 0 ) return 0;
    if ( !sp->Frdii.file ) return 0;

    // --- keep reading RDII file as need be
    while ( !feof(sp->Frdii.file) )
    {
        // --- return if date of current RDII inflow not reached yet
        if ( rd->RdiiStartDate == NO_DATE ) return 0;
        if ( aDate < rd->RdiiStartDate ) return 0;

        // --- return RDII node count if specified date falls
        //     within time interval of current RDII inflow
        if ( aDate < rd->RdiiEndDate ) return rd->NumRdiiNodes;

        // --- otherwise get next date and RDII flow values from file
        else readRdiiFlows(sp);
    }
    return 0;
}

//=============================================================================

void rdii_getRdiiFlow(SWMM_Project *sp, int i, int* j, double* q)
//
//  Input:   i = RDII node index
//           j = pointer to project node index
//           q = pointer to RDII flow rate
//  Output:  sets node index and RDII inflow for node
//  Purpose: finds index and current RDII inflow for an RDII node.
//
{
    TRdiiShared *rd = &sp->RdiiShared;

    if ( i >= 0 && i < rd->NumRdiiNodes )
    {
        *j = rd->RdiiNodeIndex[i];
        *q = rd->RdiiNodeFlow[i];
    }
}

//=============================================================================

int readRdiiFileHeader(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads header information from a binary RDII file.
//
{
    int i, j;

    TRdiiShared *rd = &sp->RdiiShared;

    // --- extract time step and number of RDII nodes
    fread(&rd->RdiiStep, sizeof(INT4), 1, sp->Frdii.file);
    if ( rd->RdiiStep <= 0 ) return ERR_RDII_FILE_FORMAT;
    fread(&rd->NumRdiiNodes, sizeof(INT4), 1, sp->Frdii.file);
    if ( rd->NumRdiiNodes <= 0 ) return ERR_RDII_FILE_FORMAT;

    // --- allocate memory for RdiiNodeIndex & RdiiNodeFlow arrays
    rd->RdiiNodeIndex = (int *) calloc(rd->NumRdiiNodes, sizeof(int));
    if ( !rd->RdiiNodeIndex ) return ERR_MEMORY;
    rd->RdiiNodeFlow = (REAL4 *) calloc(rd->NumRdiiNodes, sizeof(REAL4));              //(5.1.003)
    if ( !rd->RdiiNodeFlow ) return ERR_MEMORY;

    // --- read indexes of RDII nodes
    if ( feof(sp->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
    fread(rd->RdiiNodeIndex, sizeof(INT4), rd->NumRdiiNodes, sp->Frdii.file);
    for ( i = 0; i < rd->NumRdiiNodes; i++ )
    {
        j = rd->RdiiNodeIndex[i];
        if ( sp->Node[j].rdiiInflow == NULL ) return ERR_RDII_FILE_FORMAT;
    }
    if ( feof(sp->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
    return 0;
}

//=============================================================================

int readRdiiTextFileHeader(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads header information from a text RDII file.
//
{
    int   i;
    char  line[MAXLINE+1];             // line from RDII data file
    char  s1[MAXLINE+1];               // general string variable
    char  s2[MAXLINE+1];

    TRdiiShared *rd = &sp->RdiiShared;

    // --- check for correct file type
    fgets(line, MAXLINE, sp->Frdii.file);
    sscanf(line, "%s", s1);
    if ( strcmp(s1, "SWMM5") != 0 ) return ERR_RDII_FILE_FORMAT;

    // --- skip title line
    fgets(line, MAXLINE, sp->Frdii.file);

    // --- read RDII UH time step interval (sec)
    rd->RdiiStep = 0;
    fgets(line, MAXLINE, sp->Frdii.file);
    sscanf(line, "%d", &rd->RdiiStep);
    if ( rd->RdiiStep <= 0 ) return ERR_RDII_FILE_FORMAT;

    // --- skip over line with number of constituents (= 1 for RDII)
    fgets(line, MAXLINE, sp->Frdii.file);

    // --- read flow units
    fgets(line, MAXLINE, sp->Frdii.file);
    sscanf(line, "%s %s", s1, s2);
    rd->RdiiFlowUnits = findmatch(s2, FlowUnitWords);
    if ( rd->RdiiFlowUnits < 0 ) return ERR_RDII_FILE_FORMAT;

    // --- read number of RDII nodes
    fgets(line, MAXLINE, sp->Frdii.file);
    if ( sscanf(line, "%d", &rd->NumRdiiNodes) < 1 ) return ERR_RDII_FILE_FORMAT;

    // --- allocate memory for RdiiNodeIndex & RdiiNodeFlow arrays
    rd->RdiiNodeIndex = (int *) calloc(rd->NumRdiiNodes, sizeof(int));
    if ( !rd->RdiiNodeIndex ) return ERR_MEMORY;
    rd->RdiiNodeFlow = (REAL4 *) calloc(rd->NumRdiiNodes, sizeof(REAL4));              //(5.1.003)
    if ( !rd->RdiiNodeFlow ) return ERR_MEMORY;

    // --- read names of RDII nodes from file & save their indexes
    for ( i=0; i<rd->NumRdiiNodes; i++ )
    {
        if ( feof(sp->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
        fgets(line, MAXLINE, sp->Frdii.file);
        sscanf(line, "%s", s1);
        rd->RdiiNodeIndex[i] = project_findObject(sp, NODE, s1);
    }

    // --- skip column heading line
    if ( feof(sp->Frdii.file) ) return ERR_RDII_FILE_FORMAT;
    fgets(line, MAXLINE, sp->Frdii.file);
    return 0;
}

//=============================================================================

void readRdiiFlows(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reads date and flow values of next RDII inflows from RDII file.
//
{
    TRdiiShared *rd = &sp->RdiiShared;

    if ( rd->RdiiFileType == TEXT ) readRdiiTextFlows(sp);
    else
    {
        rd->RdiiStartDate = NO_DATE;
        rd->RdiiEndDate = NO_DATE;
        if ( feof(sp->Frdii.file) ) return;
        fread(&rd->RdiiStartDate, sizeof(DateTime), 1, sp->Frdii.file);
        if ( rd->RdiiStartDate == NO_DATE ) return;
        if ( fread(rd->RdiiNodeFlow, sizeof(REAL4), rd->NumRdiiNodes, sp->Frdii.file)      //(5.1.003)
            < (size_t)rd->NumRdiiNodes ) rd->RdiiStartDate = NO_DATE;
        else rd->RdiiEndDate = datetime_addSeconds(rd->RdiiStartDate, rd->RdiiStep);
    }
}

//=============================================================================

void readRdiiTextFlows(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reads date and flow values of next RDII inflows from RDII file.
//
{
    int    i, n;
    int    yr = 0, mon = 0, day = 0,
		   hr = 0, min = 0, sec = 0;   // year, month, day, hour, minute, second
    double x;                          // RDII flow in original units          //(5.1.003)
    char   line[MAXLINE+1];            // line from RDII data file
    char   s[MAXLINE+1];               // node ID label (not used)

    TRdiiShared *rd = &sp->RdiiShared;

    rd->RdiiStartDate = NO_DATE;
    for (i = 0; i < rd->NumRdiiNodes; i++)
    {
        if ( feof(sp->Frdii.file) ) return;
        fgets(line, MAXLINE, sp->Frdii.file);
        n = sscanf(line, "%s %d %d %d %d %d %d %f",
            s, &yr, &mon, &day, &hr, &min, &sec, &x);
        if ( n < 8 ) return;
        rd->RdiiNodeFlow[i] = (REAL4)(x / Qcf[rd->RdiiFlowUnits]);                     //(5.1.003)
    }
    rd->RdiiStartDate = datetime_encodeDate(yr, mon, day) +
                    datetime_encodeTime(hr, min, sec);
    rd->RdiiEndDate = datetime_addSeconds(rd->RdiiStartDate, rd->RdiiStep);
}


//=============================================================================
//                   Creation of a RDII Interface File
//=============================================================================

void createRdiiFile(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: computes time history of RDII inflows and saves them to file.
//
{
    int      hasRdii;                  // true when total RDII > 0
    double   elapsedTime;              // current elapsed time (sec)
    double   duration;                 // duration being analyzed (sec)
    DateTime currentDate;              // current calendar date/time

    TRdiiShared *rd = &sp->RdiiShared;

    // --- set RDII reporting time step to Runoff wet step
    rd->RdiiStep = sp->WetStep;

    // --- count nodes with RDII data
    rd->NumRdiiNodes = getNumRdiiNodes(sp);

    // --- if no RDII nodes then re-set RDII file usage to NO_FILE
    if ( rd->NumRdiiNodes == 0 )
    {
        sp->Frdii.mode = NO_FILE;
        return;
    }

    // --- otherwise set file usage to SCRATCH if originally set to NO_FILE
    else if ( sp->Frdii.mode == NO_FILE ) sp->Frdii.mode = SCRATCH_FILE;

    // --- validate RDII data
    validateRdii(sp);
    initGageData(sp);
    if ( sp->ErrorCode ) return;

    // --- open RDII processing system
    openRdiiProcessor(sp);
    if ( !sp->ErrorCode )
    {
        // --- initialize rain gage & UH processing data
        initUnitHydData(sp);

        // --- convert total simulation duration from millisec to sec
        duration = sp->TotalDuration / 1000.0;

        // --- examine rainfall record over each RdiiStep time step
        elapsedTime = 0.0;
        while ( elapsedTime <= duration && !sp->ErrorCode )
        {
            // --- compute current calendar date/time
            currentDate = sp->StartDateTime + elapsedTime / SECperDAY;

            // --- update rainfall at all rain gages
            getRainfall(sp, currentDate);

            // --- compute convolutions of past rainfall with UH's
            getUnitHydRdii(sp, currentDate);

            // --- find RDII at all nodes
            hasRdii = getNodeRdii(sp);

            // --- save RDII at all nodes to file for current date
            if ( hasRdii ) saveRdiiFlows(sp, currentDate);

            // --- advance one time step
            elapsedTime += rd->RdiiStep;
        }
    }

    // --- close RDII processing system
    closeRdiiProcessor(sp);
}

//=============================================================================

int  getNumRdiiNodes(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns node count
//  Purpose: counts number of nodes that receive RDII inflow.
//
{
    int j,                             // node index
        n;                             // node count

    n = 0;
    for (j=0; j<sp->Nobjects[NODE]; j++)
    {
        if ( sp->Node[j].rdiiInflow ) n++;
    }
    return n;
}

//=============================================================================

void validateRdii(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: validates UH and RDII inflow object data.
//
{
    int    i,                          // node index
           j,                          // UH group index
           k,                          // individual UH index
           m;                          // month index
    double rsum;                       // sum of UH r-values
//  long   gageInterval;               // rain gage time interval

    // --- check each unit hydrograph for consistency
    for (j=0; j<sp->Nobjects[UNITHYD]; j++)
    {
        for (m=0; m<12; m++)
        {
            rsum = 0.0;
            for (k=0; k<3; k++)
            {
                // --- if no base time then UH doesn't exist
                if ( sp->UnitHyd[j].tBase[m][k] == 0 ) continue;

                // --- restriction on time to peak being less than the
                //     rain gage's recording interval no longer applies

                // --- can't have negative UH parameters
                if ( sp->UnitHyd[j].tPeak[m][k] < 0.0 )
                {
                    report_writeErrorMsg(sp, ERR_UNITHYD_TIMES, sp->UnitHyd[j].ID);
                }

                // --- can't have negative UH response ratio
                if ( sp->UnitHyd[j].r[m][k] < 0.0 )
                {
                    report_writeErrorMsg(sp, ERR_UNITHYD_RATIOS, sp->UnitHyd[j].ID);
                }
                else rsum += sp->UnitHyd[j].r[m][k];
            }
            if ( rsum > 1.01 )
            {
                report_writeErrorMsg(sp, ERR_UNITHYD_RATIOS, sp->UnitHyd[j].ID);
            }
        }
    }

    // --- check each node's RDII inflow object
    for (i=0; i<sp->Nobjects[NODE]; i++)
    {
        if ( sp->Node[i].rdiiInflow )
        {
            // --- check that sewer area is non-negative
            if ( sp->Node[i].rdiiInflow->area < 0.0 )
            {
                report_writeErrorMsg(sp, ERR_RDII_AREA, sp->Node[i].ID);
            }
        }
    }
}

//=============================================================================

void openRdiiProcessor(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: opens RDII processing system.
//
{
    int j;                             // object index
    int n;                             // RDII node count

    TRdiiShared *rd = &sp->RdiiShared;

    // --- set RDII processing arrays to NULL
    rd->UHGroup = NULL;
    rd->RdiiNodeIndex = NULL;
    rd->RdiiNodeFlow = NULL;
    rd->TotalRainVol = 0.0;
    rd->TotalRdiiVol = 0.0;

    // --- allocate memory used for RDII processing
    if ( !allocRdiiMemory(sp) )
    {
        report_writeErrorMsg(sp, ERR_MEMORY, "");
        return;
    }

    // --- open & initialize RDII file
    if ( !openNewRdiiFile(sp) )
    {
        report_writeErrorMsg(sp, ERR_RDII_FILE_SCRATCH, "");
        return;
    }

    // --- identify index of each node with RDII inflow
    n = 0;
    for (j=0; j<sp->Nobjects[NODE]; j++)
    {
        if ( sp->Node[j].rdiiInflow )
        {
            rd->RdiiNodeIndex[n] = j;
            n++;
        }
    }
}

//=============================================================================

int  allocRdiiMemory(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: allocates memory used for RDII processing .
//
//
{
    int i;                             // UH group index
    int k;                             // UH index
    int n;                             // number of past rain periods

    TRdiiShared *rd = &sp->RdiiShared;

    // --- allocate memory for RDII processing data for UH groups
    rd->UHGroup = (TUHGroup *) calloc(sp->Nobjects[UNITHYD], sizeof(TUHGroup));
    if ( !rd->UHGroup ) return FALSE;

    // --- allocate memory for past rainfall data for each UH in each group
    for (i=0; i<sp->Nobjects[UNITHYD]; i++)
    {
        rd->UHGroup[i].rainInterval = getRainInterval(sp, i);
        for (k=0; k<3; k++)
        {
            rd->UHGroup[i].uh[k].pastRain = NULL;
            rd->UHGroup[i].uh[k].pastMonth = NULL;
            rd->UHGroup[i].uh[k].maxPeriods = getMaxPeriods(sp, i, k);
            n = rd->UHGroup[i].uh[k].maxPeriods;
            if ( n > 0 )
            {
                rd->UHGroup[i].uh[k].pastRain =
                    (double *) calloc(n, sizeof(double));
                if ( !rd->UHGroup[i].uh[k].pastRain ) return FALSE;
                rd->UHGroup[i].uh[k].pastMonth =
                    (char *) calloc(n, sizeof(char));
                if ( !rd->UHGroup[i].uh[k].pastMonth ) return FALSE;
            }
        }
    }

    // --- allocate memory for RDII indexes & inflow at each node w/ RDII data
    rd->RdiiNodeIndex = (int *) calloc(rd->NumRdiiNodes, sizeof(int));
    if ( !rd->RdiiNodeIndex ) return FALSE;
    rd->RdiiNodeFlow = (REAL4 *) calloc(rd->NumRdiiNodes, sizeof(REAL4));              //(5.1.003)
    if ( !rd->RdiiNodeFlow ) return FALSE;
    return TRUE;
}

//=============================================================================

int  getRainInterval(SWMM_Project *sp, int i)
//
//  Input:   i = UH group index
//  Output:  returns a time interval (sec)
//  Purpose: finds rainfall processing time interval for a unit hydrograph group.
//
{
    int ri;        // rainfal processing time interval for the UH group
    int tLimb;     // duration of a UH's rising & falling limbs
    int k, m;

    // --- begin with UH group time step equal to wet runoff step
    ri = sp->WetStep;

    // --- examine each UH in the group
    for (m=0; m<12; m++)
    {
        for (k=0; k<3; k++)
        {
            // --- make sure the UH exists
            if ( sp->UnitHyd[i].tPeak[m][k] > 0 )
            {
                // --- reduce time step if rising/falling limb is smaller
                tLimb = sp->UnitHyd[i].tPeak[m][k];
                ri = MIN(ri, tLimb);
                tLimb = sp->UnitHyd[i].tBase[m][k] - tLimb;
                if ( tLimb > 0 ) ri = MIN(ri, tLimb);
            }
        }
    }
    return ri;
}

//=============================================================================

int  getMaxPeriods(SWMM_Project *sp, int i, int k)
//
//  Input:   i = UH group index
//           k = UH index
//  Output:  returns number of past rainfall values
//  Purpose: finds number of past rainfall values to save for a UH.
//
{
    int   m,                           // month index
          n,                           // number of time periods
          nMax,                        // maximum number of time periods
          rainInterval;                // rainfall processing interval (sec)

    TRdiiShared *rd = &sp->RdiiShared;

    // --- examine each monthly set of UHs
    rainInterval = rd->UHGroup[i].rainInterval;
    nMax = 0;
    for (m=0; m<12; m++)
    {
        // --- compute number of time periods in UH base
        n = (sp->UnitHyd[i].tBase[m][k] / rainInterval) + 1;

        // --- update number of time periods to be saved
        nMax = MAX(n, nMax);
    }
    return nMax;
}

//=============================================================================

void initGageData(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes state of Unit Hydrograph rain gages.
//
{
    int i;                             // unit hyd. index
    int g;                             // rain gage index

    // --- first initialize the state of each rain gage
    for (g=0; g<sp->Nobjects[GAGE]; g++)
    {
        if ( sp->Gage[g].tSeries >= 0 )
        {
            table_tseriesInit(sp, &sp->Tseries[sp->Gage[g].tSeries]);
        }
        gage_initState(sp, g);
    }

    // --- then flag each gage that is used by a Unit Hydrograph set
    for (i=0; i<sp->Nobjects[UNITHYD]; i++)
    {
        g = sp->UnitHyd[i].rainGage;
        if ( g >= 0 )
        {
            sp->Gage[g].isUsed = TRUE;

            // --- if UH's gage uses same time series as a previous gage,
            //     then assign the latter gage to the UH
            if ( sp->Gage[g].coGage >= 0 )
            {
                sp->UnitHyd[i].rainGage = sp->Gage[g].coGage;
                sp->Gage[sp->Gage[g].coGage].isUsed = TRUE;
            }
        }
    }
}

//=============================================================================

void initUnitHydData(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes unit hydrograph processing data.
//
{
    int i,                             // UH group index
        j,                             // node index
        k,                             // UH index
        n;                             // RDII node index
//  int g,                             // rain gage index
    int month;                         // month index

    TRdiiShared *rd = &sp->RdiiShared;

    // --- initialize UHGroup entries for each Unit Hydrograph
    month = datetime_monthOfYear(sp->StartDateTime) - 1;
    for (i=0; i<sp->Nobjects[UNITHYD]; i++)
    {
        for (k=0; k<3; k++)
        {
            // --- make the first recorded rainfall begin a new RDII event
            // --- (new RDII event occurs when dry period > base of longest UH)
            rd->UHGroup[i].uh[k].drySeconds =
                (rd->UHGroup[i].uh[k].maxPeriods * rd->UHGroup[i].rainInterval) + 1;
            rd->UHGroup[i].uh[k].period = rd->UHGroup[i].uh[k].maxPeriods + 1;
            rd->UHGroup[i].uh[k].hasPastRain = FALSE;

            // --- assign initial abstraction used
            rd->UHGroup[i].uh[k].iaUsed = sp->UnitHyd[i].iaInit[month][k];
        }

        // --- initialize gage date to simulation start date
        rd->UHGroup[i].gageDate = sp->StartDateTime;
        rd->UHGroup[i].area = 0.0;
        rd->UHGroup[i].rdii = 0.0;
    }

    // --- assume each UH group is not used
    for (i=0; i<sp->Nobjects[UNITHYD]; i++) rd->UHGroup[i].isUsed = FALSE;

    // --- look at each node with RDII inflow
    for (n = 0; n < rd->NumRdiiNodes; n++)
    {
        // --- mark as used the UH group associated with the node
        j = rd->RdiiNodeIndex[n];
        i = sp->Node[j].rdiiInflow->unitHyd;
        rd->UHGroup[i].isUsed = TRUE;

        // --- add node's sewer area to UH group's area
        rd->UHGroup[i].lastDate = sp->StartDateTime;
        rd->UHGroup[i].area += sp->Node[j].rdiiInflow->area;
    }
}

//=============================================================================

int openNewRdiiFile(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns TRUE if successful, FALSE if not
//  Purpose: opens a new RDII interface file.
//
{
    int j;                             // node index

    TRdiiShared *rd = &sp->RdiiShared;

    // --- create a temporary file name if scratch file being used
    if ( sp->Frdii.mode == SCRATCH_FILE ) getTempFileName(sp, sp->Frdii.name);

    // --- open the RDII file as a formatted text file
    sp->Frdii.file = fopen(sp->Frdii.name, "w+b");
    if ( sp->Frdii.file == NULL )
    {
        return FALSE;
    }

    // --- write file stamp to RDII file
    fwrite(FileStamp, sizeof(char), strlen(FileStamp), sp->Frdii.file);

    // --- initialize the contents of the file with RDII time step (sec),
    //     number of RDII nodes, and index of each node
    fwrite(&rd->RdiiStep, sizeof(INT4), 1, sp->Frdii.file);
    fwrite(&rd->NumRdiiNodes, sizeof(INT4), 1, sp->Frdii.file);
    for (j=0; j<sp->Nobjects[NODE]; j++)
    {
        if ( sp->Node[j].rdiiInflow ) fwrite(&j, sizeof(INT4), 1, sp->Frdii.file);
    }
    return TRUE;
}

//=============================================================================

void getRainfall(SWMM_Project *sp, DateTime currentDate)
//
//  Input:   currentDate = current calendar date/time
//  Output:  none
//  Purpose: determines rainfall at current RDII processing date.
//
//
{
    int      j;                        // UH group index
    int      k;                        // UH index
    int      g;                        // rain gage index
    int      i;                        // past rainfall index
    int      month;                    // month of current date
    int      rainInterval;             // rainfall interval (sec)
    double   rainDepth;                // rainfall depth (inches or mm)
    double   excessDepth;              // excess rainfall depth (inches or mm))
    DateTime gageDate;                 // calendar date for rain gage

    TRdiiShared *rd = &sp->RdiiShared;

    // --- examine each UH group
    month = datetime_monthOfYear(currentDate) - 1;
    for (g = 0; g < sp->Nobjects[GAGE]; g++) sp->Gage[g].isCurrent = FALSE;
    for (j = 0; j < sp->Nobjects[UNITHYD]; j++)
    {
        // --- repeat until gage's date reaches or exceeds current date
        g = sp->UnitHyd[j].rainGage;
        rainInterval = rd->UHGroup[j].rainInterval;
        while ( rd->UHGroup[j].gageDate < currentDate )
        {
            // --- get rainfall volume over gage's recording interval
            //     at gage'a current date (in original depth units)
            gageDate = rd->UHGroup[j].gageDate;
            sp->Adjust.rainFactor = sp->Adjust.rain[datetime_monthOfYear(gageDate)-1]; //(5.1.007)
            if (!sp->Gage[g].isCurrent)
            {
                gage_setState(sp, g, gageDate);
                sp->Gage[g].isCurrent = TRUE;
            }
            rainDepth = sp->Gage[g].rainfall * (double)rainInterval / 3600.0;

            // --- update amount of total rainfall volume (ft3)
            rd->TotalRainVol += rainDepth / UCF(sp, RAINDEPTH) * rd->UHGroup[j].area;

            // --- compute rainfall excess for each UH in the group
            for (k=0; k<3; k++)
            {
                // --- adjust rainfall volume for any initial abstraction
                excessDepth = applyIA(sp, j, k, gageDate, rainInterval, rainDepth);

                // --- adjust extent of dry period for the UH
                updateDryPeriod(sp, j, k, excessDepth, rainInterval);

                // --- add rainfall to list of past values,
                //     wrapping array index if necessary
                i = rd->UHGroup[j].uh[k].period;
                if ( i >= rd->UHGroup[j].uh[k].maxPeriods ) i = 0;
                rd->UHGroup[j].uh[k].pastRain[i] = excessDepth;
                rd->UHGroup[j].uh[k].pastMonth[i] = (char)month;
                rd->UHGroup[j].uh[k].period = i + 1;
            }

            // --- advance rain date by gage recording interval
            rd->UHGroup[j].gageDate = datetime_addSeconds(gageDate, rainInterval);
        }
    }
}

//=============================================================================

double  applyIA(SWMM_Project *sp, int j, int k, DateTime aDate, double dt,
        double rainDepth)
//
//  Input:   j = UH group index
//           k = unit hydrograph index
//           aDate = current date/time
//           dt = time interval (sec)
//           rainDepth = unadjusted rain depth (in or mm)
//  Output:  returns rainfall adjusted for initial abstraction (IA)
//  Purpose: adjusts rainfall for any initial abstraction and updates the
//           amount of available initial abstraction actually used.
//
{
    int m;
    double ia, netRainDepth;

    TRdiiShared *rd = &sp->RdiiShared;

    // --- determine amount of unused IA
    m = datetime_monthOfYear(aDate) - 1;
    ia = sp->UnitHyd[j].iaMax[m][k] - rd->UHGroup[j].uh[k].iaUsed;
    ia = MAX(ia, 0.0);

    // --- case where there's some rainfall
    if ( rainDepth > 0.0 )
    {
        // --- reduce rain depth by unused IA
        netRainDepth = rainDepth - ia;
        netRainDepth = MAX(netRainDepth, 0.0);

        // --- update amount of IA used up
        ia = rainDepth - netRainDepth;
        rd->UHGroup[j].uh[k].iaUsed += ia;
    }

    // --- case where there's no rainfall
    else
    {
        // --- recover a portion of the IA already used
        rd->UHGroup[j].uh[k].iaUsed -= dt / 86400. * sp->UnitHyd[j].iaRecov[m][k];
        rd->UHGroup[j].uh[k].iaUsed = MAX(rd->UHGroup[j].uh[k].iaUsed, 0.0);
        netRainDepth = 0.0;
    }
    return netRainDepth;
}

//=============================================================================

void updateDryPeriod(SWMM_Project *sp, int j, int k, double rainDepth, int rainInterval)
//
//  Input:   j = UH group index
//           k = unit hydrograph index
//           rainDepth = excess rain depth (in or mm)
//           rainInterval = rainfall time interval (sec)
//  Output:  none
//  Purpose: adjusts the length of the dry period between rainfall events.
//
{
    int i;

    TRdiiShared *rd = &sp->RdiiShared;

    // --- if rainfall occurs
    if ( rainDepth > 0.0 )
    {
        // --- if previous dry period long enough then begin
        //     new RDII event with time period index set to 0
        if ( rd->UHGroup[j].uh[k].drySeconds >= rainInterval *
                rd->UHGroup[j].uh[k].maxPeriods )
        {
            for (i = 0; i < rd->UHGroup[j].uh[k].maxPeriods; i++)
            {
                rd->UHGroup[j].uh[k].pastRain[i] = 0.0;
            }
            rd->UHGroup[j].uh[k].period = 0;
        }
        rd->UHGroup[j].uh[k].drySeconds = 0;
        rd->UHGroup[j].uh[k].hasPastRain = TRUE;
    }

    // --- if no rainfall, update duration of dry period
    else
    {
        rd->UHGroup[j].uh[k].drySeconds += rainInterval;
        if ( rd->UHGroup[j].uh[k].drySeconds >=
            rainInterval * rd->UHGroup[j].uh[k].maxPeriods )
        {
            rd->UHGroup[j].uh[k].hasPastRain = FALSE;
        }
        else rd->UHGroup[j].uh[k].hasPastRain = TRUE;
    }
}

//=============================================================================

void getUnitHydRdii(SWMM_Project *sp, DateTime currentDate)
//
//  Input:   currentDate = current calendar date/time
//  Output:  none
//  Purpose: computes RDII generated by past rainfall for each UH group.
//
{
    int   j;                           // UH group index
    int   k;                           // UH index
    int   rainInterval;                // rainfall time interval (sec)

    TRdiiShared *rd = &sp->RdiiShared;

    // --- examine each UH group
    for (j=0; j<sp->Nobjects[UNITHYD]; j++)
    {
        // --- skip calculation if group not used by any RDII node or if
        //     current date hasn't reached last date RDII was computed
        if ( !rd->UHGroup[j].isUsed ) continue;
        if ( currentDate < rd->UHGroup[j].lastDate ) continue;

        // --- update date RDII last computed
        rd->UHGroup[j].lastDate = rd->UHGroup[j].gageDate;

        // --- perform convolution for each UH in the group
        rainInterval = rd->UHGroup[j].rainInterval;
        rd->UHGroup[j].rdii = 0.0;
        for (k=0; k<3; k++)
        {
            if ( rd->UHGroup[j].uh[k].hasPastRain )
            {
                rd->UHGroup[j].rdii += getUnitHydConvol(sp, j, k, rainInterval);
            }
        }
    }
}

//=============================================================================

double getUnitHydConvol(SWMM_Project *sp, int j, int k, int rainInterval)
//
//  Input:   j = UH group index
//           k = UH index
//           rainInterval = rainfall time interval (sec)
//  Output:  returns a RDII flow value
//  Purpose: computes convolution of Unit Hydrographs with past rainfall.
//
{
    int    i;                          // previous rainfall period index
    int    m;                          // month of year index
    int    p;                          // UH time period index
    int    pMax;                       // max. number of periods
    double t;                          // UH time value (sec)
    double u;                          // UH ordinate
    double v;                          // rainfall volume
    double rdii;                       // RDII flow
    TUHData* uh;                       // UH data

    TRdiiShared *rd = &sp->RdiiShared;

    // --- initialize RDII, rain period index and UH period index
    rdii = 0.0;
    uh = &rd->UHGroup[j].uh[k];
    i = uh->period - 1;
    if ( i < 0 ) i = uh->maxPeriods - 1;
    pMax = uh->maxPeriods;
    p = 1;

    // --- evaluate each time period of UH's
    while ( p < pMax )
    {
        // --- if rain period has rainfall
        v = uh->pastRain[i];
        m = uh->pastMonth[i];
        if ( v > 0.0 )
        {
            // --- find mid-point time of UH period in seconds
            t = ((double)(p) - 0.5) * (double)rainInterval;

            // --- convolute rain volume with UH ordinate
            u = getUnitHydOrd(sp, j, m, k, t) * sp->UnitHyd[j].r[m][k];
            rdii += u * v;
        }

        // --- move to next UH period & previous rainfall period
        p = p + 1;
        i = i - 1;
        if ( i < 0 ) i = uh->maxPeriods - 1;
    }
    return rdii;
}

//=============================================================================

double getUnitHydOrd(SWMM_Project *sp, int h, int m, int k, double t)
//
//  Input:   h = index of UH group
//           m = month index
//           k = individual UH index
//           t = UH time (sec)
//  Output:  returns ordinate of a unit hydrograph
//  Purpose: gets ordinate of a particular unit hydrograph at specified time.
//
{
    double qPeak;                      // peak flow of unit hydrograph
    double f;                          // fraction of time to/from peak on UH
    double t1;                         // time to peak on UH (sec)
    double t2;                         // time after peak on UH (sec)
    double tBase;                      // base time of UH (sec)

    // --- return 0 if past end of UH time base
    tBase = sp->UnitHyd[h].tBase[m][k];
    if ( t >= tBase ) return 0.0;

    // --- compute peak value of UH in original rainfall units (in/hr or mm/hr)
    qPeak = 2. / tBase * 3600.0;

    // --- break UH base into times before & after peak flow
    t1 = sp->UnitHyd[h].tPeak[m][k];
    t2 = tBase - t1;

    // --- find UH flow at time t
    if ( t <= t1 ) f = t / t1;
    else           f = 1.0 - (t - t1) / t2;
    return MAX(f, 0.0) * qPeak;
}

//=============================================================================

int getNodeRdii(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns TRUE if any node has RDII inflow, FALSE if not
//  Purpose: computes current RDII inflow at each node.
//
{
    int   hasRdii = FALSE;             // true if any node has some RDII
    int   i;                           // UH group index
    int   j;                           // node index
    int   n;                           // number of nodes w/ RDII
    double rdii;                       // RDII flow (cfs)

    TRdiiShared *rd = &sp->RdiiShared;

    // --- examine each node w/ RDII data
    for (n = 0; n < rd->NumRdiiNodes; n++)
    {
        // --- identify node's index in project's data base
        j = rd->RdiiNodeIndex[n];

        // --- apply node's sewer area to UH RDII to get node RDII in CFS
        i = sp->Node[j].rdiiInflow->unitHyd;
        rdii = rd->UHGroup[i].rdii * sp->Node[j].rdiiInflow->area / UCF(sp, RAINFALL);
        if ( rdii < ZERO_RDII ) rdii = 0.0;
        else hasRdii = TRUE;

        // --- update total RDII volume
        rd->RdiiNodeFlow[n] = (REAL4)rdii;
        if ( rdii > 0.0 )
        {
            rd->TotalRdiiVol += rdii * (double)rd->RdiiStep;
        }
    }
    return hasRdii;
}

//=============================================================================

void saveRdiiFlows(SWMM_Project *sp, DateTime currentDate)
//
//  Input:   currentDate = current calendar date/time
//  Output:  none
//  Purpose: saves current set of RDII inflows in current flow units to file.
//
{
    TRdiiShared *rd = &sp->RdiiShared;

    fwrite(&currentDate, sizeof(DateTime), 1, sp->Frdii.file);
    fwrite(rd->RdiiNodeFlow, sizeof(REAL4), rd->NumRdiiNodes, sp->Frdii.file);             //(5.1.003)
}

//=============================================================================

void  closeRdiiProcessor(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: closes RDII processing system.
//
{
    TRdiiShared *rd = &sp->RdiiShared;

    // --- write rainfall & RDII totals to report file
    if ( !sp->ErrorCode )
    {
        report_writeRdiiStats(sp, rd->TotalRainVol, rd->TotalRdiiVol);
    }

    // --- free allocated memory and close RDII file
    freeRdiiMemory(sp);
    if ( sp->Frdii.file ) fclose(sp->Frdii.file);
}

//=============================================================================

void freeRdiiMemory(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: frees memory used for RDII processing.
//
{
    int i;
    int k;

    TRdiiShared *rd = &sp->RdiiShared;

    if ( rd->UHGroup )
    {
        for (i = 0; i < sp->Nobjects[UNITHYD]; i++)
        {
            for (k=0; k<3; k++)
            {
                FREE(rd->UHGroup[i].uh[k].pastRain);
                FREE(rd->UHGroup[i].uh[k].pastMonth);
            }
        }
        FREE(rd->UHGroup);
    }
    FREE(rd->RdiiNodeIndex);
    FREE(rd->RdiiNodeFlow);
}
