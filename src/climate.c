//-----------------------------------------------------------------------------
//   climate.c
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    03/20/10 (Build 5.1.001)
//            09/15/14 (Build 5.1.007)
//            03/19/15 (Build 5.1.008)
//            08/05/15 (Build 5.1.010)
//            08/01/16 (Build 5.1.011)
//   Author:  L. Rossman
//
//   Climate related functions.
//
//   Build 5.1.007:
//   - NCDC GHCN climate file format added.
//   - Monthly adjustments for temperature, evaporation & rainfall added.
//
//   Build 5.1.008:
//   - Monthly adjustments for hyd. conductivity added.
//   - Time series evaporation rates can now vary within a day.
//   - Evaporation rates are now properly updated when only flow routing
//     is being simulated.
//
//   Build 5.1.010:
//   - Hargreaves evaporation now computed using 7-day average temperatures.
//             
//   Build 5.1.011:
//   - Monthly adjustment for hyd. conductivity <= 0 is ignored.
///-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "headers.h"
#include "climate.h"

////-----------------------------------------------------------------------------
////  Constants
////-----------------------------------------------------------------------------
enum ClimateFileFormats {UNKNOWN_FORMAT,
                         USER_PREPARED,     // SWMM 5's own user format
                         GHCND,             // NCDC GHCN Daily format          //(5.1.007)
                         TD3200,            // NCDC TD3200 format
                         DLY0204};          // Canadian DLY02 or DLY04 format

static const int    MAXCLIMATEVARS  = 4;
static const int    MAXDAYSPERMONTH = 32;

// These variables are used when processing climate files.
enum   ClimateVarType {TMIN, TMAX, EVAP, WIND};
enum   WindSpeedType  {WDMV, AWND};                                            //(5.1.007)

static const char* ClimateVarWords[] = {"TMIN", "TMAX", "EVAP", "WDMV", "AWND",      //(5.1.007)
                                  NULL};

//-----------------------------------------------------------------------------
//  External functions (defined in funcs.h)
//-----------------------------------------------------------------------------
//  climate_readParams                 // called by input_parseLine
//  climate_readEvapParams             // called by input_parseLine
//  climate_validate                   // called by project_validate
//  climate_openFile                   // called by runoff_open
//  climate_initState                  // called by project_init
//  climate_setState                   // called by runoff_execute
//  climate_getNextEvapDate            // called by runoff_getTimeStep         //(5.1.008)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int  getFileFormat(SWMM_Project *sp);
static void readFileLine(SWMM_Project *sp, int *year, int *month);
static void readUserFileLine(SWMM_Project *sp, int *year, int *month);
static void readTD3200FileLine(SWMM_Project *sp, int *year, int *month);
static void readDLY0204FileLine(SWMM_Project *sp, int *year, int *month);
static void readFileValues(SWMM_Project *sp);

static void setNextEvapDate(SWMM_Project *sp, DateTime thedate);                                 //(5.1.008)
static void setEvap(SWMM_Project *sp, DateTime theDate);
static void setTemp(SWMM_Project *sp, DateTime theDate);
static void setWind(SWMM_Project *sp, DateTime theDate);
static void updateTempTimes(SWMM_Project *sp, int day);
static void updateTempMoveAve(SWMM_Project *sp, double tmin, double tmax);                       //(5.1.010)
static double getTempEvap(SWMM_Project *sp, int day, double ta, double tr);                      //(5.1.010)

static void updateFileValues(SWMM_Project *sp, DateTime theDate);
static void parseUserFileLine(SWMM_Project *sp);
static void parseTD3200FileLine(SWMM_Project *sp);
static void parseDLY0204FileLine(SWMM_Project *sp);
static void setTD3200FileValues(SWMM_Project *sp, int param);

static int  isGhcndFormat(SWMM_Project *sp, char* line);                                         //(5.1.007)
static void readGhcndFileLine(SWMM_Project *sp, int *year, int *month);                          //(5.1.007)
static void parseGhcndFileLine(SWMM_Project *sp);                                          //(5.1.007)

//=============================================================================

int  climate_readParams(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads climate/temperature parameters from input line of data
//
//  Format of data can be
//    TIMESERIES  name
//    FILE        name
//    WINDSPEED   MONTHLY  v1  v2  ...  v12
//    WINDSPEED   FILE
//    SNOWMELT    v1  v2  ...  v6
//    ADC         IMPERV/PERV  v1  v2  ...  v10
//
{
    int      i, j, k;
    double   x[6], y;
    DateTime aDate;

    // --- identify keyword
    k = findmatch(tok[0], TempKeyWords);
    if ( k < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[0]);
    switch (k)
    {
      case 0: // Time series name
        // --- check that time series name exists
        if ( ntoks < 2 ) return error_setInpError(sp, ERR_ITEMS, "");
        i = project_findObject(sp, TSERIES, tok[1]);
        if ( i < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);

        // --- record the time series as being the data source for temperature
        sp->Temp.dataSource = TSERIES_TEMP;
        sp->Temp.tSeries = i;
        sp->Tseries[i].refersTo = TSERIES_TEMP;
        break;

      case 1: // Climate file
        // --- record file as being source of temperature data
        if ( ntoks < 2 ) return error_setInpError(sp, ERR_ITEMS, "");
        sp->Temp.dataSource = FILE_TEMP;

        // --- save name and usage mode of external climate file
        sp->Fclimate.mode = USE_FILE;
        sstrncpy(sp->Fclimate.name, tok[1], MAXFNAME);

        // --- save starting date to read from file if one is provided
        sp->Temp.fileStartDate = NO_DATE;
        if ( ntoks > 2 )
        {
            if ( *tok[2] != '*')
            {
                if ( !datetime_strToDate(sp, tok[2], &aDate) )
                    return error_setInpError(sp, ERR_DATETIME, tok[2]);
                sp->Temp.fileStartDate = aDate;
            }
        }
        break;

      case 2: // Wind speeds
        // --- check if wind speeds will be supplied from climate file
        if ( strcomp(tok[1], w_FILE) )
        {
            sp->Wind.type = FILE_WIND;
        }

        // --- otherwise read 12 monthly avg. wind speed values
        else
        {
            if ( ntoks < 14 ) return error_setInpError(sp, ERR_ITEMS, "");
            sp->Wind.type = MONTHLY_WIND;
            for (i=0; i<12; i++)
            {
                if ( !getDouble(tok[i+2], &y) )
                    return error_setInpError(sp, ERR_NUMBER, tok[i+2]);
                sp->Wind.aws[i] = y;
            }
        }
        break;

      case 3: // Snowmelt params
        if ( ntoks < 7 ) return error_setInpError(sp, ERR_ITEMS, "");
        for (i=1; i<7; i++)
        {
            if ( !getDouble(tok[i], &x[i-1]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i]);
        }
        // --- convert deg. C to deg. F for snowfall temperature
        if ( sp->UnitSystem == SI ) x[0] = 9./5.*x[0] + 32.0;
        sp->Snow.snotmp = x[0];
        sp->Snow.tipm   = x[1];
        sp->Snow.rnm    = x[2];
        sp->Temp.elev   = x[3] / UCF(sp, LENGTH);
        sp->Temp.anglat = x[4];
        sp->Temp.dtlong = x[5] / 60.0;
        break;

      case 4:  // Areal Depletion Curve data
        // --- check if data is for impervious or pervious areas
        if ( ntoks < 12 ) return error_setInpError(sp, ERR_ITEMS, "");
        if      ( match(tok[1], w_IMPERV) ) i = 0;
        else if ( match(tok[1], w_PERV)   ) i = 1;
        else return error_setInpError(sp, ERR_KEYWORD, tok[1]);

        // --- read 10 fractional values
        for (j=0; j<10; j++)
        {
            if ( !getDouble(tok[j+2], &y) || y < 0.0 || y > 1.0 )
                return error_setInpError(sp, ERR_NUMBER, tok[j+2]);
            sp->Snow.adc[i][j] = y;
        }
        break;
    }
    return 0;
}

//=============================================================================

int climate_readEvapParams(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads evaporation parameters from input line of data.
//
//  Data formats are:
//    CONSTANT  value
//    MONTHLY   v1 ... v12
//    TIMESERIES name
//    TEMPERATURE
//    FILE      (v1 ... v12)
//    RECOVERY   name
//    DRY_ONLY   YES/NO
//
{
    int i, k;
    double x;

    // --- find keyword indicating what form the evaporation data is in
    k = findmatch(tok[0], EvapTypeWords);
    if ( k < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[0]);

    // --- check for RECOVERY pattern data
    if ( k == RECOVERY )
    {
        if ( ntoks < 2 ) return error_setInpError(sp, ERR_ITEMS, "");
        i = project_findObject(sp, TIMEPATTERN, tok[1]);
        if ( i < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
        sp->Evap.recoveryPattern = i;
        return 0;
    }

    // --- check for no evaporation in wet periods
    if ( k == DRYONLY )
    {
        if ( ntoks < 2 ) return error_setInpError(sp, ERR_ITEMS, "");
        if      ( strcomp(tok[1], w_NO ) )  sp->Evap.dryOnly = FALSE;
        else if ( strcomp(tok[1], w_YES ) ) sp->Evap.dryOnly = TRUE;
        else return error_setInpError(sp, ERR_KEYWORD, tok[1]);
        return 0;
    }

    // --- process data depending on its form
    sp->Evap.type = k;
    if ( k != TEMPERATURE_EVAP && ntoks < 2 )
        return error_setInpError(sp, ERR_ITEMS, "");
    switch ( k )
    {
      case CONSTANT_EVAP:
        // --- for constant evap., fill monthly avg. values with same number
        if ( !getDouble(tok[1], &x) )
            return error_setInpError(sp, ERR_NUMBER, tok[1]);
        for (i=0; i<12; i++) sp->Evap.monthlyEvap[i] = x;
        break;

      case MONTHLY_EVAP:
        // --- for monthly evap., read a value for each month of year
        if ( ntoks < 13 ) return error_setInpError(sp, ERR_ITEMS, "");
        for ( i=0; i<12; i++)
            if ( !getDouble(tok[i+1], &sp->Evap.monthlyEvap[i]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i+1]);
        break;

      case TIMESERIES_EVAP:
        // --- for time series evap., read name of time series
        i = project_findObject(sp, TSERIES, tok[1]);
        if ( i < 0 ) return error_setInpError(sp, ERR_NAME, tok[1]);
        sp->Evap.tSeries = i;
        sp->Tseries[i].refersTo = TIMESERIES_EVAP;
        break;

      case FILE_EVAP:
        // --- for evap. from climate file, read monthly pan coeffs.
        //     if they are provided (default values are 1.0)
        if ( ntoks > 1 )
        {
            if ( ntoks < 13 ) return error_setInpError(sp, ERR_ITEMS, "");
            for (i=0; i<12; i++)
            {
                if ( !getDouble(tok[i+1], &sp->Evap.panCoeff[i]) )
                    return error_setInpError(sp, ERR_NUMBER, tok[i+1]);
            }
        }
        break;
    }
    return 0;
}

//=============================================================================

////  New function added to release 5.1.007.  ////                             //(5.1.007)

int climate_readAdjustments(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads adjustments to monthly evaporation or rainfall
//           from input line of data.
//
//  Data formats are:
//    TEMPERATURE   v1 ... v12
//    EVAPORATION   v1 ... v12
//    RAINFALL      v1 ... v12
//    CONDUCTIVITY  v1 ... v12                                                 //(5.1.008)
{
    int i;
    if (ntoks == 1) return 0;

    if ( match(tok[0], "TEMP") )
    {
        if ( ntoks < 13 )  return error_setInpError(sp, ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &sp->Adjust.temp[i-1]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i]);
        }
        return 0;
    }

    if ( match(tok[0], "EVAP") )
    {
        if ( ntoks < 13 )  return error_setInpError(sp, ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &sp->Adjust.evap[i-1]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i]);
        }
        return 0;
    }

    if ( match(tok[0], "RAIN") )
    {
        if ( ntoks < 13 )  return error_setInpError(sp, ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &sp->Adjust.rain[i-1]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i]);
        }
        return 0;
    }

////  Following code segment added for release 5.1.008.  ////                  //(5.1.008)
////
    if ( match(tok[0], "CONDUCT") )
    {
        if ( ntoks < 13 )  return error_setInpError(sp, ERR_ITEMS, "");
        for (i = 1; i < 13; i++)
        {
            if ( !getDouble(tok[i], &sp->Adjust.hydcon[i-1]) )
                return error_setInpError(sp, ERR_NUMBER, tok[i]);
            if ( sp->Adjust.hydcon[i-1] <= 0.0 ) sp->Adjust.hydcon[i-1] = 1.0;         //(5.1.011)
        }
        return 0;
    }
////
    return error_setInpError(sp, ERR_KEYWORD, tok[0]);
}

//=============================================================================

void climate_validate(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: validates climatological variables
//
{
    int       i;                                                               //(5.1.007)
    double    a, z, pa;

    // --- check if climate data comes from external data file                 //(5.1.007)
    if ( sp->Wind.type == FILE_WIND || sp->Evap.type == FILE_EVAP ||
         sp->Evap.type == TEMPERATURE_EVAP )
    {
        if ( sp->Fclimate.mode == NO_FILE )
        {
            report_writeErrorMsg(sp, ERR_NO_CLIMATE_FILE, "");
        }
    }

    // --- open the climate data file                                          //(5.1.007)
    if ( sp->Fclimate.mode == USE_FILE ) climate_openFile(sp);                       //(5.1.007)

    // --- snow melt parameters tipm & rnm must be fractions
    if ( sp->Snow.tipm < 0.0 ||
         sp->Snow.tipm > 1.0 ||
         sp->Snow.rnm  < 0.0 ||
         sp->Snow.rnm  > 1.0 ) report_writeErrorMsg(sp, ERR_SNOWMELT_PARAMS, "");

    // --- latitude should be between -90 & 90 degrees
    a = sp->Temp.anglat;
    if ( a <= -89.99 ||
         a >= 89.99  ) report_writeErrorMsg(sp, ERR_SNOWMELT_PARAMS, "");
    else sp->Temp.tanAnglat = tan(a * PI / 180.0);

    // --- compute psychrometric constant
    z = sp->Temp.elev / 1000.0;
    if ( z <= 0.0 ) pa = 29.9;
    else  pa = 29.9 - 1.02*z + 0.0032*pow(z, 2.4); // atmos. pressure
    sp->Temp.gamma = 0.000359 * pa;

    // --- convert units of monthly temperature & evap adjustments             //(5.1.007)
    for (i = 0; i < 12; i++)
    {
        if (sp->UnitSystem == SI) sp->Adjust.temp[i] *= 9.0/5.0;
        sp->Adjust.evap[i] /= UCF(sp, EVAPRATE);
    }
}

//=============================================================================

void climate_openFile(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: opens a climate file and reads in first set of values.
//
{
    int i, m, y;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- open the file
    if ( (sp->Fclimate.file = fopen(sp->Fclimate.name, "rt")) == NULL )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_FILE_OPEN, sp->Fclimate.name);
        return;
    }

    // --- initialize values of file's climate variables
    //     (sp->Temp.ta was previously initialized in project.c)
    clmt->FileValue[TMIN] = sp->Temp.ta;
    clmt->FileValue[TMAX] = sp->Temp.ta;
    clmt->FileValue[EVAP] = 0.0;
    clmt->FileValue[WIND] = 0.0;

    // --- find climate file's format
    clmt->FileFormat = getFileFormat(sp);
    if ( clmt->FileFormat == UNKNOWN_FORMAT )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_FILE_READ, sp->Fclimate.name);
        return;
    }

    // --- position file to begin reading climate file at either user-specified
    //     month/year or at start of simulation period.
    rewind(sp->Fclimate.file);
    strcpy(clmt->FileLine, "");
    if ( sp->Temp.fileStartDate == NO_DATE )
        datetime_decodeDate(sp->StartDate, &clmt->FileYear, &clmt->FileMonth, &clmt->FileDay);
    else
        datetime_decodeDate(sp->Temp.fileStartDate, &clmt->FileYear, &clmt->FileMonth, &clmt->FileDay);
    while ( !feof(sp->Fclimate.file) )
    {
        strcpy(clmt->FileLine, "");
        readFileLine(sp, &y, &m);
        if ( y == clmt->FileYear && m == clmt->FileMonth ) break;
    }
    if ( feof(sp->Fclimate.file) )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_END_OF_FILE, sp->Fclimate.name);
        return;
    }

    // --- initialize file dates and current climate variable values
    if ( !sp->ErrorCode )
    {
        clmt->FileElapsedDays = 0;
        clmt->FileLastDay = datetime_daysPerMonth(clmt->FileYear, clmt->FileMonth);
        readFileValues(sp);
        for (i=TMIN; i<=WIND; i++)
        {
            if ( clmt->FileData[i][clmt->FileDay] == MISSING ) continue;
            clmt->FileValue[i] = clmt->FileData[i][clmt->FileDay];
        }
    }
}

//=============================================================================

////  This function was re-written for release 5.1.008.  ////                  //(5.1.008)

void climate_initState(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: initializes climate state variables.
//
{
    TClimateShared *clmt = &sp->ClimateShared;

    clmt->LastDay = NO_DATE;
    sp->Temp.tmax = MISSING;
    sp->Snow.removed = 0.0;
    clmt->NextEvapDate = sp->StartDate;
    clmt->NextEvapRate = 0.0;

    // --- initialize variables for time series evaporation
    if ( sp->Evap.type == TIMESERIES_EVAP && sp->Evap.tSeries >= 0  )
    {
        // --- initialize NextEvapDate & NextEvapRate to first entry of
        //     time series whose date <= the simulation start date
        table_getFirstEntry(sp, &sp->Tseries[sp->Evap.tSeries],
                            &clmt->NextEvapDate, &clmt->NextEvapRate);
        if ( clmt->NextEvapDate < sp->StartDate )
        {  
            setNextEvapDate(sp, sp->StartDate);
        }
        sp->Evap.rate = clmt->NextEvapRate / UCF(sp, EVAPRATE);

        // --- find the next time evaporation rates change after this
        setNextEvapDate(sp, clmt->NextEvapDate);
    }

////  Following section added to release 5.1.010.  ////                        //(5.1.010)
    // --- initialize variables for temperature evaporation
    if ( sp->Evap.type == TEMPERATURE_EVAP )
    {
        clmt->Tma.maxCount = sizeof(clmt->Tma.ta) / sizeof(double);
        clmt->Tma.count = 0;
        clmt->Tma.front = 0;
        clmt->Tma.tAve = 0.0;
        clmt->Tma.tRng = 0.0;
    }
////
}

//=============================================================================

void climate_setState(SWMM_Project *sp, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: sets climate variables for current date.
//
{
    if ( sp->Fclimate.mode == USE_FILE ) updateFileValues(sp, theDate);
    if ( sp->Temp.dataSource != NO_TEMP ) setTemp(sp, theDate);
    setEvap(sp, theDate);
    setWind(sp, theDate);
    sp->Adjust.rainFactor = sp->Adjust.rain[datetime_monthOfYear(theDate)-1];          //(5.1.007)
    sp->Adjust.hydconFactor = sp->Adjust.hydcon[datetime_monthOfYear(theDate)-1];      //(5.1.008)
    setNextEvapDate(sp, theDate);                                              //(5.1.008)
}

//=============================================================================

////  New function added to release 5.1.008.  ////                             //(5.1.008)

DateTime climate_getNextEvapDate(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns the current value of NextEvapDate
//  Purpose: gets the next date when evaporation rate changes.
//
{
    TClimateShared *clmt = &sp->ClimateShared;

    return clmt->NextEvapDate;
}

//=============================================================================

////  Modified from what was previously named climate_getNextEvap.  ////       //(5.1.008)

void setNextEvapDate(SWMM_Project *sp, DateTime theDate)
//
//  Input:   theDate = current simulation date
//  Output:  sets a new value for NextEvapDate
//  Purpose: finds date for next change in evaporation after the current date.
//
{
    int    yr, mon, day, k;
    double d, e;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- do nothing if current date hasn't reached the current next date
    if ( clmt->NextEvapDate > theDate ) return;

    switch ( sp->Evap.type )
    {
      // --- for constant evaporation, use a next date far in the future
      case CONSTANT_EVAP:
         clmt->NextEvapDate = theDate + 365.;
         break;

      // --- for monthly evaporation, use the start of the next month
      case MONTHLY_EVAP:
        datetime_decodeDate(theDate, &yr, &mon, &day);
        if ( mon == 12 )
        {
            mon = 1;
            yr++;
        }
        else mon++;
        clmt->NextEvapDate = datetime_encodeDate(yr, mon, 1);
        break;

      // --- for time series evaporation, find the next entry in the
      //     series on or after the current date
      case TIMESERIES_EVAP:
        k = sp->Evap.tSeries;
        if ( k >= 0 )
        {
            clmt->NextEvapDate = theDate + 365.;
            while ( table_getNextEntry(sp, &sp->Tseries[k], &d, &e) &&
                    d <= sp->EndDateTime )
            {
                if ( d >= theDate )
                {
                    clmt->NextEvapDate = d;
                    clmt->NextEvapRate = e;
                    break;
                }
            }
        }
        break;

      // --- for climate file daily evaporation, use the next day
      case FILE_EVAP:
        clmt->NextEvapDate = floor(theDate) + 1.0;
        break;

      default: clmt->NextEvapDate = theDate + 365.;
    }
}

//=============================================================================

void updateFileValues(SWMM_Project *sp, DateTime theDate)
//
//  Input:   theDate = current simulation date
//  Output:  none
//  Purpose: updates daily climate variables for new day or reads in
//           another month worth of values if a new month begins.
//
//  NOTE:    counters FileElapsedDays, FileDay, FileMonth, FileYear and
//           FileLastDay were initialized in climate_openFile().
//
{
    int i;
    int deltaDays;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- see if a new day has begun
    deltaDays = (int)(floor(theDate) - floor(sp->StartDateTime));
    if ( deltaDays > clmt->FileElapsedDays )
    {
        // --- advance day counters
        clmt->FileElapsedDays++;
        clmt->FileDay++;

        // --- see if new month of data needs to be read from file
        if ( clmt->FileDay > clmt->FileLastDay )
        {
            clmt->FileMonth++;
            if ( clmt->FileMonth > 12 )
            {
                clmt->FileMonth = 1;
                clmt->FileYear++;
            }
            readFileValues(sp);
            clmt->FileDay = 1;
            clmt->FileLastDay = datetime_daysPerMonth(clmt->FileYear, clmt->FileMonth);
        }

        // --- set climate variables for new day
        for (i=TMIN; i<=WIND; i++)
        {
            // --- no change in current value if its missing
            if ( clmt->FileData[i][clmt->FileDay] == MISSING ) continue;
            clmt->FileValue[i] = clmt->FileData[i][clmt->FileDay];
        }
    }
}

//=============================================================================

void setTemp(SWMM_Project *sp, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: updates temperatures for new simulation date.
//
{
    int      j;                        // snow data object index
    int      k;                        // time series index
    int      mon;                      // month of year                        //(5.1.007)
    int      day;                      // day of year
    DateTime theDay;                   // calendar day
    double   hour;                     // hour of day
    double   tmp;                      // temporary temperature

    TClimateShared *clmt = &sp->ClimateShared;

    // --- see if a new day has started
    mon = datetime_monthOfYear(theDate);                                       //(5.1.007)
    theDay = floor(theDate);
    if ( theDay > clmt->LastDay )
    {
        // --- update min. & max. temps & their time of day
        day = datetime_dayOfYear(theDate);
        if ( sp->Temp.dataSource == FILE_TEMP )
        {
            clmt->Tmin = clmt->FileValue[TMIN] + sp->Adjust.temp[mon-1];                       //(5.1.007)
            clmt->Tmax = clmt->FileValue[TMAX] + sp->Adjust.temp[mon-1];                       //(5.1.007)
            if ( clmt->Tmin > clmt->Tmax )
            {
                tmp = clmt->Tmin;
                clmt->Tmin = clmt->Tmax;
                clmt->Tmax = tmp;
            }
            updateTempTimes(sp, day);
            if ( sp->Evap.type == TEMPERATURE_EVAP )
            {
                updateTempMoveAve(sp, clmt->Tmin, clmt->Tmax);                                 //(5.1.010)
                clmt->FileValue[EVAP] = getTempEvap(sp, day, clmt->Tma.tAve, clmt->Tma.tRng);        //(5.1.010)
            }
        }

        // --- compute snow melt coefficients based on day of year
        sp->Snow.season = sin(0.0172615*(day-81.0));
        for (j=0; j<sp->Nobjects[SNOWMELT]; j++)
        {
            snow_setMeltCoeffs(sp, j, sp->Snow.season);
        }

        // --- update date of last day analyzed
        clmt->LastDay = theDate;
    }

    // --- for min/max daily temps. from climate file,
    //     compute hourly temp. by sinusoidal interp.
    if ( sp->Temp.dataSource == FILE_TEMP )
    {
        hour = (theDate - theDay) * 24.0;
        if ( hour < clmt->Hrsr )
            sp->Temp.ta = clmt->Tmin + clmt->Trng1/2.0 * sin(PI/clmt->Dydif * (clmt->Hrsr - hour));
        else if ( hour >= clmt->Hrsr && hour <= clmt->Hrss )
            sp->Temp.ta = clmt->Tave + clmt->Trng * sin(PI/clmt->Dhrdy * (clmt->Hrday - hour));
        else
            sp->Temp.ta = clmt->Tmax - clmt->Trng * sin(PI/clmt->Dydif * (hour - clmt->Hrss));
    }

    // --- for user-supplied temperature time series,
    //     get temperature value from time series
    if ( sp->Temp.dataSource == TSERIES_TEMP )
    {
        k = sp->Temp.tSeries;
        if ( k >= 0)
        {
            sp->Temp.ta = table_tseriesLookup(sp, &sp->Tseries[k], theDate, TRUE);

            // --- convert from deg. C to deg. F if need be
            if ( sp->UnitSystem == SI )
            {
                sp->Temp.ta = (9./5.) * sp->Temp.ta + 32.0;
            }

            // --- apply climate change adjustment factor                      //(5.1.007)
            sp->Temp.ta += sp->Adjust.temp[mon-1];                                     //(5.1.007)
        }
    }

    // --- compute saturation vapor pressure
    sp->Temp.ea = 8.1175e6 * exp(-7701.544 / (sp->Temp.ta + 405.0265) );
}

//=============================================================================

void setEvap(SWMM_Project *sp, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: sets evaporation rate (ft/sec) for a specified date.
//
{
    int k;
    int mon = datetime_monthOfYear(theDate);                                   //(5.1.007)

    TClimateShared *clmt = &sp->ClimateShared;

    switch ( sp->Evap.type )
    {
      case CONSTANT_EVAP:
        sp->Evap.rate = sp->Evap.monthlyEvap[0] / UCF(sp, EVAPRATE);
        break;

      case MONTHLY_EVAP:
        sp->Evap.rate = sp->Evap.monthlyEvap[mon-1] / UCF(sp, EVAPRATE);
        break;

      case TIMESERIES_EVAP:
        if ( theDate >= clmt->NextEvapDate )
            sp->Evap.rate = clmt->NextEvapRate / UCF(sp, EVAPRATE);
        break;

      case FILE_EVAP:
        sp->Evap.rate = clmt->FileValue[EVAP] / UCF(sp, EVAPRATE);
        sp->Evap.rate *= sp->Evap.panCoeff[mon-1];
        break;

      case TEMPERATURE_EVAP:
        sp->Evap.rate = clmt->FileValue[EVAP] / UCF(sp, EVAPRATE);
        break;

      default: sp->Evap.rate = 0.0;
    }

    // --- apply climate change adjustment                                     //(5.1.007)
    sp->Evap.rate += sp->Adjust.evap[mon-1];                                           //(5.1.007)

    // --- set soil recovery factor
    sp->Evap.recoveryFactor = 1.0;
    k = sp->Evap.recoveryPattern;
    if ( k >= 0 && sp->Pattern[k].type == MONTHLY_PATTERN )
    {
        sp->Evap.recoveryFactor = sp->Pattern[k].factor[mon-1];                        //(5.1.007)
    }
}

//=============================================================================

void setWind(SWMM_Project *sp, DateTime theDate)
//
//  Input:   theDate = simulation date
//  Output:  none
//  Purpose: sets wind speed (mph) for a specified date.
//
{
    int yr, mon, day;

    TClimateShared *clmt = &sp->ClimateShared;

    switch ( sp->Wind.type )
    {
      case MONTHLY_WIND:
        datetime_decodeDate(theDate, &yr, &mon, &day);
        sp->Wind.ws = sp->Wind.aws[mon-1] / UCF(sp, WINDSPEED);
        break;

      case FILE_WIND:
        sp->Wind.ws = clmt->FileValue[WIND];
        break;

      default: sp->Wind.ws = 0.0;
    }
}

//=============================================================================

void updateTempTimes(SWMM_Project *sp, int day)
//
//  Input:   day = day of year
//  Output:  none
//  Purpose: computes time of day when min/max temperatures occur.
//           (min. temp occurs at sunrise, max. temp. at 3 hrs. < sunset)
//
{
    double decl;                       // earth's declination
    double hrang;                      // hour angle of sunrise/sunset
    double arg;

    TClimateShared *clmt = &sp->ClimateShared;

    decl  = 0.40928*cos(0.017202*(172.0-day));
    arg = -tan(decl)*sp->Temp.tanAnglat;
    if      ( arg <= -1.0 ) arg = PI;
    else if ( arg >= 1.0 )  arg = 0.0;
    else                    arg = acos(arg);
    hrang = 3.8197 * arg;
    clmt->Hrsr  = 12.0 - hrang + sp->Temp.dtlong;
    clmt->Hrss  = 12.0 + hrang + sp->Temp.dtlong - 3.0;
    clmt->Dhrdy = clmt->Hrsr - clmt->Hrss;
    clmt->Dydif = 24.0 + clmt->Hrsr - clmt->Hrss;
    clmt->Hrday = (clmt->Hrsr + clmt->Hrss) / 2.0;
    clmt->Tave  = (clmt->Tmin + clmt->Tmax) / 2.0;
    clmt->Trng  = (clmt->Tmax - clmt->Tmin) / 2.0;

    if ( sp->Temp.tmax == MISSING )
        clmt->Trng1 = clmt->Tmax - clmt->Tmin;
    else
        clmt->Trng1 = sp->Temp.tmax - clmt->Tmin;

    sp->Temp.tmax = clmt->Tmax;
}

//=============================================================================

////  This function was modified for release 5.1.010.  ////                    //(5.1.010)

double getTempEvap(SWMM_Project *sp, int day, double tave, double trng)
//
//  Input:   day = day of year
//           tave = 7-day average temperature (deg F)
//           trng = 7-day average daily temperature range (deg F)
//  Output:  returns evaporation rate in user's units (US:in/day, SI:mm/day)
//  Purpose: uses Hargreaves method to compute daily evaporation rate
//           from daily average temperatures and Julian day.
//
{
    double a = 2.0*PI/365.0;
    double ta = (tave - 32.0)*5.0/9.0;           //average temperature (deg C)
    double tr = trng*5.0/9.0;                    //temperature range (deg C)
    double lamda = 2.50 - 0.002361 * ta;         //latent heat of vaporization
    double dr = 1.0 + 0.033*cos(a*day);          //relative earth-sun distance
    double phi = sp->Temp.anglat*2.0*PI/360.0;       //latitude angle (rad)
    double del = 0.4093*sin(a*(284+day));        //solar declination angle (rad)
    double omega = acos(-tan(phi)*tan(del));     //sunset hour angle (rad)
    double ra = 37.6*dr*                         //extraterrestrial radiation
                (omega*sin(phi)*sin(del) +
                 cos(phi)*cos(del)*sin(omega));
    double e = 0.0023*ra/lamda*sqrt(tr)*(ta+17.8);    //evap. rate (mm/day)
    if ( e < 0.0 ) e = 0.0;
    if ( sp->UnitSystem == US ) e /= MMperINCH;           //evap rate (in/day)
    return e;
}

//=============================================================================

int  getFileFormat(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns code number of climate file's format
//  Purpose: determines what format the climate file is in.
//
{
    char recdType[4] = "";
    char elemType[4] = "";
    char filler[5] = "";
    char staID[80];
    char s[80];
    char line[MAXLINE];

    int  y, m, d, n;

    // --- read first line of file
    if ( fgets(line, MAXLINE, sp->Fclimate.file) == NULL ) return UNKNOWN_FORMAT;

    // --- check for TD3200 format
    sstrncpy(recdType, line, 3);
    sstrncpy(filler, &line[23], 4);
    if ( strcmp(recdType, "DLY") == 0 &&
         strcmp(filler, "9999")  == 0 ) return TD3200;

    // --- check for DLY0204 format
    if ( strlen(line) >= 233 )
    {
        sstrncpy(elemType, &line[13], 3);
        n = atoi(elemType);
        if ( n == 1 || n == 2 || n == 151 ) return DLY0204;
    }

    // --- check for USER_PREPARED format
    n = sscanf(line, "%s %d %d %d %s", staID, &y, &m, &d, s);
    if ( n == 5 ) return USER_PREPARED;

    // --- check for GHCND format                                              //(5.1.007)
    if ( isGhcndFormat(sp, line) ) return GHCND;                                   //(5.1.007)

    return UNKNOWN_FORMAT;
}

//=============================================================================

void readFileLine(SWMM_Project *sp, int *y, int *m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from next line of climate file.
//
{
    TClimateShared *clmt = &sp->ClimateShared;

    // --- read next line from climate data file
    while ( strlen(clmt->FileLine) == 0 )
    {
        if ( fgets(clmt->FileLine, MAXLINE, sp->Fclimate.file) == NULL ) return;
     	if ( clmt->FileLine[0] == '\n' ) clmt->FileLine[0] = '\0';
    }

    // --- parse year & month from line
    switch (clmt->FileFormat)
    {
    case  USER_PREPARED: readUserFileLine(sp, y, m);   break;
    case  TD3200:        readTD3200FileLine(sp, y, m);  break;
    case  DLY0204:       readDLY0204FileLine(sp, y, m); break;
    case  GHCND:         readGhcndFileLine(sp, y, m);   break;                      //(5.1.007)
    }
}

//=============================================================================

void readUserFileLine(SWMM_Project *sp, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of User-Prepared climate file.
//
{
    int n;
    char staID[80];

    TClimateShared *clmt = &sp->ClimateShared;

    n = sscanf(clmt->FileLine, "%s %d %d", staID, y, m);
    if ( n < 3 )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_FILE_READ, sp->Fclimate.name);
    }
}

//=============================================================================

void readTD3200FileLine(SWMM_Project *sp, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of TD-3200 climate file.
//
{
    char recdType[4] = "";
    char year[5] = "";
    char month[3] = "";
    int  len;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- check for minimum number of characters
    len = strlen(clmt->FileLine);
    if ( len < 30 )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_FILE_READ, sp->Fclimate.name);
        return;
    }

    // --- check for proper type of record
    sstrncpy(recdType, clmt->FileLine, 3);
    if ( strcmp(recdType, "DLY") != 0 )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_FILE_READ, sp->Fclimate.name);
        return;
    }

    // --- get record's date
    sstrncpy(year,  &clmt->FileLine[17], 4);
    sstrncpy(month, &clmt->FileLine[21], 2);
    *y = atoi(year);
    *m = atoi(month);
}

//=============================================================================

void readDLY0204FileLine(SWMM_Project *sp, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of DLY02 or DLY04 climate file.
//
{
    char year[5] = "";
    char month[3] = "";
    int  len;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- check for minimum number of characters
    len = strlen(clmt->FileLine);
    if ( len < 16 )
    {
        report_writeErrorMsg(sp, ERR_CLIMATE_FILE_READ, sp->Fclimate.name);
        return;
    }

    // --- get record's date
    sstrncpy(year,  &clmt->FileLine[7], 4);
    sstrncpy(month, &clmt->FileLine[11], 2);
    *y = atoi(year);
    *m = atoi(month);
}

//=============================================================================

void readFileValues(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reads next month's worth of data from climate file.
//
{
    int  i, j;
    int  y, m;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- initialize FileData array to missing values
    for ( i=0; i<MAXCLIMATEVARS; i++)
    {
        for (j=0; j<MAXDAYSPERMONTH; j++) clmt->FileData[i][j] = MISSING;
    }

    while ( !sp->ErrorCode )
    {
        // --- return when date on line is after current file date
        if ( feof(sp->Fclimate.file) ) return;
        readFileLine(sp, &y, &m);
        if ( y > clmt->FileYear || m > clmt->FileMonth ) return;

        // --- parse climate values from file line
        switch (clmt->FileFormat)
        {
        case  USER_PREPARED: parseUserFileLine(sp);   break;
        case  TD3200:        parseTD3200FileLine(sp);  break;
        case  DLY0204:       parseDLY0204FileLine(sp); break;
        case  GHCND:         parseGhcndFileLine(sp);   break;                    //(5.1.007)
        }
        strcpy(clmt->FileLine, "");
    }
}

//=============================================================================

void parseUserFileLine(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: parses climate variable values from a line of a user-prepared
//           climate file.
//
{
    int   n;
    int   y, m, d;
    char  staID[80];
    char  s0[80];
    char  s1[80];
    char  s2[80];
    char  s3[80];
    double x;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- read day, Tmax, Tmin, Evap, & Wind from file line
    n = sscanf(clmt->FileLine, "%s %d %d %d %s %s %s %s",
        staID, &y, &m, &d, s0, s1, s2, s3);
    if ( n < 4 ) return;
    if ( d < 1 || d > 31 ) return;

    // --- process TMAX
    if ( strlen(s0) > 0 && *s0 != '*' )
    {
        x = atof(s0);
        if ( sp->UnitSystem == SI ) x = 9./5.*x + 32.0;
        clmt->FileData[TMAX][d] =  x;
    }

    // --- process TMIN
    if ( strlen(s1) > 0 && *s1 != '*' )
    {
        x = atof(s1);
        if ( sp->UnitSystem == SI ) x = 9./5.*x + 32.0;
        clmt->FileData[TMIN][d] =  x;
    }

    // --- process EVAP
    if ( strlen(s2) > 0 && *s2 != '*' ) clmt->FileData[EVAP][d] = atof(s2);

    // --- process WIND
    if ( strlen(s3) > 0 && *s3 != '*' ) clmt->FileData[WIND][d] = atof(s3);
}

//=============================================================================

void parseTD3200FileLine(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: parses climate variable values from a line of a TD3200 file.
//
{
    int  i;
    char param[5] = "";

    TClimateShared *clmt = &sp->ClimateShared;

    // --- parse parameter name
    sstrncpy(param, &clmt->FileLine[11], 4);

    // --- see if parameter is temperature, evaporation or wind speed
    for (i=0; i<MAXCLIMATEVARS; i++)
    {
        if (strcmp(param, ClimateVarWords[i]) == 0 ) setTD3200FileValues(sp, i);
    }
}

//=============================================================================

void setTD3200FileValues(SWMM_Project *sp, int i)
//
//  Input:   i = climate variable code
//  Output:  none
//  Purpose: reads month worth of values for climate variable from TD-3200 file.
//
{
    char valCount[4] = "";
    char day[3] = "";
    char sign[2] = "";
    char value[6] = "";
    char flag2[2] = "";
    double x;
    int  nValues;
    int  j, k, d;
    int  lineLength;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- parse number of days with data from cols. 27-29 of file line
    sstrncpy(valCount, &clmt->FileLine[27], 3);
    nValues = atoi(valCount);
    lineLength = strlen(clmt->FileLine);

    // --- check for enough characters on line
    if ( lineLength >= 12*nValues + 30 )
    {
        // --- for each day's value
        for (j=0; j<nValues; j++)
        {
            // --- parse day, value & flag from file line
            k = 30 + j*12;
            sstrncpy(day,   &clmt->FileLine[k], 2);
            sstrncpy(sign,  &clmt->FileLine[k+4], 1);
            sstrncpy(value, &clmt->FileLine[k+5], 5);
            sstrncpy(flag2, &clmt->FileLine[k+11], 1);

            // --- if value is valid then store it in FileData array
            d = atoi(day);
            if ( strcmp(value, "99999") != 0
                 && ( flag2[0] == '0' || flag2[0] == '1')
                 &&   d > 0
                 &&   d <= 31 )
            {
                // --- convert from string value to numerical value
                x = atof(value);
                if ( sign[0] == '-' ) x = -x;

                // --- convert evaporation from hundreths of inches
                if ( i == EVAP )
                {
                    x /= 100.0;

                    // --- convert to mm if using SI units
                    if ( sp->UnitSystem == SI ) x *= MMperINCH;
                }

                // --- convert wind speed from miles/day to miles/hour
                if ( i == WIND ) x /= 24.0;

                // --- store value
                clmt->FileData[i][d] = x;
            }
        }
    }
}

//=============================================================================

void parseDLY0204FileLine(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: parses a month's worth of climate variable values from a line of
//           a DLY02 or DLY04 climate file.
//
{
    int  j, k, p;
    char param[4] = "";
    char sign[2]  = "";
    char value[6] = "";
    char code[2]  = "";
    double x;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- parse parameter name
    sstrncpy(param, &clmt->FileLine[13], 3);

    // --- see if parameter is min or max temperature
    p = atoi(param);
    if ( p == 1 ) p = TMAX;
    else if ( p == 2 ) p = TMIN;
    else if ( p == 151 ) p = EVAP;
    else return;

    // --- check for 233 characters on line
    if ( strlen(clmt->FileLine) < 233 ) return;

    // --- for each of 31 days
    k = 16;
    for (j=1; j<=31; j++)
    {
        // --- parse value & flag from file line
        sstrncpy(sign,  &clmt->FileLine[k], 1);
        sstrncpy(value, &clmt->FileLine[k+1], 5);
        sstrncpy(code,  &clmt->FileLine[k+6], 1);
        k += 7;

        // --- if value is valid then store it in FileData array

        if ( strcmp(value, "99999") != 0 && strcmp(value, "     ") != 0 )
        {
            switch (p)
            {
            case TMAX:
            case TMIN:
                // --- convert from integer tenths of a degree C to degrees F
                x = atof(value) / 10.0;
                if ( sign[0] == '-' ) x = -x;
                x = 9./5.*x + 32.0;
                break;
            case EVAP:
                // --- convert from 0.1 mm to inches or mm
                x = atof(value) / 10.0;
                if ( sp->UnitSystem == US ) x /= MMperINCH;
                break;
			default: return;
            }
            clmt->FileData[p][j] = x;
        }
    }
}

//=============================================================================

////  This function was added to release 5.1.007.  ////                        //(5.1.007)

int isGhcndFormat(SWMM_Project *sp, char* line)
//
//  Input:   line = first line of text from a climate file
//  Output:  returns TRUE if climate file is in NCDC GHCN Daily format.
//  Purpose: Checks if a climate file is in the NCDC GHCN Daily format
//           and determines the position of each climate variable field.
//
{
    int i;
    char* ptr;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- find starting position of the DATE field
    ptr = strstr(line, "DATE");
    if ( ptr == NULL ) return FALSE;
    clmt->FileDateFieldPos = ptr - line;

    // --- initialize starting position of each data field
    for ( i = TMIN; i <= WIND; i++) clmt->FileFieldPos[i] = -1;

    // --- find starting position of each climate variable's data field
    ptr = strstr(line, "TMIN");
    if ( ptr ) clmt->FileFieldPos[TMIN] = ptr - line;
    ptr = strstr(line, "TMAX");
    if ( ptr ) clmt->FileFieldPos[TMAX] = ptr - line;
    ptr = strstr(line, "EVAP");
    if ( ptr ) clmt->FileFieldPos[EVAP] = ptr - line;

    // --- WIND can either be daily movement or average speed
    clmt->FileWindType = WDMV;
    ptr = strstr(line, "WDMV");
    if ( ptr == NULL )
    {
        clmt->FileWindType = AWND;
        ptr = strstr(line, "AWND");
    }
    if ( ptr ) clmt->FileFieldPos[WIND] = ptr - line;

    // --- check if at least one climate variable was found
    for (i = TMIN; i <= WIND; i++) if (clmt->FileFieldPos[i] >= 0 ) return TRUE;
    return FALSE;
}

//=============================================================================

////  This function was added to release 5.1.007.  ////                        //(5.1.007)

void readGhcndFileLine(SWMM_Project *sp, int* y, int* m)
//
//  Input:   none
//  Output:  y = year
//           m = month
//  Purpose: reads year & month from line of a NCDC GHCN Daily climate file.
//
{
    int n;

    TClimateShared *clmt = &sp->ClimateShared;

    n = sscanf(&clmt->FileLine[clmt->FileDateFieldPos], "%4d%2d", y, m);

    if ( n != 2 )
    {
        *y = -99999;
        *m = -99999;
    }
}

//=============================================================================

////  This function was added to release 5.1.007.  ////                        //(5.1.007)

void parseGhcndFileLine(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: parses a line of a NCDC GHCN Daily file for daily
//           values of max/min temperature, pan evaporation and
//           wind speed.
//
{
    int y, m, d, n, v;
    double x;

    TClimateShared *clmt = &sp->ClimateShared;

    // --- parse day of month from date field
    n = sscanf(&clmt->FileLine[clmt->FileDateFieldPos], "%4d%2d%2d", &y, &m, &d);
    if ( n < 3 ) return;
    if ( d < 1 || d > 31 ) return;

    // --- parse temperatures (in tenths of deg. C) to deg F
    if ( clmt->FileFieldPos[TMAX] >= 0 )
    {
        if ( sscanf(&clmt->FileLine[clmt->FileFieldPos[TMAX]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
                clmt->FileData[TMAX][d] = (double)v*0.1*9.0/5.0 + 32.0;
        }
    }
    if ( clmt->FileFieldPos[TMIN] >= 0 )
    {
        if ( sscanf(&clmt->FileLine[clmt->FileFieldPos[TMIN]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
                clmt->FileData[TMIN][d] = (double)v*0.1*9.0/5.0 + 32.0;
        }
    }

    // -- parse evaporation (in tenths of mm) to user units
    if ( clmt->FileFieldPos[EVAP] >= 0 )
    {
        if ( sscanf(&clmt->FileLine[clmt->FileFieldPos[EVAP]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
            {
                x = (double)v * 0.1;
                if ( sp->UnitSystem == US ) x /= MMperINCH;
                clmt->FileData[EVAP][d] = x;
            }
        }
    }

    // --- parse wind speed (in km/day for WDMV or tenths of m/s for AWND)
    //     to miles/hr
    if ( clmt->FileFieldPos[WIND] >= 0 )
    {
        if ( sscanf(&clmt->FileLine[clmt->FileFieldPos[WIND]], "%8d", &v) > 0 )
        {
            if ( abs(v) < 9999 )
            {
                if ( clmt->FileWindType == WDMV ) x = (double)v * 0.62137 / 24.;
                else x = (double)v * 0.1 / 1000. * 0.62137 * 3600.;
                clmt->FileData[WIND][d] = x;
            }
        }
    }
}

//=============================================================================

////  New function added to release 5.1.010.  ////                             //(5.1.010)

void updateTempMoveAve(SWMM_Project *sp, double tmin, double tmax)
//
//  Input:   tmin = minimum daily temperature (deg F)
//           tmax = maximum daily temperature (deg F)
//  Output:  none
//  Purpose: updates moving averages of average daily temperature
//           and daily temperature range stored in structure Tma.
//
{
    double ta,               // new day's average temperature (deg F)
           tr;               // new day's temperature range (deg F)

    TClimateShared *clmt = &sp->ClimateShared;

    int    count = clmt->Tma.count;

    // --- find ta and tr from new day's min and max temperature
    ta = (tmin + tmax) / 2.0;
    tr = fabs(tmax - tmin);

    // --- if the array used to store previous days' temperatures is full
    if ( count == clmt->Tma.maxCount )
    {
        // --- update the moving averages with the new day's value
        clmt->Tma.tAve = (clmt->Tma.tAve * count + ta - clmt->Tma.ta[clmt->Tma.front]) / count;
        clmt->Tma.tRng = (clmt->Tma.tRng * count + tr - clmt->Tma.tr[clmt->Tma.front]) / count;

        // --- replace the values at the front of the moving average window
        clmt->Tma.ta[clmt->Tma.front] = ta;
        clmt->Tma.tr[clmt->Tma.front] = tr;

        // --- move the front one position forward
        clmt->Tma.front++;
        if ( clmt->Tma.front == count ) clmt->Tma.front = 0;
    }

    // --- array of previous day's values not full (at start of simulation)
    else
    {
        // --- find new moving averages by adding new values to previous ones
        clmt->Tma.tAve = (clmt->Tma.tAve * count + ta) / (count + 1);
        clmt->Tma.tRng = (clmt->Tma.tRng * count + tr) / (count + 1);

        // --- save new day's values
        clmt->Tma.ta[clmt->Tma.front] = ta;
        clmt->Tma.tr[clmt->Tma.front] = tr;

        // --- increment count and front of moving average window
        clmt->Tma.count++;
        clmt->Tma.front++;
        if ( clmt->Tma.count == clmt->Tma.maxCount ) clmt->Tma.front = 0;
    }
}
