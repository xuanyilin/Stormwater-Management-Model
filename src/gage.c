//-----------------------------------------------------------------------------
//   gage.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/10  (Build 5.1.001)
//             09/15/14  (Build 5.1.007)
//   Author:   L. Rossman
//
//   Rain gage functions.
//
//   Build 5.1.007:
//   - Support for monthly rainfall adjustments added.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <string.h>
#include <math.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
const double OneSecond = 1.1574074e-5;

//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  gage_readParams        (called by input_readLine)
//  gage_validate          (called by project_validate)
//  gage_initState         (called by project_init)
//  gage_setState          (called by runoff_execute & getRainfall in rdii.c)
//  gage_getPrecip         (called by subcatch_getRunoff)
//  gage_getNextRainDate   (called by runoff_getTimeStep)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    readGageSeriesFormat(SWMM_Project *sp, char* tok[], int ntoks, double x[]);
static int    readGageFileFormat(SWMM_Project *sp, char* tok[], int ntoks, double x[]);
static int    getFirstRainfall(SWMM_Project *sp, int gage);
static int    getNextRainfall(SWMM_Project *sp, int gage);
static double convertRainfall(SWMM_Project *sp, int gage, double rain);


//=============================================================================

int gage_readParams(SWMM_Project *sp, int j, char* tok[], int ntoks)
//
//  Input:   j = rain gage index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads rain gage parameters from a line of input data
//
//  Data formats are:
//    Name RainType RecdFreq SCF TIMESERIES SeriesName
//    Name RainType RecdFreq SCF FILE FileName Station Units StartDate
//
{
    int      k, err;
    char     *id;
    char     fname[MAXFNAME+1];
    char     staID[MAXMSG+1];
    double   x[7];

    // --- check that gage exists
    if ( ntoks < 2 ) return error_setInpError(sp, ERR_ITEMS, "");
    id = project_findID(sp, GAGE, tok[0]);
    if ( id == NULL ) return error_setInpError(sp, ERR_NAME, tok[0]);

    // --- assign default parameter values
    x[0] = -1.0;         // No time series index
    x[1] = 1.0;          // Rain type is volume
    x[2] = 3600.0;       // Recording freq. is 3600 sec
    x[3] = 1.0;          // Snow catch deficiency factor
    x[4] = NO_DATE;      // Default is no start/end date
    x[5] = NO_DATE;
    x[6] = 0.0;          // US units
    strcpy(fname, "");
    strcpy(staID, "");

    if ( ntoks < 5 ) return error_setInpError(sp, ERR_ITEMS, "");
    k = findmatch(tok[4], GageDataWords);
    if      ( k == RAIN_TSERIES )
    {
        err = readGageSeriesFormat(sp, tok, ntoks, x);
    }
    else if ( k == RAIN_FILE    )
    {
        if ( ntoks < 8 ) return error_setInpError(sp, ERR_ITEMS, "");
        sstrncpy(fname, tok[5], MAXFNAME);
        sstrncpy(staID, tok[6], MAXMSG);
        err = readGageFileFormat(sp, tok, ntoks, x);
    }
    else return error_setInpError(sp, ERR_KEYWORD, tok[4]);

    // --- save parameters to rain gage object
    if ( err > 0 ) return err;
    sp->Gage[j].ID = id;
    sp->Gage[j].tSeries      = (int)x[0];
    sp->Gage[j].rainType     = (int)x[1];
    sp->Gage[j].rainInterval = (int)x[2];
    sp->Gage[j].snowFactor   = x[3];
    sp->Gage[j].rainUnits    = (int)x[6];
    if ( sp->Gage[j].tSeries >= 0 ) sp->Gage[j].dataSource = RAIN_TSERIES;
    else                        sp->Gage[j].dataSource = RAIN_FILE;
    if ( sp->Gage[j].dataSource == RAIN_FILE )
    {
        sstrncpy(sp->Gage[j].fname, fname, MAXFNAME);
        sstrncpy(sp->Gage[j].staID, staID, MAXMSG);
        sp->Gage[j].startFileDate = x[4];
        sp->Gage[j].endFileDate = x[5];
    }
    sp->Gage[j].unitsFactor = 1.0;
    sp->Gage[j].coGage = -1;
    sp->Gage[j].isUsed = FALSE;
    return 0;
}

//=============================================================================

int readGageSeriesFormat(SWMM_Project *sp, char* tok[], int ntoks, double x[])
{
    int m, ts;
    DateTime aTime;

    if ( ntoks < 6 ) return error_setInpError(sp, ERR_ITEMS, "");

    // --- determine type of rain data
    m = findmatch(tok[1], RainTypeWords);
    if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[1]);
    x[1] = (double)m;

    // --- get data time interval & convert to seconds
    if ( getDouble(tok[2], &x[2]) ) x[2] = floor(x[2]*3600 + 0.5);
    else if ( datetime_strToTime(tok[2], &aTime) )
    {
        x[2] = floor(aTime*SECperDAY + 0.5);
    }
    else return error_setInpError(sp, ERR_DATETIME, tok[2]);
    if ( x[2] <= 0.0 ) return error_setInpError(sp, ERR_DATETIME, tok[2]);

    // --- get snow catch deficiency factor
    if ( !getDouble(tok[3], &x[3]) )
        return error_setInpError(sp, ERR_DATETIME, tok[3]);;

    // --- get time series index
    ts = project_findObject(sp, TSERIES, tok[5]);
    if ( ts < 0 ) return error_setInpError(sp, ERR_NAME, tok[5]);
    x[0] = (double)ts;
    strcpy(tok[2], "");
    return 0;
}

//=============================================================================

int readGageFileFormat(SWMM_Project *sp, char* tok[], int ntoks, double x[])
{
    int   m, u;
    DateTime aDate;
    DateTime aTime;

    // --- determine type of rain data
    m = findmatch(tok[1], RainTypeWords);
    if ( m < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[1]);
    x[1] = (double)m;

    // --- get data time interval & convert to seconds
    if ( getDouble(tok[2], &x[2]) ) x[2] *= 3600;
    else if ( datetime_strToTime(tok[2], &aTime) )
    {
        x[2] = floor(aTime*SECperDAY + 0.5);
    }
    else return error_setInpError(sp, ERR_DATETIME, tok[2]);
    if ( x[2] <= 0.0 ) return error_setInpError(sp, ERR_DATETIME, tok[2]);

    // --- get snow catch deficiency factor
    if ( !getDouble(tok[3], &x[3]) )
        return error_setInpError(sp, ERR_NUMBER, tok[3]);
 
    // --- get rain depth units
    u = findmatch(tok[7], RainUnitsWords);
    if ( u < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[7]);
    x[6] = (double)u;

    // --- get start date (if present)
    if ( ntoks > 8 && *tok[8] != '*')
    {
        if ( !datetime_strToDate(sp, tok[8], &aDate) )
            return error_setInpError(sp, ERR_DATETIME, tok[8]);
        x[4] = (float) aDate;
    }
    return 0;
}

//=============================================================================

void  gage_validate(SWMM_Project *sp, int j)
//
//  Input:   j = rain gage index
//  Output:  none
//  Purpose: checks for valid rain gage parameters
//
//  NOTE: assumes that any time series used by a rain gage has been
//        previously validated.
//
{
    int i, k;
    int gageInterval;

    // --- for gage with time series data:
    if ( sp->Gage[j].dataSource == RAIN_TSERIES )
    {
        // --- check gage's recording interval against that of time series
        k = sp->Gage[j].tSeries;
        if ( sp->Tseries[k].refersTo >= 0 )
        {
            report_writeErrorMsg(sp, ERR_RAIN_GAGE_TSERIES, sp->Gage[j].ID);
        }
        gageInterval = (int)(floor(sp->Tseries[k].dxMin*SECperDAY + 0.5));
        if ( gageInterval > 0 && sp->Gage[j].rainInterval > gageInterval )
        {
            report_writeErrorMsg(sp, ERR_RAIN_GAGE_INTERVAL, sp->Gage[j].ID);
        } 
        if ( sp->Gage[j].rainInterval < gageInterval )
        {
            report_writeWarningMsg(sp, WARN09, sp->Gage[j].ID);
        }
        if ( sp->Gage[j].rainInterval < sp->WetStep )
        {
            report_writeWarningMsg(sp, WARN01, sp->Gage[j].ID);
            sp->WetStep = sp->Gage[j].rainInterval;
        }

        // --- see if gage uses same time series as another gage
        for (i=0; i<j; i++)
        {
            if ( sp->Gage[i].dataSource == RAIN_TSERIES && sp->Gage[i].tSeries == k )
            {
                sp->Gage[j].coGage = i;

                // --- check that both gages record same type of data
                if ( sp->Gage[j].rainType != sp->Gage[i].rainType )
                {
                    report_writeErrorMsg(sp, ERR_RAIN_GAGE_FORMAT, sp->Gage[j].ID);
                }
                return;
            }
        }
    }
}

//=============================================================================

void  gage_initState(SWMM_Project *sp, int j)
//
//  Input:   j = rain gage index
//  Output:  none
//  Purpose: initializes state of rain gage.
//
{
    // --- assume gage not used by any subcatchment
    //     (will be updated in subcatch_initState)
    sp->Gage[j].isUsed = FALSE;
    sp->Gage[j].rainfall = 0.0;
    sp->Gage[j].reportRainfall = 0.0;
    // --- rainfall api sets external rainfall rate
    sp->Gage[j].externalRain = 0.0;
    if ( sp->IgnoreRainfall ) return;

    // --- for gage with file data:
    if ( sp->Gage[j].dataSource == RAIN_FILE)
    {
        // --- set current file position to start of period of record
        sp->Gage[j].currentFilePos = sp->Gage[j].startFilePos;

        // --- assign units conversion factor
        //     (rain depths on interface file are in inches)
        if ( sp->UnitSystem == SI ) sp->Gage[j].unitsFactor = MMperINCH;
    }

    // --- get first & next rainfall values
    if ( getFirstRainfall(sp, j) )
    {
        // --- find date at end of starting rain interval
        sp->Gage[j].endDate = datetime_addSeconds(
                          sp->Gage[j].startDate, sp->Gage[j].rainInterval);

        // --- if rainfall record begins after start of simulation,
        if ( sp->Gage[j].startDate > sp->StartDateTime )
        {
            // --- make next rainfall date the start of the rain record
            sp->Gage[j].nextDate = sp->Gage[j].startDate;
            sp->Gage[j].nextRainfall = sp->Gage[j].rainfall;

            // --- make start of current rain interval the simulation start
            sp->Gage[j].startDate = sp->StartDateTime;
            sp->Gage[j].endDate = sp->Gage[j].nextDate;
            sp->Gage[j].rainfall = 0.0;
        }

        // --- otherwise find next recorded rainfall
        else if ( !getNextRainfall(sp, j) ) sp->Gage[j].nextDate = NO_DATE;
    }
    else sp->Gage[j].startDate = NO_DATE;
}

//=============================================================================

void gage_setState(SWMM_Project *sp, int j, DateTime t)
//
//  Input:   j = rain gage index
//           t = a calendar date/time
//  Output:  none
//  Purpose: updates state of rain gage for specified date. 
//
{
    // --- return if gage not used by any subcatchment
    if ( sp->Gage[j].isUsed == FALSE ) return;

    // --- set rainfall to zero if disabled
    if ( sp->IgnoreRainfall )
    {
        sp->Gage[j].rainfall = 0.0;
        return;
    }

    // --- use rainfall from co-gage (gage with lower index that uses
    //     same rainfall time series or file) if it exists
    if ( sp->Gage[j].coGage >= 0)
    {
        sp->Gage[j].rainfall = sp->Gage[sp->Gage[j].coGage].rainfall;
        return;
    }

    // --- otherwise march through rainfall record until date t is bracketed
    t += OneSecond;
    for (;;)
    {
	// --- no rainfall if no interval start date
	if ( sp->Gage[j].startDate == NO_DATE )
	{
	    sp->Gage[j].rainfall = 0.0;
	    return;
	}

	// --- no rainfall if time is before interval start date
	if ( t < sp->Gage[j].startDate )
	{
	    sp->Gage[j].rainfall = 0.0;
	    return;
	}

	// --- use current rainfall if time is before interval end date
	if ( t < sp->Gage[j].endDate )
	{
	    return;
	}

	// --- no rainfall if t >= interval end date & no next interval exists
	if ( sp->Gage[j].nextDate == NO_DATE)
	{
	    sp->Gage[j].rainfall = 0.0;
	    return;
	}

	// --- no rainfall if t > interval end date & <  next interval date
	if ( t < sp->Gage[j].nextDate )
	{
	    sp->Gage[j].rainfall = 0.0;
	    return;
	}

	// --- otherwise update next rainfall interval date
	sp->Gage[j].startDate = sp->Gage[j].nextDate;
	sp->Gage[j].endDate = datetime_addSeconds(sp->Gage[j].startDate,
			  sp->Gage[j].rainInterval);
	sp->Gage[j].rainfall = sp->Gage[j].nextRainfall;

	if ( !getNextRainfall(sp, j) ) sp->Gage[j].nextDate = NO_DATE;
    }
}

//=============================================================================

DateTime gage_getNextRainDate(SWMM_Project *sp, int j, DateTime aDate)
//
//  Input:   j = rain gage index
//           aDate = calendar date/time
//  Output:  next date with rainfall occurring
//  Purpose: finds the next date from  specified date when rainfall occurs.
//
{
    if ( sp->Gage[j].isUsed == FALSE ) return aDate;
    aDate += OneSecond;
    if ( aDate < sp->Gage[j].startDate ) return sp->Gage[j].startDate;
    if ( aDate < sp->Gage[j].endDate   ) return sp->Gage[j].endDate;
    return sp->Gage[j].nextDate;
}

//=============================================================================

double gage_getPrecip(SWMM_Project *sp, int j, double *rainfall, double *snowfall)
//
//  Input:   j = rain gage index
//  Output:  rainfall = rainfall rate (ft/sec)
//           snowfall = snow fall rate (ft/sec)
//           returns total precipitation (ft/sec)
//  Purpose: determines whether gage's recorded rainfall is rain or snow.
//
{
    *rainfall = 0.0;
    *snowfall = 0.0;
    if ( !sp->IgnoreSnowmelt && sp->Temp.ta <= sp->Snow.snotmp )
    {
       *snowfall = sp->Gage[j].rainfall * sp->Gage[j].snowFactor / UCF(sp, RAINFALL);
    }
    else *rainfall = sp->Gage[j].rainfall / UCF(sp, RAINFALL);
    return (*rainfall) + (*snowfall);
} 

//=============================================================================

void gage_setReportRainfall(SWMM_Project *sp, int j, DateTime reportDate)
//
//  Input:   j = rain gage index
//           reportDate = date/time value of current reporting time
//  Output:  none
//  Purpose: sets the rainfall value reported at the current reporting time.
//
{
    double result;

    // --- use value from co-gage if it exists
    if ( sp->Gage[j].coGage >= 0)
    {
        sp->Gage[j].reportRainfall = sp->Gage[sp->Gage[j].coGage].reportRainfall;
        return;
    }

    // --- otherwise increase reporting time by 1 second to avoid
    //     roundoff problems
    reportDate += OneSecond;

    // --- use current rainfall if report date/time is before end
    //     of current rain interval
    if ( reportDate < sp->Gage[j].endDate ) result = sp->Gage[j].rainfall;

    // --- use 0.0 if report date/time is before start of next rain interval
    else if ( reportDate < sp->Gage[j].nextDate ) result = 0.0;

    // --- otherwise report date/time falls right on end of current rain
    //     interval and start of next interval so use next interval's rainfall
    else result = sp->Gage[j].nextRainfall;
    sp->Gage[j].reportRainfall = result;
}

//=============================================================================

int getFirstRainfall(SWMM_Project *sp, int j)
//
//  Input:   j = rain gage index
//  Output:  returns TRUE if successful
//  Purpose: positions rainfall record to date with first rainfall.
//
{
    int    k;                          // time series index
    float  vFirst;                     // first rain volume (ft or m)
    double rFirst;                     // first rain intensity (in/hr or mm/hr)

    // --- assign default values to date & rainfall
    sp->Gage[j].startDate = NO_DATE;
    sp->Gage[j].rainfall = 0.0;

    // --- initialize internal cumulative rainfall value
    sp->Gage[j].rainAccum = 0;

    // --- use rain interface file if applicable
    if ( sp->Gage[j].dataSource == RAIN_FILE )
    {
        if ( sp->Frain.file && sp->Gage[j].endFilePos > sp->Gage[j].startFilePos )
        {
            // --- retrieve 1st date & rainfall volume from file
            fseek(sp->Frain.file, sp->Gage[j].startFilePos, SEEK_SET);
            fread(&sp->Gage[j].startDate, sizeof(DateTime), 1, sp->Frain.file);
            fread(&vFirst, sizeof(float), 1, sp->Frain.file);
            sp->Gage[j].currentFilePos = ftell(sp->Frain.file);

            // --- convert rainfall to intensity
            sp->Gage[j].rainfall = convertRainfall(sp, j, (double)vFirst);
            return 1;
        }
        return 0;
    }

    // --- otherwise access user-supplied rainfall time series
    else
    {
        k = sp->Gage[j].tSeries;
        if ( k >= 0 )
        {
            // --- retrieve first rainfall value from time series
            if ( table_getFirstEntry(sp, &sp->Tseries[k], &sp->Gage[j].startDate,
                                     &rFirst) )
            {
                // --- convert rainfall to intensity
                sp->Gage[j].rainfall = convertRainfall(sp, j, rFirst);
                return 1;
            }
        }
        return 0;
    }
}

//=============================================================================

int getNextRainfall(SWMM_Project *sp, int j)
//
//  Input:   j = rain gage index
//  Output:  returns 1 if successful; 0 if not
//  Purpose: positions rainfall record to date with next non-zero rainfall
//           while updating the gage's next rain intensity value.
//
//  Note: zero rainfall values explicitly entered into a rain file or
//        time series are skipped over so that a proper accounting of
//        wet and dry periods can be maintained.
//
{
    int    k;                          // time series index
    float  vNext;                      // next rain volume (ft or m)
    double rNext;                      // next rain intensity (in/hr or mm/hr)

    sp->Gage[j].nextRainfall = 0.0;
    if (sp->Gage[j].dataSource == RAIN_API)
    {
	    rNext = sp->Gage[j].externalRain;
    }
    else
    {
	    do
	    {
		if ( sp->Gage[j].dataSource == RAIN_FILE )
		{
		    if ( sp->Frain.file && sp->Gage[j].currentFilePos < sp->Gage[j].endFilePos )
		    {
			fseek(sp->Frain.file, sp->Gage[j].currentFilePos, SEEK_SET);
			fread(&sp->Gage[j].nextDate, sizeof(DateTime), 1, sp->Frain.file);
			fread(&vNext, sizeof(float), 1, sp->Frain.file);
			sp->Gage[j].currentFilePos = ftell(sp->Frain.file);
			rNext = convertRainfall(sp, j, (double)vNext);
		    }
		    else return 0;
		}

		else if (sp->Gage[j].dataSource == RAIN_TSERIES)
		{
		    k = sp->Gage[j].tSeries;
		    if ( k >= 0 )
		    {
			if ( !table_getNextEntry(sp, &sp->Tseries[k],
				&sp->Gage[j].nextDate, &rNext) ) return 0;
			rNext = convertRainfall(sp, j, rNext);
		    }
		    else return 0;
		}


	    } while (rNext == 0.0);
    }
    sp->Gage[j].nextRainfall = rNext;
    return 1;
}

//=============================================================================

double convertRainfall(SWMM_Project *sp, int j, double r)
//
//  Input:   j = rain gage index
//           r = rainfall value (user units)
//  Output:  returns rainfall intensity (user units)
//  Purpose: converts rainfall value to an intensity (depth per hour).
//
{
    double r1;
    switch ( sp->Gage[j].rainType )
    {
      case RAINFALL_INTENSITY:
        r1 = r;
        break;

      case RAINFALL_VOLUME:
        r1 = r / sp->Gage[j].rainInterval * 3600.0;
        break;

      case CUMULATIVE_RAINFALL:
        if ( r  < sp->Gage[j].rainAccum )
             r1 = r / sp->Gage[j].rainInterval * 3600.0;
        else r1 = (r - sp->Gage[j].rainAccum) / sp->Gage[j].rainInterval * 3600.0;
        sp->Gage[j].rainAccum = r;
        break;

      default: r1 = r;
    }
    return r1 * sp->Gage[j].unitsFactor * sp->Adjust.rainFactor;                       //(5.1.007)
}

//=============================================================================
