/*
 * climate.h
 *
 *  Created on: May 1, 2018
 *      Author: mtryby
 */

#ifndef SRC_CLIMATE_H_
#define SRC_CLIMATE_H_


//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
//enum ClimateFileFormats {UNKNOWN_FORMAT,
//                         USER_PREPARED,     // SWMM 5's own user format
//                         GHCND,             // NCDC GHCN Daily format          //(5.1.007)
//                         TD3200,            // NCDC TD3200 format
//                         DLY0204};          // Canadian DLY02 or DLY04 format
//
//static const int    MAXCLIMATEVARS  = 4;
//static const int    MAXDAYSPERMONTH = 32;

// These variables are used when processing climate files.
//enum   ClimateVarType {TMIN, TMAX, EVAP, WIND};
//enum   WindSpeedType  {WDMV, AWND};                                            //(5.1.007)
//static const char* ClimateVarWords[] = {"TMIN", "TMAX", "EVAP", "WDMV", "AWND",      //(5.1.007)
//                                  NULL};

////  Added for release 5.1.010.  ////                                         //(5.1.010)
//-----------------------------------------------------------------------------
//  Data Structures
//-----------------------------------------------------------------------------
typedef struct
{
    double    tAve;          // moving avg. for daily temperature (deg F)
    double    tRng;          // moving avg. for daily temp. range (deg F)
    double    ta[7];         // data window for tAve
    double    tr[7];         // data window for tRng
    int       count;         // length of moving average window
    int       maxCount;      // maximum length of moving average window
    int       front;         // index of front of moving average window
} TMovAve;
////

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
// Temperature variables
typedef struct
{
    double    Tmin;                 // min. daily temperature (deg F)
    double    Tmax;                 // max. daily temperature (deg F)
    double    Trng;                 // 1/2 range of daily temperatures
    double    Trng1;                // prev. max - current min. temp.
    double    Tave;                 // average daily temperature (deg F)
    double    Hrsr;                 // time of min. temp. (hrs)
    double    Hrss;                 // time of max. temp (hrs)
    double    Hrday;                // avg. of min/max temp times
    double    Dhrdy;                // hrs. between min. & max. temp. times
    double    Dydif;                // hrs. between max. & min. temp. times
    DateTime  LastDay;              // date of last day with temp. data
    TMovAve   Tma;                  // moving average of daily temperatures //(5.1.010)

//// Evaporation variables
    DateTime  NextEvapDate;         // next date when evap. rate changes
    double    NextEvapRate;         // next evaporation rate (user units)
//
//// Climate file variables
    int      FileFormat;            // file format (see ClimateFileFormats)
    int      FileYear;              // current year of file data
    int      FileMonth;             // current month of year of file data
    int      FileDay;               // current day of month of file data
    int      FileLastDay;           // last day of current month of file data
    int      FileElapsedDays;       // number of days read from file
    double   FileValue[4];          // current day's values of climate data
    double   FileData[4][32];       // month's worth of daily climate data
    char     FileLine[MAXLINE+1];   // line from climate data file
//
    int      FileFieldPos[4];       // start of data fields for file record //(5.1.007)
    int      FileDateFieldPos;      // start of date field for file record  //(5.1.007)
    int      FileWindType;          // wind speed type;                     //(5.1.007)
} TClimateShared;

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
//static int  getFileFormat(SWMM_Project *sp);
//static void readFileLine(SWMM_Project *sp, int *year, int *month);
//static void readUserFileLine(SWMM_Project *sp, int *year, int *month);
//static void readTD3200FileLine(SWMM_Project *sp, int *year, int *month);
//static void readDLY0204FileLine(SWMM_Project *sp, int *year, int *month);
//static void readFileValues(SWMM_Project *sp);
//
//static void setNextEvapDate(SWMM_Project *sp, DateTime thedate);                                 //(5.1.008)
//static void setEvap(SWMM_Project *sp, DateTime theDate);
//static void setTemp(SWMM_Project *sp, DateTime theDate);
//static void setWind(SWMM_Project *sp, DateTime theDate);
//static void updateTempTimes(SWMM_Project *sp, int day);
//static void updateTempMoveAve(double tmin, double tmax);                       //(5.1.010)
//static double getTempEvap(SWMM_Project *sp, int day, double ta, double tr);                      //(5.1.010)
//
//static void updateFileValues(SWMM_Project *sp, DateTime theDate);
//static void parseUserFileLine(SWMM_Project *sp);
//static void parseTD3200FileLine(SWMM_Project *sp);
//static void parseDLY0204FileLine(SWMM_Project *sp);
//static void setTD3200FileValues(SWMM_Project *sp, int param);
//
//static int  isGhcndFormat(char* line);                                         //(5.1.007)
//static void readGhcndFileLine(SWMM_Project *sp, int *year, int *month);                          //(5.1.007)
//static void parseGhcndFileLine(SWMM_Project *sp);                                          //(5.1.007)


#endif /* SRC_CLIMATE_H_ */
