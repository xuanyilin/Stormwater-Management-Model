/*
 * climate.h
 *
 *  Created on: May 1, 2018
 *      Author: mtryby
 */

#ifndef SRC_CLIMATE_H_
#define SRC_CLIMATE_H_


typedef double DateTime;

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


#endif /* SRC_CLIMATE_H_ */
