/*
 * rain.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_RAIN_H_
#define SRC_RAIN_H_


//--------------------
// RAINFALL STATISTICS
//--------------------
typedef struct
{
   DateTime    startDate;
   DateTime    endDate;
   long        periodsRain;
   long        periodsMissing;
   long        periodsMalfunc;
}  TRainStats;


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    TRainStats RainStats;                  // see objects.h for definition
    int        Condition;                  // rainfall condition code
    int        TimeOffset;                 // time offset of rainfall reading (sec)
    int        DataOffset;                 // start of data on line of input
    int        ValueOffset;                // start of rain value on input line
    int        RainType;                   // rain measurement type code
    int        Interval;                   // rain measurement interval (sec)
    double     UnitsFactor;                // units conversion factor
    float      RainAccum;                  // rainfall depth accumulation
    char       *StationID;                 // station ID appearing in rain file
    DateTime   AccumStartDate;             // date when accumulation begins
    DateTime   PreviousDate;               // date of previous rainfall record
    int        GageIndex;                  // index of rain gage analyzed
    int        hasStationName;             // true if data contains station name
} TRainShared;


#endif /* SRC_RAIN_H_ */
