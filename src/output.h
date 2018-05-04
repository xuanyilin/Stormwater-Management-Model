/*
 * output.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_OUTPUT_H_
#define SRC_OUTPUT_H_


// Definition of 4-byte integer, 4-byte real and 8-byte real types
#define INT4  int
#define REAL4 float
#define REAL8 double

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    INT4      IDStartPos;           // starting file position of ID names
    INT4      InputStartPos;        // starting file position of input data
    INT4      OutputStartPos;       // starting file position of output data
    INT4      BytesPerPeriod;       // bytes saved per simulation time period
    INT4      NsubcatchResults;     // number of subcatchment output variables
    INT4      NnodeResults;         // number of node output variables
    INT4      NlinkResults;         // number of link output variables
    INT4      NumSubcatch;          // number of subcatchments reported on
    INT4      NumNodes;             // number of nodes reported on
    INT4      NumLinks;             // number of links reported on
    INT4      NumPolluts;           // number of pollutants reported on
    REAL4     SysResults[MAX_SYS_RESULTS];    // values of system output vars.
} TOutputShared;

//-----------------------------------------------------------------------------
//  Exportable variables (shared with report.c)
//-----------------------------------------------------------------------------
REAL4*           SubcatchResults;
REAL4*           NodeResults;
REAL4*           LinkResults;


#endif /* SRC_OUTPUT_H_ */
