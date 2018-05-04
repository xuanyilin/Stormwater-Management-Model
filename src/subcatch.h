/*
 * subcatch.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_SUBCATCH_H_
#define SRC_SUBCATCH_H_


//---------------
// SUBAREA OBJECT
//---------------
// An array of 3 subarea objects is associated with each subcatchment object.
// They describe the runoff process on 3 types of surfaces:
//   1 - impervious with no depression storage
//   2 - impervious with depression storage
//   3 - pervious
typedef struct
{
   int           routeTo;         // code indicating where outflow is sent
   double        fOutlet;         // fraction of outflow to outlet
   double        N;               // Manning's n
   double        fArea;           // fraction of total area
   double        dStore;          // depression storage (ft)
   //-----------------------------
   double        alpha;           // overland flow factor
   double        inflow;          // inflow rate (ft/sec)
   double        runoff;          // runoff rate (ft/sec)
   double        depth;           // depth of surface runoff (ft)
}  TSubarea;

//-----------------------------------------------------------------------------
// Globally shared variables
//-----------------------------------------------------------------------------
// Volumes (ft3) for a subcatchment over a time step                           //(5.1.008)
typedef struct
{
    double     Vevap;         // evaporation
    double     Vpevap;        // pervious area evaporation
    double     Vinfil;        // non-LID infiltration
    double     Vinflow;       // non-LID precip + snowmelt + runon + ponded water
    double     Voutflow;      // non-LID runoff to subcatchment's outlet
    double     VlidIn;        // impervious area flow to LID units
    double     VlidInfil;     // infiltration from LID units
    double     VlidOut;       // surface outflow from LID units
    double     VlidDrain;     // drain outflow from LID units
    double     VlidReturn;    // LID outflow returned to pervious area

//-----------------------------------------------------------------------------
// Locally shared variables
//-----------------------------------------------------------------------------
    TSubarea* theSubarea;     // subarea to which getDdDt() is applied
} TSubcatchShared;


#endif /* SRC_SUBCATCH_H_ */
