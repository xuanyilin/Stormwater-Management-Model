/*
 * iface.h
 *
 *  Created on: May 2, 2018
 *      Author: mtryby
 */

#ifndef SRC_IFACE_H_
#define SRC_IFACE_H_


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    int      IfaceFlowUnits;        // flow units for routing interface file
    int      IfaceStep;             // interface file time step (sec)
    int      NumIfacePolluts;       // number of pollutants in interface file
    int*     IfacePolluts;          // indexes of interface file pollutants
    int      NumIfaceNodes;         // number of nodes on interface file
    int*     IfaceNodes;            // indexes of nodes on interface file
    double** OldIfaceValues;        // interface flows & WQ at previous time
    double** NewIfaceValues;        // interface flows & WQ at next time
    double   IfaceFrac;             // fraction of interface file time step
    DateTime OldIfaceDate;          // previous date of interface values
    DateTime NewIfaceDate;          // next date of interface values
} TIfaceShared;


#endif /* SRC_IFACE_H_ */
