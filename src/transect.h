/*
 * transect.h
 *
 *  Created on: May 7, 2018
 *      Author: mtryby
 */

#ifndef SRC_TRANSECT_H_
#define SRC_TRANSECT_H_


//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
#define MAXSTATION 1500                // max. number of stations in a transect

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    int    Ntransects;              // total number of transects
    int    Nstations;               // number of stations in current transect
    double  Station[MAXSTATION+1];  // x-coordinate of each station
    double  Elev[MAXSTATION+1];     // elevation of each station
    double  Nleft;                  // Manning's n for left overbank
    double  Nright;                 // Manning's n for right overbank
    double  Nchannel;               // Manning's n for main channel
    double  Xleftbank;              // station where left overbank ends
    double  Xrightbank;             // station where right overbank begins
    double  Xfactor;                // multiplier for station spacing
    double  Yfactor;                // factor added to station elevations
    double  Lfactor;                // main channel/flood plain length
} TTransectShared;


#endif /* SRC_TRANSECT_H_ */
