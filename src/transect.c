//-----------------------------------------------------------------------------
//   transect.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//   Author:   L. Rossman
//
//   Geometry processing for irregular cross-section transects.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "headers.h"


//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)   
//-----------------------------------------------------------------------------
//  transect_create      (called by createObjects in project.c)
//  transect_delete      (called by deleteObjects in project.c)
//  transect_readParams  (called by parseLine in input.c)
//  transect_validate    (called by input_readData)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int    setParams(SWMM_Project *sp, int transect, char* id, double x[]);
static int    setManning(SWMM_Project *sp, double n[]);
static int    addStation(SWMM_Project *sp, double x, double y);
static double getFlow(SWMM_Project *sp, int k, double a, double wp, int findFlow);
static void   getGeometry(SWMM_Project *sp, int i, int j, double y);
static void   getSliceGeom(SWMM_Project *sp, int k, double y, double yu,
        double yd, double *w, double *a, double *wp);
static void   setMaxSectionFactor(SWMM_Project *sp, int transect);

//=============================================================================

int transect_create(SWMM_Project *sp, int n)
//
//  Input:   n = number of transect objects to create
//  Output:  returns an error code
//  Purpose: creates an array of cross-section transects.
//
{
    TTransectShared *trnsct = &sp->TransectShared;

    trnsct->Ntransects = n;
    if ( n == 0 ) return 0;
    sp->Transect = (TTransect *) calloc(trnsct->Ntransects, sizeof(TTransect));
    if ( sp->Transect == NULL ) return ERR_MEMORY;
    trnsct->Nchannel = 0.0;
    trnsct->Nleft = 0.0;
    trnsct->Nright = 0.0;
    trnsct->Nstations = 0;
    return 0;
}

//=============================================================================

void transect_delete(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: deletes memory allocated for all transects.
//
{
    TTransectShared *trnsct = &sp->TransectShared;

    if ( trnsct->Ntransects == 0 ) return;
    FREE(sp->Transect);
    trnsct->Ntransects = 0;
}

//=============================================================================

int transect_readParams(SWMM_Project *sp, int* count, char* tok[], int ntoks)
//
//  Input:   count = transect index
//           tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  updated value of count,
//           returns an error code
//  Purpose: read parameters of a transect from a tokenized line of input data.
//
//  Format of transect data follows that used for HEC-2 program:
//    NC  nLeft  nRight  nChannel
//    X1  name  nSta  xLeftBank  xRightBank  0  0  0  xFactor  yFactor
//    GR  Elevation  Station  ... 
//
{
    int    i, k;
    int    index = *count;             // transect index
    int    errcode;                    // error code
    double x[10];                      // parameter values
    char*  id;                         // transect ID name

    // --- match first token to a transect keyword
    k = findmatch(tok[0], TransectKeyWords);
    if ( k < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);

    // --- read parameters associated with keyword
    switch ( k )
    {
      // --- NC line: Manning n values
      case 0:

        // --- finish processing the previous transect
        transect_validate(sp, index - 1);

        // --- read Manning's n values
        if ( ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
        for (i = 1; i <= 3; i++)
        {
            if ( ! getDouble(tok[i], &x[i]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }
        return setManning(sp, x);

      // --- X1 line: identifies start of next transect
      case 1:

        // --- check that transect was already added to project
        //     (by input_countObjects)
        if ( ntoks < 10 ) return error_setInpError(ERR_ITEMS, "");
        id = project_findID(sp, TRANSECT, tok[1]);
        if ( id == NULL ) return error_setInpError(ERR_NAME, tok[1]);

        // --- read in rest of numerical values on data line
        for ( i = 2; i < 10; i++ )
        {
            if ( ! getDouble(tok[i], &x[i]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
        }

        // --- update total transect count
        *count = index + 1;

        // --- transfer parameter values to transect's properties
        return setParams(sp, index, id, x);

      // --- GR line: station elevation & location data
      case 2:

        // --- check that line contains pairs of data values
        if ( (ntoks - 1) % 2 > 0 ) return error_setInpError(ERR_ITEMS, "");

        // --- parse each pair of Elevation-Station values
        i = 1;
        while ( i < ntoks )
        {
            if ( ! getDouble(tok[i], &x[1]) )
                return error_setInpError(ERR_NUMBER, tok[i]);
            if ( ! getDouble(tok[i+1], &x[2]) )
                return error_setInpError(ERR_NUMBER, tok[i+1]);
            errcode = addStation(sp, x[1], x[2]);
            if ( errcode ) return errcode;
            i += 2;
        }
        return 0;
    }
    return 0;
}

//=============================================================================

void  transect_validate(SWMM_Project *sp, int j)
//
//  Input:   j = transect index
//  Output:  none
//  Purpose: validates transect data and creates its geometry tables.
//
{
    int    i, nLast;
    double dy, y, ymin, ymax;

    TTransectShared *trnsct = &sp->TransectShared;

    double oldNchannel = trnsct->Nchannel;

    // --- check for valid transect data
    if ( j < 0 || j >= trnsct->Ntransects ) return;
    if ( trnsct->Nstations < 2 )
    {
        report_writeErrorMsg(sp, ERR_TRANSECT_TOO_FEW, sp->Transect[j].ID);
        return;
    }
    if ( trnsct->Nstations >= MAXSTATION )
    {
        report_writeErrorMsg(sp, ERR_TRANSECT_TOO_MANY, sp->Transect[j].ID);
        return;
    }
    if ( trnsct->Nchannel <= 0.0 )
    {
        report_writeErrorMsg(sp, ERR_TRANSECT_MANNING, sp->Transect[j].ID);
        return;
    }
    if ( trnsct->Xleftbank > trnsct->Xrightbank )
    {
        report_writeErrorMsg(sp, ERR_TRANSECT_OVERBANK, sp->Transect[j].ID);
        return;
    }

    // --- adjust main channel's Mannings n to make its equivalent
    //     length equal to that of entire flood plain
    trnsct->Nchannel = trnsct->Nchannel * sqrt(trnsct->Lfactor);
    sp->Transect[j].lengthFactor = trnsct->Lfactor;

    // --- find max. depth across transect
    ymax = trnsct->Elev[1];
    ymin = trnsct->Elev[1];
    for (i = 2; i <= trnsct->Nstations; i++)
    {
        ymax = MAX(trnsct->Elev[i], ymax);
        ymin = MIN(trnsct->Elev[i], ymin);
    }
    if ( ymin >= ymax )
    {
        report_writeErrorMsg(sp, ERR_TRANSECT_NO_DEPTH, sp->Transect[j].ID);
        return;
    }
    sp->Transect[j].yFull = ymax - ymin;

    // --- add vertical sides to transect to reach full ht. on both ends
    trnsct->Station[0] = trnsct->Station[1];
    trnsct->Elev[0] = ymax;
    trnsct->Nstations++;
    trnsct->Station[trnsct->Nstations] = trnsct->Station[trnsct->Nstations-1];
    trnsct->Elev[trnsct->Nstations] = trnsct->Elev[0];

    // --- determine size & depth increment for geometry tables
    sp->Transect[j].nTbl = N_TRANSECT_TBL;
    dy = (ymax - ymin) / (double)(sp->Transect[j].nTbl - 1);

    // --- set 1st table entries to zero
    sp->Transect[j].areaTbl[0] = 0.0;
    sp->Transect[j].hradTbl[0] = 0.0;
    sp->Transect[j].widthTbl[0] = 0.0;

    // --- compute geometry for each depth increment
    y = ymin;
    sp->Transect[j].wMax = 0.0;
    for (i = 1; i < sp->Transect[j].nTbl; i++)
    {
        y += dy;
        sp->Transect[j].areaTbl[i] = 0.0;
        sp->Transect[j].hradTbl[i] = 0.0;
        sp->Transect[j].widthTbl[i] = 0.0;
        getGeometry(sp, i, j, y);
    }

    // --- determine max. section factor 
    setMaxSectionFactor(sp, j);

    // --- normalize geometry table entries
    //     (full cross-section values are last table entries)
    nLast = sp->Transect[j].nTbl - 1;
    sp->Transect[j].aFull = sp->Transect[j].areaTbl[nLast];
    sp->Transect[j].rFull = sp->Transect[j].hradTbl[nLast];
    sp->Transect[j].wMax = sp->Transect[j].widthTbl[nLast];

    for (i = 1; i <= nLast; i++)
    {
        sp->Transect[j].areaTbl[i] /= sp->Transect[j].aFull;
        sp->Transect[j].hradTbl[i] /= sp->Transect[j].rFull;
        sp->Transect[j].widthTbl[i] /= sp->Transect[j].wMax;
    }

    // --- set width at 0 height equal to width at 4% of max. height
    sp->Transect[j].widthTbl[0] = sp->Transect[j].widthTbl[1];

    // --- save unadjusted main channel roughness 
    sp->Transect[j].roughness = oldNchannel;
}

//=============================================================================

int  setManning(SWMM_Project *sp, double n[])
//
//  Input:   n[] = array of Manning's n values
//  Output:  returns an error code
//  Purpose: sets Manning's n for overbanks and main channel of a transect.
//
{
    int i;

    TTransectShared *trnsct = &sp->TransectShared;

    for (i=1; i<=3; i++)
    {
        if ( n[i] < 0.0 ) return ERR_NUMBER;
    }
    if ( n[1] > 0.0 ) trnsct->Nleft = n[1];
    if ( n[2] > 0.0 ) trnsct->Nright = n[2];
    if ( n[3] > 0.0 ) trnsct->Nchannel = n[3];
    if ( trnsct->Nleft == 0.0  ) trnsct->Nleft = trnsct->Nchannel;
    if ( trnsct->Nright == 0.0 ) trnsct->Nright = trnsct->Nchannel;
    return 0;
}

//=============================================================================

int  setParams(SWMM_Project *sp, int j, char* id, double x[])
//
//  Input:   j = transect index
//           id = transect ID name
//           x[] = array of parameter values
//  Output:  returns an error code
//  Purpose: assigns parameter values to current transect being processed.
//
{
    TTransectShared *trnsct = &sp->TransectShared;

    if ( j < 0 || j >= trnsct->Ntransects ) return ERR_NUMBER;
    sp->Transect[j].ID = id;                         // ID name
    trnsct->Xleftbank = x[3] / UCF(sp, LENGTH);              // left overbank location
    trnsct->Xrightbank = x[4] / UCF(sp, LENGTH);             // right overbank location
    trnsct->Lfactor = x[7];                              // channel/bank length
    if ( trnsct->Lfactor == 0.0 ) trnsct->Lfactor = 1.0;
    trnsct->Xfactor = x[8];                              // station location multiplier
    if ( trnsct->Xfactor == 0.0 ) trnsct->Xfactor = 1.0;
    trnsct->Xleftbank *= trnsct->Xfactor;                        // adjusted left bank
    trnsct->Xrightbank *= trnsct->Xfactor;                       // adjusted right bank
    trnsct->Yfactor = x[9] / UCF(sp, LENGTH);                // elevation offset
    trnsct->Nstations = 0;
    return 0;
}

//=============================================================================

int  addStation(SWMM_Project *sp, double y, double x)
//
//  Input:   y = station elevation value
//           x = station distance value
//  Output:  returns an error code
//  Purpose: adds a new station to the transect currently being processed.
//
{
    TTransectShared *trnsct = &sp->TransectShared;

    // --- check for valid number of stations
    if ( trnsct->Nstations < 0 ) return ERR_TRANSECT_UNKNOWN;
    trnsct->Nstations++;
    if ( trnsct->Nstations >= MAXSTATION ) return 0;

    // --- add station distance, modified by distance multiplier
    trnsct->Station[trnsct->Nstations] = x * trnsct->Xfactor / UCF(sp, LENGTH);

    // --- add station elevation, modified by offset elevation
    trnsct->Elev[trnsct->Nstations] = (y + trnsct->Yfactor) / UCF(sp, LENGTH);

    // --- check if station distances are non-increasing
    if ( trnsct->Nstations > 1
        && trnsct->Station[trnsct->Nstations] < trnsct->Station[trnsct->Nstations-1] )
        return ERR_TRANSECT_SEQUENCE;
    return 0;    
}

//=============================================================================

void  getGeometry(SWMM_Project *sp, int i, int j, double y)
//
//  Input:   i = index of current entry in geometry tables
//           j = transect index
//           y = depth of current entry in geometry tables
//  Output:  none
//  Purpose: computes entries in a transect's geometry tables at a given depth. 
//
{
    int    k;                // station index
    double ylo,              // lower elev. of transect slice
           yhi,              // higher elev. of transect slice
           w,                // top width of transect slice
           wp,               // wetted perimeter of transect slice
           wpSum,            // total wetted perimeter across transect
           a,                // area of transect slice
           aSum,             // total area across transect
           q,                // flow across transect slices with same roughness
           qSum;             // total flow across transect
    int   findFlow;          // true if flow thru area slice needs updating

    TTransectShared *trnsct = &sp->TransectShared;

    // --- initialize
    wpSum = 0.0;
    aSum = 0.0;
    qSum = 0.0;

    // --- examine each horizontal station from left to right
    for (k = 1; k <= trnsct->Nstations; k++)
    {
        // --- determine low & high elevations for transect sub-section
        if ( trnsct->Elev[k-1] >= trnsct->Elev[k] )
        {
            yhi = trnsct->Elev[k-1];
            ylo = trnsct->Elev[k];
        }
        else
        {
            yhi = trnsct->Elev[k];
            ylo = trnsct->Elev[k-1];
        }

        // --- skip station if its totally dry
        if ( ylo >= y ) continue;

        // --- get top width, area & wetted perimeter values for transect
        //     slice between station k and k-1
        getSliceGeom(sp, k, y, ylo, yhi, &w, &a, &wp);

        // --- update total transect values
        wpSum += wp;
        aSum += a;
        sp->Transect[j].areaTbl[i] += a;
        sp->Transect[j].widthTbl[i] += w;

        // --- must update flow if station elevation is above water level
        if ( trnsct->Elev[k] >= y ) findFlow = TRUE;
        else findFlow = FALSE;

        // --- update flow across transect if called for
        q = getFlow(sp, k, aSum, wpSum, findFlow);
        if ( q > 0.0 )
        {
            qSum += q;
            aSum = 0.0;
            wpSum = 0.0;
        }

    }   // next station k 

    // --- find hyd. radius table entry solving Manning eq. with
    //     total flow, total area, and main channel n
    aSum = sp->Transect[j].areaTbl[i];
    if ( aSum == 0.0 ) sp->Transect[j].hradTbl[i] = sp->Transect[j].hradTbl[i-1];
    else sp->Transect[j].hradTbl[i] = pow(qSum * trnsct->Nchannel / 1.49 / aSum, 1.5);
}

//=============================================================================

void getSliceGeom(SWMM_Project *sp, int k, double y, double ylo, double yhi,
        double *w, double *a, double *wp)
//
//  Input:   k = station index
//           y = water elevation
//           ylo = transect elevation on low side of slice
//           yhi = transect elevation on high side of slice
//  Output   w = width of transect slice
//           a = area of transect slice
//           wp = wetted perimeter of transect slice
//  Purpose: finds area, width & wetted perim. for slice of transect that
//           is covered by given water depth.
//
//      yhi  |           
//           |
//        y  |**********
//           |********** --> slice of transect being analyzed
//      ylo  |**********|
//           |**********|
//           |**********|
//         Station    Station
//           k-1        k
//
{
    double width, ratio;

    TTransectShared *trnsct = &sp->TransectShared;

    // --- compute width & wetted perimeter of transect slice
    width = fabs(trnsct->Station[k] - trnsct->Station[k-1]);
    (*w) = width;
    (*wp) = sqrt(width * width + (yhi - ylo) * (yhi - ylo));
    (*a)  = 0.0;

    // --- find area for completely submerged slice
    if ( y > yhi )
    {
        (*a) = width * ( (y - yhi) + (y - ylo) ) / 2.0;
    }

    // --- otherwise find area and adjust width & wetted perim. for
    //     partly submerged slice
    else if ( yhi > ylo )
    {
         ratio = (y - ylo) / (yhi - ylo);
         (*a) = width * (yhi - ylo) / 2.0 * ratio * ratio;
         (*w) *= ratio;
         (*wp) *= ratio;
     }
}

//=============================================================================

double getFlow(SWMM_Project *sp, int k, double a, double wp, int findFlow)
//
//  Input:   k = index of station at end of transect sub-section
//           a = flow area of sub-section
//           wp = wetted perimeter of flow area of sub-section
//           findFlow = TRUE if flow needs updating 
//  Output:  returns normal flow (per unit of slope)
//  Purpose: finds flow through a sub-section of a transect.
//
{
    double n;                          // Manning's n

    TTransectShared *trnsct = &sp->TransectShared;

    if ( findFlow == FALSE)
    {
        // --- flow needs updating if we are at last station
        if ( k == trnsct->Nstations - 1 ) findFlow = TRUE;

        // --- flow needs updating if we are at end of left overbank and
        //     there is a change in Manning's n and section not vertical
        else if ( trnsct->Station[k] == trnsct->Xleftbank )
        {
            if ( trnsct->Nleft != trnsct->Nchannel &&
                    trnsct->Station[k] != trnsct->Station[k-1] ) findFlow = TRUE;
        }

        // --- flow needs updating if we are at start of right overbank and
        //     there is a change in Manning's n and section not vertical
        else if ( trnsct->Station[k] == trnsct->Xrightbank )
        {
            if ( trnsct->Nright != trnsct->Nchannel &&
                    trnsct->Station[k] != trnsct->Station[k+1] ) findFlow = TRUE;
        }
    }

    // --- if flow needs updating
    if ( findFlow )
    {
        // --- find value of Manning's n to use
        n = trnsct->Nchannel;
        if ( trnsct->Station[k-1] < trnsct->Xleftbank ) n = trnsct->Nleft;
        if ( trnsct->Station[k] > trnsct->Xrightbank )  n = trnsct->Nright;

        // --- compute flow through flow area
        return PHI / n * a * pow(a/wp, 2./3.);
    }
    return 0.0;
}

//=============================================================================

void setMaxSectionFactor(SWMM_Project *sp, int j)
//
//  Input:   j = transect index
//  Output:  none
//  Purpose: determines the maximum section factor for a transect and the
//           area where this maxumum occurs.
//
{
    int    i;
    double sf;

    sp->Transect[j].aMax = 0.0;
    sp->Transect[j].sMax = 0.0;
    for (i=1; i<sp->Transect[j].nTbl; i++)
    {
        sf = sp->Transect[j].areaTbl[i] * pow(sp->Transect[j].hradTbl[i], 2./3.);
        if ( sf > sp->Transect[j].sMax )
        {
            sp->Transect[j].sMax = sf;
            sp->Transect[j].aMax = sp->Transect[j].areaTbl[i];
        }
    }
}

//=============================================================================
