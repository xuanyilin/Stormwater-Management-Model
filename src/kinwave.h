/*
 * kinwave.h
 *
 *  Created on: May 3, 2018
 *      Author: mtryby
 */

#ifndef SRC_KINWAVE_H_
#define SRC_KINWAVE_H_


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    double   Beta1;
    double   C1;
    double   C2;
    double   Afull;
    double   Qfull;
    TXsect*  pXsect;
} TKinwaveShared;


#endif /* SRC_KINWAVE_H_ */
