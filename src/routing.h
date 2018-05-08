/*
 * routing.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_ROUTING_H_
#define SRC_ROUTING_H_


//-----------------------------------------------------------------------------
// Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    int* SortedLinks;
    int  NextEvent;                                                         //(5.1.011)
    int  BetweenEvents;                                                     //(5.1.012)
} TRoutingShared;


#endif /* SRC_ROUTING_H_ */
