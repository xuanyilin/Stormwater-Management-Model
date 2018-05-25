/*
 * input.h
 *
 *  Created on: May 3, 2018
 *      Author: mtryby
 */

#ifndef SRC_INPUT_H_
#define SRC_INPUT_H_


#include "enums.h"

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    char *Tok[MAXTOKS];             // String tokens from line of input
    int  Ntokens;                   // Number of tokens in line of input
    int  Mobjects[MAX_OBJ_TYPES];   // Working number of objects of each type
    int  Mnodes[MAX_NODE_TYPES];    // Working number of node objects
    int  Mlinks[MAX_LINK_TYPES];    // Working number of link objects
    int  Mevents;                   // Working number of event periods      //(5.1.011)
} TInputShared;


#endif /* SRC_INPUT_H_ */
