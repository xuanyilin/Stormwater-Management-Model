/*
 * project.h
 *
 *  Created on: May 4, 2018
 *      Author: mtryby
 */

#ifndef SRC_PROJECT_H_
#define SRC_PROJECT_H_


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    HTtable* Htable[MAX_OBJ_TYPES]; // Hash tables for object ID names
    char     MemPoolAllocated;      // TRUE if memory pool allocated
} TProjectShared;


#endif /* SRC_PROJECT_H_ */
