/*
 * toposort.h
 *
 *  Created on: May 7, 2018
 *      Author: mtryby
 */

#ifndef SRC_TOPOSORT_H_
#define SRC_TOPOSORT_H_


//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    int* InDegree;                  // number of incoming links to each node
    int* StartPos;                  // start of a node's outlinks in AdjList
    int* AdjList;                   // list of outlink indexes for each node
    int* Stack;                     // array of nodes "reached" during sorting
    int  First;                     // position of first node in stack
    int  Last;                      // position of last node added to stack

    char* Examined;                 // TRUE if node included in spanning tree
    char* InTree;                   // state of each link in spanning tree:
                                       // 0 = unexamined,
                                       // 1 = in spanning tree,
                                       // 2 = chord of spanning tree
    int*  LoopLinks;                // list of links which forms a loop
    int   LoopLinksLast;            // number of links in a loop
} TToposortShared;


#endif /* SRC_TOPOSORT_H_ */
