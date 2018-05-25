//-----------------------------------------------------------------------------
//   toposort.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//   Author:   L. Rossman
//
//   Topological sorting of conveyance network links
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
enum AdjListType {UNDIRECTED, DIRECTED};    // type of nodal adjacency list


//-----------------------------------------------------------------------------
//  External functions (declared in funcs.h)   
//-----------------------------------------------------------------------------
//  toposort_sortLinks (called by routing_open)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void createAdjList(SWMM_Project *sp, int listType);
static void adjustAdjList(SWMM_Project *sp);
static int  topoSort(SWMM_Project *sp, int sortedLinks[]);
static void findCycles(SWMM_Project *sp);
static void findSpanningTree(SWMM_Project *sp, int startNode);
static void evalLoop(SWMM_Project *sp, int startLink);
static int  traceLoop(SWMM_Project *sp, int i1, int i2, int k);
static void checkDummyLinks(SWMM_Project *sp);
//=============================================================================

void toposort_sortLinks(SWMM_Project *sp, int sortedLinks[])
//
//  Input:   none
//  Output:  sortedLinks = array of link indexes in sorted order
//  Purpose: sorts links from upstream to downstream.
//
{
    int i, n = 0;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- no need to sort links for Dyn. Wave routing
    for ( i=0; i<sp->Nobjects[LINK]; i++) sortedLinks[i] = i;
    if ( sp->RouteModel == DW )
    {

        // --- check for nodes with both incoming and outgoing
        //     dummy links (creates ambiguous ordering)
        checkDummyLinks(sp);
        if ( sp->ErrorCode ) return;

        // --- find number of outflow links for each node
        for ( i=0; i<sp->Nobjects[NODE]; i++ ) sp->Node[i].degree = 0;
        for ( i=0; i<sp->Nobjects[LINK]; i++ )
        {
            // --- if upstream node is an outfall, then increment outflow
            //     count for downstream node, otherwise increment count
            //     for upstream node
            n = sp->Link[i].node1;
            if ( sp->Link[i].direction < 0 ) n = sp->Link[i].node2;
            if ( sp->Node[n].type == OUTFALL )
            {
                if ( sp->Link[i].direction < 0 ) n = sp->Link[i].node1;
                else n = sp->Link[i].node2;
                sp->Node[n].degree++;
            }
            else sp->Node[n].degree++;
        }
        return;
    }

    // --- allocate arrays used for topo sorting
    if ( sp->ErrorCode ) return;
    tpsrt->InDegree = (int *) calloc(sp->Nobjects[NODE], sizeof(int));
    tpsrt->StartPos = (int *) calloc(sp->Nobjects[NODE], sizeof(int));
    tpsrt->AdjList  = (int *) calloc(sp->Nobjects[LINK], sizeof(int));
    tpsrt->Stack    = (int *) calloc(sp->Nobjects[NODE], sizeof(int));
    if ( tpsrt->InDegree == NULL || tpsrt->StartPos == NULL ||
            tpsrt->AdjList == NULL || tpsrt->Stack == NULL )
    {
        report_writeErrorMsg(sp, ERR_MEMORY, "");
    }
    else
    {
        // --- create a directed adjacency list of links leaving each node
        createAdjList(sp, DIRECTED);

        // --- adjust adjacency list for DIVIDER nodes
        adjustAdjList(sp);

        // --- find number of links entering each node
        for (i = 0; i < sp->Nobjects[NODE]; i++)
            tpsrt->InDegree[i] = 0;
        for (i = 0; i < sp->Nobjects[LINK]; i++)
            tpsrt->InDegree[ sp->Link[i].node2 ]++;

        // --- topo sort the links
        n = topoSort(sp, sortedLinks);
    }   

    // --- free allocated memory
    FREE(tpsrt->InDegree);
    FREE(tpsrt->StartPos);
    FREE(tpsrt->AdjList);
    FREE(tpsrt->Stack);

    // --- check that all links are included in SortedLinks
    if ( !sp->ErrorCode &&  n != sp->Nobjects[LINK] )
    {
        report_writeErrorMsg(sp, ERR_LOOP, "");
        findCycles(sp);
    }
}

//=============================================================================

void createAdjList(SWMM_Project *sp ,int listType)
//
//  Input:   lsitType = DIRECTED or UNDIRECTED
//  Output:  none
//  Purpose: creates listing of links incident on each node.
//
{
    int i, j, k;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- determine degree of each node
    //     (for DIRECTED list only count link at its upstream node;
    //      for UNDIRECTED list count link at both end nodes)
    for (i = 0; i < sp->Nobjects[NODE]; i++) sp->Node[i].degree = 0;
    for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        sp->Node[ sp->Link[j].node1 ].degree++;
        if ( listType == UNDIRECTED )
            sp->Node[ sp->Link[j].node2 ].degree++;
    }

    // --- determine start position of each node in the adjacency list
    //     (the adjacency list, AdjList, is one long vector containing
    //      the individual node lists one after the other)
    tpsrt->StartPos[0] = 0;
    for (i = 0; i < sp->Nobjects[NODE]-1; i++)
    {
        tpsrt->StartPos[i+1] = tpsrt->StartPos[i] + sp->Node[i].degree;
        sp->Node[i].degree = 0;
    }
    sp->Node[sp->Nobjects[NODE]-1].degree = 0;

    // --- traverse the list of links once more,
    //     adding each link's index to the proper 
    //     position in the adjacency list
    for (j = 0; j < sp->Nobjects[LINK]; j++)
    {
        i = sp->Link[j].node1;
        k = tpsrt->StartPos[i] + sp->Node[i].degree;
        tpsrt->AdjList[k] = j;
        sp->Node[i].degree++;
        if ( listType == UNDIRECTED )
        {
            i = sp->Link[j].node2;
            k = tpsrt->StartPos[i] + sp->Node[i].degree;
            tpsrt->AdjList[k] = j;
            sp->Node[i].degree++;
        }
    }
}

//=============================================================================

void adjustAdjList(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: adjusts adjacency list for Divider nodes so that non-
//           diversion link appears before diversion link.
//
{
    int i, j, k, m;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- check each node
    for (i=0; i<sp->Nobjects[NODE]; i++)
    {
        // --- skip nodes that are not Dividers
        if ( sp->Node[i].type != DIVIDER ) continue;
        if ( sp->Node[i].degree != 2 ) continue;

        // --- switch position of outgoing links at the node if the
        //     diversion link appears first in the adjacency list
        k = sp->Node[i].subIndex;
        m = tpsrt->StartPos[i];
        j = tpsrt->AdjList[m];
        if ( j == sp->Divider[k].link )
        {
            tpsrt->AdjList[m] = tpsrt->AdjList[m+1];
            tpsrt->AdjList[m+1] = j;
        }
    }
}

//=============================================================================

int topoSort(SWMM_Project *sp, int sortedLinks[])
//
//  Input:   none
//  Output:  sortedLinks = array of sorted link indexes,
//           returns number of links successfully sorted
//  Purpose: performs a stack-based topo sort of the drainage network's links.
//
{
    int i, j, k, n;
    int i1, i2, k1, k2;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- initialize a stack which contains nodes with zero in-degree
    tpsrt->First = 0;
    tpsrt->Last = -1;
    for (i = 0; i < sp->Nobjects[NODE]; i++)
    {
        if ( tpsrt->InDegree[i] == 0 )
        {
            tpsrt->Last++;
            tpsrt->Stack[tpsrt->Last] = i;
        }
    }

    // --- traverse the stack, adding each node's outgoing link indexes
    //     to the SortedLinks array in the order processed
    n = 0;
    while ( tpsrt->First <= tpsrt->Last )
    {
        // --- determine range of adjacency list indexes belonging to 
        //     first node remaining on the stack
        i1 = tpsrt->Stack[tpsrt->First];
        k1 = tpsrt->StartPos[i1];
        k2 = k1 + sp->Node[i1].degree;

        // --- for each outgoing link from first node on stack
        for (k = k1; k < k2; k++)
        {
            // --- add link index to current position in SortedLinks
            j = tpsrt->AdjList[k];
            sortedLinks[n] = j;
            n++;

            // --- reduce in-degree of link's downstream node
            i2 = sp->Link[j].node2;
            tpsrt->InDegree[i2]--;

            // --- add downstream node to stack if its in-degree is zero
            if ( tpsrt->InDegree[i2] == 0 )
            {
                tpsrt->Last++;
                tpsrt->Stack[tpsrt->Last] = i2;
            }  
        }
        tpsrt->First++;
    }
    return n;
}

//=============================================================================

void  findCycles(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: finds all cycles in the drainage network (i.e., closed loops that
//           start and end at the same node).
//
{
    int i;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- allocate arrays
    tpsrt->AdjList  = (int *) calloc(2*sp->Nobjects[LINK], sizeof(int));
    tpsrt->StartPos = (int *) calloc(sp->Nobjects[NODE], sizeof(int));
    tpsrt->Stack    = (int *) calloc(sp->Nobjects[NODE], sizeof(int));
    tpsrt->Examined = (char *) calloc(sp->Nobjects[NODE], sizeof(char));
    tpsrt->InTree   = (char *) calloc(sp->Nobjects[LINK], sizeof(char));
    tpsrt->LoopLinks = (int *) calloc(sp->Nobjects[LINK], sizeof(int));
    if ( tpsrt->StartPos && tpsrt->AdjList && tpsrt->Stack && tpsrt->Examined &&
            tpsrt->InTree && tpsrt->LoopLinks )
    {
        // --- create an undirected adjacency list for the nodes
        createAdjList(sp, UNDIRECTED);

        // --- set to empty the list of nodes examined and the list
        //     of links in the spanning tree
        for ( i=0; i<sp->Nobjects[NODE]; i++)
            tpsrt->Examined[i] = 0;
        for ( i=0; i<sp->Nobjects[LINK]; i++)
            tpsrt->InTree[i] = 0;

        // --- find a spanning tree for each unexamined node
        //     (cycles are identified as tree is constructed)
        for ( i=0; i<sp->Nobjects[NODE]; i++)
        {
            if ( tpsrt->Examined[i] ) continue;
            tpsrt->Last = -1;
            findSpanningTree(sp, i);
        }
    }
    FREE(tpsrt->StartPos);
    FREE(tpsrt->AdjList);
    FREE(tpsrt->Stack);
    FREE(tpsrt->Examined);
    FREE(tpsrt->InTree);
    FREE(tpsrt->LoopLinks);
}

//=============================================================================

void  findSpanningTree(SWMM_Project *sp, int startNode)
//
//  Input:   i = index of starting node of tree
//  Output:  none
//  Purpose: finds continuation of network's spanning tree of links.
//
{
    int nextNode, j, k, m;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- examine each link connected to node i
    for ( m = tpsrt->StartPos[startNode];
          m < tpsrt->StartPos[startNode]+sp->Node[startNode].degree; m++ )
    {
        // --- find which node (j) connects link k from start node
        k = tpsrt->AdjList[m];
        if ( sp->Link[k].node1 == startNode ) j = sp->Link[k].node2;
        else j = sp->Link[k].node1;

        // --- if link k is not in the tree
        if ( tpsrt->InTree[k] == 0 )
        {
            // --- if connecting node already examined,
            //     then link k forms a loop; mark it as a chord
            //     and check if loop forms a cycle
            if ( tpsrt->Examined[j] )
            {
                tpsrt->InTree[k] = 2;
                evalLoop(sp, k);
            }

            // --- otherwise mark connected node as being examined,
            //     add it to the node stack, and mark the connecting
            //     link as being in the spanning tree
            else
            {
                tpsrt->Examined[j] = 1;
                tpsrt->Last++;
                tpsrt->Stack[tpsrt->Last] = j;
                tpsrt->InTree[k] = 1;
            }
        }
    }

    // --- continue to grow the spanning tree from
    //     the last node added to the stack
    if ( tpsrt->Last >= 0 )
    {
        nextNode = tpsrt->Stack[tpsrt->Last];
        tpsrt->Last--;
        findSpanningTree(sp, nextNode);
    }
}

//=============================================================================

void evalLoop(SWMM_Project *sp, int startLink)
//
//  Input:   startLink = index of starting link of a loop
//  Output:  none
//  Purpose: checks if loop starting with a given link forms a closed cycle.
//
{
    int i;                             // loop list index
    int i1, i2;                        // start & end node indexes
    int j;                             // index of link in loop
    int kount;                         // items per line counter
    int isCycle;                       // TRUE if loop forms a cycle

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- make startLink the first link in the loop
    tpsrt->LoopLinksLast = 0;
    tpsrt->LoopLinks[0] = startLink;

    // --- trace a path on the spanning tree that starts at the
    //     tail node of startLink and ends at its head node
    i1 = sp->Link[startLink].node1;
    i2 = sp->Link[startLink].node2;
    if ( !traceLoop(sp, i1, i2, startLink) ) return;

    // --- check if all links on the path are oriented head-to-tail
    isCycle = TRUE;
    j = tpsrt->LoopLinks[0];
    i2 = sp->Link[j].node2;
    for (i = 1; i <= tpsrt->LoopLinksLast; i++)
    {
        j = tpsrt->LoopLinks[i];
        i1 = sp->Link[j].node1;
        if ( i1 == i2 ) i2 = sp->Link[j].node2;
        else
        {
            isCycle = FALSE;
            break;
        }
    }

    // --- print cycle to report file
    if ( isCycle )
    {
        kount = 0;
        for (i = 0; i <= tpsrt->LoopLinksLast; i++)
        {
            if ( kount % 5 == 0 ) fprintf(sp->Frpt.file, "\n");
            kount++;
            fprintf(sp->Frpt.file, "  %s", sp->Link[tpsrt->LoopLinks[i]].ID);
            if ( i < tpsrt->LoopLinksLast ) fprintf(sp->Frpt.file, "  -->");
        }
    }
}

//=============================================================================

int traceLoop(SWMM_Project *sp, int i1, int i2, int k1)
//
//  Input:   i1 = index of current node reached while tracing a loop
//           i2 = index of final node on the loop
//           k1 = index of link which extends loop to node i1
//  Output:  returns TRUE if loop can be extended through node i1
//  Purpose: tries to extend closed loop through current node.
//
{
    int j, k, m;

    TToposortShared *tpsrt = &sp->ToposortShared;

    // --- if current node equals final node then return with loop found
    if ( i1 == i2 ) return TRUE;

    // --- examine each link connected to current node
    for (m = tpsrt->StartPos[i1]; m < tpsrt->StartPos[i1] + sp->Node[i1].degree; m++)
    {
        // --- ignore link if it comes from predecessor node or if
        //     it is not on the spanning tree
        k = tpsrt->AdjList[m];
        if ( k == k1 || tpsrt->InTree[k] != 1 ) continue;

        // --- identify other node that link is connected to
        if ( sp->Link[k].node1 == i1 ) j = sp->Link[k].node2;
        else                       j = sp->Link[k].node1;

        // --- try to continue tracing the loop from this node;
        //     if successful, then add link to loop and return
        if ( traceLoop(sp, j, i2, k) )
        {
            tpsrt->LoopLinksLast++;
            tpsrt->LoopLinks[tpsrt->LoopLinksLast] = k;
            return TRUE;
        }
    }

    // --- return false if loop cannot be continued from current node
    return FALSE;
}

//=============================================================================

void checkDummyLinks(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: checks for nodes that have both incoming and outgoing dummy links.
//
{
    int   i, j;
    int* marked;

    // --- create an array that marks nodes
    //     (calloc initializes the array to 0)
    marked = (int *) calloc(sp->Nobjects[NODE], sizeof(int));
    if ( marked == NULL )
    {
        report_writeErrorMsg(sp, ERR_MEMORY, "");
        return;
    }

    // --- mark nodes that whose incoming links are all
    //     either dummy links or ideal pumps
    for ( i = 0; i < sp->Nobjects[LINK]; i++ )
    {
        j = sp->Link[i].node2;
        if ( sp->Link[i].direction < 0 ) j = sp->Link[i].node1;
        if ( (sp->Link[i].type == CONDUIT && sp->Link[i].xsect.type == DUMMY) ||
             (sp->Link[i].type == PUMP &&
              sp->Pump[sp->Link[i].subIndex].type == IDEAL_PUMP) )
        {
            if ( marked[j] == 0 ) marked[j] = 1;
        }
        else marked[j] = -1;
    }

    // --- find marked nodes with outgoing dummy links or ideal pumps
    for ( i = 0; i < sp->Nobjects[LINK]; i++ )
    {
        if ( (sp->Link[i].type == CONDUIT && sp->Link[i].xsect.type == DUMMY) ||
             (sp->Link[i].type == PUMP && 
              sp->Pump[sp->Link[i].subIndex].type == IDEAL_PUMP) )
        {
            j = sp->Link[i].node1;
            if ( marked[j] > 0 )
            {
                report_writeErrorMsg(sp, ERR_DUMMY_LINK, sp->Node[j].ID);
            }
        }
    }
    FREE(marked);
}

//=============================================================================
