//-----------------------------------------------------------------------------
//   iface.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14   (Build 5.1.001)
//   Author:   L. Rossman
//
//   Routing interface file functions.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include "headers.h"

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
//#ifdef __cplusplus
//extern const double Qcf[];             // flow units conversion factors
                                       // (see swmm5.c)
//#else
extern const double Qcf[];                   // flow units conversion factors
                                       // (see swmm5.c)
//#endif

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  iface_readFileParams     (called by input_readLine)
//  iface_openRoutingFiles   (called by routing_open)
//  iface_closeRoutingFiles  (called by routing_close)
//  iface_getNumIfaceNodes   (called by addIfaceInflows in routing.c)
//  iface_getIfaceNode       (called by addIfaceInflows in routing.c)
//  iface_getIfaceFlow       (called by addIfaceInflows in routing.c)
//  iface_getIfaceQual       (called by addIfaceInflows in routing.c)
//  iface_saveOutletResults  (called by output_saveResults)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static void  openFileForOutput(SWMM_Project *sp);
static void  openFileForInput(SWMM_Project *sp);
static int   getIfaceFilePolluts(SWMM_Project *sp);
static int   getIfaceFileNodes(SWMM_Project *sp);
static void  setOldIfaceValues(SWMM_Project *sp);
static void  readNewIfaceValues(SWMM_Project *sp);
static int   isOutletNode(SWMM_Project *sp, int node);

//=============================================================================

int iface_readFileParams(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns an error code
//  Purpose: reads interface file information from a line of input data.
//
//  Data format is:
//  USE/SAVE  FileType  FileName
//
{
    char  k;
    int   j;

    // --- determine file disposition and type
    if ( ntoks < 2 ) return error_setInpError(sp, ERR_ITEMS, "");
    k = (char)findmatch(tok[0], FileModeWords);
    if ( k < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[0]);
    j = findmatch(tok[1], FileTypeWords);
    if ( j < 0 ) return error_setInpError(sp, ERR_KEYWORD, tok[1]);
    if ( ntoks < 3 ) return 0;

    // --- process file name
    switch ( j )
    {
      case RAINFALL_FILE:
        sp->Frain.mode = k;
        sstrncpy(sp->Frain.name, tok[2], MAXFNAME);
        break;

      case RUNOFF_FILE:
        sp->Frunoff.mode = k;
        sstrncpy(sp->Frunoff.name, tok[2], MAXFNAME);
        break;

      case HOTSTART_FILE:
        if ( k == USE_FILE )
        {
            sp->Fhotstart1.mode = k;
            sstrncpy(sp->Fhotstart1.name, tok[2], MAXFNAME);
        }
        else if ( k == SAVE_FILE )
        {
            sp->Fhotstart2.mode = k;
            sstrncpy(sp->Fhotstart2.name, tok[2], MAXFNAME);
        }
        break;

      case RDII_FILE:
        sp->Frdii.mode = k;
        sstrncpy(sp->Frdii.name, tok[2], MAXFNAME);
        break;

      case INFLOWS_FILE:
        if ( k != USE_FILE ) return error_setInpError(sp, ERR_ITEMS, "");
        sp->Finflows.mode = k;
        sstrncpy(sp->Finflows.name, tok[2], MAXFNAME);
        break;

      case OUTFLOWS_FILE:
        if ( k != SAVE_FILE ) return error_setInpError(sp, ERR_ITEMS, "");
        sp->Foutflows.mode = k;
        sstrncpy(sp->Foutflows.name, tok[2], MAXFNAME);
        break;
    }
    return 0;
}

//=============================================================================

void iface_openRoutingFiles(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: opens routing interface files.
//
{
    TIfaceShared *ifc = &sp->IfaceShared;

    // --- initialize shared variables
    ifc->NumIfacePolluts = 0;
    ifc->IfacePolluts = NULL;
    ifc->NumIfaceNodes = 0;
    ifc->IfaceNodes = NULL;
    ifc->OldIfaceValues = NULL;
    ifc->NewIfaceValues = NULL;

    // --- check that inflows & outflows files are not the same
    if ( sp->Foutflows.mode != NO_FILE && sp->Finflows.mode != NO_FILE )
    {
        if ( strcomp(sp->Foutflows.name, sp->Finflows.name) )
        {
            report_writeErrorMsg(sp, ERR_ROUTING_FILE_NAMES, "");
            return;
        }
    }

    // --- open the file for reading or writing
    if ( sp->Foutflows.mode == SAVE_FILE ) openFileForOutput(sp);
    if ( sp->Finflows.mode == USE_FILE ) openFileForInput(sp);
}

//=============================================================================

void iface_closeRoutingFiles(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: closes routing interface files.
//
{
    TIfaceShared *ifc = &sp->IfaceShared;

    FREE(ifc->IfacePolluts);
    FREE(ifc->IfaceNodes);
    if ( ifc->OldIfaceValues != NULL ) project_freeMatrix(ifc->OldIfaceValues);
    if ( ifc->NewIfaceValues != NULL ) project_freeMatrix(ifc->NewIfaceValues);
    if ( sp->Finflows.file )  fclose(sp->Finflows.file);
    if ( sp->Foutflows.file ) fclose(sp->Foutflows.file);
}

//=============================================================================

int iface_getNumIfaceNodes(SWMM_Project *sp, DateTime currentDate)
//
//  Input:   currentDate = current date/time
//  Output:  returns number of interface nodes if data exists or
//           0 otherwise
//  Purpose: reads inflow data from interface file at current date.
//
{
    TIfaceShared *ifc = &sp->IfaceShared;

    // --- return 0 if file begins after current date
    if ( ifc->OldIfaceDate > currentDate ) return 0;

    // --- keep updating new interface values until current date bracketed
    while ( ifc->NewIfaceDate < currentDate && ifc->NewIfaceDate != NO_DATE )
    {
        setOldIfaceValues(sp);
        readNewIfaceValues(sp);
    }

    // --- return 0 if no data available
    if ( ifc->NewIfaceDate == NO_DATE ) return 0;

    // --- find fraction current date is bewteen old & new interface dates
    ifc->IfaceFrac = (currentDate - ifc->OldIfaceDate) /
            (ifc->NewIfaceDate - ifc->OldIfaceDate);
    ifc->IfaceFrac = MAX(0.0, ifc->IfaceFrac);
    ifc->IfaceFrac = MIN(ifc->IfaceFrac, 1.0);

    // --- return number of interface nodes
    return ifc->NumIfaceNodes;
}

//=============================================================================

int iface_getIfaceNode(SWMM_Project *sp, int index)
//
//  Input:   index = interface file node index
//  Output:  returns project node index
//  Purpose: finds index of project node associated with interface node index
//
{
    TIfaceShared *ifc = &sp->IfaceShared;

    if ( index >= 0 && index < ifc->NumIfaceNodes )
        return ifc->IfaceNodes[index];
    else return -1;
}

//=============================================================================

double iface_getIfaceFlow(SWMM_Project *sp, int index)
//
//  Input:   index = interface file node index
//  Output:  returns inflow to node
//  Purpose: finds interface flow for particular node index.
//
{
    double q1, q2;

    TIfaceShared *ifc = &sp->IfaceShared;

    if ( index >= 0 && index < ifc->NumIfaceNodes )
    {
        // --- interpolate flow between old and new values
        q1 = ifc->OldIfaceValues[index][0];
        q2 = ifc->NewIfaceValues[index][0];
        return (1.0 - ifc->IfaceFrac)*q1 + ifc->IfaceFrac*q2;
    }
    else return 0.0;
}

//=============================================================================

double iface_getIfaceQual(SWMM_Project *sp, int index, int pollut)
//
//  Input:   index = index of node on interface file
//           pollut = index of pollutant on interface file
//  Output:  returns inflow pollutant concentration
//  Purpose: finds interface concentration for particular node index & pollutant.
//
{
    int    j;
    double c1, c2;

    TIfaceShared *ifc = &sp->IfaceShared;

    if ( index >= 0 && index < ifc->NumIfaceNodes )
    {
        // --- find index of pollut on interface file
        j = ifc->IfacePolluts[pollut];
        if ( j < 0 ) return 0.0;

        // --- interpolate flow between old and new values
        //     (remember that 1st col. of values matrix is for flow)
        c1 = ifc->OldIfaceValues[index][j+1];
        c2 = ifc->NewIfaceValues[index][j+1];
        return (1.0 - ifc->IfaceFrac)*c1 + ifc->IfaceFrac*c2;
    }
    else return 0.0;
}

//=============================================================================

void iface_saveOutletResults(SWMM_Project *sp, DateTime reportDate, FILE* file)
//
//  Input:   reportDate = reporting date/time
//           file = ptr. to interface file
//  Output:  none
//  Purpose: saves system outflows to routing interface file.
//
{
    int i, p, yr, mon, day, hr, min, sec;
    char theDate[25];
    datetime_decodeDate(reportDate, &yr, &mon, &day);
    datetime_decodeTime(reportDate, &hr, &min, &sec);
    sprintf(theDate, " %04d %02d  %02d  %02d  %02d  %02d ",
            yr, mon, day, hr, min, sec);
    for (i=0; i<sp->Nobjects[NODE]; i++)
    {
        // --- check that node is an outlet node
        if ( !isOutletNode(sp, i) ) continue;

        // --- write node ID, date, flow, and quality to file
        fprintf(file, "\n%-16s", sp->Node[i].ID);
        fprintf(file, "%s", theDate);
        fprintf(file, " %-10f", sp->Node[i].inflow * UCF(sp, FLOW));
        for ( p = 0; p < sp->Nobjects[POLLUT]; p++ )
        {
            fprintf(file, " %-10f", sp->Node[i].newQual[p]);
        }
    }
}

//=============================================================================

void openFileForOutput(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: opens a routing interface file for writing.
//
{
    int i, n;

    // --- open the routing file for writing text
    sp->Foutflows.file = fopen(sp->Foutflows.name, "wt");
    if ( sp->Foutflows.file == NULL )
    {
        report_writeErrorMsg(sp, ERR_ROUTING_FILE_OPEN, sp->Foutflows.name);
        return;
    }

    // --- write title & reporting time step to file
    fprintf(sp->Foutflows.file, "SWMM5 Interface File");
    fprintf(sp->Foutflows.file, "\n%s", sp->Title[0]);
    fprintf(sp->Foutflows.file, "\n%-4d - reporting time step in sec", sp->ReportStep);

    // --- write number & names of each constituent (including flow) to file
    fprintf(sp->Foutflows.file, "\n%-4d - number of constituents as listed below:",
            sp->Nobjects[POLLUT] + 1);
    fprintf(sp->Foutflows.file, "\nFLOW %s", FlowUnitWords[sp->FlowUnits]);
    for (i=0; i<sp->Nobjects[POLLUT]; i++)
    {
        fprintf(sp->Foutflows.file, "\n%s %s", sp->Pollut[i].ID,
            QualUnitsWords[sp->Pollut[i].units]);
    }

    // --- count number of outlet nodes
    n = 0;
    for (i=0; i<sp->Nobjects[NODE]; i++)
    {
        if ( isOutletNode(sp, i) ) n++;
    }

    // --- write number and names of outlet nodes to file
    fprintf(sp->Foutflows.file, "\n%-4d - number of nodes as listed below:", n);
    for (i=0; i<sp->Nobjects[NODE]; i++)
    {
          if ( isOutletNode(sp, i) )
            fprintf(sp->Foutflows.file, "\n%s", sp->Node[i].ID);
    }

    // --- write column headings
    fprintf(sp->Foutflows.file,
        "\nNode             Year Mon Day Hr  Min Sec FLOW      ");
    for (i=0; i<sp->Nobjects[POLLUT]; i++)
    {
        fprintf(sp->Foutflows.file, " %-10s", sp->Pollut[i].ID);
    }

    // --- if reporting starts immediately, save initial outlet values
    if ( sp->ReportStart == sp->StartDateTime )
    {
        iface_saveOutletResults(sp, sp->ReportStart, sp->Foutflows.file);
    }
}

//=============================================================================

void openFileForInput(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: opens a routing interface file for reading.
//
{
    int   err;                         // error code
    char  line[MAXLINE+1];             // line from Routing interface file
    char  s[MAXLINE+1];                // general string variable

    TIfaceShared *ifc = &sp->IfaceShared;

    // --- open the routing interface file for reading text
    sp->Finflows.file = fopen(sp->Finflows.name, "rt");
    if ( sp->Finflows.file == NULL )
    {
        report_writeErrorMsg(sp, ERR_ROUTING_FILE_OPEN, sp->Finflows.name);
        return;
    }

    // --- check for correct file type
    fgets(line, MAXLINE, sp->Finflows.file);
    sscanf(line, "%s", s);
    if ( !strcomp(s, "SWMM5") )
    {
        report_writeErrorMsg(sp, ERR_ROUTING_FILE_FORMAT, sp->Finflows.name);
        return;
    }

    // --- skip title line
    fgets(line, MAXLINE, sp->Finflows.file);

    // --- read reporting time step (sec)
    ifc->IfaceStep = 0;
    fgets(line, MAXLINE, sp->Finflows.file);
    sscanf(line, "%d", &ifc->IfaceStep);
    if ( ifc->IfaceStep <= 0 )
    {
        report_writeErrorMsg(sp, ERR_ROUTING_FILE_FORMAT, sp->Finflows.name);
        return;
    }

    // --- match constituents in file with those in project
    err = getIfaceFilePolluts(sp);
    if ( err > 0 )
    {
        report_writeErrorMsg(sp, err, sp->Finflows.name);
        return;
    }

    // --- match nodes in file with those in project
    err = getIfaceFileNodes(sp);
    if ( err > 0 )
    {
        report_writeErrorMsg(sp, err, sp->Finflows.name);
        return;
    }

    // --- create matrices for old & new interface flows & WQ values
    ifc->OldIfaceValues =
            project_createMatrix(ifc->NumIfaceNodes, 1 + ifc->NumIfacePolluts);
    ifc->NewIfaceValues =
            project_createMatrix(ifc->NumIfaceNodes, 1 + ifc->NumIfacePolluts);
    if ( ifc->OldIfaceValues == NULL || ifc->NewIfaceValues == NULL )
    {
        report_writeErrorMsg(sp, ERR_MEMORY, "");
        return;
    }

    // --- read in new interface flows & WQ values
    readNewIfaceValues(sp);
    ifc->OldIfaceDate = ifc->NewIfaceDate;
}

//=============================================================================

int  getIfaceFilePolluts(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: reads names of pollutants saved on the inflows interface file.
//
{
    int   i, j;
    char  line[MAXLINE+1];             // line from inflows interface file
    char  s1[MAXLINE+1];               // general string variable
    char  s2[MAXLINE+1];         

    TIfaceShared *ifc = &sp->IfaceShared;

    // --- read number of pollutants (minus FLOW)
    fgets(line, MAXLINE, sp->Finflows.file);
    sscanf(line, "%d", &ifc->NumIfacePolluts);
    ifc->NumIfacePolluts--;
    if ( ifc->NumIfacePolluts < 0 ) return ERR_ROUTING_FILE_FORMAT;

    // --- read flow units
    fgets(line, MAXLINE, sp->Finflows.file);
    sscanf(line, "%s %s", s1, s2);
    if ( !strcomp(s1, "FLOW") )
        return ERR_ROUTING_FILE_FORMAT;
    ifc->IfaceFlowUnits = findmatch(s2, FlowUnitWords);
    if ( ifc->IfaceFlowUnits < 0 )
        return ERR_ROUTING_FILE_FORMAT;

    // --- allocate memory for pollutant index array
    if ( sp->Nobjects[POLLUT] > 0 )
    {
        ifc->IfacePolluts = (int *) calloc(sp->Nobjects[POLLUT], sizeof(int));
        if ( !ifc->IfacePolluts )
            return ERR_MEMORY;
        for (i=0; i<sp->Nobjects[POLLUT]; i++)
            ifc->IfacePolluts[i] = -1;
    }

    // --- read pollutant names & units
    if ( ifc->NumIfacePolluts > 0 )
    {
        // --- check each pollutant name on file with project's pollutants
        for (i=0; i < ifc->NumIfacePolluts; i++)
        {
            if ( feof(sp->Finflows.file) ) return ERR_ROUTING_FILE_FORMAT;
            fgets(line, MAXLINE, sp->Finflows.file);
            sscanf(line, "%s %s", s1, s2);
            if ( sp->Nobjects[POLLUT] > 0 )
            {
                j = project_findObject(sp, POLLUT, s1);
                if ( j < 0 ) continue;
                if ( !strcomp(s2, QualUnitsWords[sp->Pollut[j].units]) )
                    return ERR_ROUTING_FILE_NOMATCH;
                ifc->IfacePolluts[j] = i;
            }
        }
    }
    return 0;
}

//=============================================================================

int getIfaceFileNodes(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns an error code
//  Purpose: reads names of nodes contained on inflows interface file.
//
{
    int   i;
    char  line[MAXLINE+1];             // line from inflows interface file
    char  s[MAXLINE+1];                // general string variable

    TIfaceShared *ifc = &sp->IfaceShared;

    // --- read number of interface nodes
    fgets(line, MAXLINE, sp->Finflows.file);
    sscanf(line, "%d", &ifc->NumIfaceNodes);
    if ( ifc->NumIfaceNodes <= 0 ) return ERR_ROUTING_FILE_FORMAT;

    // --- allocate memory for interface nodes index array
    ifc->IfaceNodes = (int *) calloc(ifc->NumIfaceNodes, sizeof(int));
    if ( !ifc->IfaceNodes ) return ERR_MEMORY;

    // --- read names of interface nodes from file & save their indexes
    for ( i=0; i < ifc->NumIfaceNodes; i++ )
    {
        if ( feof(sp->Finflows.file) ) return ERR_ROUTING_FILE_FORMAT;
        fgets(line, MAXLINE, sp->Finflows.file);
        sscanf(line, "%s", s);
        ifc->IfaceNodes[i] = project_findObject(sp, NODE, s);
    }

    // --- skip over column headings line
    if ( feof(sp->Finflows.file) ) return ERR_ROUTING_FILE_FORMAT;
    fgets(line, MAXLINE, sp->Finflows.file);
    return 0;
}

//=============================================================================

void readNewIfaceValues(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: reads data from inflows interface file for next date.
//
{
    int    i, j;
    char*  s;
    int    yr = 0, mon = 0, day = 0,
		   hr = 0, min = 0, sec = 0;   // year, month, day, hour, minute, second
    char   line[MAXLINE+1];            // line from interface file

    TIfaceShared *ifc = &sp->IfaceShared;

    // --- read a line for each interface node
    ifc->NewIfaceDate = NO_DATE;
    for (i = 0; i < ifc->NumIfaceNodes; i++)
    {
        if ( feof(sp->Finflows.file) )
            return;

        fgets(line, MAXLINE, sp->Finflows.file);

        // --- parse date & time from line
        if ( strtok(line, SEPSTR) == NULL ) return;
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        yr  = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        mon = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        day = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        hr  = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        min = atoi(s);
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        sec = atoi(s);

        // --- parse flow value
        s = strtok(NULL, SEPSTR);
        if ( s == NULL ) return;
        ifc->NewIfaceValues[i][0] = atof(s) / Qcf[ifc->IfaceFlowUnits];

        // --- parse pollutant values
        for (j = 1; j <= ifc->NumIfacePolluts; j++)
        {
            s = strtok(NULL, SEPSTR);
            if ( s == NULL ) return;
            ifc->NewIfaceValues[i][j] = atof(s);
        }

    }

    // --- encode date & time values
    ifc->NewIfaceDate = datetime_encodeDate(yr, mon, day) +
                   datetime_encodeTime(hr, min, sec);
}

//=============================================================================

void setOldIfaceValues(SWMM_Project *sp)
//
//  Input:   none
//  Output:  none
//  Purpose: replaces old values read from routing interface file with new ones. 
//
{
    int i, j;

    TIfaceShared *ifc = &sp->IfaceShared;

    ifc->OldIfaceDate = ifc->NewIfaceDate;
    for ( i = 0; i < ifc->NumIfaceNodes; i++)
    {
        for ( j = 0; j < ifc->NumIfacePolluts+1; j++ )
        {
            ifc->OldIfaceValues[i][j] = ifc->NewIfaceValues[i][j];
        }
    }
}

//=============================================================================

int  isOutletNode(SWMM_Project *sp, int i)
//
//  Input:   i = node index
//  Output:  returns 1 if node is an outlet, 0 if not.
//  Purpose: determines if a node is an outlet point or not.
//
{
    // --- for DW routing only outfalls are outlets
    if ( sp->RouteModel == DW )
    {
        return (sp->Node[i].type == OUTFALL);
    }

    // --- otherwise outlets are nodes with no outflow links (degree is 0)
    else return (sp->Node[i].degree == 0);
}
