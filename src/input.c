//-----------------------------------------------------------------------------
//   input.c
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14  (Build 5.1.001)
//             09/15/14  (Build 5.1.007)
//             08/01/16  (Build 5.1.011)
//   Author:   L. Rossman
//
//   Input data processing functions.
//
//   Build 5.1.007:
//   - Support added for climate adjustment input data.
//
//   Build 5.1.011:
//   - Support added for reading hydraulic event dates.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"
#include "lid.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
static const int MAXERRS = 100;        // Max. input errors reported

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
static char *Tok[MAXTOKS];             // String tokens from line of input
static int  Ntokens;                   // Number of tokens in line of input
static int  Mobjects[MAX_OBJ_TYPES];   // Working number of objects of each type
static int  Mnodes[MAX_NODE_TYPES];    // Working number of node objects
static int  Mlinks[MAX_LINK_TYPES];    // Working number of link objects
static int  Mevents;                   // Working number of event periods      //(5.1.011)

//-----------------------------------------------------------------------------
//  External Functions (declared in funcs.h)
//-----------------------------------------------------------------------------
//  input_countObjects  (called by swmm_open in swmm5.c)
//  input_readData      (called by swmm_open in swmm5.c)

//-----------------------------------------------------------------------------
//  Local functions
//-----------------------------------------------------------------------------
static int  addObject(SWMM_Project *sp, int objType, char* id);
static int  getTokens(char *s);
static int  parseLine(SWMM_Project *sp, int sect, char* line);
static int  readOption(SWMM_Project *sp, char* line);
static int  readTitle(SWMM_Project *sp, char* line);
static int  readControl(SWMM_Project *sp, char* tok[], int ntoks);
static int  readNode(SWMM_Project *sp, int type);
static int  readLink(SWMM_Project *sp, int type);
static int  readEvent(char* tok[], int ntoks);                                 //(5.1.011)

//=============================================================================

int input_countObjects(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads input file to determine number of system objects.
//
{
    char  line[MAXLINE+1];             // line from input data file     
    char  wLine[MAXLINE+1];            // working copy of input line   
    char  *tok;                        // first string token of line          
    int   sect = -1, newsect;          // input data sections          
    int   errcode = 0;                 // error code
    int   errsum = 0;                  // number of errors found                   
    int   i;
    long  lineCount = 0;

    // --- initialize number of objects & set default values
    if ( sp->ErrorCode ) return sp->ErrorCode;
    error_setInpError(0, "");
    for (i = 0; i < MAX_OBJ_TYPES; i++) sp->Nobjects[i] = 0;
    for (i = 0; i < MAX_NODE_TYPES; i++) sp->Nnodes[i] = 0;
    for (i = 0; i < MAX_LINK_TYPES; i++) sp->Nlinks[i] = 0;

    // --- make pass through data file counting number of each object
    while ( fgets(line, MAXLINE, sp->Finp.file) != NULL )
    {
        // --- skip blank lines & those beginning with a comment
        lineCount++;
        strcpy(wLine, line);           // make working copy of line
        tok = strtok(wLine, SEPSTR);   // get first text token on line
        if ( tok == NULL ) continue;
        if ( *tok == ';' ) continue;

        // --- check if line begins with a new section heading
        if ( *tok == '[' )
        {
            // --- look for heading in list of section keywords
            newsect = findmatch(tok, SectWords);
            if ( newsect >= 0 )
            {
                sect = newsect;
                continue;
            }
            else
            {
                sect = -1;
                errcode = ERR_KEYWORD;
            }
        }

        // --- if in OPTIONS section then read the option setting
        //     otherwise add object and its ID name (tok) to project
        if ( sect == s_OPTION ) errcode = readOption(sp, line);
        else if ( sect >= 0 )   errcode = addObject(sp, sect, tok);

        // --- report any error found
        if ( errcode )
        {
            report_writeInputErrorMsg(sp, errcode, sect, line, lineCount);
            errsum++;
            if (errsum >= MAXERRS ) break;
        }
    }

    // --- set global error code if input errors were found
    if ( errsum > 0 ) sp->ErrorCode = ERR_INPUT;
    return sp->ErrorCode;
}

//=============================================================================

int input_readData(SWMM_Project *sp)
//
//  Input:   none
//  Output:  returns error code
//  Purpose: reads input file to determine input parameters for each object.
//
{
    char  line[MAXLINE+1];        // line from input data file
    char  wLine[MAXLINE+1];       // working copy of input line
    char* comment;                // ptr. to start of comment in input line
    int   sect, newsect;          // data sections
    int   inperr, errsum;         // error code & total error count
    int   lineLength;             // number of characters in input line
    int   i;
    long  lineCount = 0;

    // --- initialize working item count arrays
    //     (final counts in Mobjects, Mnodes & Mlinks should
    //      match those in sp->Nobjects, sp->Nnodes and sp->Nlinks).
    if ( sp->ErrorCode ) return sp->ErrorCode;
    error_setInpError(0, "");
    for (i = 0; i < MAX_OBJ_TYPES; i++)  Mobjects[i] = 0;
    for (i = 0; i < MAX_NODE_TYPES; i++) Mnodes[i] = 0;
    for (i = 0; i < MAX_LINK_TYPES; i++) Mlinks[i] = 0;
    Mevents = 0;                                                               //(5.1.011)

    // --- initialize starting date for all time series
    for ( i = 0; i < sp->Nobjects[TSERIES]; i++ )
    {
        sp->Tseries[i].lastDate = sp->StartDate + sp->StartTime;
    }

    // --- read each line from input file
    sect = 0;
    errsum = 0;
    rewind(sp->Finp.file);
    while ( fgets(line, MAXLINE, sp->Finp.file) != NULL )
    {
        // --- make copy of line and scan for tokens
        lineCount++;
        strcpy(wLine, line);
        Ntokens = getTokens(wLine);

        // --- skip blank lines and comments
        if ( Ntokens == 0 ) continue;
        if ( *Tok[0] == ';' ) continue;

        // --- check if max. line length exceeded
        lineLength = strlen(line);
        if ( lineLength >= MAXLINE )
        {
            // --- don't count comment if present
            comment = strchr(line, ';');
            if ( comment ) lineLength = comment - line;    // Pointer math here
            if ( lineLength >= MAXLINE )
            {
                inperr = ERR_LINE_LENGTH;
                report_writeInputErrorMsg(sp, inperr, sect, line, lineCount);
                errsum++;
            }
        }

        // --- check if at start of a new input section
        if (*Tok[0] == '[')
        {
            // --- match token against list of section keywords
            newsect = findmatch(Tok[0], SectWords);
            if (newsect >= 0)
            {
                // --- SPECIAL CASE FOR TRANSECTS
                //     finish processing the last set of transect data
                if ( sect == s_TRANSECT )
                    transect_validate(sp, sp->Nobjects[TRANSECT]-1);

                // --- begin a new input section
                sect = newsect;
                continue;
            }
            else
            {
                inperr = error_setInpError(ERR_KEYWORD, Tok[0]);
                report_writeInputErrorMsg(sp, inperr, sect, line, lineCount);
                errsum++;
                break;
            }
        }

        // --- otherwise parse tokens from input line
        else
        {
            inperr = parseLine(sp, sect, line);
            if ( inperr > 0 )
            {
                errsum++;
                if ( errsum > MAXERRS ) report_writeLine(sp, FMT19);
                else report_writeInputErrorMsg(sp, inperr, sect, line, lineCount);
            }
        }

        // --- stop if reach end of file or max. error count
        if (errsum > MAXERRS) break;
    }   /* End of while */

    // --- check for errors
    if (errsum > 0)  sp->ErrorCode = ERR_INPUT;
    return sp->ErrorCode;
}

//=============================================================================

int  addObject(SWMM_Project *sp, int objType, char* id)
//
//  Input:   objType = object type index
//           id = object's ID string
//  Output:  returns an error code
//  Purpose: adds a new object to the project.
//
{
    int errcode = 0;
    switch( objType )
    {
      case s_RAINGAGE:
        if ( !project_addObject(GAGE, id, sp->Nobjects[GAGE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[GAGE]++;
        break;

      case s_SUBCATCH:
        if ( !project_addObject(SUBCATCH, id, sp->Nobjects[SUBCATCH]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[SUBCATCH]++;
        break;

      case s_AQUIFER:
        if ( !project_addObject(AQUIFER, id, sp->Nobjects[AQUIFER]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[AQUIFER]++;
        break;

      case s_UNITHYD:
        // --- the same Unit Hydrograph can span several lines
        if ( project_findObject(UNITHYD, id) < 0 )
        {
            if ( !project_addObject(UNITHYD, id, sp->Nobjects[UNITHYD]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            sp->Nobjects[UNITHYD]++;
        }
        break;

      case s_SNOWMELT:
        // --- the same Snowmelt object can appear on several lines
        if ( project_findObject(SNOWMELT, id) < 0 )
        {
            if ( !project_addObject(SNOWMELT, id, sp->Nobjects[SNOWMELT]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            sp->Nobjects[SNOWMELT]++;
        }
        break;

      case s_JUNCTION:
        if ( !project_addObject(NODE, id, sp->Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[NODE]++;
        sp->Nnodes[JUNCTION]++;
        break;

      case s_OUTFALL:
        if ( !project_addObject(NODE, id, sp->Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[NODE]++;
        sp->Nnodes[OUTFALL]++;
        break;

      case s_STORAGE:
        if ( !project_addObject(NODE, id, sp->Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[NODE]++;
        sp->Nnodes[STORAGE]++;
        break;

      case s_DIVIDER:
        if ( !project_addObject(NODE, id, sp->Nobjects[NODE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[NODE]++;
        sp->Nnodes[DIVIDER]++;
        break;

      case s_CONDUIT:
        if ( !project_addObject(LINK, id, sp->Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[LINK]++;
        sp->Nlinks[CONDUIT]++;
        break;

      case s_PUMP:
        if ( !project_addObject(LINK, id, sp->Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[LINK]++;
        sp->Nlinks[PUMP]++;
        break;

      case s_ORIFICE:
        if ( !project_addObject(LINK, id, sp->Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[LINK]++;
        sp->Nlinks[ORIFICE]++;
        break;

      case s_WEIR:
        if ( !project_addObject(LINK, id, sp->Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[LINK]++;
        sp->Nlinks[WEIR]++;
        break;

      case s_OUTLET:
        if ( !project_addObject(LINK, id, sp->Nobjects[LINK]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[LINK]++;
        sp->Nlinks[OUTLET]++;
        break;

      case s_POLLUTANT:
        if ( !project_addObject(POLLUT, id, sp->Nobjects[POLLUT]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[POLLUT]++;
        break;

      case s_LANDUSE:
        if ( !project_addObject(LANDUSE, id, sp->Nobjects[LANDUSE]) )
            errcode = error_setInpError(ERR_DUP_NAME, id);
        sp->Nobjects[LANDUSE]++;
        break;

      case s_PATTERN:
        // --- a time pattern can span several lines
        if ( project_findObject(TIMEPATTERN, id) < 0 )
        {
            if ( !project_addObject(TIMEPATTERN, id, sp->Nobjects[TIMEPATTERN]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            sp->Nobjects[TIMEPATTERN]++;
        }
        break;

      case s_CURVE:
        // --- a Curve can span several lines
        if ( project_findObject(CURVE, id) < 0 )
        {
            if ( !project_addObject(CURVE, id, sp->Nobjects[CURVE]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            sp->Nobjects[CURVE]++;

            // --- check for a conduit shape curve
            id = strtok(NULL, SEPSTR);
            if ( findmatch(id, CurveTypeWords) == SHAPE_CURVE )
                sp->Nobjects[SHAPE]++;
        }
        break;

      case s_TIMESERIES:
        // --- a Time Series can span several lines
        if ( project_findObject(TSERIES, id) < 0 )
        {
            if ( !project_addObject(TSERIES, id, sp->Nobjects[TSERIES]) )
                errcode = error_setInpError(ERR_DUP_NAME, id);
            sp->Nobjects[TSERIES]++;
        }
        break;

      case s_CONTROL:
        if ( match(id, w_RULE) ) sp->Nobjects[CONTROL]++;
        break;

      case s_TRANSECT:
        // --- for TRANSECTS, ID name appears as second entry on X1 line
        if ( match(id, "X1") )
        {
            id = strtok(NULL, SEPSTR);
            if ( id ) 
            {
                if ( !project_addObject(TRANSECT, id, sp->Nobjects[TRANSECT]) )
                    errcode = error_setInpError(ERR_DUP_NAME, id);
                sp->Nobjects[TRANSECT]++;
            }
        }
        break;

      case s_LID_CONTROL:
        // --- an LID object can span several lines
        if ( project_findObject(LID, id) < 0 )
        {
            if ( !project_addObject(LID, id, sp->Nobjects[LID]) )
            {
                errcode = error_setInpError(ERR_DUP_NAME, id);
            }
            sp->Nobjects[LID]++;
        }
        break;

      case s_EVENT: sp->NumEvents++; break;                                        //(5.1.011)
    }
    return errcode;
}

//=============================================================================

int  parseLine(SWMM_Project *sp, int sect, char *line)
//
//  Input:   sect  = current section of input file
//           *line = line of text read from input file
//  Output:  returns error code or 0 if no error found
//  Purpose: parses contents of a tokenized line of text read from input file.
//
{
    int j, err;
    switch (sect)
    {
      case s_TITLE:
        return readTitle(sp, line);

      case s_RAINGAGE:
        j = Mobjects[GAGE];
        err = gage_readParams(sp, j, Tok, Ntokens);
        Mobjects[GAGE]++;
        return err;

      case s_TEMP:
        return climate_readParams(sp, Tok, Ntokens);

      case s_EVAP:
        return climate_readEvapParams(sp, Tok, Ntokens);

      case s_ADJUST:                                                           //(5.1.007)
        return climate_readAdjustments(sp, Tok, Ntokens);                      //(5.1.007)

      case s_SUBCATCH:
        j = Mobjects[SUBCATCH];
        err = subcatch_readParams(sp, j, Tok, Ntokens);
        Mobjects[SUBCATCH]++;
        return err;

      case s_SUBAREA:
        return subcatch_readSubareaParams(sp, Tok, Ntokens);

      case s_INFIL:
        return infil_readParams(sp, sp->InfilModel, Tok, Ntokens);

      case s_AQUIFER:
        j = Mobjects[AQUIFER];
        err = gwater_readAquiferParams(sp, j, Tok, Ntokens);
        Mobjects[AQUIFER]++;
        return err;

      case s_GROUNDWATER:
        return gwater_readGroundwaterParams(sp, Tok, Ntokens);

      case s_GWF:
        return gwater_readFlowExpression(sp, Tok, Ntokens);

      case s_SNOWMELT:
        return snow_readMeltParams(sp, Tok, Ntokens);

      case s_JUNCTION:
        return readNode(sp, JUNCTION);

      case s_OUTFALL:
        return readNode(sp, OUTFALL);

      case s_STORAGE:
        return readNode(sp, STORAGE);

      case s_DIVIDER:
        return readNode(sp, DIVIDER);

      case s_CONDUIT:
        return readLink(sp, CONDUIT);

      case s_PUMP:
        return readLink(sp, PUMP);

      case s_ORIFICE:
        return readLink(sp, ORIFICE);

      case s_WEIR:
        return readLink(sp, WEIR);

      case s_OUTLET:
        return readLink(sp, OUTLET);

      case s_XSECTION:
        return link_readXsectParams(sp, Tok, Ntokens);

      case s_TRANSECT:
        return transect_readParams(sp, &Mobjects[TRANSECT], Tok, Ntokens);

      case s_LOSSES:
        return link_readLossParams(sp, Tok, Ntokens);

      case s_POLLUTANT:
        j = Mobjects[POLLUT];
        err = landuse_readPollutParams(sp, j, Tok, Ntokens);
        Mobjects[POLLUT]++;
        return err;

      case s_LANDUSE:
        j = Mobjects[LANDUSE];
        err = landuse_readParams(sp, j, Tok, Ntokens);
        Mobjects[LANDUSE]++;
        return err;

      case s_BUILDUP:
        return landuse_readBuildupParams(sp, Tok, Ntokens);

      case s_WASHOFF:
        return landuse_readWashoffParams(sp, Tok, Ntokens);

      case s_COVERAGE:
        return subcatch_readLanduseParams(sp, Tok, Ntokens);

      case s_INFLOW:
        return inflow_readExtInflow(sp, Tok, Ntokens);

      case s_DWF:
        return inflow_readDwfInflow(sp, Tok, Ntokens);

      case s_PATTERN:
        return inflow_readDwfPattern(sp, Tok, Ntokens);

      case s_RDII:
        return rdii_readRdiiInflow(sp, Tok, Ntokens);

      case s_UNITHYD:
        return rdii_readUnitHydParams(sp, Tok, Ntokens);

      case s_LOADING:
        return subcatch_readInitBuildup(sp, Tok, Ntokens);

      case s_TREATMENT:
        return treatmnt_readExpression(sp, Tok, Ntokens);

      case s_CURVE:
        return table_readCurve(sp, Tok, Ntokens);

      case s_TIMESERIES:
        return table_readTimeseries(sp, Tok, Ntokens);

      case s_CONTROL:
        return readControl(sp, Tok, Ntokens);

      case s_REPORT:
        return report_readOptions(sp, Tok, Ntokens);

      case s_FILE:
        return iface_readFileParams(sp, Tok, Ntokens);

      case s_LID_CONTROL:
        return lid_readProcParams(sp, Tok, Ntokens);

      case s_LID_USAGE:
        return lid_readGroupParams(sp, Tok, Ntokens);

      case s_EVENT:
        return readEvent(Tok, Ntokens);                                        //(5.1.011)

      default: return 0;
    }
}

//=============================================================================

int readControl(SWMM_Project *sp, char* tok[], int ntoks)
//
//  Input:   tok[] = array of string tokens
//           ntoks = number of tokens
//  Output:  returns error code
//  Purpose: reads a line of input for a control rule.
//
{
    int index;
    int keyword;

    // --- check for minimum number of tokens
    if ( ntoks < 2 ) return error_setInpError(ERR_ITEMS, "");

    // --- get index of control rule keyword
    keyword = findmatch(tok[0], RuleKeyWords);
    if ( keyword < 0 ) return error_setInpError(ERR_KEYWORD, tok[0]);

    // --- if line begins a new control rule, add rule ID to the database
    if ( keyword == 0 )
    {
        if ( !project_addObject(CONTROL, tok[1], Mobjects[CONTROL]) )
        {
            return error_setInpError(ERR_DUP_NAME, Tok[1]);
        }
        Mobjects[CONTROL]++;
    }

    // --- get index of last control rule processed
    index = Mobjects[CONTROL] - 1;
    if ( index < 0 ) return error_setInpError(ERR_RULE, "");

    // --- add current line as a new clause to the control rule
    return controls_addRuleClause(sp, index, keyword, Tok, Ntokens);
}

//=============================================================================

int readOption(SWMM_Project *sp, char* line)
//
//  Input:   line = line of input data
//  Output:  returns error code
//  Purpose: reads an input line containing a project option.
//
{
    Ntokens = getTokens(line);
    if ( Ntokens < 2 ) return 0;
    return project_readOption(sp, Tok[0], Tok[1]);
}

//=============================================================================

int readTitle(SWMM_Project *sp, char* line)
//
//  Input:   line = line from input file
//  Output:  returns error code
//  Purpose: reads project title from line of input.
//
{
    int i, n;
    for (i = 0; i < MAXTITLE; i++)
    {
        // --- find next empty Title entry
        if ( strlen(sp->Title[i]) == 0 )
        {
            // --- strip line feed character from input line
            n = strlen(line);
            if (line[n-1] == 10) line[n-1] = ' ';

            // --- copy input line into Title entry
            sstrncpy(sp->Title[i], line, MAXMSG);
            break;
        }
    }
    return 0;
}

//=============================================================================
    
int readNode(SWMM_Project *sp, int type)
//
//  Input:   type = type of node
//  Output:  returns error code
//  Purpose: reads data for a node from a line of input.
//
{
    int j = Mobjects[NODE];
    int k = Mnodes[type];
    int err = node_readParams(sp, j, type, k, Tok, Ntokens);
    Mobjects[NODE]++;
    Mnodes[type]++;
    return err;
}

//=============================================================================

int readLink(SWMM_Project *sp, int type)
//
//  Input:   type = type of link
//  Output:  returns error code
//  Purpose: reads data for a link from a line of input.
//
{
    int j = Mobjects[LINK];
    int k = Mlinks[type];
    int err = link_readParams(sp, j, type, k, Tok, Ntokens);
    Mobjects[LINK]++;
    Mlinks[type]++;
    return err;
}

//=============================================================================

////  This function was added to release 5.1.011.  ////                        //(5.1.011)

int  readEvent(char* tok[], int ntoks)
{
    DateTime x[4];

    if ( ntoks < 4 ) return error_setInpError(ERR_ITEMS, "");
    if ( !datetime_strToDate(tok[0], &x[0]) )
        return error_setInpError(ERR_DATETIME, tok[0]);
    if ( !datetime_strToTime(tok[1], &x[1]) )
        return error_setInpError(ERR_DATETIME, tok[1]);
    if ( !datetime_strToDate(tok[2], &x[2]) )
        return error_setInpError(ERR_DATETIME, tok[2]);
    if ( !datetime_strToTime(tok[3], &x[3]) )
        return error_setInpError(ERR_DATETIME, tok[3]);

    Event[Mevents].start = x[0] + x[1];
    Event[Mevents].end = x[2] + x[3];
    if ( Event[Mevents].start >= Event[Mevents].end )
       return error_setInpError(ERR_DATETIME, " - start date exceeds end date");
    Mevents++;
    return 0;
}

//=============================================================================

int  findmatch(char *s, char *keyword[])
//
//  Input:   s = character string
//           keyword = array of keyword strings
//  Output:  returns index of matching keyword or -1 if no match found  
//  Purpose: finds match between string and array of keyword strings.
//
{
   int i = 0;
   while (keyword[i] != NULL)
   {
      if (match(s, keyword[i])) return(i);
      i++;
   }
   return(-1);
}

//=============================================================================

int  match(char *str, char *substr)
//
//  Input:   str = character string being searched
//           substr = sub-string being searched for
//  Output:  returns 1 if sub-string found, 0 if not
//  Purpose: sees if a sub-string of characters appears in a string
//           (not case sensitive).
//
{
    int i,j;

    // --- fail if substring is empty
    if (!substr[0]) return(0);

    // --- skip leading blanks of str
    for (i = 0; str[i]; i++)
    {
        if (str[i] != ' ') break;
    }

    // --- check if substr matches remainder of str
    for (i = i,j = 0; substr[j]; i++,j++)
    {
        if (!str[i] || UCHAR(str[i]) != UCHAR(substr[j])) return(0);
    }
    return(1);
}

//=============================================================================

int  getInt(char *s, int *y)
//
//  Input:   s = a character string
//  Output:  y = converted value of s,
//           returns 1 if conversion successful, 0 if not
//  Purpose: converts a string to an integer number.
//
{
    double x;
    if ( getDouble(s, &x) )
    {
        if ( x < 0.0 ) x -= 0.01;
        else x += 0.01;
        *y = (int)x;
        return 1;
    }
    *y = 0;
    return 0;
}

//=============================================================================

int  getFloat(char *s, float *y)
//
//  Input:   s = a character string
//  Output:  y = converted value of s,
//           returns 1 if conversion successful, 0 if not
//  Purpose: converts a string to a single precision floating point number.
//
{
    char *endptr;
    *y = (float) strtod(s, &endptr);
    if (*endptr > 0) return(0);
    return(1);
}

//=============================================================================

int  getDouble(char *s, double *y)
//
//  Input:   s = a character string
//  Output:  y = converted value of s,
//           returns 1 if conversion successful, 0 if not
//  Purpose: converts a string to a double precision floating point number.
//
{
    char *endptr;
    *y = strtod(s, &endptr);
    if (*endptr > 0) return(0);
    return(1);
}

//=============================================================================

int  getTokens(char *s)
//
//  Input:   s = a character string
//  Output:  returns number of tokens found in s
//  Purpose: scans a string for tokens, saving pointers to them
//           in shared variable Tok[].
//
//  Notes:   Tokens can be separated by the characters listed in SEPSTR
//           (spaces, tabs, newline, carriage return) which is defined
//           in CONSTS.H. Text between quotes is treated as a single token.
//
{
    int  len, m, n;
    char *c;

    // --- begin with no tokens
    for (n = 0; n < MAXTOKS; n++) Tok[n] = NULL;
    n = 0;

    // --- truncate s at start of comment 
    c = strchr(s,';');
    if (c) *c = '\0';
    len = strlen(s);

    // --- scan s for tokens until nothing left
    while (len > 0 && n < MAXTOKS)
    {
        m = strcspn(s,SEPSTR);              // find token length 
        if (m == 0) s++;                    // no token found
        else
        {
            if (*s == '"')                  // token begins with quote
            {
                s++;                        // start token after quote
                len--;                      // reduce length of s
                m = strcspn(s,"\"\n");      // find end quote or new line
            }
            s[m] = '\0';                    // null-terminate the token
            Tok[n] = s;                     // save pointer to token 
            n++;                            // update token count
            s += m+1;                       // begin next token
        }
        len -= m+1;                         // update length of s
    }
    return(n);
}

//=============================================================================
