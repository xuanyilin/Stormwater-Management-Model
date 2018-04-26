//-----------------------------------------------------------------------------
//   toolkitAPI.c
//
//   Project: EPA SWMM5
//   Version: 5.1
//   Date:    08/30/2016
//   Author:  B. McDonnell (EmNet LLC)
//            K. Ratliff
//
//   Exportable Functions for Project Definition API.
//
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "headers.h"
#include "swmm5.h"                     // declaration of exportable functions
#include "hash.h"


// Function Declarations for API
int     massbal_getRoutingFlowTotal(SM_RoutingTotals *routingTot);
int     massbal_getRunoffTotal(SM_RunoffTotals *runoffTot);
double  massbal_getTotalArea(void);
int     massbal_getNodeTotalInflow(int index, double *value);

int  stats_getNodeStat(SWMM_Project *sp, int index, SM_NodeStats *nodeStats);
int  stats_getStorageStat(SWMM_Project *sp, int index, SM_StorageStats *storageStats);
int  stats_getOutfallStat(SWMM_Project *sp, int index, SM_OutfallStats *outfallStats);
int  stats_getLinkStat(SWMM_Project *sp, int index, SM_LinkStats *linkStats);
int  stats_getPumpStat(SWMM_Project *sp, int index, SM_PumpStats *pumpStats);
int  stats_getSubcatchStat(SWMM_Project *sp, int index, SM_SubcatchStats *subcatchStats);

//-----------------------------------------------------------------------------
//  Extended API Functions
//-----------------------------------------------------------------------------
void DLLEXPORT swmm_getAPIError(int errcode, char *s)
//
// Input:   errcode = error code
// Output:  errmessage String
// Return:  API Error
// Purpose: Get an error message
{
    char *errmsg = error_getMsg(errcode);
    strcpy(s, errmsg);
}

int DLLEXPORT swmm_getSimulationDateTime(int timetype, int *year, int *month,
        int *day, int *hours, int *minutes, int *seconds)
{
    return swmm_getSimulationDateTime_project(_defaultProject, timetype, year,
            month, day, hours, minutes, seconds);
}

int DLLEXPORT swmm_getSimulationDateTime_project(SWMM_Project *sp, int timetype,
        int *year, int *month, int *day, int *hours, int *minutes, int *seconds)
//
// Input:   timetype = time type to return
// Output:  year, month, day, hours, minutes, seconds = int
// Return:  API Error
// Purpose: Get the simulation start, end and report date times
{
    int errcode = 0;
    // Check if Open
    if (swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    else
    {
        DateTime _dtime;
        switch (timetype)
        {
        //sp->StartDateTime (globals.h)
        case SM_STARTDATE: _dtime = sp->StartDateTime; break;
        //sp->EndDateTime (globals.h)
        case SM_ENDDATE: _dtime = sp->EndDateTime;  break;
        //sp->ReportStart (globals.h)
        case SM_REPORTDATE: _dtime = sp->ReportStart;  break;
        default: return(ERR_API_OUTBOUNDS);
        }
        datetime_decodeDate(_dtime, year, month, day);
        datetime_decodeTime(_dtime, hours, minutes, seconds);
    }

    return (errcode);
}

int DLLEXPORT swmm_setSimulationDateTime(int timetype, char *dtimestr)
{
    return swmm_setSimulationDateTime_project(_defaultProject, timetype, dtimestr);
}
int DLLEXPORT swmm_setSimulationDateTime_project(SWMM_Project *sp, int timetype,
        char *dtimestr)
//
// Input:   timetype = time type to return
//          DateTime String
// Return:  API Error
// Purpose: Get the simulation start, end and report date times
{
    int errcode = 0;

    char theDate[10];
    char theTime[9];

    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if Simulation is Running
    else if(swmm_IsStartedFlag() == TRUE)
    {
        errcode = ERR_API_SIM_NRUNNING;
    }
    else
    {
        strncpy(theDate, dtimestr, 10);
        strncpy(theTime, dtimestr+11, 9);

        switch(timetype)
        {
            //sp->StartDateTime (globals.h)
            case SM_STARTDATE:
                project_readOption(sp, "START_DATE", theDate);
                project_readOption(sp, "START_TIME", theTime);
                sp->StartDateTime = sp->StartDate + sp->StartTime;
                sp->TotalDuration = floor((sp->EndDateTime - sp->StartDateTime) * SECperDAY);
                // --- convert total duration to milliseconds
                sp->TotalDuration *= 1000.0;
                break;
            //sp->EndDateTime (globals.h)
            case SM_ENDDATE:
                project_readOption(sp, "END_DATE", theDate);
                project_readOption(sp, "END_TIME", theTime);
                sp->EndDateTime = sp->EndDate + sp->EndTime;
                sp->TotalDuration = floor((sp->EndDateTime - sp->StartDateTime) * SECperDAY);
                // --- convert total duration to milliseconds
                sp->TotalDuration *= 1000.0;
                break;
            //sp->ReportStart (globals.h)
            case SM_REPORTDATE:
                project_readOption(sp, "REPORT_START_DATE", theDate);
                project_readOption(sp, "REPORT_START_TIME", theTime);
                sp->ReportStart = sp->ReportStartDate + sp->ReportStartTime;
                break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }

    return (errcode);
}

int DLLEXPORT  swmm_getSimulationUnit(int type, int *value)
{
    return swmm_getSimulationUnit_project(_defaultProject, type, value);
}

int DLLEXPORT  swmm_getSimulationUnit_project(SWMM_Project *sp, int type, int *value)
//
// Input:   type = simulation unit type
// Output:  enum representation of units
// Returns: API Error
// Purpose: get simulation unit types
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    else
    {
        switch(type)
        {
            // System Unit (enum.h UnitsType)
            case SM_SYSTEMUNIT:  *value = sp->UnitSystem; break;
            // Flow Unit (enum.h sp->FlowUnitsType)
            case SM_FLOWUNIT:  *value = sp->FlowUnits; break;
            // Concentration Unit
            //case 2:  *value = sp->UnitSystem; break;
            // Type not available
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }

    return (errcode);
}

int DLLEXPORT  swmm_getSimulationAnalysisSetting(int type, int *value)
{
    return swmm_getSimulationAnalysisSetting_project(_defaultProject, type, value);
}

int DLLEXPORT  swmm_getSimulationAnalysisSetting_project(SWMM_Project *sp,
        int type, int *value)
//
// Input:   type = analysis type
// Output:  setting True or False
// Returns: API Error
// Purpose: get simulation analysis setting
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    else
    {
        switch(type)
        {
            // No ponding at nodes (True or False)
            case SM_ALLOWPOND:  *value = sp->AllowPonding; break;
            // Do flow routing in steady state periods  (True or False)
            case SM_SKIPSTEADY:  *value = sp->SkipSteadyState; break;
            // Analyze rainfall/runoff  (True or False)
            case SM_IGNORERAIN:  *value = sp->IgnoreRainfall; break;
            // Analyze RDII (True or False)
            case SM_IGNORERDII:  *value = sp->IgnoreRDII; break;
            // Analyze snowmelt (True or False)
            case SM_IGNORESNOW:  *value = sp->IgnoreSnowmelt; break;
            // Analyze groundwater (True or False)
            case SM_IGNOREGW:  *value = sp->IgnoreGwater; break;
            // Analyze flow routing (True or False)
            case SM_IGNOREROUTE:  *value = sp->IgnoreRouting; break;
            // Analyze water quality (True or False)
            case SM_IGNORERQUAL:  *value = sp->IgnoreQuality; break;
            // Type not available
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return (errcode);
}

int DLLEXPORT  swmm_getSimulationParam(int type, double *value)
{
    return swmm_getSimulationParam_project(_defaultProject, type, value);
}

int DLLEXPORT  swmm_getSimulationParam_project(SWMM_Project *sp, int type, double *value)
//
// Input:   type = analysis type
// Output:  Simulation Parameter
// Returns: error code
// Purpose: Get simulation analysis parameter
{
    int errcode = 0;

    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Output  setting
    else
    {
        switch(type)
        {
            // Routing time step (sec)
            case SM_ROUTESTEP: *value = sp->RouteStep; break;
            // Minimum variable time step (sec)
            case SM_MINROUTESTEP: *value = sp->MinRouteStep; break;
            // Time step for lengthening (sec)
            case SM_LENGTHSTEP: *value = sp->LengtheningStep; break;
            // Antecedent dry days
            case SM_STARTDRYDAYS: *value = sp->StartDryDays; break;
            // Courant time step factor
            case SM_COURANTFACTOR: *value = sp->CourantFactor; break;
            // Minimum nodal surface area
            case SM_MINSURFAREA: *value = sp->MinSurfArea; break;
            // Minimum conduit slope
            case SM_MINSLOPE: *value = sp->MinSlope; break;
            // Runoff continuity error
            case SM_RUNOFFERROR: *value = sp->RunoffError; break;
            // Groundwater continuity error
            case SM_GWERROR: *value = sp->GwaterError; break;
            // Flow routing error
            case SM_FLOWERROR: *value = sp->FlowError; break;
            // Quality routing error
            case SM_QUALERROR: *value = sp->QualError; break;
            // DW routing head tolerance (ft)
            case SM_HEADTOL: *value = sp->HeadTol; break;
            // Tolerance for steady system flow
            case SM_SYSFLOWTOL: *value = sp->SysFlowTol; break;
            // Tolerance for steady nodal inflow
            case SM_LATFLOWTOL: *value = sp->LatFlowTol; break;
            // Type not available
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return (errcode);
}

int DLLEXPORT  swmm_countObjects(int type, int *count) {

    return swmm_countObjects_project(_defaultProject, type, count);
}

int DLLEXPORT  swmm_countObjects_project(SWMM_Project *sp, int type, int *count)
//
// Input:   type = object type (Based on SM_ObjectType enum)
// Output:  count = pointer to integer
// Returns: API Error
// Purpose: uses Object Count table to find number of elements of an object
{
    if(type >= MAX_OBJ_TYPES)return ERR_API_OUTBOUNDS;
    *count = sp->Nobjects[type];
    return (0);
}

int DLLEXPORT swmm_getObjectId(int type, int index, char *id)
{
    return swmm_getObjectId_project(_defaultProject, type, index, id);
}

int DLLEXPORT swmm_getObjectId_project(SWMM_Project *sp, int type, int index, char *id)
//
// Input:   type = object type (Based on SM_ObjectType enum)
//          index = Index of desired ID
// Output:  id = pointer to id pass by reference
// Return:  API Error
// Purpose: Gets ID for any object
{
    int errcode = 0;
    //Provide Empty Character Array
    strcpy(id,"");

    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[type])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch (type)
        {
            case SM_GAGE:
                strcpy(id,sp->Gage[index].ID); break;
            case SM_SUBCATCH:
                strcpy(id,sp->Subcatch[index].ID); break;
            case SM_NODE:
                strcpy(id,Node[index].ID); break;
            case SM_LINK:
                strcpy(id,Link[index].ID); break;
            case SM_POLLUT:
                strcpy(id,Pollut[index].ID); break;
            case SM_LANDUSE:
                strcpy(id,Landuse[index].ID); break;
            case SM_TIMEPATTERN:
                strcpy(id,Pattern[index].ID); break;
            case SM_CURVE:
                strcpy(id,Curve[index].ID); break;
            case SM_TSERIES:
                strcpy(id,Tseries[index].ID); break;
            //case SM_CONTROL:
                //strcpy(id,Rules[index].ID); break;
            case SM_TRANSECT:
                strcpy(id,Transect[index].ID); break;
            case SM_AQUIFER:
                strcpy(id,sp->Aquifer[index].ID); break;
            case SM_UNITHYD:
                strcpy(id,UnitHyd[index].ID); break;
            case SM_SNOWMELT:
                strcpy(id,sp->Snowmelt[index].ID); break;
            //case SM_SHAPE:
                //strcpy(id,Shape[index].ID); break;
            //case SM_LID:
                //strcpy(id,LidProcs[index].ID); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
   }
   return(errcode);
}

int DLLEXPORT swmm_getNodeType(int index, int *Ntype)
{
    return swmm_getNodeType_project(_defaultProject, index, Ntype);
}

int DLLEXPORT swmm_getNodeType_project(SWMM_Project *sp, int index, int *Ntype)
//
// Input:   index = Index of desired ID
//          Ntype = Node type (Based on enum SM_NodeType)
// Return:  API Error
// Purpose: Gets Node Type
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[NODE])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else *Ntype = Node[index].type;

    return(errcode);
}

int DLLEXPORT swmm_getLinkType(int index, int *Ltype)
{
    return swmm_getLinkType_project(_defaultProject, index, Ltype);
}

int DLLEXPORT swmm_getLinkType_project(SWMM_Project *sp, int index, int *Ltype)
//
// Input:   index = Index of desired ID
//          Ltype = Link type (Based on enum SM_LinkType)
// Return:  API Error
// Purpose: Gets Link Type
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else *Ltype = Link[index].type;

    return(errcode);
}

int DLLEXPORT swmm_getLinkConnections(int index, int *Node1, int *Node2)
{
    return swmm_getLinkConnections_project(_defaultProject, index, Node1, Node2);
}

int DLLEXPORT swmm_getLinkConnections_project(SWMM_Project *sp, int index, int *Node1, int *Node2)
//
// Input:   index = Index of desired ID
// Output:  Node1 and Node2 indeces
// Return:  API Error
// Purpose: Gets link Connection ID Indeces
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        *Node1 = Link[index].node1;
        *Node2 = Link[index].node2;
    }
    return(errcode);
}

int DLLEXPORT swmm_getLinkDirection(int index, signed char *value)
{
    return swmm_getLinkDirection_project(_defaultProject, index, value);
}

int DLLEXPORT swmm_getLinkDirection_project(SWMM_Project *sp, int index, signed char *value)
//
// Input:   index = Index of desired ID
// Output:  Link Direction (Only changes is slope < 0)
// Return:  API Error
// Purpose: Gets Link Direction
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        *value = Link[index].direction;
    }
    return(errcode);
}

int DLLEXPORT swmm_getNodeParam(int index, int Param, double *value)
{
    return swmm_getNodeParam_project(_defaultProject, index, Param, value);
}

int DLLEXPORT swmm_getNodeParam_project(SWMM_Project *sp, int index, int Param, double *value)
//
// Input:   index = Index of desired ID
//          param = Parameter desired (Based on enum SM_NodeProperty)
// Output:  value = value to be output
// Return:  API Error
// Purpose: Gets Node Parameter
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[NODE])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch(Param)
        {
            case SM_INVERTEL:
                *value = Node[index].invertElev * UCF(sp, LENGTH); break;
            case SM_FULLDEPTH:
                *value = Node[index].fullDepth * UCF(sp, LENGTH); break;
            case SM_SURCHDEPTH:
                *value = Node[index].surDepth * UCF(sp, LENGTH); break;
            case SM_PONDAREA:
                *value = Node[index].pondedArea * UCF(sp, LENGTH) * UCF(sp, LENGTH); break;
            case SM_INITDEPTH:
                *value = Node[index].initDepth * UCF(sp, LENGTH); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_setNodeParam(int index, int Param, double value)
{
    return swmm_setNodeParam_project(_defaultProject, index, Param, value);
}

int DLLEXPORT swmm_setNodeParam_project(SWMM_Project *sp, int index, int Param, double value)
//
// Input:   index = Index of desired ID
//          param = Parameter desired (Based on enum SM_NodeProperty)
//          value = value to be input
// Return:  API Error
// Purpose: Sets Node Parameter
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
     // Check if Simulation is Running
    else if(swmm_IsStartedFlag() == TRUE)
    {
        errcode = ERR_API_SIM_NRUNNING;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[NODE])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch(Param)
        {
            case SM_INVERTEL:
                Node[index].invertElev = value / UCF(sp, LENGTH); break;
            case SM_FULLDEPTH:
                Node[index].fullDepth = value / UCF(sp, LENGTH); break;
            case SM_SURCHDEPTH:
                Node[index].surDepth = value / UCF(sp, LENGTH); break;
            case SM_PONDAREA:
                Node[index].pondedArea = value / ( UCF(sp, LENGTH) * UCF(sp, LENGTH) ); break;
            case SM_INITDEPTH:
                Node[index].initDepth = value / UCF(sp, LENGTH); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    // Re-validated a node BEM 1/20/2017 Probably need to re-validate connecting links
    //node_validate(index)
    return(errcode);
}

int DLLEXPORT swmm_getLinkParam(int index, int Param, double *value)
{
    return swmm_getLinkParam_project(_defaultProject, index, Param, value);
}

int DLLEXPORT swmm_getLinkParam_project(SWMM_Project *sp, int index, int Param, double *value)
//
// Input:   index = Index of desired ID
//          param = Parameter desired (Based on enum SM_LinkProperty)
// Output:  value = value to be output
// Return:  API Error
// Purpose: Gets Link Parameter
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch(Param)
        {
            case SM_OFFSET1:
                *value = Link[index].offset1 * UCF(sp, LENGTH); break;
            case SM_OFFSET2:
                *value = Link[index].offset2 * UCF(sp, LENGTH); break;
            case SM_INITFLOW:
                *value = Link[index].q0 * UCF(sp, FLOW); break;
            case SM_FLOWLIMIT:
                *value = Link[index].qLimit * UCF(sp, FLOW); break;
            case SM_INLETLOSS:
                *value = Link[index].cLossInlet; break;
            case SM_OUTLETLOSS:
                *value = Link[index].cLossOutlet; break;
            case SM_AVELOSS:
                *value = Link[index].cLossAvg; break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_setLinkParam(int index, int Param, double value)
{
    return swmm_setLinkParam_project(_defaultProject, index, Param, value);
}

int DLLEXPORT swmm_setLinkParam_project(SWMM_Project *sp, int index, int Param, double value)
//
// Input:   index = Index of desired ID
//          param = Parameter desired (Based on enum SM_LinkProperty)
//          value = value to be output
// Return:  API Error
// Purpose: Sets Link Parameter
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch(Param)
        {
            // offset1
            case SM_OFFSET1:
                // Check if Simulation is Running
                if(swmm_IsStartedFlag() == TRUE)
                {
                    errcode = ERR_API_SIM_NRUNNING; break;
                }
                Link[index].offset1 = value / UCF(sp, LENGTH); break;
            case SM_OFFSET2:
                // Check if Simulation is Running
                if(swmm_IsStartedFlag() == TRUE)
                {
                    errcode = ERR_API_SIM_NRUNNING; break;
                }
                Link[index].offset2 = value / UCF(sp, LENGTH); break;
            case SM_INITFLOW:
                Link[index].q0 = value / UCF(sp, FLOW); break;
            case SM_FLOWLIMIT:
                Link[index].qLimit = value / UCF(sp, FLOW); break;
//            case SM_INLETLOSS:
//                Link[index].cLossInlet; break;
//            case SM_OUTLETLOSS:
//                Link[index].cLossOutlet; break;
//            case SM_AVELOSS:
//                Link[index].cLossAvg; break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
        // re-validated link
        //link_validate(index);
    }

    return(errcode);
}

int DLLEXPORT swmm_getSubcatchParam(int index, int Param, double *value)
{
    return swmm_getSubcatchParam_project(_defaultProject, index, Param, value);
}

int DLLEXPORT swmm_getSubcatchParam_project(SWMM_Project *sp, int index, int Param, double *value)
//
// Input:   index = Index of desired ID
//          param = Parameter desired (Based on enum SM_SubcProperty)
// Output:  value = value to be output
// Return:  API Error
// Purpose: Gets Subcatchment Parameter
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[SUBCATCH])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch(Param)
        {
            case SM_WIDTH:
                *value = sp->Subcatch[index].width * UCF(sp, LENGTH); break;
            case SM_AREA:
                *value = sp->Subcatch[index].area * UCF(sp, LANDAREA); break;
            case SM_FRACIMPERV:
                *value = sp->Subcatch[index].fracImperv; break;
            case SM_SLOPE:
                *value = sp->Subcatch[index].slope; break;
            case SM_CURBLEN:
                *value = sp->Subcatch[index].curbLength * UCF(sp, LENGTH); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_setSubcatchParam(int index, int Param, double value) {
    return swmm_setSubcatchParam_project(_defaultProject, index, Param, value);
}

int DLLEXPORT swmm_setSubcatchParam_project(SWMM_Project *sp, int index,
        int Param, double value)
//
// Input:   index = Index of desired ID
//          param = Parameter desired (Based on enum SM_SubcProperty)
//          value = value to be output
// Return:  API Error
// Purpose: Sets Subcatchment Parameter
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
     // Check if Simulation is Running
    else if(swmm_IsStartedFlag() == TRUE)
    {
        errcode = ERR_API_SIM_NRUNNING;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[SUBCATCH])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch(Param)
        {
            case SM_WIDTH:
                sp->Subcatch[index].width = value / UCF(sp, LENGTH); break;
            case SM_AREA:
                sp->Subcatch[index].area = value / UCF(sp, LANDAREA); break;
            case SM_FRACIMPERV:
                sp->Subcatch[index].fracImperv; break;
            case SM_SLOPE:
                sp->Subcatch[index].slope; break;
            case SM_CURBLEN:
                sp->Subcatch[index].curbLength = value / UCF(sp, LENGTH); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    //re-validate subcatchment
    subcatch_validate(sp, index); // incorprate callback here

    return(errcode);
}

int DLLEXPORT swmm_getSubcatchOutConnection(int index, int *type, int *ObjIndex)
{
    return swmm_getSubcatchOutConnection_project(_defaultProject, index, type, ObjIndex);
}

int DLLEXPORT swmm_getSubcatchOutConnection_project(SWMM_Project *sp, int index, int *type, int *ObjIndex )
//
// Input:   index = Index of desired ID
//         (Subcatchments can load to Node or another Subcatchment)
// Output:  Type of Object
//          Index of Object
// Return:  API Error
// Purpose: Gets Subcatchment Connection ID Indeces for either Node or Subcatchment
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[SUBCATCH])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        if (sp->Subcatch[index].outNode == -1 && sp->Subcatch[index].outSubcatch == -1)
        {
            *ObjIndex = index; // Case of self Loading subcatchment
            *type = SUBCATCH;
        }
        if (sp->Subcatch[index].outNode >= 0)
        {
            *ObjIndex = sp->Subcatch[index].outNode;
            *type = NODE;
        }
        if (sp->Subcatch[index].outSubcatch >= 0)
        {
            *ObjIndex = sp->Subcatch[index].outSubcatch;
            *type = SUBCATCH;
        }
    }
    return(errcode);
}


//-------------------------------
// Active Simulation Results API
//-------------------------------
int DLLEXPORT swmm_getCurrentDateTimeStr(char *dtimestr)
{
    return swmm_getCurrentDateTimeStr_project(_defaultProject, dtimestr);
}

int DLLEXPORT swmm_getCurrentDateTimeStr_project(SWMM_Project *sp, char *dtimestr)
//
// Output:  DateTime String
// Return:  API Error
// Purpose: Get the current simulation time
{
    //Provide Empty Character Array
    char     theDate[12];
    char     theTime[9];
    char     _DTimeStr[22];
    DateTime currentTime;

    // Check if Simulation is Running
    if(swmm_IsStartedFlag() == FALSE) return(ERR_API_SIM_NRUNNING);

    // Fetch Current Time
    currentTime = getDateTime(sp, sp->NewRoutingTime);

    // Convert To Char
    datetime_dateToStr(currentTime, theDate);
    datetime_timeToStr(currentTime, theTime);

    strcpy(_DTimeStr, theDate);
    strcat(_DTimeStr, " ");
    strcat(_DTimeStr, theTime);

    strcpy(dtimestr, _DTimeStr);
    return(0);
}

int DLLEXPORT swmm_getNodeResult(int index, int type, double *result)
{
    return swmm_getNodeResult_project(_defaultProject, index, type, result);
}

int DLLEXPORT swmm_getNodeResult_project(SWMM_Project *sp, int index, int type, double *result)
//
// Input:   index = Index of desired ID
//          type = Result Type (SM_NodeResult)
// Output:  result = result data desired (byref)
// Return:  API Error
// Purpose: Gets Node Simulated Value at Current Time
{
    int errcode = 0;
    // Check if Simulation is Running
    if(swmm_IsStartedFlag() == FALSE)
    {
        errcode = ERR_API_SIM_NRUNNING;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[NODE])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch (type)
        {
            case SM_TOTALINFLOW:
                *result = Node[index].inflow * UCF(sp, FLOW); break;
            case SM_TOTALOUTFLOW:
                *result = Node[index].outflow * UCF(sp, FLOW); break;
            case SM_LOSSES:
                *result = Node[index].losses * UCF(sp, FLOW); break;
            case SM_NODEVOL:
                *result = Node[index].newVolume * UCF(sp, VOLUME); break;
            case SM_NODEFLOOD:
                *result = Node[index].overflow * UCF(sp, FLOW); break;
            case SM_NODEDEPTH:
                *result = Node[index].newDepth * UCF(sp, LENGTH); break;
            case SM_NODEHEAD:
                *result = (Node[index].newDepth
                            + Node[index].invertElev) * UCF(sp, LENGTH); break;
            case SM_LATINFLOW:
                *result = Node[index].newLatFlow * UCF(sp, FLOW); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_getLinkResult(int index, int type, double *result)
{
    return swmm_getLinkResult_project(_defaultProject, index, type, result);
}

int DLLEXPORT swmm_getLinkResult_project(SWMM_Project *sp, int index, int type, double *result)
//
// Input:   index = Index of desired ID
//          type = Result Type (SM_LinkResult)
// Output:  result = result data desired (byref)
// Return:  API Error
// Purpose: Gets Link Simulated Value at Current Time
{
    int errcode = 0;
    // Check if Simulation is Running
    if(swmm_IsStartedFlag() == FALSE)
    {
        errcode = ERR_API_SIM_NRUNNING;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch (type)
        {
            case SM_LINKFLOW:
                *result = Link[index].newFlow * UCF(sp, FLOW) ; break;
            case SM_LINKDEPTH:
                *result = Link[index].newDepth * UCF(sp, LENGTH); break;
            case SM_LINKVOL:
                *result = Link[index].newVolume * UCF(sp, VOLUME); break;
            case SM_USSURFAREA:
                *result = Link[index].surfArea1 * UCF(sp, LENGTH) * UCF(sp, LENGTH); break;
            case SM_DSSURFAREA:
                *result = Link[index].surfArea2 * UCF(sp, LENGTH) * UCF(sp, LENGTH); break;
            case SM_SETTING:
                *result = Link[index].setting; break;
            case SM_TARGETSETTING:
                *result = Link[index].targetSetting; break;
            case SM_FROUDE:
                *result = Link[index].froude; break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_getSubcatchResult(int index, int type, double *result)
{
    return swmm_getSubcatchResult_project(_defaultProject, index, type, result);
}

int DLLEXPORT swmm_getSubcatchResult_project(SWMM_Project *sp, int index, int type, double *result)
//
// Input:   index = Index of desired ID
//          type = Result Type (SM_SubcResult)
// Output:  result = result data desired (byref)
// Return:  API Error
// Purpose: Gets Subcatchment Simulated Value at Current Time
{
    int errcode = 0;
    // Check if Simulation is Running
    if(swmm_IsStartedFlag() == FALSE)
    {
        errcode = ERR_API_SIM_NRUNNING;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[SUBCATCH])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        switch (type)
        {
            case SM_SUBCRAIN:
                *result = sp->Subcatch[index].rainfall * UCF(sp, RAINFALL); break;
            case SM_SUBCEVAP:
                *result = sp->Subcatch[index].evapLoss * UCF(sp, EVAPRATE); break;
            case SM_SUBCINFIL:
                *result = sp->Subcatch[index].infilLoss * UCF(sp, RAINFALL); break;
            case SM_SUBCRUNON:
                *result = sp->Subcatch[index].runon * UCF(sp, FLOW); break;
            case SM_SUBCRUNOFF:
                *result = sp->Subcatch[index].newRunoff * UCF(sp, FLOW); break;
            case SM_SUBCSNOW:
                *result = sp->Subcatch[index].newSnowDepth * UCF(sp, RAINDEPTH); break;
            default: errcode = ERR_API_OUTBOUNDS; break;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_getNodeStats(int index, SM_NodeStats *nodeStats){
    return swmm_getNodeStats_project(_defaultProject, index, nodeStats);
}
int DLLEXPORT swmm_getNodeStats_project(SWMM_Project *sp, int index,
        SM_NodeStats *nodeStats)
//
// Output:  Node Stats Structure (SM_NodeStats)
// Return:  API Error
// Purpose: Gets Node Stats and Converts Units
{
    int errorcode = stats_getNodeStat(sp, index, nodeStats);

    if (errorcode == 0)
    {
        // Current Average Depth
        nodeStats->avgDepth *= (UCF(sp, LENGTH) / (double)sp->StepCount);
        // Current Maximum Depth
        nodeStats->maxDepth *= UCF(sp, LENGTH);
        // Current Maximum Lateral Inflow
        nodeStats->maxLatFlow *= UCF(sp, FLOW);
        // Current Maximum Inflow
        nodeStats->maxInflow *= UCF(sp, FLOW);
        // Cumulative Lateral Inflow
        nodeStats->totLatFlow *= UCF(sp, VOLUME);
        // Time Courant Critical (hrs)
        nodeStats->timeCourantCritical /= 3600.0;
        // Cumulative Flooded Volume
        nodeStats->volFlooded *= UCF(sp, VOLUME);
        // Time Flooded (hrs)
        nodeStats->timeFlooded /= 3600.0;
        // Current Maximum Overflow
        nodeStats->maxOverflow *= UCF(sp, FLOW);
        // Current Maximum Ponding Volume
        nodeStats->maxPondedVol *= UCF(sp, VOLUME);
        // Time Surcharged
        nodeStats->timeSurcharged /= 3600.0;
    }
    return (errorcode);
}

int DLLEXPORT swmm_getNodeTotalInflow(int index, double *value)
{
    return swmm_getNodeTotalInflow_project(_defaultProject, index, value);
}

int DLLEXPORT swmm_getNodeTotalInflow_project(SWMM_Project *sp, int index, double *value)
//
// Input:   Node Index
// Output:  Node Total inflow Volume.
// Return:  API Error
// Purpose: Get Node Total Inflow Volume.
{

    int errorcode = massbal_getNodeTotalInflow(index, value);

    if (errorcode == 0)
    {
        *value *= UCF(sp, VOLUME);
    }

    return(errorcode);
}

int DLLEXPORT swmm_getStorageStats(int index, SM_StorageStats *storageStats){
    return swmm_getStorageStats_project(_defaultProject, index, storageStats);
}
int DLLEXPORT swmm_getStorageStats_project(SWMM_Project *sp, int index,
        SM_StorageStats *storageStats)
//
// Output:  Storage Node Stats Structure (SM_StorageStats)
// Return:  API Error
// Purpose: Gets Storage Node Stats and Converts Units
{
    int errorcode = stats_getStorageStat(sp, index, storageStats);

    if (errorcode == 0)
    {
        // Initial Volume
        storageStats->initVol *= UCF(sp, VOLUME);
        // Current Average Volume
        storageStats->avgVol *= (UCF(sp, VOLUME) / (double)sp->StepCount);
        // Current Maximum Volume
        storageStats->maxVol *= UCF(sp, VOLUME);
        // Current Maximum Flow
        storageStats->maxFlow *= UCF(sp, FLOW);
        // Current Evaporation Volume
        storageStats->evapLosses *= UCF(sp, VOLUME);
        // Current Exfiltration Volume
        storageStats->exfilLosses *= UCF(sp, VOLUME);
    }

    return (errorcode);
}

int DLLEXPORT swmm_getOutfallStats(int index, SM_OutfallStats *outfallStats)
{
    return swmm_getOutfallStats_project(_defaultProject, index,
        outfallStats);
}

int DLLEXPORT swmm_getOutfallStats_project(SWMM_Project *sp, int index,
        SM_OutfallStats *outfallStats)
//
// Output:  Outfall Stats Structure (SM_OutfallStats)
// Return:  API Error
// Purpose: Gets Outfall Node Stats and Converts Units
// Note:    Caller is responsible for calling swmm_freeOutfallStats
//          to free the pollutants array.
{
    int p;
    int errorcode = stats_getOutfallStat(sp, index, outfallStats);

    if (errorcode == 0)
    {
        // Current Average Flow
        if ( outfallStats->totalPeriods > 0 )
        {
            outfallStats->avgFlow *= (UCF(sp, FLOW) / (double)outfallStats->totalPeriods);
        }
        else
        {
            outfallStats->avgFlow *= 0.0;
        }
        // Current Maximum Flow
        outfallStats->maxFlow *= UCF(sp, FLOW);
        // Convert Mass Units
        if (sp->Nobjects[POLLUT] > 0)
        {
            for (p = 0; p < sp->Nobjects[POLLUT]; p++)
                outfallStats->totalLoad[p] *= (LperFT3 * Pollut[p].mcf);

            if (Pollut[p].units == COUNT)
            {
                outfallStats->totalLoad[p] = LOG10(outfallStats->totalLoad[p]);
            }
        }
    }

    return (errorcode);
}

void DLLEXPORT swmm_freeOutfallStats(SM_OutfallStats *outfallStats)
//
// Return:  API Error
// Purpose: Frees Outfall Node Stats and Converts Units
// Note:    API user is responsible for calling swmm_freeOutfallStats
//          since this function performs a memory allocation.
{
    FREE(outfallStats->totalLoad);
}

int DLLEXPORT swmm_getLinkStats(int index, SM_LinkStats *linkStats)
{
    return swmm_getLinkStats_project(_defaultProject, index, linkStats);
}

int DLLEXPORT swmm_getLinkStats_project(SWMM_Project *sp, int index,
        SM_LinkStats *linkStats)
//
// Output:  Link Stats Structure (SM_LinkStats)
// Return:  API Error
// Purpose: Gets Link Stats and Converts Units
{
    int errorcode = stats_getLinkStat(sp, index, linkStats);

    if (errorcode == 0)
    {
        // Cumulative Maximum Flowrate
        linkStats->maxFlow *= UCF(sp, FLOW);
        // Cumulative Maximum Velocity
        linkStats->maxVeloc *= UCF(sp, LENGTH);
        // Cumulative Maximum Depth
        linkStats->maxDepth *= UCF(sp, LENGTH);
        // Cumulative Time Normal Flow
        linkStats->timeNormalFlow /= 3600.0;
        // Cumulative Time Inlet Control
        linkStats->timeInletControl /= 3600.0;
        // Cumulative Time Surcharged
        linkStats->timeSurcharged /= 3600.0;
        // Cumulative Time Upstream Full
        linkStats->timeFullUpstream /= 3600.0;
        // Cumulative Time Downstream Full
        linkStats->timeFullDnstream /= 3600.0;
        // Cumulative Time Full Flow
        linkStats->timeFullFlow /= 3600.0;
        // Cumulative Time Capacity limited
        linkStats->timeCapacityLimited /= 3600.0;
        // Cumulative Time Courant Critical Flow
        linkStats->timeCourantCritical /= 3600.0;
    }

    return (errorcode);
}

int DLLEXPORT swmm_getPumpStats(int index, SM_PumpStats *pumpStats)
{
    return swmm_getPumpStats_project(_defaultProject, index, pumpStats);
}

int DLLEXPORT swmm_getPumpStats_project(SWMM_Project *sp, int index,
        SM_PumpStats *pumpStats)
//
// Output:  Pump Link Stats Structure (SM_PumpStats)
// Return:  API Error
// Purpose: Gets Pump Link Stats and Converts Units
{
    int errorcode = stats_getPumpStat(sp, index, pumpStats);

    if (errorcode == 0)
    {
        // Cumulative Minimum Flow
        pumpStats->minFlow *= UCF(sp, FLOW);
        // Cumulative Average Flow
        if (pumpStats->totalPeriods > 0)
        {
            pumpStats->avgFlow *= (UCF(sp, FLOW) / (double)pumpStats->totalPeriods);
        }
        else
        {
            pumpStats->avgFlow *= 0.0;
        }
        // Cumulative Maximum Flow
        pumpStats->maxFlow *= UCF(sp, FLOW);
        // Cumulative Pumping Volume
        pumpStats->volume *= UCF(sp, VOLUME);
    }

    return (errorcode);
}

int DLLEXPORT swmm_getSubcatchStats(int index, SM_SubcatchStats *subcatchStats)
{
    return swmm_getSubcatchStats_project(_defaultProject, index, subcatchStats);
}

int DLLEXPORT swmm_getSubcatchStats_project(SWMM_Project *sp, int index,
        SM_SubcatchStats *subcatchStats)
//
// Output:  Subcatchment Stats Structure (SM_SubcatchStats)
// Return:  API Error
// Purpose: Gets Subcatchment Stats and Converts Units
// Note: Caller is responsible for calling swmm_freeSubcatchStats
//       to free the pollutants array.
{
    int p;
    int errorcode = stats_getSubcatchStat(sp, index, subcatchStats);

    if (errorcode == 0)
    {
        double a = sp->Subcatch[index].area;

        // Cumulative Runon Volume
        subcatchStats->runon *= (UCF(sp, RAINDEPTH) / a);
        // Cumulative Infiltration Volume
        subcatchStats->infil *= (UCF(sp, RAINDEPTH) / a);
        // Cumulative Runoff Volume
        subcatchStats->runoff *= (UCF(sp, RAINDEPTH) / a);
        // Maximum Runoff Rate
        subcatchStats->maxFlow *= UCF(sp, FLOW);
        // Cumulative Rainfall Depth
        subcatchStats->precip *= (UCF(sp, RAINDEPTH) / a);
        // Cumulative Evaporation Volume
        subcatchStats->evap *= (UCF(sp, RAINDEPTH) / a);

        if (sp->Nobjects[POLLUT] > 0)
        {
            for (p = 0; p < sp->Nobjects[POLLUT]; p++)
                subcatchStats->surfaceBuildup[p] /= (a * UCF(sp, LANDAREA));

            if (Pollut[p].units == COUNT)
            {
                subcatchStats->surfaceBuildup[p] =
                        LOG10(subcatchStats->surfaceBuildup[p]);
            }
        }
    }

    return (errorcode);
}


void DLLEXPORT swmm_freeSubcatchStats(SM_SubcatchStats *subcatchStats)
//
// Return:  API Error
// Purpose: Frees Subcatchment Stats
// Note:    API user is responsible for calling swmm_freeSubcatchStats
//          since this function performs a memory allocation.
{
    FREE(subcatchStats->surfaceBuildup);
}

int DLLEXPORT swmm_getSystemRoutingStats(SM_RoutingTotals *routingTot)
{
    return swmm_getSystemRoutingStats_project(_defaultProject, routingTot);
}

int DLLEXPORT swmm_getSystemRoutingStats_project(SWMM_Project *sp,
        SM_RoutingTotals *routingTot)
//
// Output:  System Routing Totals Structure (SM_RoutingTotals)
// Return:  API Error
// Purpose: Gets System Flow Routing Totals and Converts Units
{
    int errorcode = massbal_getRoutingFlowTotal(routingTot);

    if (errorcode == 0)
    {
        // Cumulative Dry Weather Inflow Volume
        routingTot->dwInflow *= UCF(sp, VOLUME);
        // Cumulative Wet Weather Inflow Volume
        routingTot->wwInflow *= UCF(sp, VOLUME);
        // Cumulative Groundwater Inflow Volume
        routingTot->gwInflow *= UCF(sp, VOLUME);
        // Cumulative I&I Inflow Volume
        routingTot->iiInflow *= UCF(sp, VOLUME);
        // Cumulative External Inflow Volume
        routingTot->exInflow *= UCF(sp, VOLUME);
        // Cumulative Flooding Volume
        routingTot->flooding *= UCF(sp, VOLUME);
        // Cumulative Outflow Volume
        routingTot->outflow  *= UCF(sp, VOLUME);
        // Cumulative Evaporation Loss
        routingTot->evapLoss *= UCF(sp, VOLUME);
        // Cumulative Seepage Loss
        routingTot->seepLoss *= UCF(sp, VOLUME);
        // Continuity Error
        routingTot->pctError *= 100;
    }

    return(errorcode);
}

int DLLEXPORT swmm_getSystemRunoffStats(SM_RunoffTotals *runoffTot)
{
    return swmm_getSystemRunoffStats_project(_defaultProject, runoffTot);
}

int DLLEXPORT swmm_getSystemRunoffStats_project(SWMM_Project *sp, SM_RunoffTotals *runoffTot)
//
// Output:  System Runoff Totals Structure (SM_RunoffTotals)
// Return:  API Error
// Purpose: Gets System Runoff Totals and Converts Units
{
    int errorcode =  massbal_getRunoffTotal(runoffTot);

    if (errorcode == 0)
    {
        double TotalArea = massbal_getTotalArea();
        // Cumulative Rainfall Volume
        runoffTot->rainfall *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Cumulative Evaporation Volume
        runoffTot->evap *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Cumulative Infiltration Volume
        runoffTot->infil *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Cumulative Runoff Volume
        runoffTot->runoff *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Cumulative Runon Volume
        runoffTot->runon *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Cumulative Drain Volume
        runoffTot->drains *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Cumulative Snow Removed Volume
        runoffTot->snowRemoved *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Initial Storage Volume
        runoffTot->initStorage *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Initial Snow Cover Volume
        runoffTot->initSnowCover *= (UCF(sp, RAINDEPTH) / TotalArea);
        // Continuity Error
        runoffTot->pctError *= 100;
    }

    return(errorcode);
}

int DLLEXPORT swmm_getGagePrecip(int index, double *rainfall, double *snowfall,
            double *total)
{
    return swmm_getGagePrecip_project(_defaultProject, index, rainfall, snowfall,
            total);
}

int DLLEXPORT swmm_getGagePrecip_project(SWMM_Project *sp, int index,
        double *rainfall, double *snowfall, double *total)
//
// Input:   index = Index of desired ID
// Output:  Rainfall intensity and snow for the gage
// Return:  API Error
// Purpose: Gets the precipitaion value in the gage. 
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
	    errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[GAGE])
    {
	    errcode = ERR_API_OBJECT_INDEX;
    }
    // Read the rainfall value
    else
    {
        *total = gage_getPrecip(sp, index, rainfall, snowfall);
    }
    return(errcode);
}
int DLLEXPORT swmm_setGagePrecip(int index, double value)
{
    return swmm_setGagePrecip_project(_defaultProject, index, value);
}

int DLLEXPORT swmm_setGagePrecip_project(SWMM_Project *sp, int index, double value)
//
// Input:   index = Index of desired ID
//          value = rainfall intensity to be set
// Return:  API Error
// Purpose: Sets the precipitation in from the external database
{
    int errcode = 0;
    // Check if Open
    if(swmm_IsOpenFlag() == FALSE)
    {
	    errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[GAGE])
    {
	    errcode = ERR_API_OBJECT_INDEX;
    }
    // Read the rainfall value
    else
    {
        if (sp->Gage[index].dataSource != RAIN_API)
        {
            sp->Gage[index].dataSource = RAIN_API;
        }
	    sp->Gage[index].externalRain = value * UCF(sp, RAINFALL);
    }
    return(errcode);
}


//-------------------------------
// Setters API
//-------------------------------
int DLLEXPORT swmm_setLinkSetting(int index, double targetSetting) {
    return swmm_setLinkSetting_project(_defaultProject, index, targetSetting);
}

int DLLEXPORT swmm_setLinkSetting_project(SWMM_Project *sp, int index,
        double targetSetting)
//
// Input:   index = Index of desired ID
//          value = New Target Setting
// Output:  returns API Error
// Purpose: Sets Link open fraction (Weir, Orifice, Pump, and Outlet)
{
    DateTime currentTime;
    int errcode = 0;
    char _rule_[11] = "ToolkitAPI";

    // Check if Open
    if (swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[LINK])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        // --- check that new setting lies within feasible limits
        if (targetSetting < 0.0) targetSetting = 0.0;
        if (Link[index].type != PUMP && targetSetting > 1.0) targetSetting = 1.0;

        Link[index].targetSetting = targetSetting;

        // Use internal function to apply the new setting
        link_setSetting(sp, index, 0.0);

        // Add control action to RPT file if desired flagged
        if (sp->RptFlags.controls)
        {
            currentTime = getDateTime(sp, sp->NewRoutingTime);
            report_writeControlAction(sp, currentTime, Link[index].ID,
                    targetSetting, _rule_);
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_setNodeInflow(int index, double flowrate)
{
    return swmm_setNodeInflow_project(_defaultProject, index, flowrate);
}

int DLLEXPORT swmm_setNodeInflow_project(SWMM_Project *sp, int index, double flowrate)
//
// Input:   index = Index of desired ID
//          value = New Inflow Rate
// Output:  returns API Error
// Purpose: Sets new node inflow rate and holds until set again
{
    int errcode = 0;

    // Check if Open
    if (swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if (index < 0 || index >= sp->Nobjects[NODE])
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else
    {
        // Check to see if node has an assigned inflow object
        TExtInflow* inflow;

        // --- check if an external inflow object for this constituent already exists
        inflow = Node[index].extInflow;
        while (inflow)
        {
            if (inflow->param == -1) break;
            inflow = inflow->next;
        }

        if (!inflow)
        {
            int param = -1;        // FLOW (-1) or Pollutant Index
            int type = FLOW_INFLOW;// Type of inflow (FLOW)
            int tSeries = -1;      // No Time Series
            int basePat = -1;      // No Base Pattern
            double cf = 1.0;       // Unit Convert (Converted during validation)
            double sf = 1.0;       // Scaling Factor
            double baseline = 0.0; // Baseline Inflow Rate

            // Initializes Inflow Object
            errcode = inflow_setExtInflow(sp, index, param, type, tSeries,
                basePat, cf, baseline, sf);

            // Get The Inflow Object
            if ( errcode == 0 )
            {
                inflow = Node[index].extInflow;
            }
        }
        // Assign new flow rate
        if ( errcode == 0 )
        {
            inflow -> extIfaceInflow = flowrate;
        }
    }
    return(errcode);
}

int DLLEXPORT swmm_setOutfallStage(int index, double stage)
{
    return swmm_setOutfallStage_project(_defaultProject, index, stage);
}

int DLLEXPORT swmm_setOutfallStage_project(SWMM_Project *sp, int index, double stage)
//
// Input:   index = Index of desired outfall
//          stage = New outfall stage (head)
// Output:  returns API Error
// Purpose: Sets new outfall stage and holds until set again.
{
    int errcode = 0;
    // Check if Open
    if (swmm_IsOpenFlag() == FALSE)
    {
        errcode = ERR_API_INPUTNOTOPEN;
    }
    // Check if object index is within bounds
    else if ( index < 0 || index >= sp->Nobjects[NODE] )
    {
        errcode = ERR_API_OBJECT_INDEX;
    }
    else if ( Node[index].type != OUTFALL )
    {
        errcode = ERR_API_WRONG_TYPE;
    }
    else
    {
        int k = Node[index].subIndex;
        if ( Outfall[k].type != STAGED_OUTFALL )
        {
            // Change Boundary Conditions Setting Type
            Outfall[k].type = STAGED_OUTFALL;
        }
        Outfall[k].outfallStage = stage / UCF(sp, LENGTH);
    }
    return(errcode);
}
