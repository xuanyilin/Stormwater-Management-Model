//-----------------------------------------------------------------------------
//   funcs.h
//
//   Project:  EPA SWMM5
//   Version:  5.1
//   Date:     03/20/14  (Build 5.1.000)
//             09/15/14  (Build 5.1.007)
//             04/02/15  (Build 5.1.008)
//             08/05/15  (Build 5.1.010)
//   Author:   L. Rossman (EPA)
//             M. Tryby (EPA)
//
//   Global interfacing functions.
//
//   Build 5.1.007:
//   - climate_readAdjustments() added.
//
//   Build 5.1.008:
//   - Function list was re-ordered and blank lines added for readability.
//   - Pollutant buildup/washoff functions for the new surfqual.c module added.
//   - Several other functions added, re-named or have modified arguments.
//
//   Build 5.1.010:
//   - New roadway_getInflow() function added.
//
//-----------------------------------------------------------------------------

#ifndef FUNCS_H
#define FUNCS_H


//-----------------------------------------------------------------------------
//   Project Manager Methods
//-----------------------------------------------------------------------------
// --- define WINDOWS

#undef WINDOWS
#ifdef _WIN32
#define WINDOWS
#endif
#ifdef __WIN32__
#define WINDOWS
#endif

// --- define DLLEXPORT

#ifdef WINDOWS
	#ifdef __MINGW32__
		// Seems to be more wrapper friendly
		#define DLLEXPORT __declspec(dllexport) __cdecl 
	#else
		#define DLLEXPORT __declspec(dllexport) __stdcall
	#endif
#else
	#define DLLEXPORT
#endif


void     project_open(SWMM_Project *sp, char *f1, char *f2, char *f3);
void     project_close(SWMM_Project *sp);

void     project_readInput(SWMM_Project *sp);
int      project_readOption(SWMM_Project *sp, char* s1, char* s2);
void     project_validate(SWMM_Project *sp);
int      project_init(SWMM_Project *sp);

int      project_addObject(int type, char* id, int n);

#ifdef __cplusplus
extern "C" {		// --- use "C" linkage for C++ programs
#endif 
	int   DLLEXPORT   project_findObject(int type, char* id);
#ifdef __cplusplus 
}   // matches the linkage specification from above */ 
#endif

char*    project_findID(int type, char* id);

double** project_createMatrix(int nrows, int ncols);
void     project_freeMatrix(double** m);

//-----------------------------------------------------------------------------
//   Input Reader Methods
//-----------------------------------------------------------------------------
int     input_countObjects(SWMM_Project *p);
int     input_readData(SWMM_Project *p);

//-----------------------------------------------------------------------------
//   Report Writer Methods
//-----------------------------------------------------------------------------
int     report_readOptions(SWMM_Project *sp, char* tok[], int ntoks);

void    report_writeLine(SWMM_Project *sp, char* line);
void    report_writeSysTime(SWMM_Project *sp);
void    report_writeLogo(SWMM_Project *sp);
void    report_writeTitle(SWMM_Project *sp);
void    report_writeOptions(SWMM_Project *sp);
void    report_writeReport(SWMM_Project *sp);

void    report_writeRainStats(SWMM_Project *sp, int gage, TRainStats* rainStats);
void    report_writeRdiiStats(SWMM_Project *sp, double totalRain, double totalRdii);

void    report_writeControlActionsHeading(SWMM_Project *sp);
void    report_writeControlAction(SWMM_Project *sp, DateTime aDate, char* linkID, double value,
        char* ruleID);

void    report_writeRunoffError(SWMM_Project *sp, TRunoffTotals* totals,
        double area);
void    report_writeLoadingError(SWMM_Project *sp, TLoadingTotals* totals);
void    report_writeGwaterError(SWMM_Project *sp, TGwaterTotals* totals,
        double area);
void    report_writeFlowError(SWMM_Project *sp, TRoutingTotals* totals);
void    report_writeQualError(SWMM_Project *sp, TRoutingTotals* totals);

void    report_writeMaxStats(SWMM_Project *sp, TMaxStats massBalErrs[],
        TMaxStats CourantCrit[], int nMaxStats);
void    report_writeMaxFlowTurns(SWMM_Project *sp, TMaxStats flowTurns[],
        int nMaxStats);
void    report_writeSysStats(SWMM_Project *sp, TSysStats* sysStats);

void    report_writeErrorMsg(SWMM_Project *sp, int code, char* msg);
void    report_writeErrorCode(SWMM_Project *sp);
void    report_writeInputErrorMsg(SWMM_Project *sp, int k, int sect,
        char* line, long lineCount);
void    report_writeWarningMsg(SWMM_Project *sp, char* msg, char* id);
void    report_writeTseriesErrorMsg(SWMM_Project *sp, int code, TTable *tseries);

void    inputrpt_writeInput(SWMM_Project *sp);
void    statsrpt_writeReport(SWMM_Project *sp);

//-----------------------------------------------------------------------------
//   Temperature/Evaporation Methods
//-----------------------------------------------------------------------------
int      climate_readParams(SWMM_Project *sp, char* tok[], int ntoks);
int      climate_readEvapParams(SWMM_Project *sp, char* tok[], int ntoks);
int      climate_readAdjustments(SWMM_Project *sp, char* tok[], int ntoks);                      //(5.1.007)
void     climate_validate(SWMM_Project *sp);
void     climate_openFile(SWMM_Project *sp);
void     climate_initState(SWMM_Project *sp);
void     climate_setState(SWMM_Project *sp, DateTime aDate);
DateTime climate_getNextEvapDate(SWMM_Project *sp);                                        //(5.1.008)

//-----------------------------------------------------------------------------
//   Rainfall Processing Methods
//-----------------------------------------------------------------------------
void    rain_open(SWMM_Project *sp);
void    rain_close(SWMM_Project *sp);

//-----------------------------------------------------------------------------
//   Snowmelt Processing Methods
//-----------------------------------------------------------------------------
int     snow_readMeltParams(SWMM_Project *sp, char* tok[], int ntoks);
int     snow_createSnowpack(SWMM_Project *sp, int subcacth, int snowIndex);

void    snow_validateSnowmelt(SWMM_Project *sp, int snowIndex);
void    snow_initSnowpack(SWMM_Project *sp, int subcatch);
void    snow_initSnowmelt(SWMM_Project *sp, int snowIndex);

void    snow_getState(SWMM_Project *sp, int subcatch, int subArea, double x[]);
void    snow_setState(SWMM_Project *sp, int subcatch, int subArea, double x[]);

void    snow_setMeltCoeffs(SWMM_Project *sp, int snowIndex, double season);
void    snow_plowSnow(SWMM_Project *sp, int subcatch, double tStep);
double  snow_getSnowMelt(SWMM_Project *sp, int subcatch, double rainfall,
        double snowfall, double tStep, double netPrecip[]);
double  snow_getSnowCover(SWMM_Project *sp, int subcatch);

//-----------------------------------------------------------------------------
//   Runoff Analyzer Methods
//-----------------------------------------------------------------------------
int     runoff_open(SWMM_Project *sp);
void    runoff_execute(SWMM_Project *sp);
void    runoff_close(SWMM_Project *sp);

//-----------------------------------------------------------------------------
//   Conveyance System Routing Methods
//-----------------------------------------------------------------------------
int     routing_open(SWMM_Project *sp);
double  routing_getRoutingStep(SWMM_Project *sp, int routingModel, double fixedStep);
void    routing_execute(SWMM_Project *sp, int routingModel, double routingStep);
void    routing_close(SWMM_Project *sp, int routingModel);

//-----------------------------------------------------------------------------
//   Output Filer Methods
//-----------------------------------------------------------------------------
int     output_open(SWMM_Project *sp);
void    output_end(SWMM_Project *sp);
void    output_close(void);
void    output_checkFileSize(SWMM_Project *sp);
void    output_saveResults(SWMM_Project *sp, double reportTime);
void    output_readDateTime(SWMM_Project *sp, int period, DateTime *aDate);
void    output_readSubcatchResults(SWMM_Project *sp, int period, int area);
void    output_readNodeResults(SWMM_Project *sp, int period, int node);
void    output_readLinkResults(SWMM_Project *sp, int period, int link);

//-----------------------------------------------------------------------------
//   Groundwater Methods
//-----------------------------------------------------------------------------
int     gwater_readAquiferParams(SWMM_Project *sp, int aquifer, char* tok[], int ntoks);
int     gwater_readGroundwaterParams(SWMM_Project *sp, char* tok[], int ntoks);
int     gwater_readFlowExpression(SWMM_Project *sp, char* tok[], int ntoks);
void    gwater_deleteFlowExpression(SWMM_Project *sp, int subcatch);

void    gwater_validateAquifer(SWMM_Project *sp, int aquifer);
void    gwater_validate(SWMM_Project *sp, int subcatch);

void    gwater_initState(SWMM_Project *sp, int subcatch);
void    gwater_getState(SWMM_Project *sp, int subcatch, double x[]);
void    gwater_setState(SWMM_Project *sp, int subcatch, double x[]);

void    gwater_getGroundwater(SWMM_Project *sp, int subcatch, double evap,
        double infil, double tStep);
double  gwater_getVolume(SWMM_Project *sp, int subcatch);

//-----------------------------------------------------------------------------
//   RDII Methods
//-----------------------------------------------------------------------------
int     rdii_readRdiiInflow(SWMM_Project *sp, char* tok[], int ntoks);
void    rdii_deleteRdiiInflow(SWMM_Project *sp, int node);
void    rdii_initUnitHyd(SWMM_Project *sp, int unitHyd);
int     rdii_readUnitHydParams(SWMM_Project *sp, char* tok[], int ntoks);
void    rdii_openRdii(SWMM_Project *sp);
void    rdii_closeRdii(SWMM_Project *sp);
int     rdii_getNumRdiiFlows(SWMM_Project *sp, DateTime aDate);
void    rdii_getRdiiFlow(int index, int* node, double* q);

//-----------------------------------------------------------------------------
//   Landuse Methods
//-----------------------------------------------------------------------------
int     landuse_readParams(SWMM_Project *sp, int landuse, char* tok[], int ntoks);
int     landuse_readPollutParams(SWMM_Project *sp, int pollut, char* tok[], int ntoks);
int     landuse_readBuildupParams(SWMM_Project *sp, char* tok[], int ntoks);
int     landuse_readWashoffParams(SWMM_Project *sp, char* tok[], int ntoks);

void    landuse_getInitBuildup(SWMM_Project *sp, TLandFactor* landFactor,
        double* initBuildup, double area, double curb);
double  landuse_getBuildup(SWMM_Project *sp, int landuse, int pollut,
        double area, double curb, double buildup, double tStep);

double  landuse_getWashoffLoad(SWMM_Project *sp, int landuse, int p, double area,                //(5.1.008)
        TLandFactor landFactor[], double runoff, double vOutflow);             //(5.1.008)
double  landuse_getAvgBmpEffic(SWMM_Project *sp, int j, int p);
double  landuse_getCoPollutLoad(SWMM_Project *sp, int p, double washoff[]);

//-----------------------------------------------------------------------------
//   Flow/Quality Routing Methods
//-----------------------------------------------------------------------------
void    flowrout_init(SWMM_Project *sp, int routingModel);
void    flowrout_close(int routingModel);
double  flowrout_getRoutingStep(SWMM_Project *sp, int routingModel, double fixedStep);
int     flowrout_execute(SWMM_Project *sp, int links[], int routingModel,
        double tStep);

void    toposort_sortLinks(SWMM_Project *sp, int links[]);
int     kinwave_execute(SWMM_Project *sp, int link, double* qin, double* qout,
        double tStep);

void    dynwave_validate(SWMM_Project *sp);                                                //(5.1.008)
void    dynwave_init(SWMM_Project *sp);
void    dynwave_close(void);
double  dynwave_getRoutingStep(SWMM_Project *sp, double fixedStep);
int     dynwave_execute(SWMM_Project *sp, double tStep);
void    dwflow_findConduitFlow(SWMM_Project *sp, int j, int steps, double omega,
        double dt);
void    qualrout_init(SWMM_Project *sp);
void    qualrout_execute(SWMM_Project *sp, double tStep);

//-----------------------------------------------------------------------------
//   Treatment Methods
//-----------------------------------------------------------------------------
int     treatmnt_open(SWMM_Project *sp);
void    treatmnt_close(void);
int     treatmnt_readExpression(SWMM_Project *sp, char* tok[], int ntoks);
void    treatmnt_delete(SWMM_Project *sp, int node);
void    treatmnt_treat(SWMM_Project *sp, int node, double q, double v,
        double tStep);
void    treatmnt_setInflow(SWMM_Project *sp, double qIn, double wIn[]);

//-----------------------------------------------------------------------------
//   Mass Balance Methods
//-----------------------------------------------------------------------------
int     massbal_open(SWMM_Project *sp);
void    massbal_close(void);
void    massbal_report(SWMM_Project *sp);

void    massbal_updateRunoffTotals(int type, double v);                        //(5.1.008)
void    massbal_updateLoadingTotals(int type, int pollut, double w);
void    massbal_updateGwaterTotals(double vInfil, double vUpperEvap,
        double vLowerEvap, double vLowerPerc, double vGwater);
void    massbal_updateRoutingTotals(SWMM_Project *sp, double tStep);

void    massbal_initTimeStepTotals(SWMM_Project *sp);
void    massbal_addInflowFlow(int type, double q);
void    massbal_addInflowQual(SWMM_Project *sp, int type, int pollut, double w);
void    massbal_addOutflowFlow(double q, int isFlooded);
void    massbal_addOutflowQual(SWMM_Project *sp, int pollut, double mass, int isFlooded);
void    massbal_addNodeLosses(double evapLoss, double infilLoss);
void    massbal_addLinkLosses(double evapLoss, double infilLoss);
void    massbal_addReactedMass(SWMM_Project *sp, int pollut, double mass);
void    massbal_addSeepageLoss(SWMM_Project *sp, int pollut, double seepLoss);                   //(5.1.008)
void    massbal_addToFinalStorage(SWMM_Project *sp, int pollut, double mass);                    //(5.1.008)
double  massbal_getStepFlowError(void);
double  massbal_getRunoffError(SWMM_Project *sp);
double  massbal_getFlowError(SWMM_Project *sp);

//-----------------------------------------------------------------------------
//   Simulation Statistics Methods
//-----------------------------------------------------------------------------
int     stats_open(SWMM_Project *sp);
void    stats_close(SWMM_Project *sp);
void    stats_report(SWMM_Project *sp);

void    stats_updateCriticalTimeCount(int node, int link);
void    stats_updateFlowStats(SWMM_Project *sp, double tStep, DateTime aDate,
        int stepCount, int steadyState);
void    stats_updateSubcatchStats(SWMM_Project *sp, int subcatch, double rainVol, double runonVol,
        double evapVol, double infilVol, double runoffVol, double runoff);
void    stats_updateGwaterStats(SWMM_Project *sp, int j, double infil, double evap,              //(5.1.008)
        double latFlow, double deepFlow, double theta, double waterTable,      //(5.1.008)
        double tStep);                                                         //(5.1.008)
void    stats_updateMaxRunoff(SWMM_Project *sp);
void    stats_updateMaxNodeDepth(int node, double depth);                      //(5.1.008)

//-----------------------------------------------------------------------------
//   Raingage Methods
//-----------------------------------------------------------------------------
int      gage_readParams(SWMM_Project *sp, int gage, char* tok[], int ntoks);
void     gage_validate(SWMM_Project *sp, int gage);
void     gage_initState(SWMM_Project *sp, int gage);
void     gage_setState(SWMM_Project *sp, int gage, DateTime aDate);
double   gage_getPrecip(SWMM_Project *sp, int gage, double *rainfall,
        double *snowfall);
void     gage_setReportRainfall(SWMM_Project *sp, int gage, DateTime aDate);
DateTime gage_getNextRainDate(SWMM_Project *sp, int gage, DateTime aDate);

//-----------------------------------------------------------------------------
//   Subcatchment Methods
//-----------------------------------------------------------------------------
int     subcatch_readParams(SWMM_Project *sp, int subcatch, char* tok[], int ntoks);
int     subcatch_readSubareaParams(SWMM_Project *sp, char* tok[], int ntoks);
int     subcatch_readLanduseParams(SWMM_Project *sp, char* tok[], int ntoks);
int     subcatch_readInitBuildup(SWMM_Project *sp, char* tok[], int ntoks);

void    subcatch_validate(SWMM_Project *sp, int subcatch);
void    subcatch_initState(SWMM_Project *sp, int subcatch);
void    subcatch_setOldState(SWMM_Project *sp, int subcatch);

double  subcatch_getFracPerv(SWMM_Project *sp, int subcatch);
double  subcatch_getStorage(SWMM_Project *sp, int subcatch);
double  subcatch_getDepth(SWMM_Project *sp, int subcatch);
double  subcatch_getBuildup(SWMM_Project *sp, int subcatch, int pollut);

void    subcatch_getRunon(SWMM_Project *sp, int subcatch);
void    subcatch_addRunonFlow(SWMM_Project *sp, int subcatch, double flow);                      //(5.1.008)
double  subcatch_getRunoff(SWMM_Project *sp, int subcatch, double tStep);

double  subcatch_getWtdOutflow(SWMM_Project *sp, int subcatch, double wt);
void    subcatch_getResults(SWMM_Project *sp, int subcatch, double wt, float x[]);

////  New functions added to release 5.1.008.  ////                            //(5.1.008)
//-----------------------------------------------------------------------------
//  Surface Pollutant Buildup/Washoff Methods
//-----------------------------------------------------------------------------
void    surfqual_initState(SWMM_Project *sp, int subcatch);
void    surfqual_getWashoff(SWMM_Project *sp, int subcatch, double runoff, double tStep);
void    surfqual_getBuildup(SWMM_Project *sp, int subcatch, double tStep);
void    surfqual_sweepBuildup(SWMM_Project *sp, int subcatch, DateTime aDate);
double  surfqual_getWtdWashoff(SWMM_Project *sp, int subcatch, int pollut, double wt);

//-----------------------------------------------------------------------------
//   Conveyance System Node Methods
//-----------------------------------------------------------------------------
int     node_readParams(SWMM_Project *sp, int node, int type, int subIndex,
        char* tok[], int ntoks);
void    node_validate(SWMM_Project *sp, int node);

void    node_initState(SWMM_Project *sp, int node);
void    node_initInflow(SWMM_Project *sp, int node, double tStep);
void    node_setOldHydState(SWMM_Project *sp, int node);
void    node_setOldQualState(SWMM_Project *sp, int node);

void    node_setOutletDepth(SWMM_Project *sp, int node, double yNorm, double yCrit, double z);
void    node_setDividerCutoff(int node, int link);

double  node_getSurfArea(SWMM_Project *sp, int node, double depth);
double  node_getDepth(SWMM_Project *sp, int node, double volume);
double  node_getVolume(SWMM_Project *sp, int node, double depth);
//double  node_getPondedDepth(int node, double volume); removed                //(5.1.008)
double  node_getPondedArea(SWMM_Project *sp, int node, double depth);

double  node_getOutflow(SWMM_Project *sp, int node, int link);
double  node_getLosses(SWMM_Project *sp, int node, double tStep);
double  node_getMaxOutflow(SWMM_Project *sp, int node, double q, double tStep);
double  node_getSystemOutflow(SWMM_Project *sp, int node, int *isFlooded);
void    node_getResults(SWMM_Project *sp, int node, double wt, float x[]);

//-----------------------------------------------------------------------------
//   Conveyance System Inflow Methods
//-----------------------------------------------------------------------------
int     inflow_readExtInflow(SWMM_Project *sp, char* tok[], int ntoks);
int     inflow_readDwfInflow(SWMM_Project *sp, char* tok[], int ntoks);
int     inflow_readDwfPattern(SWMM_Project *sp, char* tok[], int ntoks);
int     inflow_setExtInflow(SWMM_Project *sp, int j, int param, int type,
						int tSeries, int basePat, double cf, 
						double baseline, double sf);
int     inflow_validate(SWMM_Project *sp, int param, int type, int tSeries,
						int basePat, double *cf);					
						
void    inflow_initDwfInflow(SWMM_Project *sp, TDwfInflow* inflow);
void    inflow_initDwfPattern(SWMM_Project *sp, int pattern);

double  inflow_getExtInflow(SWMM_Project *sp, TExtInflow* inflow, DateTime aDate);
double  inflow_getDwfInflow(SWMM_Project *sp, TDwfInflow* inflow, int m, int d, int h);
double  inflow_getPatternFactor(SWMM_Project *sp, int p, int month, int day, int hour);

void    inflow_deleteExtInflows(SWMM_Project *sp, int node);
void    inflow_deleteDwfInflows(SWMM_Project *sp, int node);

//-----------------------------------------------------------------------------
//   Routing Interface File Methods
//-----------------------------------------------------------------------------
int     iface_readFileParams(SWMM_Project *sp, char* tok[], int ntoks);
void    iface_openRoutingFiles(SWMM_Project *sp);
void    iface_closeRoutingFiles(SWMM_Project *sp);
int     iface_getNumIfaceNodes(SWMM_Project *sp, DateTime aDate);
int     iface_getIfaceNode(int index);
double  iface_getIfaceFlow(int index);
double  iface_getIfaceQual(int index, int pollut);
void    iface_saveOutletResults(SWMM_Project *sp, DateTime reportDate, FILE* file);

//-----------------------------------------------------------------------------
//   Hot Start File Methods
//-----------------------------------------------------------------------------
int     hotstart_open(SWMM_Project *sp);
void    hotstart_close(SWMM_Project *sp);

//-----------------------------------------------------------------------------
//   Conveyance System Link Methods
//-----------------------------------------------------------------------------
int     link_readParams(SWMM_Project *sp, int link, int type, int subIndex,
            char* tok[], int ntoks);
int     link_readXsectParams(SWMM_Project *sp, char* tok[], int ntoks);
int     link_readLossParams(SWMM_Project *sp, char* tok[], int ntoks);

void    link_validate(SWMM_Project *sp, int link);
void    link_initState(SWMM_Project *sp, int link);
void    link_setOldHydState(SWMM_Project *sp, int link);
void    link_setOldQualState(SWMM_Project *sp, int link);

void    link_setTargetSetting(SWMM_Project *sp, int j);
void    link_setSetting(SWMM_Project *sp, int j, double tstep);
int     link_setFlapGate(SWMM_Project *sp, int link, int n1, int n2, double q);

double  link_getInflow(SWMM_Project *sp, int link);
void    link_setOutfallDepth(SWMM_Project *sp, int link);
double  link_getLength(SWMM_Project *sp, int link);
double  link_getYcrit(SWMM_Project *sp, int link, double q);
double  link_getYnorm(SWMM_Project *sp, int link, double q);
double  link_getVelocity(SWMM_Project *sp, int link, double q, double y);
double  link_getFroude(SWMM_Project *sp, int link, double v, double y);
double  link_getPower(SWMM_Project *sp, int link);
double  link_getLossRate(SWMM_Project *sp, int link, double q, double tStep);                    //(5.1.008)
char    link_getFullState(double a1, double a2, double aFull);                 //(5.1.008)

void    link_getResults(SWMM_Project *sp, int link, double wt, float x[]);

//-----------------------------------------------------------------------------
//   Link Cross-Section Methods
//-----------------------------------------------------------------------------
int     xsect_isOpen(int type);
int     xsect_setParams(SWMM_Project *sp, TXsect *xsect, int type, double p[], double ucf);
void    xsect_setIrregXsectParams(SWMM_Project *sp, TXsect *xsect);
void    xsect_setCustomXsectParams(SWMM_Project *sp, TXsect *xsect);
double  xsect_getAmax(SWMM_Project *sp, TXsect* xsect);

double  xsect_getSofA(SWMM_Project *sp, TXsect* xsect, double area);
double  xsect_getYofA(SWMM_Project *sp, TXsect* xsect, double area);
double  xsect_getRofA(SWMM_Project *sp, TXsect* xsect, double area);
double  xsect_getAofS(SWMM_Project *sp, TXsect* xsect, double sFactor);
double  xsect_getdSdA(SWMM_Project *sp, TXsect* xsect, double area);
double  xsect_getAofY(SWMM_Project *sp, TXsect* xsect, double y);
double  xsect_getRofY(SWMM_Project *sp, TXsect* xsect, double y);
double  xsect_getWofY(SWMM_Project *sp, TXsect* xsect, double y);
double  xsect_getYcrit(SWMM_Project *sp, TXsect* xsect, double q);

//-----------------------------------------------------------------------------
//   Culvert/Roadway Methods                                                   //(5.1.010)
//-----------------------------------------------------------------------------
double  culvert_getInflow(SWMM_Project *sp, int link, double q, double h);
double  roadway_getInflow(SWMM_Project *sp, int link, double dir, double hcrest,
        double h1, double h2);                                                 //(5.1.010)

//-----------------------------------------------------------------------------
//   Force Main Methods
//-----------------------------------------------------------------------------
double  forcemain_getEquivN(SWMM_Project *sp, int j, int k);
double  forcemain_getRoughFactor(SWMM_Project *sp, int j, double lengthFactor);
double  forcemain_getFricSlope(SWMM_Project *sp, int j, double v, double hrad);

//-----------------------------------------------------------------------------
//   Cross-Section Transect Methods
//-----------------------------------------------------------------------------
int     transect_create(SWMM_Project *sp, int n);
void    transect_delete(SWMM_Project *sp);
int     transect_readParams(SWMM_Project *sp, int* count, char* tok[], int ntoks);
void    transect_validate(SWMM_Project *sp, int j);

//-----------------------------------------------------------------------------
//   Custom Shape Cross-Section Methods
//-----------------------------------------------------------------------------
int     shape_validate(TShape *shape, TTable *curve);

//-----------------------------------------------------------------------------
//   Control Rule Methods
//-----------------------------------------------------------------------------
int     controls_create(SWMM_Project *sp, int n);
void    controls_delete(SWMM_Project *sp);
int     controls_addRuleClause(SWMM_Project *sp, int rule, int keyword,
        char* Tok[], int nTokens);
int     controls_evaluate(SWMM_Project *sp, DateTime currentTime,
        DateTime elapsedTime, double tStep);

//-----------------------------------------------------------------------------
//   Table & Time Series Methods
//-----------------------------------------------------------------------------
int     table_readCurve(SWMM_Project *sp, char* tok[], int ntoks);
int     table_readTimeseries(SWMM_Project *sp, char* tok[], int ntoks);

int     table_addEntry(TTable* table, double x, double y);
int     table_getFirstEntry(TTable* table, double* x, double* y);
int     table_getNextEntry(TTable* table, double* x, double* y);
void    table_deleteEntries(TTable* table);

void    table_init(TTable* table);
int     table_validate(TTable* table);
//      table_interpolate now defined in table.c                               //(5.1.008)

double  table_lookup(TTable* table, double x);
double  table_lookupEx(TTable* table, double x);
double  table_intervalLookup(TTable* table, double x);
double  table_inverseLookup(TTable* table, double y);

double  table_getSlope(TTable *table, double x);
double  table_getMaxY(TTable *table, double x);
double  table_getArea(TTable* table, double x);
double  table_getInverseArea(TTable* table, double a);

void    table_tseriesInit(TTable *table);
double  table_tseriesLookup(TTable* table, double t, char extend);

//-----------------------------------------------------------------------------
//   Utility Methods
//-----------------------------------------------------------------------------
double   UCF(SWMM_Project *sp, int quantity);                   // units conversion factor
int      getInt(char *s, int *y);             // get integer from string
int      getFloat(char *s, float *y);         // get float from string
int      getDouble(char *s, double *y);       // get double from string
char*    getTempFileName(SWMM_Project *sp, char *s); // get temporary file name
int      findmatch(char *s, char *keyword[]); // search for matching keyword
int      match(char *str, char *substr);      // true if substr matches part of str
int      strcomp(char *s1, char *s2);         // case insensitive string compare
char*    sstrncpy(char *dest, const char *src,
         size_t maxlen);                      // safe string copy
void     writecon(char *s);                   // writes string to console
DateTime getDateTime(SWMM_Project *sp, double elapsedMsec);     // convert elapsed time to date
void     getElapsedTime(SWMM_Project *sp, DateTime aDate,       // convert elapsed date
         int* days, int* hrs, int* mins);
void     getSemVersion(char* semver);         // get semantic version


#endif //FUNCS_H
