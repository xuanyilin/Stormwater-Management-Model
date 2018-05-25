/*
 * controls.h
 *
 *  Created on: May 2, 2018
 *      Author: mtryby
 */

#ifndef SRC_CONTROLS_H_
#define SRC_CONTROLS_H_


//-----------------------------------------------------------------------------
// Data Structures
//-----------------------------------------------------------------------------
// Rule Premise Variable
struct TVariable
{
   int      node;            // index of a node (-1 if N/A)
   int      link;            // index of a link (-1 if N/A)
   int      attribute;       // type of attribute for node/link
};

// Rule Premise Clause
struct  TPremise
{
    int     type;                 // clause type (IF/AND/OR)
    struct  TVariable lhsVar;     // left hand side variable                   //(5.1.008)
    struct  TVariable rhsVar;     // right hand side variable                  //(5.1.008)
    int     relation;             // relational operator (>, <, =, etc)
    double  value;                // right hand side value
    struct  TPremise *next;       // next premise clause of rule
};

// Rule Action Clause
struct  TAction
{
   int     rule;             // index of rule that action belongs to
   int     link;             // index of link being controlled
   int     attribute;        // attribute of link being controlled
   int     curve;            // index of curve for modulated control
   int     tseries;          // index of time series for modulated control
   double  value;            // control setting for link attribute
   double  kp, ki, kd;       // coeffs. for PID modulated control
   double  e1, e2;           // PID set point error from previous time steps
   struct  TAction *next;    // next action clause of rule
};

// List of Control Actions
struct  TActionList
{
   struct  TAction* action;
   struct  TActionList* next;
};

// Control Rule
struct  TRule
{
   char*    ID;                        // rule ID
   double   priority;                  // priority level
   struct   TPremise* firstPremise;    // pointer to first premise of rule
   struct   TPremise* lastPremise;     // pointer to last premise of rule
   struct   TAction*  thenActions;     // linked list of actions if true
   struct   TAction*  elseActions;     // linked list of actions if false
};

//-----------------------------------------------------------------------------
//  Shared variables
//-----------------------------------------------------------------------------
typedef struct
{
    struct   TRule*       Rules;           // array of control rules
    struct   TActionList* ActionList;      // linked list of control actions
    int      InputState;                   // state of rule interpreter
    int      RuleCount;                    // total number of rules
    double   ControlValue;                 // value of controller variable
    double   SetPoint;                     // value of controller setpoint
    DateTime CurrentDate;                  // current date in whole days
    DateTime CurrentTime;                  // current time of day (decimal)
} TControlsShared;

// Avoid Duplicate symbol error on C++
//DateTime ElapsedTime;                  // elasped simulation time (decimal days)


#endif /* SRC_CONTROLS_H_ */
