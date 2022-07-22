#ifndef CALLBACKS_HPP
#define CALLBACKS_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "Timer.h"
#include "Random.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <list>
#include <set>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <limits>
#include <filesystem>
#include <ilcplex/ilocplex.h>
namespace fs = std::filesystem;

/**************** MIP INFO CALLBACK ********************

Informational Callback: Does not affect performance 
(reasonable use of informational callbacks makes use 
of only the functions explicitly documented in the 
documentation)

Calls the user-written callback regularly during the
branch-and-cut search.

********************************************************/
ILOMIPINFOCALLBACK0(mipInfoCallback)
{
    IloNum bestObjValue = getBestObjValue();
    IloNum cutOff = getCutoff();
    IloNum incumbentObjValue = getIncumbentObjValue();
    IloNum mipRelativeGap = getMIPRelativeGap();
    IloInt nIterations = getNiterations();
    IloInt nNodes = getNnodes();
    IloInt nRemainingNodes = getNremainingNodes();
    IloBool hasIncumbnt = hasIncumbent();

    cout << "Best objective function value = " << bestObjValue << "\n";
    cout << "Cutoff = " << cutOff << "\n";
    if(hasIncumbnt) cout << "Incumbent solution = " << incumbentObjValue << "\n";
    cout << "Relative gap = "  << mipRelativeGap << "\n";
    cout << "Number of iterations = " << nIterations << "\n";
    cout << "Processed nodes = " << nNodes << "\n";
    cout << "Remaining nodes = " << nRemainingNodes << "\n\n";
}

/***************** TUNING CALLBACK **********************

Informational Callback: Does not affect performance 
(reasonable use of informational callbacks makes use 
of only the functions explicitly documented in the 
documentation)

This class enables you to access information on the 
progress of tuning.

*********************************************************/

ILOTUNINGCALLBACK0(tuningCallback)
{
    // Fraction of completion of the tuning process.
    IloNum progress = getProgress();
}

/*********** DISJUNCTIVE CUT INFO CALLBACK **************

Informational Callback: Does not affect performance 
(reasonable use of informational callbacks makes use 
of only the functions explicitly documented in the 
documentation)

This class offers a method to check on the progress of 
the generation of disjunctive cuts.

*********************************************************/

/*********** FLOW MIR CUT INFO CALLBACK *****************

Informational Callback: Does not affect performance 
(reasonable use of informational callbacks makes use 
of only the functions explicitly documented in the 
documentation)

This class offers a member function to check on the 
progress of the generation of Flow MIR cuts.

*********************************************************/

/*********** FRACTIONAL CUT INFO CALLBACK ***************

Informational Callback: Does not affect performance 
(reasonable use of informational callbacks makes use 
of only the functions explicitly documented in the 
documentation)

This class offers a method to check on the progress of
the generation of fractional cuts.

*********************************************************/

/************** PROBING INFO CALLBACK *******************

Informational Callback: Does not affect performance 
(reasonable use of informational callbacks makes use 
of only the functions explicitly documented in the 
documentation)

This class offers a method to check on the progress of
a probing operation.

*********************************************************/

/**************** CONTINUOUS CALLBACK ******************

Query/Diagnostic Callback: May slow the process

Calls the user-written callback after each iteration 
during an optimization of a problem solved at a node.

*******************************************************/
ILOCONTINUOUSCALLBACK0(continuousCallback)
{
    /*IloNum dualInfeasibility = getDualInfeasibility();
    IloNum infeasibility = getInfeasibility();
    IloInt nIterations = getNiterations();
    IloNum objValue = getObjValue();
    IloBool isDualFeasibl = isDualFeasible();
    IloBool isFeasibl = isFeasible();

    cout << "Dual infeasibility = " << dualInfeasibility << "\n";
    cout << "Infeasibility = " << infeasibility << "\n";
    cout << "Number of iterations = " << nIterations << "\n";
    cout << "Objective function = " << objValue << "\n";
    if(isDualFeasibl) cout << "Is dual feasible\n";
    else cout << "Is not dual feasible\n";
    if(isFeasibl) cout << "Is primal feasible\n\n";
    else cout << "Is not primal feasible\n\n";*/

    /*IloNum start = getStartDetTime();
    IloNum tiempo = getDetTime();

    cout << "tiempo: " << tiempo - start  << "\n";*/

    //if(nIterations > 5) abort();
}

/***************** SIMPLEX CALLBACK *******************

Query/Diagnostic Callback: May slow the process

Calls the user-written callback after each iteration 
during optimization with the simplex algorithm.

*******************************************************/

ILOSIMPLEXCALLBACK3(simplexCallback, Timer&, timer, double, start, bool*, abortt)
{
    IloNum dualInfeasibility = getDualInfeasibility();
    IloNum infeasibility = getInfeasibility();
    IloInt nIterations = getNiterations();
    IloNum objValue = getObjValue();
    IloBool isDualFeasibl = isDualFeasible();
    IloBool isFeasibl = isFeasible();

    cout << "Dual infeasibility = " << dualInfeasibility << "\n";
    cout << "Infeasibility = " << infeasibility << "\n";
    cout << "Number of iterations = " << nIterations << "\n";
    cout << "Objective function = " << objValue << "\n";
    if(isDualFeasibl) cout << "Is dual feasible\n";
    else cout << "Is not dual feasible\n";
    if(isFeasibl) cout << "Is primal feasible\n\n";
    else cout << "Is not primal feasible\n\n";

    /*double time = timer.elapsed_time(Timer::VIRTUAL);
    if(time - start > 0.5 ){
        *abortt = true;
        abort();
    }*/

   // cout << "tiempo: " << time - start << "\n";
}

/***************** BARRIER CALLBACK *******************

Query/Diagnostic Callback: May slow the process

Calls the user-written callback after each iteration
during optimization with the barrier method.

*******************************************************/

#endif
