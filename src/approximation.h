/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2023 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================

                       Definition of class 'approximation'

  ==============================================================================*/

#pragma once

#include "OFELI.h"
#include "io/Tabulation.h"
#include "rita.h"
#include "cmd.h"
#include "configure.h"
#include <map>

namespace RITA {
 
class approximation
{

 public:

    enum ApproxType {
       LAGRANGE,
       PIECEWISE_LAGRANGE,
       HERMITE,
       FITTING,
       BSPLINE,
       BEZIER,
       NURBS
    };

    enum FitType {
       POLYNOMIAL,
       EXPONENTIAL,
       DEFINED
    };

    approximation(rita* r, cmd* command, configure* config);
    ~approximation();
    int set();
    int run();
    int go();

 private:

    rita *_rita;
    configure *_configure;
    cmd *_cmd;
    data *_data;
    OFELI::Tabulation *_tab;
    ApproxType _method;
    FitType ft;
    vector<double> _x, _y;
    int eval(double x);
    int _verb, _lagrange_degree, _hermite_degree;
    int lagrange();
    int piecewise_lagrange();
    int hermite();
    int fitting();
    int bspline();
    int bezier();
    int nurbs();
    int lagrange_eval(double x);
    const vector<string> _kw {"help","?","set","file","lagrange","fitting","bspline","bezier","nurbs"
                              "end","<","quit","exit"};
    map<ApproxType,string> rApp = {{LAGRANGE,"lagrange"},
                                   {FITTING,"fitting"},
                                   {BSPLINE,"bspline"},
                                   {BEZIER,"bezier"},
                                   {NURBS,"nurbs"}};
};

} /* namespace RITA */
