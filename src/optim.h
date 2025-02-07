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

                          Definition of class 'optim'

  ==============================================================================*/

#pragma once

#include "solvers/OptSolver.h"
#include "rita.h"
#include "solve.h"
#include "io/Fct.h"
#include <map>

namespace RITA {
 
#define NO_OPT  log = true; cout << "No optimization problem defined." << endl;

class optim
{
 public:

    optim(rita *r, cmd* command, configure* config);
    ~optim();
    int set();
    int run();
    OFELI::OptSolver::OptMethod Alg;
    int size, nb_eqc, nb_lec, nb_gec, igrad, ihess, iincons, ieqcons, verbose;
    OFELI::Fct *J_Fct;
    bool G_ok, H_ok, log, solved, lp;
    double penal, b, obj;
    vector<OFELI::Fct *> G_Fct, H_Fct, inC_Fct, eqC_Fct;
    vector<double> init;
    vector<OFELI::Vect<double> *> a_le, a_ge, a_eq;
    OFELI::Vect<double> lb, ub, a, b_eq, b_ge, b_le;
    void print(ostream& s) const;

 private:
    
    rita *_rita;
    int _ret, _verb;
    string _fn, _alg;
    configure *_configure;
    cmd *_cmd;
    data *_data;

    map<string,OFELI::OptSolver::OptMethod> Nopt = {{"gradient",OFELI::OptSolver::GRADIENT},
                                                    {"truncated-newton",OFELI::OptSolver::TRUNCATED_NEWTON},
                                                    {"simulated-annealing",OFELI::OptSolver::SIMULATED_ANNEALING},
                                                    {"nelder-mead",OFELI::OptSolver::NELDER_MEAD},
                                                    {"newton",OFELI::OptSolver::NEWTON}};
};

ostream& operator<<(ostream& s, const optim& o);

} /* namespace RITA */
