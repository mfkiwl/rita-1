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

                        Definition of class 'transient'

  ==============================================================================*/

#pragma once

#include "mesh/Mesh.h"
#include "solvers/LinearSolver.h"
#include "solvers/TimeStepping.h"
#include "OFELI.h"
#include "rita.h"
#include "solve.h"
#include <map>

namespace RITA {

class equa;
struct odae;
 
class transient
{

 public:

    transient(rita *r);
    ~transient();
    void setLinearSolver(OFELI::Iteration ls, OFELI::Preconditioner prec);
    void setSave(vector<int>& isave, vector<int>& fformat, vector<string>& save_file,
                 bool phase, vector<string>& phase_file);
    int run();

 private:

    bool _phase;
    rita *_rita;
    data *_data;
    double _init_time, _final_time, _time_step;
    int _nb_vectors, _nb_ae, _nb_ode, _nb_pde, _rs;
    vector<int> *_fformat, *_isave;
    vector<string> *_save_file, *_phase_file;
    odae *_ae_eq, *_ode_eq;
    equa *_pde_eq;
    int setPDE(OFELI::TimeStepping& ts, int e);
};

} /* namespace RITA */
