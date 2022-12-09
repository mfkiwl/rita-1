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

                        Definition of class 'eigen'

  ==============================================================================*/

#pragma once


#include "solvers/EigenProblemSolver.h"
#include "rita.h"
#include "solve.h"
#include "linear_algebra/Matrix.h"
#include <map>

namespace RITA {
 
#define NO_EIGEN  log = true; cout << "No eigen problem defined." << endl;

class eigen
{
 public:

    eigen(rita *r, cmd* command, configure* config);
    ~eigen();
    int set();
    int run();
    OFELI::EigenMethod Alg;
    int size, nb_eigv, verbose;
    bool eig_vec, symm, log, solved;
    string evect, eval;
    OFELI::Matrix<double> *M;
    void print(ostream& s) const;

 private:

    rita *_rita;
    int _ret, _verb;
    string _alg;
    configure *_configure;
    cmd *_cmd;
    data *_data;

    map<string,OFELI::EigenMethod> meth = {{"subspace",OFELI::SUBSPACE},
                                           {"qr",OFELI::QR}};

};

ostream& operator<<(ostream& s, const eigen& es);

} /* namespace RITA */
