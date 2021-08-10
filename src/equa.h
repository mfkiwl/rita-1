/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 Rachid Touzani

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

                           Definition of class 'equa'

  ==============================================================================*/

#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <map>
using std::map;

#include "data.h"

#include "Laplace.h"
#include "Therm.h"
#include "Solid.h"
#include "Fluid.h"
#include "Electromagnetics.h"
#include "io/Fct.h"

using namespace OFELI;

namespace RITA {

class cmd;
class rita;

#define MAX_NB_FIELDS 6

class equa
{

 public:

    struct Log {
      bool pde, field, spd, ls, nl, mesh;
      Log() { pde = field = spd = ls = nl = mesh = false; }
      bool fail() { return (pde || field || spd || ls || nl || mesh); }
    };

    struct PdeData {
       map<int,string> cexp;
       string exp;
       string in_file, out_file;
       int size;
       PdeData() { exp=""; in_file=""; out_file=""; size=0; };
    };

    struct FieldData {
       int            field, nb_dof;
       string         fn;
       data::dataSize ds; 
    };

    equa(rita *r);
    ~equa();

    enum pde_eq {
       LAPLACE,
       HEAT,
       WAVE,
       TRANSPORT,
       LINEAR_ELASTICITY,
       TRUSS,
       BEAM,
       INCOMPRESSIBLE_NAVIER_STOKES,
       COMPRESSIBLE_EULER,
       COMPRESSIBLE_NAVIER_STOKES,
       INCOMPRESSIBLE_POROUS_1PHASE,
       EDDY_CURRENTS,
       MAXWELL,
       HELMHOLTZ,
       EIKONAL
    };

    enum sdm {
       FD,
       FE_P1,
       FE_P2,
       FE_Q1,
       FV,
       DG
    };

    int nb_fields, pde, ieq;
    sdm Sdm;
    string eq, nls;
    bool axi;
    PdeData in_data, bc_data, bf_data, sf_data;
    Iteration ls;
    Preconditioner prec;
    vector<string> analytic;
    FieldData fd[MAX_NB_FIELDS];
    Equa *theEquation;
    void setFields();
    void set(string e);
    int setSpD(string spd);
    int set(Grid* gr);
    int setEq();
    void set();
    void set(data *d);
    int setCoef();
    int getIn();
    int getBC();
    int getBF();
    int getSF();
    int setIn();
    int setBC();
    int setBF();
    int setSF();
    void check();
    void set(cmd* cmd) { _cmd = cmd; }
    void setNodeBC(int code, string exp, double t, Vect<double>& v);
    void setSize(Vect<double>& v, data::dataSize s);
    Log log;
    bool set_u, set_bc, set_bf, set_sf, set_in, set_coef;
    Vect<double> u, b, bc, bf, sf, *theSolution[5];
    string regex_u;
    map<string,sdm> pde_sdm {{"fd",FD},{"feP1",FE_P1},{"feP2",FE_P2},{"feQ1",FE_Q1},{"fv",FV},{"dg",DG}};

 private:

    rita *_rita;
    int _verb, _dim, _ret, _nb_dof, _nb_fields;
    cmd *_cmd;
    bool _rho_set, _Cp_set, _kappa_set, _mu_set, _sigma_set, _Mu_set, _epsilon_set, _omega_set;
    bool _beta_set, _v_set, _young_set, _poisson_set;
    map<string,int> pde_map = {{"laplace",LAPLACE},
                               {"heat",HEAT},
                               {"wave",WAVE},
                               {"transport",TRANSPORT},
                               {"linear-elasticity",LINEAR_ELASTICITY},
                               {"truss",TRUSS},
                               {"beam",BEAM},
                               {"incompressible-navier-stokes",INCOMPRESSIBLE_NAVIER_STOKES}, 
                               {"compressible-euler",COMPRESSIBLE_EULER},
                               {"incompressible-porous-1phase",INCOMPRESSIBLE_POROUS_1PHASE}};
    const vector<string> _kw = {"expression","value","file","save"};
    data *_data;
    Mesh *_theMesh;
    Grid *_theGrid;
    string _rho_exp, _Cp_exp, _kappa_exp, _mu_exp,_sigma_exp, _Mu_exp, _epsilon_exp, _omega_exp;
    string _beta_exp, _v_exp, _young_exp, _poisson_exp;
    OFELI::Fct _theFct;
};

} /* namespace RITA */
