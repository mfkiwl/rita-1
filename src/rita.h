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

                            Definition of class 'rita'

  ==============================================================================*/

#pragma once

#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>

#include <iostream>
using std::ostream;
using std::cin;
using std::cout;
using std::endl;
using std::cerr;
using std::string;

#include "../muparserx/mpParser.h"
#include "../muparserx/mpDefines.h"
#include "../muparserx/mpTest.h"

#include "help.h"
#include "ritaException.h"
#include "data.h"
#include "OFELI.h"

const string sPrompt = "rita>";

namespace RITA {
/*!
 *  \addtogroup RITA
 *  @{
 */

class cmd;
//class data;
class calc;
class configure;
class transient;
class stationary;
class optim;
class integration;
class approximation;
class eigen;
class help;
class solve;
class mesh;
class equa;

#define CATCH                                                   \
   catch(ritaException &e) {                                    \
      std::cout << "RITA error: " << e.what() << endl;          \
      return 1;			                                        \
   }                                                            \
   catch(runtime_error &e) {                                    \
      std::cout << "RITA Runtime error: " << e.what() << endl;  \
      return 1;                                                 \
   }                                                            \
   catch( ... ) {                                               \
      std::cout << "RITA Unexpected error: " << endl;           \
      return 1;                                                 \
   }

#define NO_AE   cout << "No algebraic equation defined." << endl;
#define NO_ODE  cout << "No ordinary differential equation defined." << endl;
#define NO_PDE  cout << "No partial differential equation defined." << endl;

enum type {
   NONE            = 0,
   APPROXIMATION   = 1,
   INTEGRATION     = 2,
   STEADY_STATE    = 3,
   TRANSIENT       = 4,
   EIGEN           = 5,
   OPTIMIZATION    = 6,
};

enum objective_type {
   ANALYTIC_FUNCTION,
   PDE_BASED
};


class rita
{

 public:

    rita();
    ~rita();
    int run();
    void setVerbose(int verb) { _verb = verb; }
    void setInput(string file, int opt=1);
    void initConfig();
    bool meshOK, solveOK, dataOK;
    ofstream *ofh, *ofl, ocf;
    void finish();

    friend class configure;
    friend class mesh;
    friend class solve;
    friend class transient;
    friend class stationary;
    friend class optim;
    friend class eigen;
    friend class integration;
    friend class approximation;
    friend class data;
    friend class calc;
    friend class equa;

 private:

   bool _load, _obj_analytic;
   data *_data;
   calc *_calc;
   odae *_ae;
   equa *_pde;
   string _script_file, _scheme, _sLine;
   ifstream _icf, *_in;
   cmd *_cmd;
   int _verb, _ret, _opt;
   mesh *_mesh;
   help *_help;
   configure *_configure;
   solve *_solve;
   transient *_transient;
   stationary *_stationary;
   optim *_optim;
   eigen *_eigen;
   approximation *_approx;
   integration *_integration;
   double _init_time, _time_step, _final_time;
   int _adapted_time_step, _nb_eigv, _nb_args;
   bool _eigen_vectors, _default_vector;
   bool _analysis_ok;
   int _dim, _analysis_type;

   void set(cmd* com) { _cmd = com; }
   void setDim(int dim) { _dim = dim; }
   void set(data *d) { _data = d; }
   int runAE();
   int runODE();
   int runPDE();

   void getLicense();
   void Load();
   void setData();
   int setMesh();
   void setStationary();
   void setTransient();
   void setEigen();
   void setOptim();
   void setApproximation();
   void setIntegration();
   void setPDE();
   void setClear();
   void setMesh(OFELI::Mesh* ms);
   void setSolve();
   int set_ls(string ls, string prec);
   int set_nls(string nls);
   int findVector(const string& s);
   int parse();
   void setParse();
   int CheckKeywords();
   void ListConst();
   void ListVar();
   void ListExprVar();
   void msg(const string& loc, const string& m1, const string& m2="", int c=0);

   const vector<string> _rita_kw {"load","unload","stat$ionary","trans$ient","eigen","optim",
                                  "approx$imation","integ$ration","algebraic","ode","pde","solve"};
   const vector<string> _gkw {"?","help","lic$ense","set","end","<"};
   const vector<string> _data_kw {"grid","mesh","vect$or","tab$ulation","func$tion","matr$ix","save",
                                  "remove","desc$ription","hist$ory","data","list","print","="};
   map<string,OFELI::Iteration> Ls = {{"direct",OFELI::DIRECT_SOLVER},
                                      {"cg",OFELI::CG_SOLVER},
                                      {"cgs",OFELI::CGS_SOLVER},
                                      {"bicg",OFELI::BICG_SOLVER},
                                      {"bicg-stab",OFELI::BICG_STAB_SOLVER},
                                      {"gmres",OFELI::GMRES_SOLVER}};
   map<OFELI::Iteration,string> rLs = {{OFELI::DIRECT_SOLVER,"direct"},
                                       {OFELI::CG_SOLVER,"cg"},
                                       {OFELI::CGS_SOLVER,"cgs"},
                                       {OFELI::BICG_SOLVER,"bicg"},
                                       {OFELI::BICG_STAB_SOLVER,"bicg-stab"},
                                       {OFELI::GMRES_SOLVER,"gmres"}};
   map<string,NonLinearIter> NLs = {{"bisection",BISECTION},
                                    {"regula-falsi",REGULA_FALSI},
                                    {"picard",PICARD},
                                    {"secant",SECANT},
                                    {"newton",NEWTON}};
   map<string,OFELI::TimeScheme> _sch = {{"forward-euler",OFELI::FORWARD_EULER},
                                         {"backward-euler",OFELI::BACKWARD_EULER},
                                         {"crank-nicolson",OFELI::CRANK_NICOLSON},
                                         {"heun",OFELI::HEUN},
                                         {"newmark",OFELI::NEWMARK},
                                         {"leap-frog",OFELI::LEAP_FROG},
                                         {"AB2",OFELI::ADAMS_BASHFORTH},
                                         {"RK4",OFELI::RK4},
                                         {"RK3-TVD",OFELI::RK3_TVD},
                                         {"BDF2",OFELI::BDF2},
                                         {"builtin",OFELI::BUILTIN}};
   map<string,OFELI::Preconditioner> Prec = {{"ident",OFELI::IDENT_PREC},
                                             {"diag",OFELI::DIAG_PREC},
                                             {"dilu",OFELI::DILU_PREC},
                                             {"ilu",OFELI::ILU_PREC},
                                             {"ssor",OFELI::SSOR_PREC}};
   map<OFELI::Preconditioner,string> rPrec = {{OFELI::IDENT_PREC,"ident"},
                                              {OFELI::DIAG_PREC,"diag"},
                                              {OFELI::DILU_PREC,"dilu"},
                                              {OFELI::ILU_PREC,"ilu"},
                                              {OFELI::SSOR_PREC,"ssor"}};
};

} /* namespace RITA */
