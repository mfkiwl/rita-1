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

                            Main calling program and
                          Implementation of class 'rita'

  ==============================================================================*/

#include <fstream>
#include <ctime>

#include "rita.h"
#include "data.h"
#include "calc.h"
#include "cmd.h"
#include "mesh.h"
#include "equa.h"
#include "solve.h"
#include "optim.h"
#include "eigen.h"
#include "integration.h"
#include "approximation.h"
#include "configure.h"

using std::cout;
using std::exception;

const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf,sizeof(buf),"%Y-%m-%d  %X",&tstruct);
    return buf;
}


namespace RITA {

rita::rita()
     : meshOK(false), solveOK(false), dataOK(false), _load(false), _ae(nullptr),
       _script_file(""), _in(nullptr), _verb(1), _ret(0),
       _default_vector(true), _analysis_type(NONE)
{
   _cmd = new cmd(this);
   _configure = new configure(this,_cmd);
   _calc = new calc(this,_cmd);
   _data = new data(this,_cmd,_configure);
   _mesh = new mesh(this,_cmd,_configure);
   _help = new help;
   _solve = new solve(this,_cmd,_configure);
   ofl = _configure->getOStreamLog();
   ofh = _configure->getOStreamHistory();
   _solve->setSave(_configure->getSaveResults());
   _optim = new optim(this,_cmd,_configure);
   _integration = new integration(this,_cmd,_configure);
   _approx = new approximation(this,_cmd,_configure);
   _eigen = new eigen(this,_cmd,_configure);
   _init_time = 0.;
   _final_time = 1.;
   _time_step = 0.1;
   _adapted_time_step = 0;
   _scheme = "backward-euler";
}


rita::~rita()
{
   if (_in!=nullptr)
      delete _in;
   delete _cmd;
   delete _mesh;
   delete _data;
   delete _calc;
   delete _solve;
   delete _optim;
   delete _integration;
   delete _help;
   delete _configure;
   if (_ae!=nullptr)
      delete _ae;
   if (_pde!=nullptr)
      delete _pde;
}


void rita::setInput(string file,
                    int    opt)
{
   _script_file = file;
   _opt = opt;
   Load();
}


int rita::run()
{
   int key=0;
   string fn="", td="";
   for (;;) {

      if (_cmd->readline(sPrompt+" ")<0)
         continue;
      _nb_args = _cmd->getNbArgs();
      key = _cmd->getKW(_rita_kw,_gkw,_data_kw);
      if (key>=200) {
         _data->setDataExt(key);
         continue;
      }
      switch (key) {

         case   0:
            Load();
            break;

         case   1:
            if (_in!=nullptr)
               delete _in;
            _in = nullptr;
            break;

         case 201:
            setMesh();
            break;

         case   2:
            setStationary();
            break;

         case   3:
            setTransient();
            break;

         case   4:
            setEigen();
            break;

         case   5:
            setOptim();
            break;

         case   6:
            setApproximation();
            break;

         case   7:
            setIntegration();
            break;

         case   8:
            _ret = runAE();
            if (_ret && _ae!=nullptr) {
               delete _ae;
              _ae = nullptr;
            }
            break;

         case   9:
            _ret = runODE();
            if (_ret && _ae!=nullptr) {
               delete _ae;
              _ae = nullptr;
            }
            break;

         case  10:
            _ret = runPDE();
            if (_ret && _pde!=nullptr) {
               delete _pde;
              _pde = nullptr;
            }
            break;

         case  11:
            setSolve();
            break;

         case  12:
            setClear();
            break;

         case 100:
            _help->run(0);
            break;

         case 101:
            _help->run(1);
            break;

         case 102:
            getLicense();
            break;

         case 103:
            _ret = _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 104:
         case 105:
            msg("","No higher level available");
            break;

         default:
            _ret = _calc->run();
            break;
      }

   }
   return 0;
}


void rita::finish()
{
   *ofh << "exit" << endl;
   exit(0);
}


void rita::getLicense()
{
   _cmd->setNbArg(0);
   cout << "Copyright (C) 2022, Rachid Touzani\n";
   cout << "rita is free software: you can redistribute it and/or modify\n";
   cout << "it under the terms of the GNU General Public License as published by\n";
   cout << "the Free Software Foundation, either version 3 of the License, or\n";
   cout << "(at your option) any later version.\n\n";
   cout << "rita is distributed in the hope that it will be useful,\n";
   cout << "but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
   cout << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n";
   cout << "GNU General Public License for more details." << endl;
}


void rita::Load()
{
   if (_script_file.size()==0) {
      if (_cmd->setNbArg(1,"Give input script file.")) {
         msg("load>","Missing script file to load.","",1);
         return;
      }
      _cmd->get(_script_file);
   }
   if (_in!=nullptr)
      delete _in, _in = nullptr;
   _in = new ifstream(_script_file);
   if (_in->is_open())
      _cmd->setIFStream(_in);
   else {
      msg("load>","Unable to open file: "+_script_file);
      if (_opt)
         _ret = 1;
      else
         finish();
      return;
   }
   _script_file = "";
   _load = true;
   run();
}

/*
void rita::setData()
{
   if (_verb>1)
      cout << "Entering module 'data' ..." << endl;
   _cmd->setNbArg(0);
   _data->setVerbose(_verb);
   _ret = _data->run();
   if (_ret==0)
      dataOK = true;
   if (_ret==10) {
      dataOK = true;
      setMesh();
   }
   if (_verb>1)
      cout << "Leaving module 'data' ..." << endl;
}*/


int rita::setMesh()
{
   if (_verb>1)
      cout << "Entering module 'mesh' ..." << endl;
   _cmd->setNbArg(0);
   *ofh << "mesh" << endl;
   _mesh->set(_cmd);
   _mesh->setVerbose(_verb);
   _ret = _mesh->run();
   if (_verb>1)
      cout << "Leaving module 'mesh' ..." << endl;
   return _ret;
}


void rita::setStationary()
{
   _analysis_type = STEADY_STATE;
}


void rita::setTransient()
{
   string fn="";
   int ret = 0;
   _analysis_type = TRANSIENT;
   _data->save_freq = 0;
   _scheme = "backward-euler";
   _adapted_time_step = 0;
   double it=_init_time, ft=_final_time, ts=_time_step;
   _ret = 0;
   static const string H = "transient [initial-time=it] [final-time=ft] [time-step=ts]  [scheme=s] [adapted]\n\n"
                           "it: Initial value of time. Default value is 0.\n"
                           "ft: Final (maximal) value of time. Default value is 1.\n"
                           "ts: Time step value. Default value is 0.1.\n"
                           "s: Time integration scheme. This is a string to choose among the values: forward-euler,\n"
                           "   backward-euler, crank-nicolson, heun, newmark, leap-frog, AB2 (Adams-Bashforth, 2nd Order),\n"
                           "   RK4 (Runge-Kutta, 4th Order), RK3-TVD (Runge-Kutta, 3rd order, TVD), BDF2 (Backward Difference\n"
                           "   Formula, 2nd Order), builtin (Any scheme built in the chosen PDE). The default value for this\n"
                           "   argument is backward-euler\n"
                           "adapted: Toggle to choose (or not) adaptive time stepping.\n";
   static const vector<string> kw_scheme {"forward-euler","backward-euler","crank-nicolson","heun","newmark",
                                          "leap-frog","AB2","RK4","RK3-TVD","BDF2","builtin"};
   static const vector<string> kw {"initial$-time","final$-time","time$-step","adapted","scheme","save-every"};
   _cmd->set(kw,_gkw,_data_kw);
   _nb_args = _cmd->getNbArgs();
   if (_nb_args==0) {
      msg("transient>","No argument for command.",H);
      _ret = 1;
      return;
   }
   if (_nb_args<1) {
      msg("transient>","No valid argument for command.");
      _ret = 1;
   }
   for (int k=0; k<_nb_args; ++k) {
      int key = _cmd->getArg();
      if (key>=200) {
         _data->setDataExt(key);
         continue;
      }
      switch (key) {

         case 100:
         case 101:
            cout << H << endl;
            _ret = 0;

         case   0:
            it = _cmd->double_token();
            break;

         case   1:
            ft = _cmd->double_token();
            break;

         case   2:
            ts = _cmd->double_token();
            break;

         case   3:
            _adapted_time_step = _cmd->int_token();
            break;

         case   4:
            _scheme = _cmd->string_token();
            break;

         case   5:
            _data->save_freq = _cmd->int_token();
            break;

         default:
            msg("transient>","Unknown argument: "+_cmd->Arg());
            _ret = 1;
            return;
      }
   }

   if (_nb_args>0) {
      _analysis_type = TRANSIENT;
      _time_step = ts;
      _init_time = it;
      _final_time = ft;
      if (ts<0)
         _time_step = -_time_step, _adapted_time_step = 1;
      if (find(kw_scheme.begin(),kw_scheme.end(),_scheme)==kw_scheme.end()) {
         msg("transient>","The scheme "+_scheme+" is unknown or unimplemented.");
         _ret = 1;
         return;
      }
      *ofh << "transient  initial-time=" << _init_time << "  final-time=" << _final_time
            << "  time-step=" << _time_step << "  adapted=" << _adapted_time_step
            << "  scheme=" << _scheme << endl;
   }
   else {
      if (_verb) {
         cout << "Default values:\n";
         cout << "   Initial time:            " << _init_time << endl;
         cout << "   Final time:              " << _final_time << endl;
         cout << "   Time step:               " << _time_step << endl;
         cout << "   Adapted time step ?      " << _adapted_time_step << endl;
         cout << "   Time integration scheme: " << _scheme << endl;
      }
      _cmd->setNbArg(0);
      *ofh << "transient" << endl;
      while (1) {
         int nb = _cmd->readline(sPrompt+"transient> ");
         if (nb<0)
            continue;
         int key = _cmd->getKW(kw,_gkw,_data_kw);
         if (key>=200) {
            _data->setDataExt(key);
            continue;
         }
         switch (key) {

            case   0:
               if (_cmd->setNbArg(1,"Initial time to be given.")) {
                  msg("transient>initial-time>","Missing initial time value.","",1);
                  break;
               }
               ret = _cmd->get(it);
               if (!ret)
                  *ofh << "  initial-time " << it << endl;
               break;

            case   1:
               if (_cmd->setNbArg(1,"Final time to be given.")) {
                  msg("transient>final-time>","Missing final time value.","",1);
                  break;
               }
               ret = _cmd->get(ft);
               if (!ret)
                  *ofh << "  final-time " << ft << endl;
               break;

            case   2:
               if (_cmd->setNbArg(1,"Time step to be given.")) {
                  msg("transient>time-step>","Missing time step value.","",1);
                  break;
               }
               ret = _cmd->get(ts);
               if (!ret)
                  *ofh << "  time-step " << ts << endl;
               break;

            case   3:
               if (_cmd->setNbArg(1,"Time integration scheme.")) {
                  msg("transient>scheme>","Missing time integration scheme.","",1);
                  break;
               }
               ret = _cmd->get(kw_scheme,_scheme);
               if (ret<0) {
                  msg("transient>scheme>","Unknown time integration scheme.",
                      "Unknown time integration scheme.\n"
                      "Available values: forward-euler, backward-euler, crank-nicolson\n"
                      "                  heun, newmark, leap-frog, AB2, RK4, RK3-TVD, BDF2, builtin");
                  break;
               }
               *ofh << "  scheme " << _scheme << endl;
               break;

            case 100:
            case 101:
               cout << "\nAvailable Commands:\n";
               cout << "initial-time: Initial time value (Default = 0.)\n";
               cout << "final-time:   Final time value (Default = 1.)\n";
               cout << "time-step:    Time step (Default = 0.1, < 0: Adapted)\n";
               cout << "scheme:       Time integration scheme (Default = backward-euler)\n";
               cout << "              Available schemes: forward-euler, backward-euler, crank-nicolson, heun,\n";
               cout << "                                 newmark, leap-frog, AB2, RK4, RK3-TVD, BDF2, builtin\n";
               cout << "end or <:     back to higher level" << endl;
               break;

            case 102:
               getLicense();
               break;

            case 103:
               _ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 108:
            case 109:
               *ofh << "  end" << endl;
               _analysis_type = TRANSIENT;
               _time_step = ts;
               _init_time = it;
               _final_time = ft;
               _adapted_time_step = 0;
               if (ts<0)
                  _time_step = -_time_step, _adapted_time_step = 1;
               _ret = 0;
               return;

            case -2:
            case -3:
            case -4:
               break;

            default:
               msg("transient>","Unknown command "+_cmd->token(),
                   "Available commands: initial-time, final-time, time-step, scheme, end, <\n"
                   "Global commands:    help, ?, license, end, <");
               break;
         }
      }
   }
}


void rita::setOptim()
{
   if (_verb>1)
      cout << "Entering module 'optimization' ..." << endl;
   _optim->run();
   if (_verb>1)
      cout << "Leaving module 'optimization' ..." << endl;
}


void rita::setEigen()
{
   if (_verb>1)
      cout << "Entering module 'eigen' ..." << endl;
   _analysis_type = EIGEN;
   _eigen->run();
   if (_verb>1)
      cout << "Leaving module 'eigen' ..." << endl;
}


void rita::setPDE()
{
/*   if (!_mesh->MeshIsOK()) {
      msg("pde>","No valid mesh created.");
      return;
   }*/
   string pde_name;
   int ret=0;
   if (_cmd->setNbArg(1,"Give PDE name.")) {
      msg("pde>","Missing pde name.","",1);
      _ret = 1;
      return;
   }
   static const vector<string> kw {"laplace","heat","wave","transport","linear-elasticity",
                                   "truss","beam","incompressible-navier-stokes","clear"};
   ret = _cmd->get(kw,pde_name);
   if (ret<0) {
      msg("pde>","Unknown pde "+pde_name,
          "Unknown pde name. Available pde's are:\n"
          "laplace, heat, wave, transport, linear-elasticity");
      _ret = 1;
      return;
   }
   _pde = new equa(this);
   _pde->set(pde_name);
   *ofh << "pde " << pde_name << endl;
   _ret = runPDE();
   if (_ret && _pde!=nullptr) {
      delete _pde; 
      _pde = nullptr;
   }
   else {
      _pde->set();
   }
}


void rita::setApproximation()
{
   if (!_approx->run())
      _approx->go();
   _ret = 0;
   return;
}


void rita::setIntegration()
{
   if (!_integration->run())
      _integration->go();
   _ret = 0;
   return;
}


int rita::set_nls(string nls)
{
   auto it = NLs.find(nls);
   if (it==NLs.end()) {
      msg("equation>pde>nls>","Unknown nonlinear iterative solver: "+nls);
      return 1;
   }
   return 0;
}


int rita::set_ls(string ls,
                 string prec)
{
   auto it1 = Ls.find(ls);
   if (it1==Ls.end()) {
      msg("equation>pde>ls>","Unknown linear solver: "+ls);
      return 1;
   }
   auto it2 = Prec.find(prec);
   if (it2==Prec.end()) {
      msg("equation>pde>ls>","Unknown linear preconditioner: "+prec);
      return 1;
   }
   return 0;
}


void rita::setClear()
{
}


void rita::setSolve()
{
   if (_verb>1)
      cout << "Entering module 'solve' ..." << endl;
   _cmd->setNbArg(0);
   *ofh << "solve" << endl;
   _solve->setVerbose(_verb);
   _ret = _solve->run();
   if (_ret==0)
      solveOK = true;
   if (_ret>=500)
      return;
   if (_verb>1)
      cout << "Leaving module 'solve' ..." << endl;
}


void rita::msg(const string& loc, const string& m1, const string& m2, int c)
{
   if (c==0) {
      cout << "Error: " << m1 << endl;
      if (m2!="")
         cout << m2 << endl;
   }
   *ofl << "In " + sPrompt+loc + ": " << m1 << endl;
}

} /* namespace RITA */
