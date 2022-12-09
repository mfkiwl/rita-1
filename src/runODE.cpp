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

                         Implementation of function runODE

  ==============================================================================*/

#include "rita.h"
#include "data.h"
#include "calc.h"
#include "cmd.h"
#include "configure.h"

using std::cout;
using std::exception;

namespace RITA {

int rita::runODE()
{
   _analysis_type = TRANSIENT;
   bool vector_ok=false;
   string ode_name="", phase="";
   int size=1, ret=0, count_fct=0, count_vector=0, count_def=0, count_init=0, ind=-1;
   _init_time = 0., _time_step=0.1, _final_time=1.;

   vector<string> def, name, var;
   vector<double> init;
   odae *ode = new odae;
   ode->isSet = false;
   _analysis_type = TRANSIENT;
   string str="", var_name="y", scheme="forward-euler", fn="";
   static const vector<string> kw {"size","func$tion","def$inition","var$iable","vect$or","init$ial",
                                   "final$-time","time-step","scheme","phase","summary","clear"};
   _cmd->set(kw,_gkw);
   for (int k=0; k<_nb_args; ++k) {

      switch (_cmd->getArg("=")) {

         case 0:
            size = _cmd->int_token();
            break;

         case 1:
            name.push_back(_cmd->string_token());
            count_fct++, vector_ok = true;
            break;

         case 2:
            def.push_back(_cmd->string_token());
            name.push_back(" ");
            count_def++;
            break;

         case 3:
         case 4:
            var_name = _cmd->string_token();
            vector_ok = true;
            count_vector++;
            break;

         case 5:
            init.push_back(_cmd->double_token());
            count_init++;
            break;

         case 6:
            _final_time = _cmd->double_token();
            break;

         case 7:
            _time_step = _cmd->double_token();
            break;

         case 8:
            scheme = _cmd->string_token();
            break;

         case 9:
            phase = _cmd->string_token();
            break;

         default:
            msg("ode>","Unknown argument: "+_cmd->Arg());
            return 1;
      }
   }

   if (_nb_args>0) {
      if (size<=0) {
         msg("ode>","Illegal size value.");
         return 1;
      }
      if (count_fct && count_vector) {
         msg("ode>","Function already defined in data module.");
         return 1;
      }
      if (count_vector>1) {
         msg("ode>","Only one variable must be defined for an ode system.");
         return 1;
      }
      if (count_fct && count_def) {
         msg("ode>","Function already defined.");
         return 1;
      }
      if (count_fct>size || count_def>size) {
         msg("ode>","Number of function names is larger than system size.");
         return 1;
      }
      if (_data->nb_ode>0 && count_vector==0) {
         msg("ode>","No variable defined as unknown for new ode system.");
         return 1;
      }
      if (!vector_ok) {
         msg("ode>","Missing a variable (or vector) name.");
         return 1;
      }
      if (count_init>size) {
         msg("ode>","Number of initial conditions is larger than system size.");
         return 1;
      }
      if (count_fct>0 && count_def<size-1) {
         msg("ode>","Number of function definitions is larger than system size.");
         return 1;
      }
      if (count_fct && size>1) {
         msg("ode>","The option 'function' is not available with an ODE system.");
         return 1;
      }
      if (count_init<size) {
         for (int i=count_init; i<size; ++i)
            init.push_back(0.);
      }
      *ofh << "ode";
      ode->theFct.resize(size);
      if (count_fct>0) {
         ode->isFct = true;
         if (count_fct<=size-1) {
            msg("ode>","Number of function names is lower than system size.");
            return 1;
         }
         for (int k=0; k<size; ++k) {
            int n = _data->FctName[name[k]];
            if (n==0) {
               msg("ode>","Non defined function "+name[k]);
               return 1;
            }
            if (ode->theFct[k].set(name[k],_data->theFct[n]->expr,_data->theFct[n]->var,1)) {
               msg("ode>","Error in function evaluation: "+ode->theFct[n].getErrorMessage());
               return 1;
            }
            *ofh << " function=" << name[k];
         }
      }
      else {
         ode->isFct = false;
         *ofh << " var=" << var_name;
         var.resize(size+1);
         var[0] = "t";
         var[1] = var_name;
         if (size>1) {
            for (int i=1; i<=size; ++i)
               var[i] = var_name + to_string(i+1);
         }
         for (int i=0; i<size; ++i) {
            if (ode->theFct[i].set(name[i],def[i],var,1)) {
               msg("ode>","Error in function evaluation: "+ode->theFct[i].getErrorMessage());
               return 1;
            }
         }
         ode->ind_fct = ind;
         for (int j=0; j<size; ++j)
            *ofh << " definition=" << def[j];
      }
      ode->size = size;
      ode->isSet = true;
      ode->log = false;
      _data->addVector(var_name,0.,size);
      _data->VectorEquation[_data->iVector] = _data->iEq;
      ode->vect = _data->iVector;
      ode->fn = var_name;
      ode->phase = phase;
      if (phase!="")
         _data->addVector(phase,0,size);
      ode->log = false;
      _data->VectorType.push_back(data::eqType::ODE);
      for (const auto& v: init) {
         *ofh << " init=" << v;
         ode->y.push_back(v);
      }
      ode->scheme = _sch[scheme];
      *ofh << " scheme=" << scheme;
      *ofh << " time-step=" << _time_step << " final-time=" << _final_time << endl;
      ode->type = DataType::ODE;
      _data->addODE(ode,ode_name);
   }

   else {
      *ofh << "ode" << endl;
      for (;;) {
         if (_cmd->readline(sPrompt+"ode> ")<0)
            continue;
         _nb_args = _cmd->getNbArgs();
         int key = _cmd->getKW(kw,_gkw,_data_kw);
         if (key>=200) {
            _data->setDataExt(key);
            continue;
         }
         switch (key) {

            case   0:
               if (_cmd->setNbArg(1,"Size of differential system to be given.")) {
                  msg("ode>size>","Missing system size.","",1);
                  break;
               }
               if (!_cmd->get(size))
                  *ofh << "  size " << size << endl;
               _ret = 0;
               break;

            case   1:
               if (count_def>0) {
                  msg("ode>function>","Function already defined by an expression.");
                  break;
               }
               if (count_fct==size)
                  msg("ode>function>","Too many functions defining ODE.");
               if (_cmd->setNbArg(1,"Function F to define equation y'(t) = F(t,y(t)) to be given.")) {
                  msg("ode>function>","Missing function expression.","",1);
                  break;
               }
               for (int i=0; i<_data->nb_fcts; ++i) {
                  if (_data->theFct[i]->name==name[count_fct])
                     ind = i;
               }
               if (ind==-1) {
                  msg("ode>function>","Non defined function "+name[count_fct]);
                  ret = 1;
                  break;
               }
               ret = _cmd->get(name[count_fct]); 
               if (!ret) {
                  *ofh << "    function " << name[count_fct++] << endl;
                  vector_ok = true;
                  ode->theFct[count_fct].name = name[count_fct-1];
                  ode->theFct[count_fct].set(_data->theFct[ind]->expr,_data->theFct[ind]->var);
                  count_fct++;
               }
               _ret = 0;
               break;

            case   2:
               if (count_fct>0) {
                  msg("ode>definition>","Function already defined by its name.");
                  break;
               }
               if (count_def==size)
                  msg("ode>definition>","Too many functions defining ODE.");
               if (_cmd->setNbArg(1,"Function F to define equation y'(t) = F(t,y(t)) to be given.")) {
                  msg("ode>definition>","Missing function expression.","",1);
                  break;
               }
               if (!_cmd->get(str)) {
                  def.push_back(str);
                  *ofh << "    definition " << str << endl;
                  count_def++;
               }
               _ret = 0;
               break;

            case   3:
            case   4:
               if (ode->log) {
                  msg("ode>variable>","pde must be set first.");
                  break;
               }
               if (_cmd->setNbArg(1,"Give name of associated variable (vector).")) {
                  msg("ode>variable>:","Missing name of associated variable (vector).","",1);
                  break;
               }
               if (!_cmd->get(var_name)) {
                  vector_ok = true, count_vector++;
                  ode->fn = var_name;
                  *ofh << "  variable " << var_name << endl;
               }
               else
                  msg("ode>variable>","Unknown variable (vector) "+var_name);
               break;

            case   5:
               if (_cmd->setNbArg(size,"Initial conditions to be given.")) {
                  msg("ode>initial>","Missing initial conditions.","",1);
                  break;
               }
               ret = 0;
               init.resize(size);
               count_init++;
               for (int i=0; i<size; ++i)
                  ret += _cmd->get(init[i]);
               if (!ret) {
                  *ofh << "  initial  ";
                  for (const auto& v: init)
                     *ofh << v << "  ";
                  *ofh << endl;
               }
               else
                  msg("ode>initial>","Error in initial data.");
               _ret = 0;
               break;

            case   6:
               if (_cmd->setNbArg(1,"Final time to be given.")) {
                  msg("ode>final-time>","Missing final time value.","",1);
                  break;
               }
               if (!_cmd->get(_final_time))
                  *ofh << "    final-time " << _final_time << endl;
               _ret = 0;
               break;

            case   7:
               if (_cmd->setNbArg(1,"Missing time step.")) {
                  msg("ode>time-step>","Missing time step value.","",1);
                  break;
               }
               if (!_cmd->get(_time_step))
                  *ofh << "    time-step " << _time_step << endl;
               _ret = 0;
               break;

            case   8:
               if (_cmd->setNbArg(1,"Time integration scheme to be given.")) {
                  msg("ode>scheme>","Missing time integration scheme.","",1);
                  break;
               }
               if (!_cmd->get(scheme))
                  *ofh << "    scheme " << scheme << endl;
               _ret = 0;
               break;

            case   9:
               if (_cmd->setNbArg(1,"Name of vector to contain phase.")) {
                  msg("ode>phase>","Missing phase variable name.","",1);
                  break;
               }
               if (!_cmd->get(phase))
                  *ofh << "    phase " << phase << endl;
               _ret = 0;
               break;

            case  10:
               cout << "Summary of ODE attributes:\n";
               *ofh << "    summary" << endl;
               if (size>1)
                  cout << "Size:           " << size << endl;
               if (ode->isFct) {
                  if (size==1)
                     cout << "Equation defined by function: " << _ae->theFct[0].name << endl;
                  else {
                     for (int i=0; i<size; ++i)
                        cout << "Equation: " << i+1 << ", defined by function: " << ode->theFct[i].name << endl;
                  }
               }
               else {
                  if (size==1) {
                     cout << "Equation defined by: " << ode->theFct[0].expr << endl;
                     cout << "Variable is          " << ode->theFct[0].var[0] << endl;
                  }
                  else {
                     for (int i=0; i<size; ++i)
                        cout << "Equation: " << i+1 << ", defined by: " << ode->theFct[i].expr << endl;
                     for (int i=0; i<size; ++i)
                        cout << "Variable " << i+1 << " is " << ode->theFct[0].var[i] << endl;
                  }
               }
               _ret = 0;
               break;

            case  11:
               ode->log = false;
               cout << "ODE equation removed from model." << endl;
               *ofh << "    clear" << endl;
               _ret = 10;
               return _ret;

            case 100:
            case 101:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands\n";
               cout << "size:       Size of differential system: (Number of equations)\n";
               cout << "function:   Give already defined function defining ode\n";
               cout << "definition: Give expression defining ode\n";
               cout << "variable:   Variable (or vector) name as unknown of the equation\n";
               cout << "initial:    Give an initial condition\n";
               cout << "final-time: Give final time\n";
               cout << "time-step:  Give time step\n";
               cout << "scheme:     Time integration scheme\n";
               cout << "phase:      Give name of a vector that will contain phase\n";
               cout << "summary:    Summary of ODE attributes\n";
               cout << "clear:      Remove ODE from model\n" << endl;
               break;

            case 102:
               getLicense();
               break;

            case 103:
               _ret = _configure->run();
               break;

            case 104:
            case 105:
               if (count_fct==0 && count_def==0) {
                  NO_ODE
                  *ofh << "  end" << endl;
                  return 0;
               }
               if ((count_fct>0 && count_fct<size) || (count_fct==0 && count_def<size)) {
                  msg("ode>end>","Insufficient number of functions defining system.");
                  *ofh << "  end" << endl;
                  break;
               }
               if (!vector_ok) {
                  msg("ode>end>","No variable defined for ode system.");
                  *ofh << "  end" << endl;
                  break;
               }
               if (count_fct && count_vector) {
                  msg("ode>end>","Function already defined in data module.");
                  return 1;
               }
               if (count_vector>1) {
                  msg("ode>end>","Only one variable must be defined for an ode system.");
                  return 1;
               }
               *ofh << "  end" << endl;
               var.clear();
               var.push_back("t");
               if (size==1)
                  var.push_back(var_name);
               else {
                  for (int i=1; i<=size; ++i)
                     var.push_back(var_name+to_string(i));
               }
               if (!count_init) {
                  init.resize(size);
                  for (int i=0; i<size; ++i)
                     init[i] = 0.;
               }
               ode->theFct.resize(size);
               ode->size = size;
               ode->ind_fct = ind;
               ode->isFct = count_fct;
               ode->y.resize(size);
               ode->scheme = _sch[scheme];
               if (phase!="") {
                  ode->phase = phase;
                  _data->addVector(phase,0,size);
               }
               for (int j=0; j<size; ++j) {
                  ode->y[j] = init[j];
                  if (!count_fct) {
                     if (ode->theFct[j].set(def[j],var,1)) {
                        msg("ode>end>","Error in function evaluation: "+ode->theFct[j].getErrorMessage());
                        break;
                     }
                  }
               }
               _data->addVector(var_name,0.,size);
               ode->vect = _data->iVector;
               _data->VectorType.push_back(data::eqType::ODE);
               ode->isSet = true;
               ode->log = false;
               _data->VectorEquation[_data->iVector] = _data->iEq;
               ode->isFct = false;
               if (count_fct)
                  ode->isFct = true;
               _data->addODE(ode,ode_name);
               _ret = 0;
               return _ret;

            case -2:
            case -3:
            case -4:
               break;

            default:
               _ret = _calc->run();
               break;
         }
      }
   }
   _ret = 0;
   return _ret;
}

} /* namespace RITA */
