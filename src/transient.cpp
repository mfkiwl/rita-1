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

                      Implementation of class 'transient'

  ==============================================================================*/

#include "OFELI_Config.h"
#include "util/util.h"
#include "io/IOField.h"
#include "io/saveField.h"
#include <iostream>

#include "transient.h"
#include "equa.h"
#include "solvers/ODESolver.h"
#include "solvers/NLASSolver.h"

namespace RITA {

transient::transient(rita *r,
                     int  eq)
          : _phase(false), _rita(r), _rs(1)
{
   _eq = eq;
   _data = _rita->_data;
   _nb_fields = _data->nb_fields;
   _nb_eq = _data->nb_eq;
   theTimeStep = _time_step = _rita->_time_step;
   _init_time = _rita->_init_time;
   theFinalTime = _final_time = _rita->_final_time;
   if (_data->eq_type[_eq]==ALGEBRAIC_EQ)
      _ae_eq = _data->theAE[_eq];
   else if (_data->eq_type[_eq]==PDE_EQ)
      _pde_eq = _data->thePDE[_eq];
}


transient::~transient()
{
}


int transient::setPDE(OFELI::TimeStepping& ts)
{
   try {
      ts.set(_rita->_sch[_rita->_scheme],_time_step,_final_time);
      ts.setPDE(*(_pde_eq->theEquation));
      ts.setLinearSolver(_pde_eq->ls,_pde_eq->prec);
      ts.setInitial(*_data->u[_pde_eq->fd[0].field]);
      if (_pde_eq->eq=="incompressible-navier-stokes" && _pde_eq->Sdm==equa::FE_P1)
         _pde_eq->theEquation->setInput(PRESSURE_FIELD,*_data->u[_pde_eq->fd[1].field]);
   } CATCH
   return 0;
}


void transient::setLinearSolver(OFELI::Iteration      ls,
                                OFELI::Preconditioner prec)
{
  //   _ls = ls;
  //   _prec = prec;
}


void transient::setSave(vector<int>&    isave,
                        vector<int>&    fformat,
                        vector<string>& save_file,
                        bool            phase,
                        vector<string>& phase_file)
{
   _isave = &isave;
   _fformat = &fformat;
   _save_file = &save_file;
   _phase = phase;
   _phase_file = &phase_file;
}


int transient::run()
{
   OFELI::Verbosity = 1;
   vector<ofstream> fs(_nb_fields), ffs(_nb_fields), pfs(_nb_fields);
   vector<OFELI::IOField> ff(_nb_fields);
   vector<string> fn(_nb_fields);
   OFELI::ODESolver ode;
   OFELI::NLASSolver nlas;
   OFELI::TimeStepping ts;
   for (int e=0; e<_nb_eq; ++e) {
      if (e==_eq) {
         if (_data->eq_type[e]==ODE_EQ) {
            _ode_eq = _data->theODE[e];
            ode.set(_ode_eq->scheme,_time_step,_final_time);
            ode.setNbEq(_ode_eq->size);
            int f = _ode_eq->field;
            fn[f] = "rita-" + to_string(10*(e+1)+f+1) + ".sol";
         }
         else if (_data->eq_type[e]==ALGEBRAIC_EQ) {
            nlas.set(_ae_eq->nls);
            nlas.setNbEq(_ae_eq->size);
            int f = _data->theAE[e]->field;
            fn[f] = "rita-" + to_string(10*(e+1)+f+1) + ".sol";
         }
         else if (_data->eq_type[e]==PDE_EQ) {
            _pde_eq = _data->thePDE[e];
            setPDE(ts);
            for (int i=0; i<_pde_eq->nb_fields; ++i) {
               int f = _pde_eq->fd[i].field;
               fn[f] = "rita-" + to_string(10*(e+1)+f+1) + ".sol";
            }
         }
      }
   }
//   if (_nb_eq==1) {
      if (_data->eq_type[_eq]==ODE_EQ) {
         _ode_eq = _data->theODE[_eq];
         int f = _ode_eq->field;
         _data->u[f]->resize(_ode_eq->size);
         *_data->u[f] = _ode_eq->y;
         if (_rs) {
            fs[f].open(fn[f].c_str(),std::fstream::out);
            fs[f] << "# Saved by rita: Solution of ODE, equation: 1" << endl;
            fs[f] << 0.;
            for (int i=0; i<_ode_eq->size; ++i)
               fs[f] << "  " << _ode_eq->y[i];
            fs[f] << endl;
         }
         if ((*_isave)[_eq]) {
            ffs[f].open((*_save_file)[f].c_str());
            ffs[f] << "# Saved by rita: Solution of ODE, equation: 1" << endl;
            ffs[f] << 0.;
            for (int i=0; i<_ode_eq->size; ++i)
               ffs[f] << "  " << _ode_eq->y[i];
            ffs[f] << endl;
            if (_phase) {
               pfs[f].open((*_phase_file)[f].c_str());
               pfs[f] << "# Saved by rita: Phase portrait of ODE, equation: 1" << endl;
            }
         }
         if (_ode_eq->size==1)
            ode.setInitial(_ode_eq->y[0]);
         else
            ode.setInitial(_ode_eq->y);
      }
      else if (_data->eq_type[_eq]==PDE_EQ && _rs) {
         _pde_eq = _data->thePDE[_eq];
         for (int i=0; i<_pde_eq->nb_fields; ++i) {
            int f = _pde_eq->fd[i].field;
            ff[f].open(fn[f],OFELI::IOField::OUT);
         }
      }
//   }
/*   else {
      for (int e=0; e<_nb_eq; ++e) {
         int f = _ode_eq[e]->field;
         _data->u[f]->resize(_ode_eq[e]->size);
         *_data->u[f] = _ode_eq[e]->y;
         if ((*_isave)[e]) {
            ffs[f].open((*_save_file)[e].c_str(),std::fstream::out);
            ffs[f] << "  " << _ode_eq[e]->y;
            if ((*_isave)[f]) {
               ffs[f].open((*_save_file)[f].c_str(),std::fstream::out);
               ffs[f] << "# Saved by rita: Solution of ODE, equation: " << e+1 << endl;
               ffs[f] << 0.;
               for (int i=0; i<_rita->_ode[e].size; ++i)
                  ffs[f] << "  " << _ode_eq[e]->y[i];
               ffs[f] << endl;
               if (_phase) {
                  pfs[f].open((*_phase_file)[f].c_str());
                  pfs[f] << "# Saved by rita: Phase portrait of ODE, equation: " << e+1 << endl;
               }
            }
         }
         if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
            if (_algebraic_eq[e]->size==1)
               nlas.setInitial(_ode_eq[e]->y[0]);
            else
               nlas.setInitial(_ode_eq[e]->y);
         }
         if ((*_eq_type)[e]==ODE_EQ) {
            if (_ode_eq[e]->size==1)
               _ode->setInitial(_ode_eq[e]->y[0]);
            else
               _ode->setInitial(_ode_eq[e]->y);
         }
         else if ((*_eq_type)[e]==PDE_EQ && _rs) {
            for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
               int f = _pde_eq[e]->fd[i].field;
               ff[_pde_eq[e]->fd[i].field].open(fn[f],OFELI::IOField::OUT);
            }
         }
      }
   }*/
   theStep = 1;
   if (_data->eq_type[_eq]==ALGEBRAIC_EQ) {
      _ae_eq = _data->theAE[_eq];
      for (int i=0; i<_ae_eq->size; ++i)
         nlas.setf(_ae_eq->theFct[i]);
   }
   else if (_data->eq_type[_eq]==ODE_EQ) {
      _ode_eq = _data->theODE[_eq];
      for (int i=0; i<_ode_eq->size; ++i)
         ode.setF(_ode_eq->theFct[i]);
   }
   else if (_data->eq_type[_eq]==PDE_EQ) {
      _pde_eq = _data->thePDE[_eq];
      ts.setLinearSolver(_pde_eq->ls,_pde_eq->prec);
   }

// Loop on time steps
   try {
      TimeLoop {

         if (_rita->_verb)
            cout << "Performing time step " << theStep <<", Time = " << theTime << endl;

//       Loop on equations to solve
         for (int e=0; e<_nb_eq; ++e) {
            if (e==_eq) {

//             Case of an algebraic equation
               if (_data->eq_type[e]==ALGEBRAIC_EQ) {
//                  int f = _ae_eq->field;
                  cout << "No algebraic equation solver implemented." << endl;
               }

//             Case of an ODE
               else if (_data->eq_type[e]==ODE_EQ) {
                  _ode_eq = _data->theODE[e];
                  int f = _ode_eq->field;
                  OFELI::Vect<double> z(_ode_eq->size);
                  ode.runOneTimeStep();
                  if (_ode_eq->size==1)
                     _ode_eq->y[0] = ode.get();
                  *_data->u[f] = _ode_eq->y;
                  if (_rs) {
                     fs[f] << theTime;
                     for (int i=0; i<_ode_eq->size; ++i)
                        fs[f] << "  " << _ode_eq->y[i];
                     fs[f] << endl;
                  }
                  if ((*_isave)[f]) {
                     if (theStep%(*_isave)[f]==0) {
                        ffs[f] << theTime;
                        for (int i=0; i<_ode_eq->size; ++i)
                           ffs[f] << "  " << _ode_eq->y[i];
                        ffs[f] << endl;
                        if (_phase) {
                           ode.getTimeDerivative(z);
                           pfs[f] << _ode_eq->y[0] << "  ";
                           for (int i=0; i<_ode_eq->size; ++i)
                              pfs[f] << ode.getTimeDerivative(i+1) << "  ";
                           pfs[f] << endl;
                        }
                     }
                  }
               }

//             Case of a PDE
               else if (_data->eq_type[e]==PDE_EQ) {
                  _pde_eq = _data->thePDE[e];

//                Body force
                  if (_pde_eq->set_bf) {
                     _pde_eq->bf.setTime(theTime);
                     if (_pde_eq->bf.withRegex(1))
                        _pde_eq->bf.set(_pde_eq->bf_data.exp);
                     ts.setRHS(_pde_eq->bf);
                     _pde_eq->theEquation->setInput(BODY_FORCE,_pde_eq->bf);
                  }

//                Boundary condition
                  if (_pde_eq->set_bc) {
                     _pde_eq->bc.setTime(theTime);
                     if (_pde_eq->bc.withRegex(1)) {
                        for (auto const& v: _pde_eq->bc_data.cexp)
                           _pde_eq->setNodeBC(v.first,v.second,theTime,_pde_eq->bc);
                     }
                     ts.setBC(_pde_eq->bc);
                  }

//                Boundary force
                  if (_pde_eq->set_sf) {
                     _pde_eq->sf.setTime(theTime);
                     if (_pde_eq->sf.withRegex(1)) {
                        for (auto const& v: _pde_eq->sf_data.cexp)
                           _pde_eq->setNodeBC(v.first,v.second,theTime,_pde_eq->sf);
                     }
            //            ts.setSF(_data->sf[i]);
                  }

//                Run
                  ts.runOneTimeStep();

//                Save in native OFELI format file
                  for (int i=0; i<_pde_eq->nb_fields; ++i) {
                     int f = _pde_eq->fd[i].field;
                     _data->u[f]->setTime(theTime);
                     _data->u[f]->setName(_data->Field[f]);
                     if (_rs)
                        ff[f].put(*_data->u[f]);
                  }
               }
            }
         }
      }
   } CATCH

   theTime -= theTimeStep;
   for (int e=0; e<_nb_eq; ++e) {
      if (e==_eq) {
         if (_data->eq_type[e]==PDE_EQ) {
            _pde_eq = _data->thePDE[e];
            for (int i=0; i<_pde_eq->nb_fields; ++i) {
               int f = _pde_eq->fd[i].field;
               if (_rs) {
                  ff[f].close();
                  OFELI::saveFormat(*(_data->theMesh[0]),fn[f],(*_save_file)[f],(*_fformat)[f],
                                    (*_isave)[f]);
               }
            }
         }
      }
   }
   return 0;
}

} /* namespace RITA */
