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

transient::transient(rita *r)
          : _phase(false), _rita(r), _rs(1)
{
   _data = _rita->_data;
   _nb_ae = _data->nb_ae;
   _nb_ode = _data->nb_ode;
   _nb_pde = _data->nb_pde;
   _nb_vectors = _data->nb_vectors;
   theTimeStep = _time_step = _rita->_time_step;
   _init_time = _rita->_init_time;
   theFinalTime = _final_time = _rita->_final_time;
}


transient::~transient()
{
}


int transient::setPDE(OFELI::TimeStepping& ts,
                      int                  e)
{
   try {
      _pde_eq = _data->thePDE[e];
      ts.set(_rita->_sch[_rita->_scheme],_time_step,_final_time);
      ts.setPDE(*(_pde_eq->theEquation));
      ts.setLinearSolver(_pde_eq->ls,_pde_eq->prec);
      ts.setInitial(*_data->theVector[_pde_eq->fd[0].vect]);
      if (_pde_eq->eq=="incompressible-navier-stokes" && _pde_eq->Sdm==equa::FE_P1)
         _pde_eq->theEquation->setInput(PRESSURE_FIELD,*_data->theVector[_pde_eq->fd[1].vect]);
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
   vector<ofstream> fs(_nb_vectors), ffs(_nb_vectors), pfs(_nb_vectors);
   vector<OFELI::IOField> ff(_nb_vectors);
   vector<string> fn(_nb_vectors);
   OFELI::ODESolver ode;
   OFELI::NLASSolver nlas;
   OFELI::TimeStepping ts;

   for (int e=1; e<=_nb_ode; ++e) {
      _ode_eq = _data->theODE[e];
      ode.set(_ode_eq->scheme,_time_step,_final_time);
      ode.setNbEq(_ode_eq->size);
//      int f = _ode_eq->vect;
//      fn[f-1] = "rita-" + to_string(10*e+f) + ".sol";
   }

  for (int e=1; e<=_nb_ae; ++e) {
      _ae_eq = _data->theAE[e];
      nlas.set(_ae_eq->nls);
      nlas.setNbEq(_ae_eq->size);
//      int f = _data->theAE[e]->vect;
//      fn[f-1] = "rita-" + to_string(10*e+f) + ".sol";
   }

   for (int e=1; e<=_nb_pde; ++e) {
      _pde_eq = _data->thePDE[e];
      setPDE(ts,e);
/*      for (int i=0; i<_pde_eq->nb_vectors; ++i) {
         int f = _pde_eq->fd[i].vect;
         fn[f-1] = "rita-" + to_string(10*e+f) + ".sol";
      }*/
   }
//   if (_nb_eq==1) {
      for (int e=1; e<=_nb_ode; ++e) {
         _ode_eq = _data->theODE[e];
         int f = _ode_eq->vect;
         _data->theVector[f]->resize(_ode_eq->size);
         *_data->theVector[f] = _ode_eq->y;
/*         if (_rs) {
            fs[f-1].open(fn[f-1].c_str(),std::fstream::out);
            fs[f-1] << "# Saved by rita: Solution of ODE, equation: 1" << endl;
            fs[f-1] << 0.;
            for (int i=0; i<_ode_eq->size; ++i)
               fs[f-1] << "  " << _ode_eq->y[i];
            fs[f-1] << endl;
         }
         if ((*_isave)[f]) {
            ffs[f-1].open((*_save_file)[f].c_str());
            ffs[f-1] << "# Saved by rita: Solution of ODE, equation: 1" << endl;
            ffs[f-1] << 0.;
            for (int i=0; i<_ode_eq->size; ++i)
               ffs[f-1] << "  " << _ode_eq->y[i];
            ffs[f-1] << endl;
            if (_phase) {
               pfs[f-1].open((*_phase_file)[f].c_str());
               pfs[f-1] << "# Saved by rita: Phase portrait of ODE, equation: 1" << endl;
            }
         }*/
         if (_ode_eq->size==1)
            ode.setInitial(_ode_eq->y[0]);
         else
            ode.setInitial(_ode_eq->y);
         string fh = _data->vect_hist[_ode_eq->fn];
         if (fh!="%$§&")
            _data->theHVector[_data->checkName(fh,DataType::HVECTOR)]->set(_ode_eq->y,theTime);
      }
      for (int e=1; e<=_nb_pde && _rs; ++e) {
         _pde_eq = _data->thePDE[e];
         for (int i=0; i<_pde_eq->nb_vectors; ++i) {
            int f = _pde_eq->fd[i].vect;
            string fh = _data->vect_hist[_pde_eq->fd[i].fn];
            _data->theVector[f]->setTime(theTime);
            _data->theVector[f]->setName(_data->Vector[f]);
            if (fh!="%$§&")
               _data->theHVector[_data->checkName(fh,DataType::HVECTOR)]->set(*(_data->theVector[f]),theTime);
         }
/*         for (int i=0; i<_pde_eq->nb_vectors; ++i) {
            int f = _pde_eq->fd[i].vect;
            ff[f-1].open(fn[f-1],OFELI::IOField::OUT);
         }*/
      }
//   }
/*   else {
      for (int e=0; e<_nb_eq; ++e) {
         int f = _ode_eq[e]->field;
         _data->theField[f]->resize(_ode_eq[e]->size);
         *_data->theField[f] = _ode_eq[e]->y;
         if ((*_isave)[e]) {
            ffs[f].open((*_save_file)[e].c_str(),std::fstream::out);
            ffs[f] << "  " << _ode_eq[e]->y;
            if ((*_isave)[f]) {
               ffs[f].open((*_save_file)[f].c_str(),std::fstream::out);
               ffs[f] << "# Saved by rita: Solution of ODE, equation: " << e << endl;
               ffs[f] << 0.;
               for (int i=0; i<_rita->_ode[e].size; ++i)
                  ffs[f] << "  " << _ode_eq[e]->y[i];
               ffs[f] << endl;
               if (_phase) {
                  pfs[f].open((*_phase_file)[f].c_str());
                  pfs[f] << "# Saved by rita: Phase portrait of ODE, equation: " << e << endl;
               }
            }
         }
         if ((*_eq_type)[e]==data::eqType::AE) {
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
         else if ((*_eq_type)[e]==data::eqType::PDE && _rs) {
            for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
               int f = _pde_eq[e]->fd[i].field;
               ff[_pde_eq[e]->fd[i].field].open(fn[f],OFELI::IOField::OUT);
            }
         }
      }
   }*/
   theStep = 1;
   for (int e=1; e<=_nb_ae; ++e) {
      _ae_eq = _data->theAE[e];
      for (int i=0; i<_ae_eq->size; ++i)
         nlas.setf(_ae_eq->theFct[i]);
   }
   for (int e=1; e<=_nb_ode; ++e) {
      _ode_eq = _data->theODE[e];
      for (int i=0; i<_ode_eq->size; ++i)
         ode.setF(_ode_eq->theFct[i]);
   }
   for (int e=1; e<=_nb_pde; ++e) {
      _pde_eq = _data->thePDE[e];
      ts.setLinearSolver(_pde_eq->ls,_pde_eq->prec);
   }

// Loop on time steps
   try {
      TimeLoop {

         if (_rita->_verb)
            cout << "Performing time step " << theStep <<", Time = " << theTime << endl;

         for (int e=1; e<=_nb_ae; ++e) {
            cout << "No algebraic equation solver implemented." << endl;
         }

         for (int e=1; e<=_nb_ode; ++e) {
            _ode_eq = _data->theODE[e];
            int f = _ode_eq->vect;
            ode.runOneTimeStep();
            if (_ode_eq->size==1)
               _ode_eq->y[0] = ode.get();
            *_data->theVector[f] = _ode_eq->y;
            string fh = _data->vect_hist[_ode_eq->fn];
            if (fh!="%$§&")
               _data->theHVector[_data->checkName(fh,DataType::HVECTOR)]->set(_ode_eq->y,theTime);
            if (_ode_eq->phase!="") {
               _ode_eq->ph.setSize(_ode_eq->size);
               ode.getTimeDerivative(_ode_eq->ph);
               string fh = _data->vect_hist[_ode_eq->phase];
               if (fh!="%$§&")
                  _data->theHVector[_data->checkName(fh,DataType::HVECTOR)]->set(_ode_eq->ph,theTime);
            }
         }

         for (int e=1; e<=_nb_pde; ++e) {
            _pde_eq = _data->thePDE[e];

            if (_pde_eq->set_bf) {
               _pde_eq->bf.setTime(theTime);
               if (_pde_eq->bf.withRegex(1))
                  _pde_eq->bf.set(_pde_eq->bf_data.exp);
               ts.setRHS(_pde_eq->bf);
               _pde_eq->theEquation->setInput(BODY_FORCE,_pde_eq->bf);
            }

            if (_pde_eq->set_bc) {
               _pde_eq->bc.setTime(theTime);
               if (_pde_eq->bc.withRegex(1)) {
                  for (auto const& v: _pde_eq->bc_data.cexp)
                     _pde_eq->setNodeBC(v.first,v.second,theTime,_pde_eq->bc);
               }
               ts.setBC(_pde_eq->bc);
            }

            if (_pde_eq->set_sf) {
               _pde_eq->sf.setTime(theTime);
               if (_pde_eq->sf.withRegex(1)) {
                  for (auto const& v: _pde_eq->sf_data.cexp)
                     _pde_eq->setNodeBC(v.first,v.second,theTime,_pde_eq->sf);
               }
            //            ts.setSF(_data->sf[i]);
            }

            ts.runOneTimeStep();

            for (int i=0; i<_pde_eq->nb_vectors; ++i) {
               int f = _pde_eq->fd[i].vect;
               string fh = _data->vect_hist[_pde_eq->fd[i].fn];
               _data->theVector[f]->setTime(theTime);
               _data->theVector[f]->setName(_data->Vector[f]);
               if (fh!="%$§&")
                  _data->theHVector[_data->checkName(fh,DataType::HVECTOR)]->set(*(_data->theVector[f]),theTime);
            }
         }
      }
   } CATCH

   return 0;
}

} /* namespace RITA */
