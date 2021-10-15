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

                      Implementation of class 'stationary'

  ==============================================================================*/

#include "stationary.h"
#include "io/IOField.h"
#include "io/saveField.h"
#include "equa.h"
#include <iostream>

using std::map;

namespace RITA {

stationary::stationary(rita *r)
           : _rita(r), _rs(1)
{
   _data = _rita->_data;
   _nb_fields = _data->nb_fields;
}


stationary::~stationary()
{
}


void stationary::setSave(vector<int>&    isave,
                         vector<int>&    fformat,
                         vector<string>& save_file)
{
   _isave = &isave;
   _fformat = &fformat;
   _save_file = &save_file;
}


int stationary::run()
{
   int ret = 0;
   vector<ofstream> fs(_nb_fields+1), ffs(_nb_fields+1), pfs(_nb_fields+1);
   vector<OFELI::IOField> ff(_nb_fields+1);
   vector<string> fn(_nb_fields+1);

   try {
      for (int e=1; e<=_data->nb_ae; ++e) {
         _ae_eq = _data->theAE[e];
         int f = _ae_eq->field;
         fn[f] = "rita-" + to_string(10*e+f) + ".sol";
      }
      for (int e=1; e<=_data->nb_pde; ++e) {
         _pde_eq = _data->thePDE[e];
         for (int i=0; i<_pde_eq->nb_fields; ++i) {
            int f = _pde_eq->fd[i].field;
            fn[f] = "rita-" + to_string(10*e+f) + ".sol";
         }
      }

//      if (_nb_eq==1) {
         for (int e=1; e<=_data->nb_ae && _rs; ++e) {
            _ae_eq = _data->theAE[e];
            int f = _ae_eq->field;
            _data->u[f]->resize(_ae_eq->size);
            *_data->u[f] = _ae_eq->y;
            fs[f].open(fn[f].c_str(),std::fstream::out);
            fs[f] << "# Saved by rita: Solution of Algebraic Equation, equation: 1" << endl;
            fs[f] << 0.;
            for (int i=0; i<_ae_eq->size; ++i)
               fs[f] << "  " << _ae_eq->y[i];
            fs[f] << endl;
            if ((*_isave)[e]) {
               ffs[f].open((*_save_file)[f].c_str());
               ffs[f] << "# Saved by rita: Solution of Algebraic Equation, equation: 1" << endl;
               ffs[f] << 0.;
               for (int i=0; i<_ae_eq->size; ++i)
                  ffs[f] << "  " << _ae_eq->y[i];
               ffs[f] << endl;
            }
         }

         for (int e=1; e<=_data->nb_pde && _rs; ++e) {
            _pde_eq = _data->thePDE[e];
            for (int i=0; i<_pde_eq->nb_fields; ++i) {
               int f = _pde_eq->fd[i].field;
               ff[f].open(fn[f],OFELI::IOField::OUT);
            }
         }
//      }
/*      else {
         for (int e=0; e<_nb_eq; ++e) {
            if ((*_eq_type)[e]==ALGEBRAIC_EQ) {
               int f = _alg_eq[e]->field;
               _data->u[f]->resize(_alg_eq[e]->size);
               *_data->u[f] = _alg_eq[e]->y;
               if ((*_isave)[f] && _rs) {
                  ffs[f].open((*_save_file)[e].c_str(),std::fstream::out);
                  ffs[f] << "  " << _alg_eq[e]->y;
                  if ((*_isave)[f]) {
                     ffs[f].open((*_save_file)[f].c_str(),std::fstream::out);
                     ffs[f] << "# Saved by rita: Solution of ODE, equation: " << e << endl;
                     ffs[f] << 0.;
                     for (int i=0; i<_alg_eq[e]->size; ++i)
                        ffs[f] << "  " << _alg_eq[e]->y[i];
                  }
               }
            }
            else if ((*_eq_type)[e]==PDE_EQ && _rs) {
               for (int i=0; i<_pde_eq[e]->nb_fields; ++i) {
                  int f = _pde_eq[e]->fd[i].field;
                  ff[_pde_eq[e]->fd[i].field].open(fn[f],OFELI::IOField::OUT);
               }
            }
         }
      }*/

      for (int e=1; e<=_data->nb_ae; ++e) {
         _ae_eq = _data->theAE[e];
         NLASSolver nls(_ae_eq->nls,_ae_eq->size);
         if (_ae_eq->size==1)
            nls.setInitial(_ae_eq->y[0]);
         else
            nls.setInitial(_ae_eq->y);
         for (int i=0; i<_ae_eq->size; ++i)
            nls.setf(_ae_eq->theFct[i]);
         nls.run();
         *_data->u[_ae_eq->field] = _ae_eq->y;
      }

      for (int e=1; e<=_data->nb_pde; ++e) {
         _pde_eq = _data->thePDE[e];
         _pde_eq->theEquation->setInput(SOLUTION,*_data->u[_pde_eq->fd[0].field]);
         _pde_eq->theEquation->setSolver(_pde_eq->ls,_pde_eq->prec);
         if (_pde_eq->set_bc) {
            if (_pde_eq->bc.withRegex(1)) {
               for (auto const& v: _pde_eq->bc_data.cexp)
                  _pde_eq->bc.setNodeBC(v.first,v.second);
            }
            _pde_eq->theEquation->setInput(BOUNDARY_CONDITION,_pde_eq->bc);
         }
         if (_pde_eq->set_bf) {
            if (_pde_eq->bf.withRegex(e))
               _pde_eq->bf.set(_pde_eq->bf_data.exp);
            _pde_eq->theEquation->setInput(BODY_FORCE,_pde_eq->bf);
         }
         if (_pde_eq->set_sf) {
            if (_pde_eq->sf.withRegex(e)) {
               for (auto const& v: _pde_eq->sf_data.cexp)
                  _pde_eq->sf.setSideBC(v.first,v.second);
            }
            _pde_eq->theEquation->setInput(BOUNDARY_FORCE,_pde_eq->sf);
         }

         ret = _pde_eq->theEquation->run();

         for (int i=0; i<_pde_eq->nb_fields; ++i) {
            int f = _pde_eq->fd[i].field;
	    //            _data->u[f]->setName(_data->Field[f]);
            if (_rs) {
               ff[f].open(fn[f],OFELI::IOField::OUT);
               ff[f].put(*_data->u[f]);
               ff[f].close();
               if ((*_isave)[f])
                  OFELI::saveFormat(*_data->theMesh[_data->iMesh],fn[f],(*_save_file)[f],
                                    (*_fformat)[f],(*_isave)[f]);
            }
         }
      }
   } CATCH
   return ret;
}

} /* namespace RITA */
