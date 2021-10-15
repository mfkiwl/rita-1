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

                       Implementation of function runPDE

  ==============================================================================*/

#include "rita.h"
#include "data.h"
#include "equa.h"
#include "cmd.h"
#include "configure.h"

using std::cout;

namespace RITA {

int rita::runPDE()
{
   _pde = new equa(this);
   string pde_name, fn="";
   int ret=0;
   if (_cmd->setNbArg(1,"Give PDE name.")) {
      msg("pde>","Missing pde name.","",1);
      return 1;
   }
   static const vector<string> pn {"laplace","heat","wave","transport","linear-elasticity",
                                   "truss","beam","incompressible-navier-stokes","clear"};
   ret = _cmd->get(pn,pde_name);
   if (ret<0) {
      msg("pde>","Unknown pde "+pde_name,
          "Unknown pde name. Available pde's are:\n"
          "laplace, heat, wave, transport, linear-elasticity");
      _ret = 1;
      return 1;
   }
   _pde->set(pde_name);
   *ofh << "pde " << pde_name << endl;

   int nb_args = 0, nb=0, nb_fields=0;
   bool field_ok = false;
   vector<string> field_name;
   string str = "", str1 = "", ff="";
   _pde->set(_cmd);
   _pde->log.field = true;
   if (_analysis_type==NONE)
      _analysis_type = STEADY_STATE;
   if (_data->nb_meshes==0) {
      msg("pde>","No mesh created");
      _pde->log.mesh = true;
      return 1;
   }
   else if (_data->theMesh[_data->iMesh]->getNbNodes()==0) {
      msg("pde>","Empty mesh");
      _pde->log.mesh = true;
      return 1;
   }
   _pde->ls = OFELI::CG_SOLVER;
   _pde->prec = OFELI::DILU_PREC;
   string spd = "feP1";
   const static vector<string> kw {"field","coef","axi","in$it","bc","bf","source","sf",
                                   "traction","space","ls","nls","clear"};
   while (1) {
      if ((nb_args=_cmd->readline("rita>pde> "))<0)
         continue;
      int key = _cmd->getKW(kw,_gkw);
      switch (key) {

         case 100:
         case 101:
            cout << "\nAvailable Commands:\n";
            cout << "field:  Field name of an unknown of the equation\n";
            cout << "coef:   PDE coefficients\n";
            cout << "axi:    Choose axisymmetric geometry\n";
            cout << "init:   Set initial condition or guess for pde\n";
            cout << "bc:     Set boundary conditions\n";
            cout << "source: Set sources or body forces\n";
            cout << "sf:     Set side (boundary) forces\n";
            cout << "space:  Space discretization method\n";
            cout << "ls:     Set linear system solver\n";
            cout << "nls:    Set nonlinear system iteration procedure\n";
            cout << "clear:  Remove pde from model\n" << endl;
            break;

         case 102:
            _ret = _configure->run();
            break;

         case   0:
            if (_cmd->setNbArg(1,"Give name of an associated field.")) {
               msg("pde>field>","Missing name of an associated field.","",1);
               break;
            }
            if (!_cmd->get(str)) {
               field_name.push_back(str);
               nb_fields++;
               field_ok = true;
               *ofh << "  field " << str << endl;
            }
            break;

         case   1:
            _pde->set_coef = true;
            _ret = _pde->setCoef();
            break;

         case   2:
            _pde->axi = true;
            break;

         case   3:
            _pde->getIn();
            break;

         case   4:
            _pde->getBC();
           break;

         case   5:
         case   6:
            _pde->getBF();
            break;

         case   7:
         case   8:
            _pde->getSF();
            break;

         case   9:
            if (_cmd->setNbArg(1,"Give space discretization method.")) {
               msg("pde>space>","Missing space discretization method.","",1);
               cout << "Available Commands\n";
               cout << "fd:   Finite Differences\n";
               cout << "feP1: P1 finite elements\n";
               cout << "feP2: P2 finite elements\n";
               cout << "feQ1: Q1 finite elements\n";
               cout << "fv:   Finite volumes\n";
               cout << "dg:   Discontinuous Galerkin" << endl;
               break;
            }
            _cmd->get(spd);
            if (_pde->pde_sdm.find(spd)==_pde->pde_sdm.end()) {
               msg("pde>spd>","Unknown or unimplemented space discretization method: "+spd,"",1);
               break;
            }
            break;

         case  10:
            if (_cmd->setNbArg(1,"Linear solver and optional preconditioner to be supplied.",1)) {
               msg("pde>ls>","Missing linear solver data.","",1);
               break;
            }
            nb = _cmd->getNbArgs();
            if (nb==0)
               msg("pde>ls>","Missing linear solver data.");
            _ret = _cmd->get(str);
            str1 = "ident";
            if (nb>1)
               _ret += _cmd->get(str1);
            if (!_ret) {
               *ofh << "  ls " << str << " " << str1 << endl;
               if (!set_ls(str,str1)) {
                  _pde->ls = Ls[str];
                  _pde->prec = Prec[str1];
               }
            }
            else
               _pde->log.ls = true;
            break;

         case  11:
            if (_cmd->setNbArg(1,"Nonlinear solver to be supplied.",1)) {
               msg("pde>nls>","Missing nonlinear solver data.","",1);
               break;
            }
            _ret = _cmd->get(str);
            if (!_ret) {
               *ofh << "  nls " << str << endl;
               _ret = set_nls(str);
               if (!_ret)
                  _pde->nls = str;
               else
                  _pde->log.nl = true;
            }
            break;

         case  12:
            cout << "PDE removed from model." << endl;
            *ofh << "  clear" << endl;
            _ret = 10;
            return _ret;

         case 104:
         case 105:
            setParam();
            break;

         case 106:
            if (_cmd->setNbArg(1,"Data name to be given.",1)) {
               msg("print>","Missing data name.","",1);
               break;
            }
            if (!_cmd->get(fn))
               _data->print(fn);
            break;

         case 107:
         case 108:
           _cmd->setNbArg(0);
            if (_ret) {
               msg("pde>end>","No PDE data created.");
               *ofh << "end" << endl;
            }
            if (!field_ok) {
               msg("pde>end>","No field(s) defined for PDE.");
               break;
            }
            _pde->setSpD(spd);
            _pde->setEq();
            *ofh << "  space " << spd << endl;
            if (nb_fields>_pde->nb_fields) {
               msg("pde>end>","Too many fields for defined PDE.");
               break;
            }
            if (nb_fields<_pde->nb_fields) {
               msg("pde>end>","Not enough fields for defined PDE.");
               break;
            }
            ff = spd.substr(0,2);
            _data->addPDE(_pde);
            for (int i=0; i<nb_fields; ++i) {
               if (ff=="fd")
                  _data->addGridField(field_name[i],_pde->fd[i].nb_dof);
               else if (ff=="fe" || ff=="fv" || ff=="dg")
                  _data->addMeshField(field_name[i],data::DataSize::NODES,_pde->fd[i].nb_dof);
               _data->FieldEquation.push_back(_data->getNbEq());
               _pde->fd[i].fn = field_name[i];
               _pde->fd[i].field = _data->iField;
            }
            _pde->log.field = false;
            _pde->b.setSize(_data->theMesh[_data->iMesh]->getNbEq());
            if (_pde->set_in)
               _pde->setIn();
            if (_pde->set_bc)
               _pde->setBC();
            if (_pde->set_sf)
               _pde->setSF();
            if (_pde->set_bf)
               _pde->setBF();
            if (_verb) {
               cout << "Summary of PDE settings:" << endl;
               cout << "   PDE name: " << _pde->eq << endl;
               cout << "   PDE unknown field(s): ";
               for (int i=0; i<_pde->nb_fields-1; ++i)
                  cout << _pde->fd[i].fn << ", ";
               cout << _pde->fd[nb_fields-1].fn << endl;
               cout << "   PDE space discretization: " << spd << endl;
               cout << "   PDE linear solver: " << rLs[_pde->ls] << endl;
            }
            _ret = 0;
            *ofh << "  end" << endl;
            _pde->log.pde = false;
            _pde->set();
            return 0;

         default:
            msg("pde>","Unknown Command "+_cmd->token(),
                "Available commands: field, coef, axi, in, bc, bf, sf, space, ls, nls, clear");
            break;
      }
   }
   _ret = 0;
   return 0;
}

} /* namespace RITA */
