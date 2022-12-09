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

                       Implementation of function runPDE

  ==============================================================================*/

#include "rita.h"
#include "data.h"
#include "calc.h"
#include "equa.h"
#include "cmd.h"
#include "configure.h"

using std::cout;

namespace RITA {

int rita::runPDE()
{
   _pde = new equa(this);
   string pde_id, name="", fn="", file="", lin_solv="";
   int ret=0, every=0;
   if (_cmd->setNbArg(1,"Give PDE name.")) {
      msg("pde>","Missing pde name.","",1);
      return 1;
   }
   static const vector<string> pn {"laplace","heat","wave","transport","linear-elasticity",
                                   "truss","beam","incompressible-navier-stokes","clear"};
   ret = _cmd->get(pn,pde_id);
   if (ret<0) {
      msg("pde>","Unknown pde "+pde_id,
          "Unknown pde name. Available pde's are:\n"
          "laplace, heat, wave, transport, linear-elasticity");
      _ret = 1;
      return 1;
   }
   _pde->set(pde_id);
   *ofh << "pde " << pde_id << endl;

   int nb_args = 0, nb=0, nb_vectors=0;
   bool vector_ok = false;
   vector<string> vector_name;
   string str = "", str1 = "", ff="";
   _pde->set(_cmd);
   _pde->log.vect = true;
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
   _pde->lsolv = "cg";
   _pde->lprec = "dilu";
   _pde->prec = OFELI::DILU_PREC;
   string spd = "feP1";
   const static vector<string> kw {"var$iable","vect$or","coef","axi","in$it","bc","bf","source","sf",
                                   "traction","space","ls","nls","clear","name","save-every","save-file"};
   while (1) {
      if ((nb_args=_cmd->readline(sPrompt+"pde> "))<0)
         continue;
      int key = _cmd->getKW(kw,_gkw);
      if (key>=200) {
         _data->setDataExt(key);
         continue;
      }
      switch (key) {

         case   0:
         case   1:
            if (_cmd->setNbArg(1,"Give name of an associated variable/vector.")) {
               msg("pde>vector>","Missing name of an associated vector.","",1);
               break;
            }
            if (!_cmd->get(str)) {
               vector_name.push_back(str);
               nb_vectors++;
               vector_ok = true;
               *ofh << "  vector " << str << endl;
            }
            break;

         case   2:
            _pde->set_coef = true;
            _ret = _pde->setCoef();
            break;

         case   3:
            _pde->axi = true;
            break;

         case   4:
            _pde->getIn();
            break;

         case   5:
            _pde->getBC();
           break;

         case   6:
         case   7:
            _pde->getBF();
            break;

         case   8:
         case   9:
            _pde->getSF();
            break;

         case  10:
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

         case  11:
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
                  _pde->lsolv = str;
                  _pde->lprec = str1;
                  _pde->ls = Ls[str];
                  _pde->prec = Prec[str1];
               }
            }
            else
               _pde->log.ls = true;
            break;

         case  12:
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

         case  13:
            cout << "PDE removed from model." << endl;
            *ofh << "  clear" << endl;
            _ret = 10;
            return _ret;

         case  14:
            if (!_cmd->get(name))
               *ofh << "  name " << name << endl;
            break;

         case  15:
            if (!_cmd->get(every)) {
               *ofh << "  save-every " << every << endl;
            }
            break;

         case  16:
            if (!_cmd->get(file)) {
               *ofh << "  save-file " << file << endl;
            }
            break;

         case 100:
         case 101:
            cout << "\nAvailable Commands:\n";
            cout << "vector:     Vector name of an unknown of the equation\n";
            cout << "coef:       PDE coefficients\n";
            cout << "axi:        Choose axisymmetric geometry\n";
            cout << "init:       Set initial condition or guess for pde\n";
            cout << "bc:         Set boundary conditions\n";
            cout << "source:     Set sources or body forces\n";
            cout << "sf:         Set side (boundary) forces\n";
            cout << "space:      Space discretization method\n";
            cout << "ls:         Set linear system solver\n";
            cout << "nls:        Set nonlinear system iteration procedure\n";
            cout << "name:       Set name of the defined PDE\n";
            cout << "save-every: Save PDE solution every n time steps\n";
            cout << "save-file:  Save file (if available) in given file\n";
            cout << "clear:  Remove pde from model\n" << endl;
            break;

         case 102:
            getLicense();
            break;

         case 103:
            _ret = _configure->run();
            break;

         case 104:
         case 105:
           _cmd->setNbArg(0);
            if (_ret) {
               msg("pde>end>","No PDE data created.");
               *ofh << "end" << endl;
            }
            if (!vector_ok) {
               msg("pde>end>","No vector(s) defined for PDE.");
               break;
            }
            _pde->setSpD(spd);
            _pde->setEq();
            *ofh << "  space " << spd << endl;
            if (nb_vectors>_pde->nb_vectors) {
               msg("pde>end>","Too many vectors for defined PDE.");
               break;
            }
            if (nb_vectors<_pde->nb_vectors) {
               msg("pde>end>","Not enough vectors for defined PDE.");
               break;
            }
            ff = spd.substr(0,2);
            if (name=="")
               name = "PDE-" + to_string(_data->theMesh.size()-1);
            _data->addPDE(_pde,name);
            for (int i=0; i<nb_vectors; ++i) {
               if (ff=="fd")
                  _data->addGridVector(vector_name[i],_pde->fd[i].nb_dof);
               else if (ff=="fe" || ff=="fv" || ff=="dg")
                  _data->addMeshVector(vector_name[i],data::DataSize::NODES,_pde->fd[i].nb_dof);
               _data->VectorEquation.push_back(_data->getNbEq());
               _pde->fd[i].fn = vector_name[i];
               _pde->fd[i].vect = _data->iVector;
            }
            _pde->log.vect = false;
            _pde->b.setSize(_data->theMesh[_data->iMesh]->getNbEq());
            if (_pde->set_in)
               _pde->setIn();
            if (_pde->set_bc)
               _pde->setBC();
            if (_pde->set_sf)
               _pde->setSF();
            if (_pde->set_bf)
               _pde->setBF();
            _pde->every = every;
            _pde->file = file;
            _pde->name = name;
            _ret = 0;
            *ofh << "  end" << endl;
            _pde->log.pde = false;
            _pde->set();
            return 0;

         default:
            _ret = _calc->run();
            break;
      }
   }
   _ret = 0;
   return 0;
}

} /* namespace RITA */
