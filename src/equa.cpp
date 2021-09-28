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

                       Implementation of class 'equa'

  ==============================================================================*/

#include "equa.h"
#include "cmd.h"
#include "rita.h"

namespace RITA {

equa::equa(rita *r)
     : eq("laplace"), nls(""), axi(false),
       ls(CG_SOLVER), prec(DILU_PREC), _nb_fields(0), _data(r->_data), _theMesh(nullptr),
       _theGrid(nullptr)
{
   _rita = r;
   _verb = 0;
   ieq = pde_map["laplace"];
   _rho_set = _Cp_set = _kappa_set = _mu_set = _sigma_set = _Mu_set = false;
   _epsilon_set = _omega_set = _beta_set = _v_set = _young_set = _poisson_set = false;
   set_u = set_bc = set_bf = set_sf = set_in = set_coef = false;
}


equa::~equa()
{
   if (theEquation!=nullptr)
      delete theEquation;
}


void equa::set(string e)
{
   ieq = pde_map[e];
   eq = e;
}


int equa::setSpD(string spd)
{
   Sdm = pde_sdm[spd];
   if (Sdm==FE_P1 || Sdm==FE_P2 || Sdm==FE_Q1 || Sdm==FV || Sdm==DG) {
      if (_data->nb_meshes==0) {
         _rita->msg("pde>","No defined mesh.");
         return 1;
      }
      _theMesh = _data->theMesh[_data->iMesh];
      _dim = _theMesh->getDim();
   }
   else if (Sdm==FD) {
      if (_data->nb_grids==0) {
         _rita->msg("pde>","No defined grid.");
         return 1;
      }
      _theGrid = _data->theGrid[_data->iGrid];
      _dim = _theGrid->getDim();
   }
   return 0;
}


void equa::setFields()
{
   nb_fields = 1;
   fd[0].fn = "u";
   fd[0].nb_dof = 1;
   fd[0].ds = data::DataSize::NODES;

   switch (ieq) {

      case LAPLACE:
         break;

      case HEAT:
         break;

      case WAVE:
         break;

      case TRANSPORT:
         break;

      case LINEAR_ELASTICITY:
         fd[0].nb_dof = _dim;
         break;

      case TRUSS:
         fd[0].nb_dof = 2;
         break;

      case BEAM:
         if (_dim==2)
            fd[0].nb_dof = 2;
         else if (_dim==3)
            fd[0].nb_dof = 6;
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         nb_fields = 2;
         fd[0].fn = "v", fd[1].fn = "p";
         fd[0].nb_dof = _dim, fd[1].nb_dof = 1;
         break;
         
      case COMPRESSIBLE_EULER:
         nb_fields = 4;
         fd[0].fn = "u", fd[1].fn = "rho", fd[2].fn = "p", fd[3].fn = "e";
         fd[0].nb_dof = _dim, fd[1].nb_dof = fd[2].nb_dof = fd[3].nb_dof = 1;
         break;
         
      case COMPRESSIBLE_NAVIER_STOKES:
         nb_fields = 4;
         fd[0].fn = "u", fd[1].fn = "rho", fd[2].fn = "p", fd[3].fn = "e";
         fd[0].nb_dof = _dim, fd[1].nb_dof = fd[2].nb_dof = fd[3].nb_dof = 1;
         break;

      case INCOMPRESSIBLE_POROUS_1PHASE:
         fd[0].fn = "p";
         break;

      case EDDY_CURRENTS:
         fd[0].fn = "A";
         if (_dim==3)
            fd[0].nb_dof = _dim;
         break;

      case MAXWELL:
         fd[0].nb_dof = _dim;
         break;

      case HELMHOLTZ:
         fd[0].fn = "p";
         break;
   }
}


void equa::setSize(Vect<double>& v, data::DataSize s)
{
   _nb_dof = 1;
   if (_theMesh!=nullptr) {
      if (_theMesh->getDOFSupport()==NODE_DOF)
         _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbNodes();
      else if (_theMesh->getDOFSupport()==SIDE_DOF)
         _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbSides();
      else if (_theMesh->getDOFSupport()==ELEMENT_DOF)
         _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbElements();
      else if (_theMesh->getDOFSupport()==EDGE_DOF)
         _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbEdges();

      if (s==data::DataSize::NODES) {
         if (_theMesh->getNbNodes()==0) {
            _rita->msg("pde>","Mesh has no nodes");
            _ret = 1;
            return;
         }
         v.setMesh(*_theMesh,NODE_DOF,_nb_dof);
      }
      else if (s==data::DataSize::ELEMENTS) {
         if (_theMesh->getNbElements()==0) {
            _rita->msg("pde>:","Mesh has no elements");
            _ret = 1;
            return;
         }
         v.setMesh(*_theMesh,ELEMENT_DOF,_nb_dof);
      }
      else if (s==data::DataSize::SIDES) {
         if (_theMesh->getNbSides()==0) {
            _rita->msg("pde>:","Mesh has no sides");
            _ret = 1;
            return;
         }
         v.setMesh(*_theMesh,SIDE_DOF,_nb_dof);
      }
      else if (s==data::DataSize::EDGES) {
         if (_theMesh->getNbEdges()==0) {
            _rita->msg("pde>","Mesh has no edges");
            _ret = 1;
            return;
         }
         v.setMesh(*_theMesh,EDGE_DOF,_nb_dof);
      }
   }
   else if (_theGrid!=nullptr) {
      _nb_dof = _theGrid->getNbDOF()/_theGrid->getNbNodes();
      v.setGrid(*_theGrid);
   }
   else {
      _rita->msg("pde>","No mesh or grid data available.");
      _ret = 1;
      return;
   }
}


int equa::getIn()
{
   bool val_ok=false, file_ok=false, save_ok=false;
   _ret = 0;
   static const vector<string> kw {"val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("pde>initial>","No argument given");
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            in_data.exp = _cmd->string_token();
            val_ok = true;
            break;

         case 1:
            in_data.in_file = _cmd->string_token();
            file_ok = true;
            break;

         case 2:
            in_data.out_file = _cmd->string_token();
            save_ok = true;
            break;

         default:
            _rita->msg("pde>initial>","Unknown argument: "+kw[n]);
            return _ret;
      }
   }
   if (!val_ok) {
      _rita->msg("pde>initial>","No value or expression given for initial condition.");
      return 1;
   }
   *_rita->ofh << "  in  value=" << in_data.exp;
   if (file_ok)
      *_rita->ofh << "  file=" << in_data.in_file;
   if (save_ok) {
      if (_verb)
         cout << "Initial condition saved in file: " << in_data.out_file << endl;
      *_rita->ofh << "  save=" << in_data.out_file;
   }
   *_rita->ofh << endl;
   in_data.size++;
   set_in = true;
   return 0;
}


int equa::setIn()
{
   if (in_data.size==0) {
      _rita->msg("pde>initial>","No defined initial data");
      return 1;
   }
   if (_theMesh!=nullptr) {
      _nb_dof = _theMesh->getNbDOF()/_theMesh->getNbNodes();
      u.setMesh(*_theMesh,NODE_DOF,_nb_dof);
   }
   theSolution[0] = &u;
   u.setRegex(1);
   if (in_data.in_file.size()) {
      OFELI::IOField ffi(in_data.in_file,OFELI::IOField::IN);
      ffi.get(u);
   }
   if (in_data.out_file.size()) {
      OFELI::IOField ffo(in_data.out_file,OFELI::IOField::OUT);
      ffo.put(u);
   }
   set_u = true;
   return 0;
}


int equa::getBC()
{
   _ret = 0;
   string val="";
   bool code_ok=false, val_ok=false, file_ok=false, save_ok=false;
   int code=0;
   const static vector<string> kw {"code","val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("pde>bc>","No arguments");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            code = _cmd->int_token();
            code_ok = true;
            break;

         case 1:
            val = _cmd->string_token();
            val_ok = true;
            break;

         case 2:
            bc_data.in_file = _cmd->string_token();
            file_ok = true;
            break;

         case 3:
            bc_data.out_file = _cmd->string_token();
            save_ok = true;
            break;

         default:
            _rita->msg("pde>bc>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (!code_ok) {
      _rita->msg("pde>bc>","No code given for boundary condition.");
      return 1;
   }
   if (!val_ok) {
      _rita->msg("pde>bc>","No value or expression given for boundary condition.");
      return 1;
   }
   if (code<=0) {
      _rita->msg("pde>bc>","Illegal boundary condition for a nonpositive code "+to_string(code));
      return 1;
   }
   bc_data.cexp[code] = val;
   bc_data.size++;
   if (_verb)
      cout << "Nodes with code " << code << " have prescribed value by the expression: " << val << endl;
   *_rita->ofh << "  bc  code=" << code << "  value=" << val;
   if (file_ok)
      *_rita->ofh << "  file=" << bc_data.in_file;
   if (save_ok) {
      *_rita->ofh << "  save=" << bc_data.out_file;
      if (_verb)
         cout << "Boundary condition saved in file: " << bc_data.out_file << endl;
   }
   *_rita->ofh << endl;
   set_bc = true;
   return 0;
}


int equa::setBC()
{
   if (bc_data.size==0) {
      _rita->msg("pde>bc>","No defined boundary condition");
      return 1;
   }
   _theMesh = _data->theMesh[_data->iMesh];
   setSize(bc,data::DataSize::NODES);
   bc.setRegex(1);
   if (bc_data.in_file.size()) {
      OFELI::IOField ffi(bc_data.in_file,OFELI::IOField::IN);
      ffi.get(bc);
   }
   if (bc_data.out_file.size()) {
      OFELI::IOField ffo(bc_data.out_file,OFELI::IOField::OUT);
      ffo.put(bc);
   }
   return 0;
}


int equa::getSF()
{
   int code=0;
   _ret = 0;
   string val="";
   bool code_ok=false, val_ok=false, file_ok=false, save_ok=false;
   const static vector<string> kw {"code","val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("pde>sf>","No arguments");
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            code = _cmd->int_token();
            code_ok = true;
            break;

         case 1:
            val = _cmd->string_token();
            val_ok = true;
            break;

         case 2:
            sf_data.in_file = _cmd->string_token();
            file_ok = true;
            break;

         case 3:
            sf_data.out_file = _cmd->string_token();
            save_ok = true;
            break;

         default:
            _rita->msg("pde>sf>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (!code_ok) {
      _rita->msg("pde>sf>","No code given.");
      return 1;
   }
   if (!val_ok) {
      _rita->msg("pde>sf>","No value or expression given for surface force.");
      return 1;
   }
   if (code<=0) {
      _rita->msg("pde>sf>","Illegal value of code: "+to_string(code));
      return 1;
   }
   sf_data.cexp[code] = val;
   sf_data.size++;
   *_rita->ofh << "  sf  code=" << code << "  value=" << val;
   if (file_ok)
      *_rita->ofh << "  file=" << sf_data.in_file;
   if (save_ok)
      *_rita->ofh << "  save=" << sf_data.out_file;
   *_rita->ofh << endl;
   set_sf = true;
   return 0;
}


int equa::setSF()
{
   if (sf_data.size==0) {
      _rita->msg("pde>sf>","No defined boundary force");
      return 1;
   }
   _theMesh = _data->theMesh[_data->iMesh];
   setSize(sf,data::DataSize::SIDES);
   sf.setRegex(1);
   if (sf_data.in_file.size()) {
      OFELI::IOField ffi(sf_data.in_file,OFELI::IOField::IN);
      ffi.get(sf);
   }
   if (sf_data.out_file.size()) {
      OFELI::IOField ffo(sf_data.out_file,OFELI::IOField::OUT);
      ffo.put(sf);
   }
   return 0;
}


int equa::getBF()
{
   _ret = 0;
   bool val_ok=false, file_ok=false, save_ok=false;
   static const vector<string> kw {"val$ue","file","save"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("pde>bf>","No arguments");
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            bf_data.exp = _cmd->string_token();
            val_ok = true;
            break;

         case 1:
            bf_data.in_file = _cmd->string_token();
            file_ok = true;
            break;

         case 2:
            bf_data.out_file = _cmd->string_token();
            save_ok = true;
            break;

         default:
            _rita->msg("pde>source>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (!val_ok)
      _rita->msg("pde>source>","No value or expression given for source.");
   *_rita->ofh << "  source  value=" << bf_data.exp;
   set_bf = true;
   if (file_ok)
      *_rita->ofh << "  file=" << bf_data.in_file;
   if (save_ok)
      *_rita->ofh << "  save=" << bf_data.out_file;
   *_rita->ofh << endl;
   bf_data.size++;
   set_bf = true;
   return 0;
}


int equa::setBF()
{
   if (bf_data.size==0) {
      _rita->msg("pde>bf>","No defined body force");
      return 1;
   }
   _theMesh = _data->theMesh[_data->iMesh];
   setSize(bf,data::DataSize::NODES);
   _ret = 0;
   bf.setRegex(1);
   if (bf_data.in_file.size()) {
      OFELI::IOField ffi(bf_data.in_file,OFELI::IOField::IN);
      ffi.get(bf);
   }
   if (bf_data.out_file.size()) {
      OFELI::IOField ffo(bf_data.out_file,OFELI::IOField::OUT);
      ffo.put(bf);
   }
   return 0;
}


void equa::setNodeBC(int code, string exp, double t, Vect<double>& v)
{
   const static vector<string> var {"x","y","z","t"};
   _theFct.set(exp,var);
   for (size_t n=1; n<=_theMesh->getNbNodes(); ++n) {
      Node *nd = (*_theMesh)[n];
      for (size_t i=1; i<=nd->getNbDOF(); ++i) {
         if (nd->getCode(i)==code)
            v(nd->n(),i) = _theFct(nd->getCoord());
      }
   }
}


int equa::setEq()
{
   int ret = 0;
   setFields();

   switch (ieq) {

//    Laplace equation
      case LAPLACE:
         switch (_dim) {
  
            case 1:
               if (Sdm!=FE_P1 && Sdm!=FE_P2) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.pde = true;
               }
               break;

            case 2:
               if (Sdm!=FE_P1 && !axi) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 3:
               if (Sdm!=FE_P1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;
         }
         break;

//    Heat equation
      case HEAT:
         switch (_dim) {

            case 1:
               if (Sdm!=FE_P1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.pde = true;
               }
               break;

            case 2:
               if (Sdm!=FE_P1 && Sdm!=FE_P2) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 3:
               if (Sdm!=FE_P1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;
         }
         break;

//    Wave equation
      case WAVE:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 2:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 3:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;
         }
         break;

//    Linear transport equation
      case TRANSPORT:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 2:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 3:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;
         }
         break;

      case LINEAR_ELASTICITY:
         switch (_dim) {

            case 1:
                if (Sdm==FE_P1) {
                 _rita->msg("pde>","Equation not implemented in rita.");
                  log.pde = true;
               }
               else {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;

            case 2:
               if (Sdm!=FE_P1 && Sdm!=FE_Q1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.pde = true;
               }
               break;

            case 3:
               if (Sdm!=FE_P1 && Sdm!=FE_Q1) {
                  _rita->msg("pde>","Approximation of equation not implemented in rita.");
                  log.spd = true;
               }
               break;
         }
         break;

      case TRUSS:
         switch (_dim) {

            case 1:
               _rita->msg("pde>","No implementation for 1-D available for bar equation.");
               log.pde = true;
               break;

            case 2:
               if (Sdm!=FE_P1) {
                  _rita->msg("pde>","Only 2-D P1 finite element is available for bar equation.");
                  log.pde = true;
               }
               break;

            case 3:
               _rita->msg("pde>","No implementation for 3-D available for bar equation.");
               log.pde = true;
               break;
         }
         break;
         
      case BEAM:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","No implementation for 1-D P1 finite element available for beam equation.");
                  log.pde = true;
               }
               break;

            case 2:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","No implementation for 2-D P1 finite element available for beam equation.");
                  log.pde = true;
               }
               else if (Sdm==FE_P2) {
                  _rita->msg("pde>","No implementation for 2-D P2 finite element available for beam equation.");
                   log.pde = true;
               }
               break;

            case 3:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","No implementation for 3-D P1 finite element available for beam equation.");
                  log.pde = true;
               }
               break;
         }
         break;

      case INCOMPRESSIBLE_NAVIER_STOKES:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","No implementation for 1-D available for incompressible Navier-Stokes equations.");
                  log.pde = true;
               }
               break;

            case 2:
               if (Sdm!=FE_P1) {
                  _rita->msg("pde>","Only P1 finite element is implemented is available for incompressible Navier-Stokes equations.");
                  log.pde = true;
               }
               break;

            case 3:
               if (Sdm==FE_P1)
                  log.pde = true;
               break;
         }
         break;
         
      case COMPRESSIBLE_EULER:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                   log.pde = true;
               }
               break;

            case 2:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               else if (Sdm==FE_P2) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               break;

            case 3:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               break;
         }
         break;

      case INCOMPRESSIBLE_POROUS_1PHASE:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               break;

            case 2:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               else if (Sdm==FE_P2) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               break;

            case 3:
               if (Sdm==FE_P1) {
                  _rita->msg("pde>","Space discretization method not implemented for this PDE.");
                  log.pde = true;
               }
               break;
         }
         break;

      default:
         _rita->msg("pde>","Equation not implemented in rita.");
         log.pde = true;
         break;

   }
   return ret;
}


int equa::setCoef()
{
   static const string H = "Command: coef [rho=x] [Cp=x] [kappa=x] [Mu=x] [sigma=x] [mu=x] [epsilon=x] [omega=x]\n"
                           "              [beta=x] [v=x] [young=x] [poisson=x]\n\n";
   const static vector<string> kw {"rho","density","Cp","specific-heat","kappa","thermal-conductivity",
                                   "Mu","magnetic-permeability","sigma","electric-conductivity","mu","viscosity","epsilon",
                                   "electric-permittivity","omega","angular-frequency","beta","thermal-dilatation","v",
                                   "velocity","young","poisson"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("pde>coef>","No argument");
      _ret = 1;
      return _ret;
   }
   if (nb_args<1) {
      _rita->msg("pde>coef>","No argument.");
      _ret = 1;
      return _ret;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case  0:
         case  1:
            _rho_exp = _cmd->string_token();
            _rho_set = true;
            break;

         case  2:
         case  3:
            _Cp_exp = _cmd->string_token();
            _Cp_set = true;
            break;

         case  4:
         case  5:
            _kappa_exp = _cmd->string_token();
            _kappa_set = true;
            break;

         case  6:
         case  7:
            _Mu_exp = _cmd->string_token();
            _Mu_set = true;
            break;

         case  8:
         case  9:
            _sigma_exp = _cmd->string_token();
            _sigma_set = true;
            break;

         case 10:
         case 11:
            _mu_exp = _cmd->string_token();
            _mu_set = true;
            break;

         case 12:
         case 13:
            _epsilon_exp = _cmd->string_token();
            _epsilon_set = true;
            break;

         case 14:
         case 15:
            _omega_exp = _cmd->string_token();
            _omega_set = true;
            break;

         case 16:
         case 17:
            _beta_exp = _cmd->string_token();
            _beta_set = true;
            break;

         case 18:
         case 19:
            _v_exp = _cmd->string_token();
            _v_set = true;
            break;

         case 20:
            _young_exp = _cmd->string_token();
            _young_set = true;
            break;

         case 21:
            _poisson_exp = _cmd->string_token();
            _poisson_set = true;
            break;

         default:
            _rita->msg("pde>coef>","Unknown argument: "+kw[n]);
            _ret = 1;
            return _ret;
      }
   }
   if (nb_args>0) {
      if (_rho_set) {
         if (ieq!=HEAT && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            _rita->msg("pde>coef>","This PDE doesn't need density input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_Cp_set) {
         if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            _rita->msg("pde>coef>","This PDE doesn't need specific heat input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_kappa_set) {
         if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            _rita->msg("pde>coef>","This PDE doesn't need thermal conductivity input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_mu_set) {
         if (ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            _rita->msg("pde>coef>","This PDE doesn't need viscosity input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_sigma_set) {
         if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            _rita->msg("pde>coef>","This PDE doesn't need electric conductivity input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_Mu_set) {
         if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            _rita->msg("pde>coef>","This PDE doesn't need viscosity input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_epsilon_set) {
         if (ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            _rita->msg("pde>coef>","This PDE doesn't need electric permittivity input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_omega_set) {
         if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
            _rita->msg("pde>coef>","This PDE doesn't need angular frequency input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_beta_set) {
         if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
            _rita->msg("pde>coef>","This PDE doesn't need thermal expansion coefficient input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_v_set) {
         if (ieq!=WAVE && ieq!=TRANSPORT) {
            _rita->msg("pde>coef>","This PDE doesn't need velocity input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_young_set) {
         if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
            _rita->msg("pde>coef>","This PDE doesn't need Young's modulus input.");
            _ret = 1;
            return _ret;
         }
      }
      if (_poisson_set) {
         if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
            _rita->msg("pde>coef>","This PDE doesn't need Poisson ratio input");
            _ret = 1;
            return _ret;
         }
      }
   }
   /*
   else {
      while (1) {
         if (_cmd->readline("rita>pde>coef> ")<0)
            continue;
         switch (key=_cmd->getKW(kw)) {

            case  0:
            case  1:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands:\n";
               cout << "rho:      Density\n";
               cout << "Cp:       Specific heat at constant pressure\n";
               cout << "kappa:    Thermal conductivity\n";
               cout << "mu:       Viscosity\n";
               cout << "sigma:    Electric conductivity\n";
               cout << "Mu:       Magnetic permeability\n";
               cout << "epsilon:  Electric permittivity\n";
               cout << "omega:    Angular frequency\n";
               cout << "beta:     Thermal dilatation coefficient\n";
               cout << "v:        Velocity\n";
               cout << "young:    Young modulus\n";
               cout << "poisson:  Poisson ratio\n";
               cout << "end or <: go back to higher level" << endl;
               break;

            case  2:
               _rita->setConfigure();
               break;

            case  3:
            case  4:
               if (ieq!=HEAT && ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need density input" << endl;
                  *_ofl << "In rita>pde>coef>rho>: This PDE doesn't need density input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for density as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>rho>: Missing regular expression for density" << endl;
                  break;
               }
               ret = _cmd->get(_rho_exp);
               if (!ret) {
                  **_rita->ofh << "    rho " << _rho_exp << endl;
                  _rho_set = true;
	       }
               break;

            case  5:
            case  6:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need specific heat input" << endl;
                  *_ofl << "In rita>pde>coef>Cp>: This PDE doesn't need specific heat input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for specific heat as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Cp>: Missing regular expression for specific heat" << endl;
                  break;
               }
               ret = _cmd->get(_Cp_exp);
               if (!ret)
                  **_rita->ofh << "    Cp " << _Cp_exp << endl;
               break;

            case  7:
            case  8:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need thermal conductivity input" << endl;
                  *_ofl << "In rita>pde>coef>kappa>: This PDE doesn't need thermal conductivity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for thermal conductivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Cp>: Missing regular expression for thermal conductivity" << endl;
                  break;
               }
               ret = _cmd->get(_kappa_exp);
               if (!ret) {
                  **_rita->ofh << "    kappa " << _kappa_exp << endl;
                  _kappa_set = true;
               }
               break;

            case  9:
            case 10:
               if (ieq!=INCOMPRESSIBLE_NAVIER_STOKES && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need viscosity input" << endl;
                  *_ofl << "In rita>pde>coef>mu>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for viscosity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>mu>: Missing regular expression for viscosity" << endl;
                  break;
               }
               ret = _cmd->get(_mu_exp);
               if (!ret) {
                  **_rita->ofh << "    mu " << _mu_exp << endl;
                  _mu_set = true;
               }
               break;

            case 11:
            case 12:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need electric conductivity input" << endl;
                  *_ofl << "In rita>pde>coef>sigma>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for electric conductivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>sigma>: Missing regular expression for electric conductivity" << endl;
                  break;
               }
               ret = _cmd->get(_sigma_exp);
               if (!ret) {
                  **_rita->ofh << "    sigma " << _sigma_exp << endl;
                  _sigma_set = true;
               }
               break;

            case 13:
            case 14:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need magnetic permeability input" << endl;
                  *_ofl << "In rita>pde>coef>Mu>: This PDE doesn't need viscosity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for magnetic permeability as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>Mu>: Missing regular expression for magnetic permeability" << endl;
                  break;
               }
               ret = _cmd->get(_Mu_exp);
               if (!ret) {
                  **_rita->ofh << "    Mu " << _Mu_exp << endl;
                  _Mu_set = true;
               }
               break;

            case 15:
            case 16:
               if (ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need electric permittivity input" << endl;
                  *_ofl << "In rita>pde>coef>epsilon>: This PDE doesn't need electric permittivity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for electric permittivity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>epsilon>: Missing regular expression for electric permittivity" << endl;
                  break;
               }
               ret = _cmd->get(_epsilon_exp);
               if (!ret) {
                  **_rita->ofh << "    epsilon " << _epsilon_exp << endl;
                  _epsilon_set = true;
               }
               break;

            case 17:
            case 18:
               if (ieq!=EDDY_CURRENTS && ieq!=MAXWELL && ieq!=HELMHOLTZ) {
                  cout << "Error: This PDE doesn't need angular frequency input" << endl;
                  *_ofl << "In rita>pde>coef>omega>: This PDE doesn't need angular frequency input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give value of angular frequency.")) {
                  *_ofl << "In rita>pde>coef>omega>: Missing value of angular frequency" << endl;
                  break;
               }
               ret = _cmd->get(_omega_exp);
               if (!ret) {
                  **_rita->ofh << "    omega " << _omega_exp << endl;
                  _omega_set = true;
               }
               break;

            case 19:
            case 20:
               if (ieq!=HEAT && ieq!=COMPRESSIBLE_NAVIER_STOKES) {
                  cout << "Error: This PDE doesn't need thermal expansion coefficient input" << endl;
                  *_ofl << "In rita>pde>coef>beta>: This PDE doesn't need thermal expansion coefficient input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for thermal expansion coefficient as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>beta>: Missing regular expression for thermal expansion coefficient" << endl;
                  break;
               }
               ret = _cmd->get(_beta_exp);
               if (!ret) {
                  **_rita->ofh << "    beta " << _beta_exp << endl;
                  _beta_set = true;
               }
               break;

            case 21:
            case 22:
               if (ieq!=WAVE && ieq!=TRANSPORT) {
                  cout << "Error: This PDE doesn't need velocity input" << endl;
                  *_ofl << "In rita>pde>coef>v>: This PDE doesn't need velocity input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for velocity as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>v>: Missing regular expression for velocity" << endl;
                  break;
               }
               ret = _cmd->get(_v_exp);
               if (!ret) {
                  **_rita->ofh << "    v " << _v_exp << endl;
                  _v_set = true;
               }
               break;

            case 23:
               if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
                  cout << "Error: This PDE doesn't need Young modulus input" << endl;
                  *_ofl << "In rita>pde>coef>young>: This PDE doesn't need Young's modulus input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for Young's modulus as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>young>: Missing regular expression for Young's modulus" << endl;
                  break;
               }
               ret = _cmd->get(_young_exp);
               if (!ret) {
                  **_rita->ofh << "    young " << _young_exp << endl;
                  _young_set = true;
               }
               break;

            case 24:
               if (ieq!=LINEAR_ELASTICITY && ieq!=BEAM) {
                  cout << "Error: This PDE doesn't need Poisson ratio input" << endl;
                  *_ofl << "In rita>pde>coef>poisson>: This PDE doesn't need Poisson ratio input" << endl;
                  break;
               }
               if (_cmd->setNbArg(1,"Give regular expression for Poisson ratio as function of x,y,z,t.")) {
                  *_ofl << "In rita>pde>coef>poisson>: Missing regular expression for Poisson ratio" << endl;
                  break;
               }
               ret = _cmd->get(_poisson_exp);
               if (!ret) {
                  **_rita->ofh << "    poisson " << _poisson_exp << endl;
                  _poisson_set = true;
               }
               break;

            case 25:
            case 26:
               _ret = 0;
               **_rita->ofh << "    end" << endl;
               return;

            case 27:
            case 28:
               _ret = 100;
               return;

            case 29:
               _ret = 200;
               return;

            case -2:
            case -3:
            case -4:
               break;

            default:
               cout << "Unknown Command: " << _cmd->token() << endl;
               cout << "Available commands: rho, Cp, kappa, mu, sigma, Mu, epsilon" << endl;
	       cout << "                    omega, beta, v, young, poisson, end, <" << endl;
               cout << "Global commands:    help, ?, set, quit, exit" << endl;
               *_ofl << "In rita>pde>coef>: Unknown PDE Coefficient " << _cmd->token() << endl;
               break;
         }
      }
      }*/
   return 0;
}


void equa::set()
{
   int ret = 0;
   switch (ieq) {

      case LAPLACE:
         switch (_dim) {
  
            case 1:
               if (Sdm==FE_P1)
                  theEquation = new Laplace1DL2(*_theMesh);
               else if (Sdm==FE_P2)
                  theEquation = new Laplace1DL3(*_theMesh);
               break;

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new Laplace2DT3(*_theMesh);
               break;

            case 3:
               if (Sdm==FE_P1)
                  theEquation = new Laplace3DT4(*_theMesh);
               break;
         }
         break;

      case HEAT:
         switch (_dim) {

            case 1:
               if (Sdm==FE_P1) {
                  theEquation = new DC1DL2(*_theMesh);
                  theEquation->setTerms(LUMPED_CAPACITY|DIFFUSION);
               }
               break;

            case 2:
               if (Sdm==FE_P1 && !axi)
                  theEquation = new DC2DT3(*_theMesh);
               if (Sdm==FE_P1 && axi)
                  theEquation = new DC3DAT3(*_theMesh);
               else if (Sdm==FE_P2 && !axi)
                  theEquation = new DC2DT6(*_theMesh);
               theEquation->setTerms(LUMPED_CAPACITY|DIFFUSION);
               break;

            case 3:
               if (Sdm==FE_P1) {
                  theEquation = new DC3DT4(*_theMesh);
                  theEquation->setTerms(LUMPED_CAPACITY|DIFFUSION);
               }
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_Cp_set)
            theEquation->set_Cp(_Cp_exp);
         if (_kappa_set)
            theEquation->set_kappa(_kappa_exp);
         break;


      case LINEAR_ELASTICITY:
         switch (_dim) {

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new Elas2DT3(*_theMesh);
               else if (Sdm==FE_Q1)
                  theEquation = new Elas2DQ4(*_theMesh);
               break;

            case 3:
               if (Sdm==FE_P1)
                  theEquation = new Elas3DT4(*_theMesh);
               else if (Sdm==FE_Q1)
                  theEquation = new Elas3DH8(*_theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_young_set)
            theEquation->set_young(_young_exp);
         if (_poisson_set)
            theEquation->set_poisson(_poisson_exp);
         break;

      case TRUSS:
         switch (_dim) {

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new Bar2DL2(*_theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_young_set)
            theEquation->set_young(_young_exp);
         break;

   case INCOMPRESSIBLE_NAVIER_STOKES:
         switch (_dim) {

            case 2:
               if (Sdm==FE_P1)
                  theEquation = new TINS2DT3S(*_theMesh);
               break;

            case 3:
               if (Sdm==FE_P1)
                  theEquation = new TINS3DT4S(*_theMesh);
               break;
         }
         if (_rho_set)
            theEquation->set_rho(_rho_exp);
         if (_mu_set)
            theEquation->set_mu(_mu_exp);
         if (_beta_set)
            theEquation->set_beta(_beta_exp);
         break;
   }
   if (ret==0 && theEquation->SolverIsSet()==false)
      ls = CG_SOLVER, prec = DILU_PREC;
}

} /* namespace RITA */
