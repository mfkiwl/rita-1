/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2022 Rachid Touzani

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

                         Implementation of class 'solve'

  ==============================================================================*/

#include "rita.h"
#include "solve.h"
#include "configure.h"
#include "data.h"
#include "equa.h"
#include "stationary.h"
#include "transient.h"
#include "optim.h"
#include "eigen.h"
#include "util/macros.h"

namespace RITA {

solve::solve()
{
   _ret = 0;
   _key = 0;
}


solve::solve(rita *r, cmd *command, configure *config)
      : _rita(r), _set_analytic(false), _solved(false), _phase(false),
        _verb(1), _configure(config), _cmd(command)
{
   _ret = 0;
   _key = 0;
   _data = _rita->_data;
}


int solve::run()
{
   _nb_eq = _data->getNbEq();
   string fn="";
   int eq=1, i=0;
   static const vector<string> kw {"select","run","save","disp$lay","plot","analytic","error","post"};
   if ((_rita->_analysis_type==STEADY_STATE||_rita->_analysis_type==TRANSIENT) && _nb_eq==0) {
      _rita->msg("solve>","No equation defined.");
      _ret = 1;
      return _ret;
   }
   _nb_fields = _data->nb_fields;
   if (_nb_fields==0) {
      _rita->msg("solve>","No field defined or no problem to solve.");
      _ret = 1;
      return _ret;
   }
   _fformat.resize(_nb_fields+1);
   _isave.resize(_nb_fields+1,0);
   _save_file.resize(_nb_fields+1);
   _phase_file.resize(_nb_fields+1);
   for (int e=1; e<=_data->nb_ae; ++e)
      _data->theAE[e]->analytic.resize(_data->theAE[e]->size);
   for (int e=1; e<=_data->nb_ode; ++e)
      _data->theODE[e]->analytic.resize(_data->theODE[e]->size);
   for (int e=1; e<=_data->nb_pde; ++e) {
      if (_data->thePDE[e]->log.fail()) {
         if (_data->thePDE[e]->log.pde)
            _rita->msg("solve>","No PDE defined for equation "+to_string(e)+". Solution procedure aborted.");
         else if (_data->thePDE[e]->log.field)
            _rita->msg("solve>","Field definition incorrect for equation "+to_string(e)+". Solution procedure aborted.");
         else if (_data->thePDE[e]->log.spd)
            _rita->msg("solve>","No space discretization method available for PDE "+to_string(e)+". Solution procedure aborted.");
         else if (_data->thePDE[e]->log.ls)
            _rita->msg("solve>","No defined linear solver for PDE "+to_string(e)+". Solution procedure aborted.");
         else if (_data->thePDE[e]->log.nl)
            _rita->msg("solve>","No defined nonlinear solver for PDE "+to_string(e)+". Solution procedure aborted.");
         if (_verb>1)
            cout << "Getting back to higher level ..." << endl;
         _ret = 0;
         return _ret;
      }
      else
         _data->thePDE[e]->analytic.resize(1);
   }

   while (1) {
      int nb = _cmd->readline("rita>solve> ");
      if (nb<0)
         continue;
      switch (_key=_cmd->getKW(kw,_rita->_gkw)) {

         case   0:
            if (_cmd->setNbArg(1)==0)
               _ret = _cmd->get(eq);
            break;

         case   1:
            if (_rita->_analysis_type==STEADY_STATE) {
               if (_verb)
                  cout << "Running stationary solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_steady();
            }
            else if (_rita->_analysis_type==TRANSIENT) {
               if (_verb)
                  cout << "Running transient solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_transient();
            }
            else if (_rita->_analysis_type==OPTIMIZATION) {
               if (_verb)
                  cout << "Running optimization problem solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_optim();
            }
            else if (_rita->_analysis_type==EIGEN) {
               if (_verb)
                  cout << "Running eigen problem solver ..." << endl;
               *_rita->ofh << "  run" << endl;
               _ret = run_eigen();
            }
            else
               ;
            break;

         case   2:
            save();
            break;

         case   3:
            *_rita->ofh << "  display" << endl;
            display();
            break;

         case   4:
            _ret = plot();
            break;

         case   5:
            if (_verb)
               cout << "Setting analytical solution ..." << endl;
            *_rita->ofh << "  analytic" << endl;
            setAnalytic();
            _set_analytic = true;
            break;

         case   6:
            if (_verb)
               cout << "Computing error ..." << endl;
            eq = 1, i = 0;
            if (nb>1)
               _ret = _cmd->get(eq);
            if (nb>2)
               _ret = _cmd->get(i);
            *_rita->ofh << "  error  " << eq << "  " << i << endl;
            get_error(eq,i);
            break;

         case   7:
            if (_verb)
               cout << "Setting post calculations ..." << endl;
            *_rita->ofh << "  post" << endl;
            break;

         case 100:
         case 101:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands:\n";
            cout << "select:   Select equation (or system of equations) to solve (Default: last defined equation)\n";
            cout << "run:      Run the model\n";
            cout << "save:     Save output, can be executed before run\n";
            cout << "display:  Print solution and related data\n";
            cout << "plot:     Plot solution\n";
            cout << "analytic: Give analytic solution to test accuracy\n";
            cout << "error:    Compute error in various norms\n";
            cout << "post:     Post calculations" << endl;
            break;

         case 102:
            _rita->getLicense();
            break;

         case 103:
            if (_verb)
               cout << "Setting configuration parameter(s) ..." << endl;
            _configure->setVerbose(_verb);
            _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 104:
         case 105:
            _rita->setParam();
            break;

         case 106:
            if (_cmd->setNbArg(1,"Data name to be given.",1)) {
               _rita->msg("print>","Missing data name.","",1);
               break;
            }
            if (!_cmd->get(fn))
               _data->print(fn);
            break;

         case 107:
            _data->Summary();
            break;

         case 108:
         case 109:
            if (_verb>1)
               cout << "Getting back to higher level ..." << endl;
            *_rita->ofh << "  end" << endl;
            _ret = 0;
            return _ret;

         case -2:
            break;

         default:
            _rita->msg("solve>","Unknown command: "+_cmd->token(),
                       "Available commands for this mode:\n"
                       "run, save, display, plot, analytic, error, post");
            break;
      }
   }
}


int solve::run_steady()
{
   stationary st(_rita);
   st.setSave(_isave,_fformat,_save_file);
   int ret = 0;
   try {
      ret = st.run();
      _solved = true;
   } CATCH
   return ret;
}


int solve::run_transient()
{
   transient *ts = new transient(_rita);
   ts->setSave(_isave,_fformat,_save_file,_phase,_phase_file);
   for (int e=1; e<=_data->nb_ae; ++e) {
   }
   for (int e=1; e<=_data->nb_ode; ++e) {
   }
   for (int e=1; e<=_data->nb_pde; ++e)
      ts->setLinearSolver(_data->thePDE[e]->ls,_data->thePDE[e]->prec);
   int ret = ts->run();
/*   if (_isave[1] && theStep%_isave[1]==0) {
      if (_data->eq_type[1]==ALGEBRAIC_EQ) {
      }
      else if (_data->eq_type[1]==ODE_EQ) {
      }
      else if (_data->eq_type[1]==PDE_EQ)
         ;
   }*/
   delete ts;
   _solved = true;
   return ret;
}


int solve::run_eigen()
{
   _eigen = _rita->_eigen;
   if (!_eigen->log) {
      _rita->msg("solve>run_eigen>","Eigenproblem undefined or improperly defined.");
      return 1;
   }
   EigenProblemSolver es;
   es.setMatrix(_eigen->M);
   if (_eigen->symm) {
      es.set(OFELI::SUBSPACE,true);
      es.setNbEigv(_eigen->nb_eigv);
   }
   else {
      es.set(OFELI::QR,false);
      if (_eigen->eig_vec)
         es.setEigenVectors();
   }
   es.run();
   for (int i=1; i<=_eigen->nb_eigv; ++i) {
      (*_data->u[_data->FieldName[_eigen->eval+"-r"]])(i) = es.getEigenValue(i,1);      
      (*_data->u[_data->FieldName[_eigen->eval+"-i"]])(i) = es.getEigenValue(i,2);      
      if (_eigen->eig_vec) {
         es.getEigenVector(i,*(_data->u[_data->FieldName[_eigen->evect+"-"+to_string(i)+"r"]]),
                             *(_data->u[_data->FieldName[_eigen->evect+"-"+to_string(i)+"i"]]));
      }
   }
   cout << "Eigenvalues saved in vectors: " << _eigen->eval+"-r, " << _eigen->eval+"-i" << endl;
   cout << "Eigenvectors saved in vectors: " << _eigen->evect+"-*r*, " << _eigen->evect+"-*i*, " << endl;
   return 0;
}


int solve::run_optim()
{
   _optim = _rita->_optim;
   if (_optim->log) {
      _rita->msg("solve>run_optim>","Optimization problem undefined or improperly defined.");
      return 1;
   }
   int size=_optim->size;
   try {
      if (_optim->lp) {
         OFELI::LPSolver s;
         s.setSize(size,_optim->nb_lec,_optim->nb_gec,_optim->nb_eqc);
         s.set(*_data->u[_data->iField]);
         s.set(OFELI::LPSolver::OBJECTIVE,_optim->a,_optim->b);
         for (int i=0; i<_optim->nb_eqc; ++i)
            s.set(OFELI::LPSolver::EQ_CONSTRAINT,*_optim->a_eq[i],_optim->b_eq[i]);
         for (int i=0; i<_optim->nb_lec; ++i)
            s.set(OFELI::LPSolver::LE_CONSTRAINT,*_optim->a_le[i],_optim->b_le[i]);
         for (int i=0; i<_optim->nb_gec; ++i)
            s.set(OFELI::LPSolver::GE_CONSTRAINT,*_optim->a_ge[i],_optim->b_ge[i]);
         int ret = s.run();
         _optim->obj = s.getObjective();
         if (ret==0)
            _optim->solved = true;
      }
      else {
         OFELI::OptSolver s(*_data->u[_data->iField]);
         s.setOptMethod(_optim->Alg);
         s.setObjective(*_optim->J_Fct);
         if (_optim->G_ok) {
            for (int i=0; i<size; ++i)
               s.setGradient(*_data->theFct[_optim->igrad+i],i+1);
         }
         if (_optim->H_ok) {
            for (int i=0; i<size*size; ++i)
               s.setHessian(*_data->theFct[_optim->ihess+i],i+1);
         }
         for (int i=0; i<_optim->nb_lec; ++i)
            s.setIneqConstraint(*_optim->inC_Fct[i],_optim->penal);
         for (int i=0; i<_optim->nb_eqc; ++i)
            s.setEqConstraint(*_optim->eqC_Fct[i],_optim->penal);
         s.setLowerBounds(_optim->lb);
         s.setUpperBounds(_optim->ub);
         if (!s.run()) {
            _optim->obj = s.getObjective();
            _optim->solved = true;
         }
      }
   } CATCH
   return 0;
}


void solve::save()
{
   int k=0, freq=1, eq=0, field_ok=0;
   string file="rita-1.pos", fformat="gmsh", ext="", fd="u", phase_f="";
   _phase = false;
   if (_nb_fields==1) {
      fd = _data->Field[1];
      eq = _data->FieldEquation[1];
      _fformat[1] = GNUPLOT;
      _save_file[1] = "u.dat";
      _isave[1] = 1;
      field_ok = 1;
/*      if (_data->eq_type[0]==data::eqType::PDE) {
         _fformat[_data->iField] = GMSH;
         _save_file[_data->iField] = "u.pos";
      }*/
   }
   _ret = 0;

   static const vector<string> kw {"field","format","freq$uency","phase","file"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            fd = _cmd->string_token();
            field_ok = 2;
            break;

         case 1:
            fformat = _cmd->string_token();
            break;

         case 2:
            freq = _cmd->int_token();
            break;

         case 3:
            phase_f = _cmd->string_token();
            _phase = true;
            break;

         case 4:
            file = _cmd->string_token();
            break;

         default:
            _rita->msg("solve>save>","Unknown argument: "+kw[n]);
            return;
      }
   }
   if (nb_args>0) {
      if (_data->nb_fields==0) {
         _rita->msg("solve>","No field to save.");
         _ret = 1;
         return;
      }
      if (!field_ok) {
         _rita->msg("solve>:","No field name given.");
         _ret = 1;
         return;
      }
      k = _data->checkField(fd);
      if (k<0) {
         _rita->msg("solve>:","Unknown field: "+fd);
         _ret = 1;
         return;
      }
      eq = _data->FieldEquation[k];
      if (_data->eq_type[eq]==data::eqType::PDE) {
         _isave[k] = freq;
         _fformat[k] = GMSH;
         _save_file[k] = fd + ".pos";
      }
      else if (_data->eq_type[eq]==data::eqType::ODE || (_data->eq_type[eq]==data::eqType::AE)) {
         _isave[k] = freq;
         _fformat[k] = GNUPLOT;
         _save_file[k] = fd + ".dat";
      }
      eq = _data->FieldEquation[k];
      if (_data->eq_type[eq]!=data::eqType::PDE && fformat != "gnuplot") {
         _rita->msg("solve>save>","Only Gnuplot format is available for this type of equation.");
         _ret = 1;
         return;
      }
      if (!_ret)
         _fformat[k] = _ff[fformat];
      else {
         _rita->msg("solve>save>","Unknown file format: "+fformat);
         _ret = 1;
         return;
      }
      _save_file[k] = file;
      *_rita->ofh << "  save field=" << fd << " file=" << file << " format=" << fformat
                  << " frequency=" << freq;
      if (_phase) {
         _phase_file[k] = phase_f;
         *_rita->ofh << "  phase=" << phase_f;
      }
      *_rita->ofh << endl;
   }
   _ret = 0;
   return;
}


void solve::setAnalytic()
{
   int eq=1, comp=1;
   bool exp_ok=false;
   string exp="";
   _ret = 0;

   static const vector<string> kw {"eq$uation","comp$onent","def$inition"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case  0:
            eq = _cmd->int_token();
            break;

         case 1:
            comp = _cmd->int_token();
            break;

         case 2:
            exp = _cmd->string_token();
            exp_ok = true;
            break;

         default:
            _rita->msg("solve>analytic>","Unknown argument: "+kw[n]);
            _ret = 1;
            return;
      }
   }
   if (nb_args>0) {
      if (!exp_ok) {
         _rita->msg("solve>analytic>","No expression of analytical solution provided.");
         _ret = 1;
         return;
      }
      if (_data->eq_type[eq]==data::eqType::AE)
         _data->theAE[eq]->analytic[comp-1] = exp;
      else if (_data->eq_type[eq]==data::eqType::ODE)
         _data->theODE[eq]->analytic[comp-1] = exp;
      else if (_data->eq_type[eq]==data::eqType::PDE)
         _data->thePDE[eq]->analytic[comp-1] = exp;
      *_rita->ofh << "analytic equation=" << eq << " component=" << comp << " expression=" << exp << endl;
   }
   _ret = 0;
   return;
}


void solve::display(int f)
{
   if (f==-1)
      f = _data->iField;
   if (_data->FieldType[f]==data::eqType::AE) {
      if (_data->u[f]->size()==1)
         cout << "Solution of algebraic equation: " << (*_data->u[f])[0] << endl;
      else {
         cout << "Solution of algebraic equation: " << endl;
         for (size_t i=0; i<_data->u[f]->size(); ++i)
            cout << (*_data->u[f])[i] << endl;
      }
   }
   else if (_data->FieldType[f]==data::eqType::ODE) {
      if (_data->u[f]->size()==1)
         cout << "Solution of ordinary differential equation: "
              << (*_data->u[f])[0] << endl;
      else {
         cout << "Solution of ordinary differential equation: " << endl;
         for (size_t i=0; i<_data->u[f]->size(); ++i)
            cout << (*_data->u[f])[i] << endl;
      }
   }
   else if (_data->FieldType[f]==data::eqType::PDE) {
      cout << "Solution of partial differential equation: " << endl;
      cout << "Solution vector:\n" << *_data->u[f];
   }
   else if (_data->FieldType[f]==data::eqType::OPT) {
      if (!_optim->solved) {
         _rita->msg("solve>display>","Optimization problem ot properly solved.");
         _ret = 1;
         return;
      }
      if (_data->u[f]->size()==1)
         cout << "Solution of optimization problem: " << (*_data->u[f])[0] << endl;
      else {
         cout << "Solution of optimization problem: " << endl;
         for (size_t i=0; i<_data->u[f]->size(); ++i)
            cout << (*_data->u[f])[i] << endl;
      }
      cout << "Optimal objective: " << _rita->_optim->obj << endl;
   }
   else if (_data->FieldType[f]==data::eqType::EIGEN) {
      _eigen->verbose = 2;
      cout << *_eigen;
   }
}


int solve::plot()
{
   OFELI::saveField(*_data->u[1],"rita.pos",GMSH);
   string com = "gmsh rita.pos";
   int ret = system(com.c_str());
   remove("rita.pos");
   return ret;
}


void solve::get_error(int eq, int i)
{
   double err2=0., errI=0.;
   vector<double> y;
   if (!_set_analytic) {
      _rita->msg("solve>error>","No analytical solution given.");
      return;
   }
   if (!_solved) {
      _rita->msg("solve>error>","Problem has not been solved yet.");
      return;
   }

// Case of an algebraic equation
   if (_data->eq_type[eq]==data::eqType::AE) {
      int e = _data->eqq[eq];
      _ae = _data->theAE[e];
      _var.resize(_ae->size);
      for (int i=0; i<_ae->size; ++i) {
         y.resize(_ae->size);
         for (int j=0; j<_ae->size; ++j) {
            _var[j] = _ae->fn+to_string(j+1);
            y[j] = _ae->y[j];
         }
         _theFct.set(_ae->analytic[i],_var);
         double u=_ae->y[i], v=_theFct(y);
         errI += (u-v)*(u-v);
      }
      cout << "Error: " << sqrt(errI) << endl;
   }

// Case of an ODE
   else if (_data->eq_type[eq]==data::eqType::ODE) {
      int e = _data->eqq[eq];
       _ode = _data->theODE[e];
       y.resize(_ode->size+1);
      _var.resize(_ode->size+1);
      _var[0] = "t";
      y[0] = _rita->_final_time;
      for (int i=0; i<_ode->size; ++i) {
         for (int j=0; j<_ode->size; ++j) {
            _var[j+1] = _ode->fn+to_string(j+1);
            y[j+1] = _ode->y[j];
         }
         _theFct.set(_ode->analytic[i],_var);
         double u=_ode->y[i], v=_theFct(y);
         errI += (u-v)*(u-v);
      }
      cout << "Error: " << sqrt(errI) << endl;
   }

// Case of a PDE
   else if (_data->eq_type[eq]==data::eqType::PDE) {
      int e = _data->eqq[eq];
      _pde = _data->thePDE[e];
      int nb_dof = 0, neq = 0;
      OFELI::Mesh *theMesh = _data->theMesh[_data->iMesh];
      OFELI::Grid *theGrid = _data->theGrid[_data->iGrid];
      if (theMesh!=nullptr) {
         neq = theMesh->getNbDOF();
         nb_dof = neq/theMesh->getNbNodes();
         for (int i=0; i<nb_dof; ++i) {
            _theFct.set(_pde->analytic[i]);
            for (int n=1; n<=int(theMesh->getNbNodes()); ++n) {
               double u = (*_data->u[eq])(n,i+1);
               double v = _theFct((*theMesh)[n]->getCoord());
               err2 += (u-v)*(u-v);
               errI  = std::max(fabs(u-v),errI);
            }
         }
         err2 = sqrt(err2/theMesh->getNbNodes());
      }
      else if (theGrid!=nullptr) {
         neq = theGrid->getNbDOF();
         nb_dof = neq/theGrid->getNbNodes();
      }
      cout << "Error in L2-Norm:         " << err2 << endl;
      cout << "Error in L-Infinity-Norm: " << errI << endl;
   }

// Case of an Optimization problem
   else if (_data->eq_type[eq]==data::eqType::OPT)
      _rita->msg("solve>error>","Error computation is not available for optimization problems.");
}

} /* namespace RITA */
