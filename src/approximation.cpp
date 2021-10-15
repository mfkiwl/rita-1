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

                      Implementation of class 'approx'

  ==============================================================================*/


#include "approximation.h"
#include "data.h"

namespace RITA {

approximation::approximation(rita* r, cmd* command, configure* config)
              : _rita(r), _configure(config), _cmd(command)
{
   _data = _rita->_data;
   _tab = nullptr;
   _verb = _configure->getVerbose();
}


approximation::~approximation()
{
   if (_tab==nullptr)
      delete _tab;
}


int approximation::run()
{
   _rita->_analysis_type = APPROXIMATION;
   string file="", name="";
   int nb=0, nb_tab=0;
   int file_count=0, lagrange_count=0, bspline_count=0, fitting_count=0, bezier_count=0;
   int nurbs_count=0, approx_count=0;
   static const vector<string> kw {"file","name","lagrange","piecewise$_lagrange","hermite",
                                   "fitting","bspline","bezier","nurbs"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<1) {
      _rita->msg("approximation>","No argument for command.");
      return 1;
   }
   for (int k=0; k<nb_args; ++k) {

      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            file = _cmd->string_token(0);
            nb_tab = 1;
            file_count++;
            break;

         case 1:
            name = _cmd->string_token(0);
            break;

         case 2:
            _lagrange_degree = _cmd->int_token(0);
            lagrange_count++, approx_count++;
            _method = LAGRANGE;
            break;

         case 3:
            approx_count++;
            _method = PIECEWISE_LAGRANGE;
            break;

         case 4:
            approx_count++;
            _method = HERMITE;
            _rita->msg("approximation>","This approximation is not implemented yet.");
            return 1;

         case 5:
            fitting_count++, approx_count++;
            _method = FITTING;
            _rita->msg("approximation>","This approximation is not implemented yet.");
            return 1;

         case 6:
            bspline_count++, approx_count++;
            _method = BSPLINE;
            _rita->msg("approximation>","This approximation is not implemented yet.");
            return 1;

         case 7:
            bezier_count++, approx_count++;
            _method = BEZIER;
            _rita->msg("approximation>","This approximation is not implemented yet.");
            return 1;

         case 8:
            nurbs_count++, approx_count++;
            _method = NURBS;
            _rita->msg("approximation>","This approximation is not implemented yet.");
            return 1;

         default:
            _rita->msg("approximation>","Unknown argument: "+_cmd->Arg());
            return 1;
      }
   }

   if (nb_args>0) {
     if (file_count==0) {
         _rita->msg("approximation>","No data file given.");
         return 1;
      }
      if (approx_count>1) {
         _rita->msg("approximation>","More than one approximation method given.");
         return 1;
      }
      *_rita->ofh << "approximation";
      *_rita->ofh << " file=" << file;
      if (lagrange_count++)
         *_rita->ofh << " lagrange=" << _lagrange_degree;
      *_rita->ofh << endl;
      _tab = new OFELI::Tabulation(file);
      _data->addTab(_tab,name);
   }
   /*   else {
      *_ofh << "approximation " << endl;
      while (1) {
      if (_cmd->readline("approximation> ")<0)
         continue;
      switch (_key=_cmd->getKW(_kw_approx)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands:\n";
            cout << "file:         File containing tabulated data to approximate\n";
            cout << "name:         "
            cout << "lagrange:     \n";
            cout << "fitting:     Least Square fitting\n";
            cout << "bspline";
            cout << "bezier";
            cout << "nurbs";
            cout << "summary:      Summary of approximation problem attributes\n";
            cout << "clear:        Remove problem\n";
            cout << "end or <:     go back to higher level" << endl;
            break;

         case 2:
            setConfigure();
            break;

         case 3:
            if (!_cmd->get(size))
               *_ofh << "    size " << size << endl;
            _ret = 0;
            break;

         case 4:
            _ret = 0;
            break;

         case 5:
            _ret = 0;
            break;

         case 7:
            break;

         case 8:
            _ret = 0;
            break;

         case 9:
            break;

         case 11:
            cout << "Summary of optimization problem attributes:\n";
            *_ofh << "      summary" << endl;
            _ret = 0;
            break;

         case 13:
         case 14:
            _cmd->setNbArg(0);
            if (fct_ok<size) {
               cout << "Error: Insufficient number of functions defining system." << endl;
               *_ofl << "In optimization>end>: Insufficient number of functions defining system." << endl;
               break;
            }
            if (field_ok==0) {
               cout << "Error: No field defined for optimization problem." << endl;
               *_ofl << "In optimization>end>: No field defined for optimization problem." << endl;
               break;
            }
            *_ofh << "      end" << endl;
            _ret = 0;
            return;

         case 15:
         case 16:
            _ret = 100;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            cout << "Unknown Command: " << _cmd->token() << endl;
            cout << "Available commands: size, objective, gradient, hessian, constraint, init, field, algorithm" << endl;
	    cout << "                    summary, clear, end, <" << endl;
            cout << "Global commands:    help, ?, set, quit, exit" << endl;
            *_ofl << "In optimization>: Unknown Command " << _cmd->token() << endl;
         break;
      }
   }
   _ret = 0;*/
   return 0;
}


int approximation::go()
{
   switch (_method) {

      case LAGRANGE:
         lagrange();
         break;

      case PIECEWISE_LAGRANGE:
         piecewise_lagrange();
         break;

      case HERMITE:
         hermite();
         break;

      case FITTING:
         fitting();
         break;

      case BSPLINE:
         bspline();
         break;

      case BEZIER:
         bezier();
         break;

      case NURBS:
         nurbs();
         break;
   }
   return 0;
}


int approximation::lagrange()
{
   if (_lagrange_degree>4)
      _rita->msg("approximation>","Polynomial degree is too large.");
   if (_tab->getNbVar(1)>1) { 
      _rita->msg("approximation>","This approximation method is available for one-variable cases only.");
      return 1;
   }
   int np = _tab->getSize(1,1);
   if (np != _lagrange_degree+1) {
      _rita->msg("approximation>","Required interpolation degree is incompatible with number of points in tabulation.");
      return 1;
   }
   _x.resize(np);
   _y.resize(np);
   double h = (_tab->getMaxVar(1,1)-_tab->getMinVar(1,1))/(np-1);
   _x[0] = _tab->getMinVar(1,1);
   _y[0] = _tab->Funct[0].Val(1);
   for (int i=1; i<np; ++i) {
      _x[i] = _x[i-1] + h;
      _y[i] = _tab->Funct[0].Val(i+1);
   }
   
// Lagrange polynomial
   string p="";
   for (int i=0; i<np; ++i) {
      double d = 1;
      for (int j=0; j<np; ++j) {
         if (i!=j) {
            d *= (_x[i]-_x[j]);
            p = p + "(x-" + to_string(_x[j]) + ")*";
         }
      }
      p = p + "(" + to_string(_y[i]/d) + ")";
      if (i<np-1)
         p = p + " + ";
   }
   vector<string> var;
   var.push_back("x");
   _data->addFunction(p,var,"");
   _data->setTab2Grid(_tab);
   _data->setTab2Field(_tab);
   if (_verb)
      cout << "Lagrange interpolation created. Lagrange polynomial is: " 
           << _data->fct_name[_data->iFct] << endl;
   return 0;
}


int approximation::piecewise_lagrange()
{
   if (_tab->getNbVar(1)>1) { 
      _rita->msg("approximation>","This approximation method is available for one-variable cases only.");
      return 1;
   }
   int np = _tab->getSize(1,1);
   _x.resize(np);
   _y.resize(np);
   double h = (_tab->getMaxVar(1,1)-_tab->getMinVar(1,1))/(np-1);
   _x[0] = _tab->getMinVar(1,1);
   _y[0] = _tab->Funct[0].Val(1);
   for (int i=1; i<np; ++i) {
      _x[i] = _x[i-1] + h;
      _y[i] = _tab->Funct[0].Val(i+1);
   }
   _data->setTab2Grid(_tab);
   _data->setTab2Field(_tab);
   if (_verb)
      cout << "Piecewise Lagrange interpolation created. Interpolated field is: " 
           << _data->Field[_data->iField] << endl;
   return 0;
}


int approximation::hermite()
{
   return 0;
}


int approximation::fitting()
{
   return 0;
}


int approximation::bspline()
{
   return 0;
}


int approximation::bezier()
{
   return 0;
}


int approximation::nurbs()
{
   return 0;
}


int approximation::eval(double x)
{
   switch (_method) {

      case LAGRANGE:
         break;

      case PIECEWISE_LAGRANGE:
         break;

      case HERMITE:
         break;

      case FITTING:
         break;

      case BSPLINE:
         break;

      case BEZIER:
         break;

      case NURBS:
         break;
   }
   return 0;

}

} /* namespace RITA */
