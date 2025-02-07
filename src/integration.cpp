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

                      Implementation of class 'integration'

  ==============================================================================*/


#include "integration.h"
#include "data.h"

namespace RITA {

integration::integration(rita*      r,
                         cmd*       command,
                         configure* config)
            : _rita(r), _configure(config), _cmd(command)
{
}


integration::~integration()
{
}


int integration::run()
{
   int ret = 0;
   _rita->_analysis_type = INTEGRATION;
   string fct="", def="", var_name="", form="trapezoidal";
   int count_fct=0, count_def=0, count_vector=0;
   int nb=0, ind=0;
   nx = ny = nz = 10;
   xmin = ymin = zmin = 0.;
   xmax = ymax = zmax = 1.;
   dim = 1;
   unif = 1;
   ng = 1;
   var.clear();
   data *theData = _rita->_data;
   static const string H = "Command: integration [function=f] [definition=exp] [interval=min,max] [variable=x] [ne=nx]\n"
                           "                     [formula=f] [display]\n\n"
                           "f: Name of already defined function to integrate\n"
                           "exp: Expression defining function to integrate\n"
                           "x: name of variable\n"
                           "min, max: Integration is made on the interval (min,max)\n"
                           "ne: Number of subdivisions of the interval\n"
                           "m: Numerical integration formula. To choose among the values: left-rectangle,\n"
                           "   right-rectangle, mid-point, trapezoidal, simpson, gauss-legendre, gauss-lobatto.\n"
	                        "   For the Gauss formulae, the number of points can be specified by typing \n"
 	                        "   gauss-legendre,2 for instance.";
   static const vector<string> kw {"func$tion","def$inition","var$iable","vect$or","interval","domain",
                                   "ne","form$ula"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<1) {
      _rita->msg("integration>","No argument for command.");
      return 1;
   }
   for (int k=0; k<nb_args; ++k) {

      int n = _cmd->getArgs(nb);
      if (nb<1 && n!=11) {
         _rita->msg("integration>","Insufficient number of arguments.");
         return 1;
      }

      switch (n) {

         case 100:
         case 101:
            cout << H << endl;
            return 0;

         case   0:
            fct = _cmd->string_token(0);
            count_fct++;
            break;

         case   1:
            def = _cmd->string_token(0);
            count_def++;
            break;

         case   2:
         case   3:
            var_name = _cmd->string_token(0);
            count_vector++;
            break;

         case   4:
         case   5:
            ret  = theData->getPar(0,"integration>",xmin);
            ret += theData->getPar(1,"integration>",xmax);
            break;

         case   6:
            if (nb==1) {
               ret  = theData->getPar(0,"integration>",nx);
               ny = nz = nx;
            }
            else {
               ret  = theData->getPar(0,"integration>",nx);
               ret += theData->getPar(1,"integration>",ny);
               ret += theData->getPar(2,"integration>",nz);
            }
            break;

         case   7:
            form = _cmd->string_token(0);
            if (nb>1)
               ng = _cmd->int_token(1);
            break;

         default:
            _rita->msg("integration>","Unknown argument: "+_cmd->Arg());
            return 1;
      }
   }

   if (nb_args>0) {
      if (count_fct && count_vector) {
         _rita->msg("integration>","Function already defined in data module.");
         return 1;
      }
      if (count_fct && count_def) {
         _rita->msg("integration>","Function already defined.");
         return 1;
      }
      if (count_fct>1 || count_def>1) {
         _rita->msg("integration>","Too many functions defined.");
         return 1;
      }
      if (count_def && !count_vector) {
         _rita->msg("integration>","Missing a variable name.");
         return 1;
      }
      *_rita->ofh << "integration";
      if (dim==1)
         *_rita->ofh << " interval=" << xmin << "," << xmax << " ne=" << nx;
      if (count_fct) {
         ind = theData->checkName(fct,DataType::FCT);
         if (ind==-1) {
            _rita->msg("integration>","Non defined function "+fct);
            return 1;
         }
         IFct = theData->theFct[ind];
         *_rita->ofh << " function=" << fct;
      }
      else {
         theData->addVector(var_name,0.,dim);
         *_rita->ofh << " var=" << var_name << " definition=" << def;
         if (dim==1)
            var.push_back(var_name);
         else {
            for (int i=0; i<dim; ++i)
               var.push_back(var_name+to_string(i+1));
         }
         theData->addFunction(def,var);
         if (theData->iVector<theData->nb_vectors-1)
            theData->VectorType[theData->iVector] = data::eqType::INTEGR;
         else
            theData->VectorType.push_back(data::eqType::INTEGR);
         IFct = theData->theFct[theData->iFct];
      }
      nim = Nint[form];
      *_rita->ofh << " formula=" << form;
      if (nim==GAUSS_LEGENDRE || nim==GAUSS_LOBATTO)
         *_rita->ofh << "," << ng;
      *_rita->ofh << endl;
   }
   return 0;
}


int integration::go()
{
   res = 0.;
   double dx=0., dy=0., dz=0.;
   if (unif==1)
      dx = (xmax-xmin)/nx, dy = (ymax-ymin)/ny, dz = (zmax-zmin)/nz;
   static const vector<double> xg {0.,-0.5773502691896257,0.5773502691896257,0.,-0.7745966692414834,
                                   0.7745966692414834,-0.3399810435848563,0.3399810435848563,
                                   -0.8611363115940526,0.8611363115940526,0.,-0.5384693101056831,
                                   0.5384693101056831,-0.9061798459386640,0.9061798459386640,
                                   -0.6612093864662645,0.6612093864662645,-0.2386191860831969,
                                   0.2386191860831969,-0.9324695142031521,0.9324695142031521};
    static const vector<double> wg {2.,1.,1.,0.8888888888888889,0.5555555555555556,0.5555555555555556,
                                    0.6521451548625461,0.6521451548625461,0.3478548451374538,0.3478548451374538,
                                    0.5688888888888889,0.4786286704993665,0.4786286704993665,0.2369268850561891,
                                    0.2369268850561891,0.3607615730481386,0.3607615730481386,0.4679139345726910,
                                    0.4679139345726910,0.1713244923791704,0.1713244923791704};
    static const vector<double> xl {0.,-1.0,1.0,-0.447213595499958,0.447213595499958,-1.0,1.0,0.,-0.475481063909541,
                                    0.475481063909541,-1.0,1.0};
    static const vector<double> wl {1.333333333333333,1.333333333333333,0.333333333333333,0.833333333333333,
                                    0.833333333333333,0.166666666666667,0.166666666666667,0.711111111111111,
                                    0.544444444444444,0.544444444444444,0.1,0.1};

   _x.resize(nx+1), _y.resize(ny+1), _z.resize(nz+1);
   _x[0] = xmin, _y[0] = ymin, _z[0] = zmin;
   for (int i=1; i<=nx; ++i)
      _x[i] = _x[i-1] + dx;
   for (int i=1; i<=ny; ++i)
      _y[i] = _y[i-1] + dy;
   for (int i=1; i<=nz; ++i)
      _z[i] = _z[i-1] + dz;
   for (int ii=0; ii<nx; ++ii) {
      double x1=_x[ii], x2=_x[ii+1], x12=0.5*(_x[ii]+_x[ii+1]); 
      if (nim==LRECTANGLE)
         res += dx*(*IFct)(x1);
      else if (nim==RRECTANGLE)
         res += dx*(*IFct)(x2);
      else if (nim==MIDPOINT)
         res += dx*(*IFct)(x12);
      else if (nim==TRAPEZOIDAL)
         res += 0.5*dx*((*IFct)(x1)+(*IFct)(x2));
      else if (nim==SIMPSON)
         res += ((*IFct)(x1)+4*(*IFct)(x12)+(*IFct)(x2))*dx/6.0;
      else if (nim==GAUSS_LEGENDRE) {
         if (ng<1 || ng>5) {
            _rita->msg("integration>","For Gauss-Legendre formula, number of points must be between 1 and 6.");
            return 1;
         }
         for (int i=0; i<ng; ++i) {
            int j = ng*(ng-1)/2 + i;
            res += 0.5*dx*wg[j]*(*IFct)(x12+0.5*dx*xg[j]);
         }
      }
      else if (nim==GAUSS_LOBATTO) {
         if (ng<3 || ng>5) {
            _rita->msg("integration>","For Gauss-Lobatto formula, number of points must be between 3 and 5.");
            return 1;
         }
         for (int i=0; i<ng; ++i) {
            int j = ng*(ng-1)/2 + i - 3;
            res += 0.5*dx*wl[j]*(*IFct)(x12+0.5*dx*xl[j]);
         }
      }
   }
   cout << "Approximate Integral: " << res << endl;
   return 0;
}

} /* namespace RITA */
