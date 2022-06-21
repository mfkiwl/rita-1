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

                      Implementation of class 'eigen'

  ==============================================================================*/

#include "data.h"
#include "configure.h"
#include "eigen.h"
#include "linear_algebra/Matrix_impl.h"

namespace RITA {

eigen::eigen(rita*      r,
             cmd*       command,
             configure* config)
      : _rita(r), _configure(config), _cmd(command)
{
   _verb = _rita->_verb;
   nb_eigv = 0;
   eig_vec = false;
   verbose = 0;
   _data = _rita->_data;
   solved = false;
   log = true;
}


eigen::~eigen()
{
}


int eigen::run()
{
   _rita->_analysis_type = EIGEN;
   string mat_name="", method="qr";
   evect = "";
   symm = eig_vec = false;
   nb_eigv = 0;

   static const vector<string> kw {"matrix","symm$etric","method","nb","eigv","evect","summary"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("eigen>","No argument to command");
      NO_EIGEN
      return 0;
   }
   for (int k=0; k<nb_args; ++k) {

      int n = _cmd->getArg("=");
      switch (n) {

         case 0:
            mat_name = _cmd->string_token();
            break;

         case 1:
            symm = true;
            break;

         case 2:
            method = _cmd->string_token();
            break;

         case 3:
            nb_eigv = _cmd->int_token();
            break;

         case 4:
            eig_vec = _cmd->int_token();
            break;

         case 5:
            evect = _cmd->string_token();
            break;

         case 100:
         case 101:
            cout << "\nAvailable Commands:\n";
            cout << "matrix, symmetric, method, nb, eigv, evect, summary\n";
            break;

         default:
            _rita->msg("eigen>","Unknown argument: "+_cmd->Arg());
            log = true;
            return 1;
      }
   }

   if (nb_args>0) {
      if (mat_name=="") {
         _rita->msg("eigen>","No matrix given.");
         NO_EIGEN
         return 1;
      }
      int k = _data->checkMatrix(mat_name);
      if (k<0) {
         _rita->msg("eigen>","Matrix "+mat_name+" not defined.");
         NO_EIGEN
         return 1;
      }
      M = _data->theMatrix[k];
      int nr=M->getNbRows(), nc=M->getNbColumns();
      if (nr!=nc) {
         _rita->msg("eigen>","Matrix "+mat_name+" must be a square matrix.");
         NO_EIGEN
         return 1;
      }
      size = nr;
      if (nb_eigv<0) {
         _rita->msg("eigen>","Illegal number of eigen values: "+to_string(nb_eigv));
         NO_EIGEN
         return 1;
      }
      if (nb_eigv==0)
         nb_eigv = size;
      Alg = meth[method];
      if (Alg==OFELI::SUBSPACE)
         eig_vec = true;
      if (Alg!=OFELI::SUBSPACE && Alg!=OFELI::QR) {
         _rita->msg("eigen>","Method "+to_string(Alg)+" not available");
         NO_EIGEN
         return 1;
      }
      *_rita->ofh << "eigen matrix=" << mat_name;
      if (symm)
         *_rita->ofh << " symmetric";
      if (nb_eigv<size)
         *_rita->ofh << " nb=" << nb_eigv;
      *_rita->ofh << " method=" << method << endl;
      if (evect=="")
         evect = mat_name+"-ev";
      if (eig_vec) {
         for (int i=1; i<=nb_eigv; ++i) {
            _data->addField(evect+"-"+to_string(i)+"r");
            _data->addField(evect+"-"+to_string(i)+"i");
         }
      }
      eval = mat_name + "-ev";
      _data->addField(eval+"-r");
      _data->addField(eval+"-i");
      _data->u[_data->FieldName[eval+"-r"]]->setSize(nb_eigv);
      _data->u[_data->FieldName[eval+"-i"]]->setSize(nb_eigv);
   }

   else {
   }
   _rita->_ret = 0;
   return 0;
}


void eigen::print(ostream& s) const
{
   if (log) {
      cout << "Eigen problem improperly or not defined !" << endl;
      return;
   }
   if (!solved) {
      cout << "Nothing to output: Eigen problem unsolved !" << endl;
      return;
   }
   s << "Eigen problem solver" << endl;
   s << "Matrix: " << " " << ", size: " << endl;
   if (verbose==0)
      return;
   s << "Number of computed eigenvalues: " << nb_eigv << endl;
   for (int i=1; i<=nb_eigv; ++i) {
      s << "Eigenvalue #" << i << ": ";
      if (symm)
         s << (*_data->u[_data->FieldName[eval+"-r"]])(i) << endl;
      else
      s << (*_data->u[_data->FieldName[eval+"-r"]])(i)
        << " + " << (*_data->u[_data->FieldName[eval+"-i"]])(i) << "I" << endl;
   }
   if (!eig_vec || verbose==1)
      return;
   s << "Eigenvectors" << endl;
   for (int i=1; i<=nb_eigv; ++i) {
      s << "Eigenvector " << i << ": ";
      if (symm)
         for (int j=1; j<=size; ++i)
            s << (*_data->u[_data->FieldName[eval+"-r"]])(i) << endl;
      else
         for (int j=1; j<=size; ++i)
            s << (*_data->u[_data->FieldName[eval+"-r"]])(i) << " + " 
              <<  (*_data->u[_data->FieldName[eval+"-i"]])(i) << "I" << endl;
   }
}


ostream& operator<<(ostream& s, const eigen& es)
{
   es.print(s);
   return s;
}

} /* namespace RITA */
