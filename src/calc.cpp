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

                          Implementation of class 'calc'

  ==============================================================================*/


#define _USE_MATH_DEFINES 
#undef __STRICT_ANSI__
#if defined(_WIN32) 
  // Memory leak dumping
  #if defined(_DEBUG)
    #define _CRTDBG_MAP_ALLOC
    #include <stdlib.h>
    #include <crtdbg.h>
    #define CREATE_LEAKAGE_REPORT
  #endif

  // Needed for windows console UTF-8 support
  #include <fcntl.h>
  #include <io.h>
#endif

#include "calc.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <locale>
#include <limits>
#include <ios> 
#include <iomanip>
#include <numeric>

#include "rita.h"
#include "cmd.h"
#include "data.h"

namespace RITA {

calc::calc(rita* r,
           cmd*  command)
     : _rita(r), _cmd(command), _data(r->_data)
{
   parser = new ParserX(pckALL_NON_COMPLEX);
   parser->EnableAutoCreateVar(true);
}


calc::~calc()
{
   delete parser;
}


Value calc::Help()
{
   cout << "-----------------------------------------------------------\n";
   cout << "Commands:\n\n";
   cout << "  list var     - list parser variables\n";
   cout << "  list exprvar - list expression variables\n";
   cout << "  list const   - list all numeric parser constants\n";
   cout << "  end (or <)   - exits the parser\n";
   cout << "\nConstants:\n\n";
   cout << "  \"e\"   2.718281828459045235360287\n";
   cout << "  \"pi\"  3.141592653589793238462643\n";
   cout << "-----------------------------------------------------------\n";
   return 0;
}


void calc::ListVar(const ParserXBase* p)
{
   var_maptype variables = p->GetVar();
   if (!variables.size())
      return;
   cout << "\nParser variables:\n";
   cout << "-----------------\n";
   cout << "Number: " << int(variables.size()) << "\n";
   for (const auto& v: variables)
      cout << "Name: " << v.first << " = " << *v.second << std::endl;
}


void calc::ListConst(const ParserXBase* p)
{
   cout << "\nParser constants:\n";
   cout << "-----------------\n";
   var_maptype cmap = p->GetConst();
   if (!cmap.size())
      cout << "Expression does not contain constants\n";
	else {
      for (const auto& v: cmap)
         cout << "  " << v.first << " = " << *v.second << std::endl;
   }
}


void calc::setData(ParserXBase* p)
{
   var_maptype variables = p->GetVar();
   if (!variables.size())
      return;
   for (const auto& v: variables) {
      Variable &w = (Variable&)(*v.second);
      switch (w.GetType()) {

         case 'i': 
         case 'f': 
            {
               _data->addParam(v.first,w.GetFloat());
               break;
            }

         case 'm':
            {
               int nr=w.GetRows(), nc=w.GetCols();
               if (nr==1 || nc==1) {
                  _u = new OFELI::Vect<double>(nr,nc);
                  _data->addField(v.first,std::max(nr,nc));
                  for (int i=1; i<=std::max(nr,nc); ++i)
                     (*_u)(i,1) = w.GetArray().At(i-1,0).GetFloat();
                  _data->u[_data->iField] = _u;
               }
               else {
                  _data->addMatrix(v.first,nr,nc);
                  _M = _data->theMatrix[_data->iMatrix];
                  for (int i=1; i<=nr; ++i)
                     for (int j=1; j<=nc; ++j)
                        (*_M)(i,j) = w.GetArray().At(i-1,j-1).GetFloat();
               }
               break;
            }

         default:
            break;
      }
   }
}


int calc::getVar(string_type& s)
{
   parser->Eval();
   var_maptype vmap = parser->GetVar();
   if (vmap.size()==_vmap.size())
      return 1;
   for (const auto& v: vmap) {
      s = v.first;
      if (_vmap.find(s)==_vmap.end()) {
         _vmap[s] = v.second;
         return 0;
      }
      _vmap[s] = v.second;
   }
   return 0;
}


void calc::addVar()
{
   string s;
   if (getVar(s)==0) {
      switch (parser->Eval().GetType()) {

         case 'i': 
         case 'f': 
            {
//               _data->addParam(s,parser->Eval().GetFloat());
               break;
            }

         case 'm': 
            {
/*               int nr=parser->Eval().GetRows(), nc=parser->Eval().GetCols();
               if (nr==1 || nc==1) {
                  _u = new OFELI::Vect<double>(nr,nc);
                  _data->addField(s,std::max(nr,nc));
                  for (int i=1; i<=std::max(nr,nc); ++i)
                     (*_u)(i,1) = parser->Eval().GetArray().At(i-1,0).GetFloat();
                  _data->u[_data->iField] = _u;
               }
               else {
                  _data->addMatrix(s,nr,nc);
                  _M = _data->theMatrix[_data->iMatrix];
                  for (int i=1; i<=nr; ++i)
                     for (int j=1; j<=nc; ++j)
                        (*_M)(i,j) = parser->Eval().GetArray().At(i-1,j-1).GetFloat();
               }*/
               break;
            }

         default:
            break;
      }
   }
}


void calc::ListExprVar(const ParserXBase* p)
{
   string_type sExpr = p->GetExpr();
   if (sExpr.length() == 0) {
      cout << _T("Expression string is empty\n");
      return;
   }

// Query the used variables (must be done after calc)
   cout << _T("\nExpression variables:\n");
   cout << _T("---------------------\n");
   cout << _T("Expression: ") << p->GetExpr() << "\n";

   var_maptype variables = p->GetVar();
   if (!variables.size())
      cout << "Expression does not contain variables\n";
   else {
      cout << "Number: " << (int)variables.size() << "\n";
      for (const auto& v: variables)
         cout << "Name: " << v.first << ", Value: " << *v.second << endl;
   }
}


int calc::CheckKeywords(const string& ln, ParserXBase* p)
{
   if (ln=="end" || ln=="<") {
      setData(p);
      return -1;
   }
   else if (ln.substr(0,5)=="print") {
      _rita->msg("calc>","Command print is not allowed in this mode","");
      return 0;
   }
   else if (ln=="list") {
      ListConst(p);
      ListVar(p);
      ListExprVar(p);
      return 1;
   }
   return 0;
}


int calc::run()
{
   for (;;) {
      try {
         int verb = 1;
         int ret = _cmd->readline("rita>calc> ");
         string ln = _cmd->buffer();
         if (ln=="?" || ln=="help") {
            cout << "In this module variables can be defined by regular algebraic expressions.\n";
            cout << "Other available commands:\n";
            cout << "list: To list all defined variables in the module\n";
            cout << "end or <: To go back to main module\n";
            cout << "exit: To exit rita" << endl;
            continue;
         }
         if (ret==-3 || ret==-4)
            continue;
         *_rita->ofh << ln << endl;
         if (ln[ln.size()-1]==';') {
            verb = 0;
            ln.pop_back();
         }

         switch (CheckKeywords(ln,parser))
         {
            case  0: break;
            case  1: continue;
            case -1: return -1;
         }
         parser->SetExpr(ln);
//         addVar();
         if (verb)
            cout << std::setprecision(12) << "= " << parser->Eval() << endl;
      }
      catch(ParserError &e) {
         if (e.GetPos()!=-1) {
            string_type sMarker;
            sMarker.insert(0,e.GetPos()+10,' ');
            sMarker += _T("^\n");
            cout << sMarker;
         }
         _rita->msg("calc>",e.GetMsg()+" (Errc: "+to_string(e.GetCode())+")","",1);
      }
   }
   return 0;
}

} /* namespace RITA */
