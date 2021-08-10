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
     : _cmd(command), _data(r->_data)
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
   if (getVar(s)==0)
      _data->addParam(s,parser->Eval().GetFloat());
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
   if (ln=="end" || ln=="<")
      return -1;
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
   int verb = 1;
   for (;;) {
      try {
         int ret = _cmd->readline("rita>calc> ");
         string ln = _cmd->buffer();
         if (ln=="?" || ln=="help") {
            cout << "In this module variables can be defined by regular algebraic expressions.\n";
            cout << "Other available commands:\n";
            cout << "list: To list all defined variables in the module\n";
            cout << "end or <: To go back to main module\n";
            cout << "exit: To exit rita" << endl;
         }
         if (ret==-3 || ret==-4)
            continue;
         if (ln=="end" || ln=="<")
            return 100;
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
         addVar();
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
         cout << e.GetMsg() << _T(" (Errc: ") << std::dec << e.GetCode() << _T(")") << _T("\n\n");
      }
   }
   return 0;
}

} /* namespace RITA */
