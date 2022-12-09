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


void FctMatrix::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc)
{
   if (argc < 1 || argc>2) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }
   int_type m = a_pArg[0]->GetInteger(),
   n = (argc == 1) ? m : a_pArg[1]->GetInteger();
   if (m==n && n==1)
      *ret = 0.0;
        else
      *ret = matrix_type(m, n, 0.0);
}


const char_type* FctMatrix::GetDesc() const
{
   return _T("matrix(x [, y]) - Returns a matrix whose elements are all 0.");
}


IToken* FctMatrix::Clone() const
{
   return new FctMatrix(*this);
}


void FunMatrix::Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/)
{
   cout << a_pArg[0].Get()->ToString() << "\n";
   *ret = (float_type)0.0;
}


const char_type* FunMatrix::GetDesc() const
{
   return _T("");
}


IToken* FunMatrix::Clone() const
{
   return new FunMatrix(*this);
}


void FunPrint::Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/)
{
   console() << a_pArg[0].Get()->ToString() << _T("\n");
   *ret = (float_type)0.0;
}


const char_type* FunPrint::GetDesc() const
{
   return _T("");
}


IToken* FunPrint::Clone() const
{
   return new FunPrint(*this);
}


void FunListVar::Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/)
{
   ParserXBase& parser = *GetParent();
   cout << "\nParser variables:\n-----------------\n";

// Query the used variables (must be done after calc)
   var_maptype vmap = parser.GetVar();
   if (!vmap.size())
      cout << "Expression does not contain variables\n";
   else {
      for (const auto& w: vmap) {
         cout << "  " << w.first << " =  " << *(w.second)/* << "\n"*/;

//       If you need more specific information cast the token to a variable object
         Variable& v = (Variable&)(*(w.second));
         cout << "  (type=\"" << v.GetType() << "\"; ptr=0x" << hex << v.GetPtr() << ")\n";
      }
   }
   *ret = (float_type)vmap.size();
}


const char_type* FunListVar::GetDesc() const
{
   return "list_var() - List all variables of the parser bound to this function and returns the number of defined variables.";
}


IToken* FunListVar::Clone() const
{
   return new FunListVar(*this);
}


void FunListConst::Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/)
{
   ParserXBase& parser = *GetParent();
   cout << _T("\nParser constants:\n-----------------\n");
   val_maptype cmap = parser.GetConst();
   if (!cmap.size())
      cout << "No constants defined\n";
   else {
      for (const auto& v: cmap)
         cout << "  " << v.first << " = " << (Value&)(*(v.second)) << endl;
   }
   *ret = (float_type)cmap.size();
}


const char_type* FunListConst::GetDesc() const
{
   return "list_const() - List all constants of the parser bound to this function and returns the number of defined constants.";
}


IToken* FunListConst::Clone() const
{
   return new FunListConst(*this);
}


void FunListFunctions::Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/)
{
   ParserXBase& parser = *GetParent();
   cout << "\nParser functions:\n----------------\n";
   fun_maptype fmap = parser.GetFunDef();
   if (!fmap.size())
      cout << "No functions defined\n";
   else {
      for (const auto& v: fmap) {
         ICallback* pFun = (ICallback*)v.second.Get();
         cout << pFun->GetDesc() << endl;
      }
   }
   *ret = (float_type)fmap.size();
}

 
const char_type* FunListFunctions::GetDesc() const
{
   return "list_fun() - List all parser functions and returns the total number of defined functions.";
}


IToken* FunListFunctions::Clone() const
{
   return new FunListFunctions(*this);
}


void FunEnableOptimizer::Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/)
{
   ParserXBase& parser = *GetParent();
   parser.EnableOptimizer(a_pArg[0]->GetBool());
   *ret = a_pArg[0]->GetBool();
}


const char_type* FunEnableOptimizer::GetDesc() const
{
   return _T("enable_optimizer(bool) - Enables the parsers built in expression optimizer.");
}


IToken* FunEnableOptimizer::Clone() const
{
   return new FunEnableOptimizer(*this);
}


void FunSelfTest::Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/)
{
   ParserXBase::EnableDebugDump(0, 0);
   ParserTester pt;
   pt.Run();
   *ret = (float_type)0.0;
}


const char_type* FunSelfTest::GetDesc() const
{
   return _T("test() - Runs the unit test of muparserx.");
}


IToken* FunSelfTest::Clone() const
{
   return new FunSelfTest(*this);
}


void FunEnableDebugDump::Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/)
{
   ParserXBase::EnableDebugDump(a_pArg[0]->GetBool(), a_pArg[1]->GetBool());
   if (a_pArg[0]->GetBool())
      cout << "Bytecode output activated.\n";
   else
      cout << "Bytecode output deactivated.\n";
   if (a_pArg[1]->GetBool())
      cout << "Stack output activated.\n";
   else
      cout << _T("Stack output deactivated.\n");
   *ret = (float_type)0.0;
}


const char_type* FunEnableDebugDump::GetDesc() const
{
   return "debug(bDumpRPN, bDumpStack) - Enable dumping of RPN and stack content.";
}


IToken* FunEnableDebugDump::Clone() const
{
   return new FunEnableDebugDump(*this);
}


FunGeneric::FunGeneric(string_type sIdent, string_type sFunction)
           : ICallback(cmFUNC, sIdent.c_str())
            ,m_parser()
            ,m_vars()
            ,m_val()
{
   m_parser.SetExpr(sFunction);
   m_vars = m_parser.GetExprVar();
   SetArgc(m_vars.size());

// Create values for the undefined variables and bind them
// to the variables
   for (const auto& v: m_vars) {
      ptr_val_type val(new Value());
      m_val.push_back(val);
      m_parser.DefineVar(v.second->GetIdent(), Variable(val.Get()));
   }
}


void FunGeneric::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc)
{
// Set the variables
   for (std::size_t i=0; i<(std::size_t)a_iArgc; ++i)
      *m_val[i] = *a_pArg[i];
   *ret = m_parser.Eval();
}


const char_type* FunGeneric::GetDesc() const
{
   return _T("xxx(...) - Dynamically defined function");
}


IToken* FunGeneric::Clone() const
{
   return new FunGeneric(*this);
}


void FunDefine::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc)
{
   string_type sFun = a_pArg[0]->GetString();
   string_type sDef = a_pArg[1]->GetString();
   ParserXBase &parser = *GetParent();
   parser.DefineFun(new FunGeneric(sFun, sDef));
}


const char_type* FunDefine::GetDesc() const
{
   return _T("define(Function, Definition) - Define a new parser function.");
}


IToken* FunDefine::Clone() const
{
   return new FunDefine(*this);
}


void calc::ListExprVar()
{
   cout << "\nVariables found in : \"" << _parser.GetExpr() << "\"\n";
   cout << "-----------------------------\n";
   var_maptype vmap = _parser.GetExprVar();
   if (!vmap.size())
      cout << "Expression does not contain variables" << endl;
	else {
      for (const auto& v: vmap)
         cout << "  " << v.first << " = " << (Variable&)(*(v.second)) << endl;
	}
}


void calc::ListVar()
{
   var_maptype variables = _parser.GetVar();
   if (!variables.size())
      return;
   cout << "\nParser variables:\n";
   cout << "-----------------\n";
   cout << "Number: " << variables.size() << "\n";
   for (const auto& v: variables)
      cout << "Name: " << v.first << " = " << *v.second << std::endl;
}


void calc::ListConst()
{
   cout << "\nParser constants:\n";
   cout << "-----------------\n";
   var_maptype cmap = _parser.GetConst();
   if (!cmap.size())
      cout << "Expression does not contain constants\n";
	else {
      for (const auto& v: cmap)
         cout << "  " << v.first << " = " << *v.second << std::endl;
   }
}


FctVector::FctVector() : ICallback(cmFUNC, _T("vector"), -1)
{ }


FctVector::~FctVector()
{ }


void FctVector::Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc)
{
   if (argc != 1) {
      ErrorContext err;
      err.Errc = ecINVALID_NUMBER_OF_PARAMETERS;
      err.Arg = argc;
      err.Ident = GetIdent();
      throw ParserError(err);
   }

   int_type n = a_pArg[0]->GetInteger();
   if (n==1)
      *ret = 0.0;
   else
      *ret = matrix_type(n, 1, 0.0);
}


const char_type* FctVector::GetDesc() const
{
   return _T("vector(x) - Returns a vector whose elements are all 0.");
}


IToken* FctVector::Clone() const
{
   return new FctVector(*this);
}


calc::calc(rita* r,
           cmd*  command)
     : _rita(r), _cmd(command)
{
   _parser.EnableAutoCreateVar(true);
   _parser.DefineFun(new FctMatrix);
   _parser.DefineFun(new FctVector);
}


calc::~calc()
{
}


Value calc::Help()
{
   cout << "-----------------------------------------------------------\n";
   cout << "Commands:\n\n";
   cout << "  list var     - list parser variables\n";
   cout << "  list exprvar - list expression variables\n";
   cout << "  list const   - list all numeric parser constants\n";
   cout << "  end (or <)   - exits the parser\n";
   cout << "-----------------------------------------------------------\n";
   return 0;
}


void calc::setData()
{
   var_maptype variables = _parser.GetVar();
   if (!variables.size())
      return;
   for (const auto& v: variables) {
      Variable &w = (Variable&)(*v.second);
      switch (w.GetType()) {

         case 'i': 
         case 'f': 
            {
               _rita->_data->addParam(v.first,w.GetFloat());
               break;
            }

         case 'm':
            {
               int nr=w.GetRows(), nc=w.GetCols();
               if (nr==1 || nc==1) {
                  OFELI::Vect<double> *u = _rita->_data->theVector[_rita->_data->iVector];
                  _rita->_data->addVector(v.first,0.,std::max(nr,nc),"",true);
                  for (int i=1; i<=std::max(nr,nc); ++i)
                     (*u)(i,1) = w.GetArray().At(i-1,0).GetFloat();
                  _rita->_data->theVector[_rita->_data->iVector] = u;
               }
               else {
                  _rita->_data->addMatrix(v.first,nr,nc,"","dense",true);
                  OFELI::Matrix<double> *M = _rita->_data->theMatrix[_rita->_data->iMatrix];
                  for (int i=1; i<=nr; ++i)
                     for (int j=1; j<=nc; ++j)
                        (*M)(i,j) = w.GetArray().At(i-1,j-1).GetFloat();
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
   _parser.Eval();
   var_maptype vmap = _parser.GetVar();
   if (vmap.size()==0) {
      _rita->msg("","Expression does not contain variables.");
      return -1;
   }
   if (vmap.size()==_vmap.size())
      return 1;
   for (const auto& v: vmap) {
      s = v.first;
      if (_rita->_data->checkName(s,DataType::PARAM,1)<0)
         return 1;
      if (_vmap.find(s)==_vmap.end()) {
         _vmap[s] = v.second;
         return 0;
      }
      _vmap[s] = v.second;
   }
   return 0;
}


int calc::setVector(const string&        name,
                    OFELI::Vect<double>* u)
{
   int n = u->size();
   _v = new Value(n,0.0);
   for (int i=0; i<n; ++i)
      _v->At(i) = (*u)[i];
   _parser.DefineVar(name,Variable(_v));
   return 0;
}


int calc::setMatrix(const string&          name,
                    OFELI::Matrix<double>* M)
{
   int nr=M->getNbRows(), nc=M->getNbColumns();
   _v = new Value(nr,nc,0.0);
   _theV.push_back(_v);
   _parser.DefineVar(name,Variable(_v));
   for (int i=1; i<=nr; ++i)
      for (int j=1; j<=nc; ++j)
         _v->At(i-1,j-1) = (*M)(i,j);
   return 0;
}


void calc::addVar()
{
   string s;
   if (getVar(s)==0) {
      switch (_parser.Eval().GetType()) {

         case 'i': 
         case 'f': 
            {
//               _rita->_data->addParam(s,parser.Eval().GetFloat());
               break;
            }

         case 'm': 
            {
               int nr=_parser.Eval().GetRows(), nc=_parser.Eval().GetCols();
               if (nr==1 || nc==1) {
                  OFELI::Vect<double> *u = _rita->_data->theVector[_rita->_data->iVector];
                  _rita->_data->addVector(s,0.,std::max(nr,nc));
                  for (int i=1; i<=std::max(nr,nc); ++i)
                     (*u)(i,1) = _parser.Eval().GetArray().At(i-1,0).GetFloat();
                  _rita->_data->theVector[_rita->_data->iVector] = u;
               }
               else {
                  _rita->_data->addMatrix(s,nr,nc);
                  OFELI::Matrix<double> *M = _rita->_data->theMatrix[_rita->_data->iMatrix];
                  for (int i=1; i<=nr; ++i)
                     for (int j=1; j<=nc; ++j)
                        (*M)(i,j) = _parser.Eval().GetArray().At(i-1,j-1).GetFloat();
               }
               break;
            }

         default:
            break;
      }
   }
}


int calc::CheckKeywords()
{
   if (_sLine.substr(0,5)=="print") {
      _rita->msg("","Command print is not allowed in this mode","");
      return 0;
   }
   else if (_sLine=="list") {
      ListConst();
      ListVar();
      ListExprVar();
      return 1;
   }
   return 0;
}


int calc::run()
{
   try {
      int verb=1, ret=0;
      _sLine = _cmd->buffer();
      *_rita->ofh << _sLine << endl;
      if (_sLine=="?" || _sLine=="help") {
         cout << "In this module variables can be defined by regular algebraic expressions.\n";
         cout << "Other available commands:\n";
         cout << "list: To list all defined variables in the module\n";
         cout << "end or <: To go back to main module\n";
         cout << "exit: To exit rita" << endl;
      }
      if (ret==-3 || ret==-4)
         return ret;
      if (_sLine[_sLine.size()-1]==';') {
         verb = 0;
         _sLine.pop_back();
      }

      switch (CheckKeywords())
      {
         case  0: break;
         case  1: break;
         case -1: return -1;
      }
      parse();
//      if (verb)
//         cout << std::setprecision(12) << "= " << _parser.Eval() << endl;
   }
   catch(ParserError &e) {
      if (e.GetPos()!=-1) {
         string_type sMarker;
         sMarker.insert(0,e.GetPos()+10,' ');
         sMarker += "^\n";
         cout << sMarker;
      }
      cout << "Error: " << e.GetMsg() << std::dec << endl;
      *_rita->ofl << "In " << sPrompt << e.GetMsg() << std::dec << endl;
   }
   return 0;
}


void calc::parse()
{
   _parser.SetExpr(_sLine);
   _parser.Eval();
   var_maptype var = _parser.GetVar();
   for (auto const& v: var) {

      Variable &w = (Variable&)(*(v.second));
      switch (w.GetType()) {

         case 'i': 
         case 'f': 
            _rita->_data->addParam(v.first,w.GetFloat());
            break;

         case 'v': 
            _rita->_data->print(v.first);
            break;

         case 'm':
            {
               int nr=w.GetRows(), nc=w.GetCols();
               if (nr==1 || nc==1) {
                  _rita->_data->addVector(v.first,0.,std::max(nr,nc),"",true);
                  OFELI::Vect<double> *u = _rita->_data->theVector[_rita->_data->iVector];
                  for (int i=1; i<=std::max(nr,nc); ++i)
                     (*u)(i) = w.GetArray().At(i-1,0).GetFloat();
               }
               else {
                  _rita->_data->addMatrix(v.first,nr,nc,"","dense",true);
                  OFELI::Matrix<double> *M = _rita->_data->theMatrix[_rita->_data->iMatrix];
                  for (int i=1; i<=nr; ++i)
                     for (int j=1; j<=nc; ++j)
                        (*M)(i,j) = w.GetArray().At(i-1,j-1).GetFloat();
               }
               break;
            }

         default:
            break;
      }
   }
}

} /* namespace RITA */
