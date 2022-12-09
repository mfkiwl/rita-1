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

                            Definition of class 'calc'

  ==============================================================================*/

#pragma once

#include "../muparserx/mpParser.h"
#include "../muparserx/mpDefines.h"
#include "../muparserx/mpTest.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"

namespace RITA {

using namespace mup;

class rita;
class cmd;
class data;


class FunMatrix : public ICallback
{
 public:
    FunMatrix() : ICallback(cmFUNC, _T(""), 1) {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FunPrint : public ICallback
{
 public:
    FunPrint() : ICallback(cmFUNC, _T("print"), 1) {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FunListVar : public ICallback
{
public:

	FunListVar() : ICallback(cmFUNC, _T("list_var"), 0) {}
	virtual void Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/);
	virtual const char_type* GetDesc() const;
	virtual IToken* Clone() const;
};


class FunListConst : public ICallback
{
 public:

   FunListConst() : ICallback(cmFUNC, _T("list_const"), 0) {}
   virtual void Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/);
	virtual const char_type* GetDesc() const;
	virtual IToken* Clone() const;
};


class FunListFunctions : public ICallback
{
 public:
    FunListFunctions() : ICallback(cmFUNC, _T("list_fun"), 0) {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FunEnableOptimizer : public ICallback
{
 public:
    FunEnableOptimizer() : ICallback(cmFUNC, _T("enable_optimizer"), 1) {}
    void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FctMatrix : public ICallback
{
  public:
     FctMatrix() : ICallback(cmFUNC, _T("matrix"), -1) {}
     virtual ~FctMatrix() {}
     virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc) override;
     virtual const char_type* GetDesc() const override;
     virtual IToken* Clone() const override;
};


class FunSelfTest : public ICallback
{
 public:
    FunSelfTest() : ICallback(cmFUNC, _T("test"), 0) {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* /*a_pArg*/, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FunEnableDebugDump : public ICallback
{
 public:
    FunEnableDebugDump() : ICallback(cmFUNC, _T("debug"), 2) {}
    virtual void Eval(ptr_val_type& ret, const ptr_val_type* a_pArg, int /*a_iArgc*/);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FunGeneric : public ICallback
{
 public:
    FunGeneric(string_type sIdent, string_type sFunction);
    virtual ~FunGeneric() {}
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;

 private:
    ParserX m_parser;
    mup::var_maptype m_vars;
    val_vec_type m_val;
};


class FunDefine : public ICallback
{
 public:
    FunDefine() : ICallback(cmFUNC, _T("define"), 2) {}
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int a_iArgc);
    virtual const char_type* GetDesc() const;
    virtual IToken* Clone() const;
};


class FctVector : public ICallback
{
 public:
    FctVector();
    virtual ~FctVector();
    virtual void Eval(ptr_val_type &ret, const ptr_val_type *a_pArg, int argc) override;
    virtual const char_type* GetDesc() const override;
    virtual IToken* Clone() const override;
};


class calc {

 public:

    calc(rita* r, cmd* command);
    ~calc();
    int run();
    int CheckKeywords();
    int getVar(string_type& s);
    int setVector(const string& name, OFELI::Vect<double> *u);
    int setMatrix(const string& name, OFELI::Matrix<double> *M);

 private:
 
    rita *_rita;
    string _sLine;
    cmd *_cmd;
    ParserX _parser;
    var_maptype _vmap;
    vector<Value *> _theV;
    Value *_v;
    map<string,Value> _sv;

    Value Help();
    void ListVar();
    void ListConst();
    void ListExprVar();
    void setData();
    void parse();

    void addVar();

};

} /* namespace RITA */
