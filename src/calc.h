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

class calc {

 public:

    calc(rita* r, cmd* command);
    ~calc();
    int run();
    int CheckKeywords(const std::string& ln, ParserXBase* p);
    int getVar(string_type& s);
    ParserX *parser;

 private:
 
    rita *_rita;
    cmd *_cmd;
    data *_data;
    var_maptype _vmap;
    OFELI::Vect<double> *_u;
    OFELI::Matrix<double> *_M;

    static Value Help();
    static void ListVar(const ParserXBase* p);
    static void ListConst(const ParserXBase* p);
    static void ListExprVar(const ParserXBase* p);
    void setData(ParserXBase* p);

// Operator callback functions
/*   static value_type Mega(value_type a_fVal) { return a_fVal * 1e6; }
   static value_type Milli(value_type a_fVal) { return a_fVal / (value_type)1e3; }
   static value_type Rnd(value_type v) { return v * std::rand() / (value_type)(RAND_MAX + 1.0); }
   static value_type Not(value_type v) { return v == 0; }
   static value_type Add(value_type v1, value_type v2) { return v1 + v2; }
   static value_type Mul(value_type v1, value_type v2) { return v1 * v2; }
   static value_type Arg2Of2(value_type v1, value_type v2) { return v2; }
   static value_type Arg1Of2(value_type v1, value_type v2) { return v1; }*/
   void addVar();

};

} /* namespace RITA */
