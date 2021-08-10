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

                        Definition of class 'stationary'

  ==============================================================================*/

#pragma once

#include <map>
#include "OFELI.h"
#include "rita.h"

namespace RITA {

class stationary
{

 public:

    stationary(rita *r, int eq);
    ~stationary();
    void setSave(vector<int>& isave, vector<int>& format, vector<string>& file);
    int run();

 private:

    rita *_rita;
    data *_data;
    int _nb_fields, _nb_eq, _rs, _eq;
    vector<int> *_fformat, *_isave;
    vector<string> *_save_file;
    equa *_pde_eq;
    odae *_ae_eq;
};

} /* namespace RITA */
