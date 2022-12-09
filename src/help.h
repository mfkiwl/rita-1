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

                          Definition of class 'help'

  ==============================================================================*/

#pragma once

#include <string>
#include <iostream>

namespace RITA {

const string H0 = "rita is an interactive application for numerical analysis and simulation.\n"
                  "rita mainly uses the OFELI llibrary.\n"
                   "General purpose commands: help or ?, license, set, load, unload, end or <, exit or quit.\n"
                   "Data commands: grid, mesh, vector, matrix, tabulation, function, data, =\n"
                   "Modelling commands: approximation, integration, algebraic, ode, pde, stationary,\n"
                   "                    transient, optim, eigen, solve\n";
const string H1 =  "rita is an interactive application to solve main numerical analysis problems and in particular\n"
                   "scientific computing problems governed by partial differential equations.\n"
                   "The numerical solution is based on the OFELI library.\n"
                   "rita can be either run interactively by typing a series of commands or by using a shell script file.\n"
                   "The latter option can be chosen either by executing\n   'rita script-file'\n"
                   "or once in the rita execution (without argument) by typing\n   'load script-file'\n"
                   "The file script-file must contain all commands to be executed, one on each line.\n\n"
                   "rita execution consists in defining specific modelling problems to solve (partial differential equations,\n" 
                   "optimization, eigen problems, ...) and defining data these problems have to manipulate (vectors, matrices,\n"
                   "meshes, ...)\n\n"
                   "The file 'script-file' must contain all commands to be executed. Some characters in this file\n"
                   "(or in the input) have special meaning:\n"
                   "- A line starting with # is considered as a comment line\n"
                   "- If the first character is a ! the line is considered as a system command to be executed\n\n"
                   "Commands in rita are classified in 4 categories:\n\n"
                   "1. General purpose commands: \n"
                   "help or ?: Display this text or an abriged version\n"
                   "license: Display rita license information\n"
                   "set: Set configuration parameters\n"
                   "load: Read commands in a given script file\n"
                   "unload: Redirect input to console\n"
                   "end or <: Go back to a higher level of commands, if available\n"
                   "exit or quit: End execution\n\n"
                   "2. Calculator: rita can be considered as a console calculator\n\n"
                   "3. Commands that enable manipulating various data:\n"
                   "grid: Define a grid\n"
                   "mesh: Define a finite element mesh\n"
                   "vector: Define a vector\n"
                   "matrix: Define a matrix\n"
                   "tabulation: Define a tabulation: (xi,yi) array\n"
                   "function: Define analytically a function\n"
                   "data: List all data\n"
                   "=: Display specific data\n\n"
                   "4. Commands to define problem to solve:\n"
                   "approximation: Define and solvea data approximation problem\n"
                   "integration: Define and solve a numerical integration problem\n"
                   "algebraic: Define an algebraic equation (or set of equations)\n"
                   "ode: Define an ordinary differential equation (or set of equations)\n"
                   "pde: Define a partial differential equation (or set of equations)\n"
                   "stationary: Set problem to solve as stationary\n"
                   "transient: Set problem to solve as transien (time-dependent)\n"
                   "optim: Define an optimization problem\n"
                   "eigen: Define an eigenvalue problem\n"
                   "solve: Solve defined problem\n";

class help
{

 public:

    help()
    {
       _H0 = "--------------------------------------------------------------------------------------------------\n"
             + H0 +
             "--------------------------------------------------------------------------------------------------\n";
       _H1 = "--------------------------------------------------------------------------------------------------\n"
             + H1 +
             "--------------------------------------------------------------------------------------------------\n";
    }

    ~help() { }
    

    int run(int opt)
    {
       if (opt==0)
          std::cout << _H0 << std::endl;
       else if (opt==1)
          std::cout << _H1 << std::endl;
       return 0;
    }

 private:
   std::string _H0, _H1;

};

} /* namespace RITA */
