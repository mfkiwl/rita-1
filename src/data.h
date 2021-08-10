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

                        Definition of class 'data'

  ==============================================================================*/

#pragma once

#include <string>
#include <iostream>
using std::string;

#include "mesh/Mesh.h"
#include "linear_algebra/Vect.h"
#include "linear_algebra/Matrix.h"
#include "io/Fct.h"
#include "io/Tabulation.h"

namespace RITA {

class rita;
class configure;
class cmd;
class equa;
struct odae;

enum eqType {
   ALGEBRAIC_EQ,
   ODE_EQ,
   PDE_EQ,
   OPT,
   INTEGR
};

class data
{

 public:

    enum dataSize {
      GIVEN_SIZE,
      GRID,
      NODES,
      ELEMENTS,
      SIDES,
      EDGES
    };

    data(rita *r, cmd *command, configure *config);
    ~data();
    void addField(const string& name);
    void addField(const string& name, int n);
    int addMeshField(const string& name, dataSize s, int nb_dof=1);
    int addGridField(const string& name, int nb_dof=1);
    int addParam(const string& name, double value);
    int addAE(const string& name, odae *ae);
    int addODE(const string& name, odae *ode);
    int addPDE(const string& name, equa *pde);
    int checkField(const string& name) const;
    int checkFct(const string& name) const;
    int checkParam(const string& name, double& value);
    int checkParam(const string& name, int& value);
    int checkVect(const string& name) const;
    int checkMesh(const string& name) const;
    int checkGrid(const string& name) const;
    int run();
    void setVerbose(int verb) { _verb = verb; }
    void setSave(int s) { _sr = s; }
    int getVerbose() const { return _verb; }
    void set(cmd* command) { _cmd = command; }
    int ret() const { return _ret; }
    int addFunction(const string &name, const string &def, const vector<string> &var);
    int addMesh(OFELI::Mesh* ms, const string& name);
    int getPar(int n, const string& msg, double& v);
    int getPar(int n, const string& msg, int& v);
    void print(const string &s);

    int nb_fields, nb_fcts, nb_tabs, nb_meshes, nb_grids, nb_params, nb_vectors, nb_matrices;
    int nb_eq, nb_pde, nb_ode, nb_ae, nb_opt, nb_int, nb_eigen;
    vector<int> nb_dof;
    vector<OFELI::Vect<double> *> u;
    vector<double *> theParam;
    vector<string> Field, Param;
    map<string,int> FieldName, ParamName;
    vector<eqType> FieldType;
    vector<dataSize> FieldSizeType;
    void setNodeBC(int code, string exp, double t, OFELI::Vect<double>& v);
    vector<int> FieldEquation;
    vector<OFELI::Grid *> theGrid;
    vector<OFELI::Mesh *> theMesh;
    vector<OFELI::Tabulation *> theTab;
    vector<OFELI::Fct *> theFct;
    vector<OFELI::Matrix<double> *> theMatrix;
    vector<OFELI::Vect<double> *> theVector;
    vector<odae *> theAE, theODE;
    vector<equa *> thePDE;
    vector<string> grid_name, tab_name, mesh_name, fct_name, vector_name, param_name, matrix_name;
    vector<string> ae_name, ode_name, pde_name;
    double obj, integral;
    bool ok;
    int iFct, iMesh, iGrid, iField, iTab, iParam, iAE, iODE, iPDE, iEq;
    vector<int> eq_type;

 private:

    rita *_rita;
    int _size, _default_field;
    int _nb_args, _verb, _ret, _nb_dof, _nb, _sr;
    configure *_configure;
    cmd *_cmd;
    int _theMesh_alloc, _theTab_alloc, _theGrid_alloc, _theFct_alloc, _theVector_alloc;
    int _theMatrix_alloc, _u_alloc, _theParam_alloc;


    OFELI::Mesh *_theMesh;
    OFELI::Tabulation *_theTab;
    OFELI::Grid *_theGrid;
    OFELI::Fct *_theFct;
    OFELI::Vect<double> *_theVector;
    OFELI::Matrix<double> *_theMatrix;
    OFELI::Vect<double> *_u;
    double *_theParam;
    void getHelp();
    int setNbDOF();

    int setConfigure();
    int setParam();
    int setVector();
    int setMatrix();
    int setGrid();
    int setField();
    int setTab();
    int setFunction();
    int setDerivative();
    //    void Clear();
    void ListParams(int opt);
    void ListGrids(int opt);
    void ListMeshes(int opt);
    void ListVectors(int opt);
    void ListFields(int opt);
    void ListFunctions(int opt);
    void ListTabs(int opt);
    void ListMatrices(int opt);
    void ListAE(int opt);
    void ListODE(int opt);
    void ListPDE(int opt);
    void Summary();
    vector<double> _xv;
};

} /* namespace RITA */
