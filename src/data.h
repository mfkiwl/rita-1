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

class data
{

 public:

    enum class eqType { AE, ODE, PDE, OPT, INTEGR };
    enum class DataType { PARAM, FIELD, MATRIX, GRID, MESH, TAB, FCT, AE, ODE, PDE };
    enum class DataSize { GIVEN_SIZE, GRID, NODES, ELEMENTS, SIDES, EDGES };
    enum class Storage { DENSE, SPARSE, SKYLINE, BAND, TRIDIAGONAL, DIAGONAL };
    struct Dat { int i; DataType dt; };

    data(rita *r, cmd *command, configure *config);
    ~data();
    int addField(const string& name, int n=1);
    int addMatrix(const string& name, int nr=1, int nc=1, Storage s=Storage::DENSE);
    int addMeshField(const string& name, DataSize s, int nb_dof=1);
    int addGridField(const string& name, int nb_dof=1);
    int addParam(const string& name, double value);
    int addAE(odae *ae, string name="");
    int addODE(odae *ode, string name="");
    int addPDE(equa *pde, string name="");
    int checkName(const string& name, const DataType& dt, int opt=0);
    int checkField(const string& name);
    int checkMatrix(const string& name);
    int checkFct(const string& name);
    int checkParam(const string& name, double& value);
    int checkParam(const string& name, int& value);
    int checkVect(const string& name);
    int checkMesh(const string& name);
    int checkGrid(const string& name);
    int run();
    void setVerbose(int verb) { _verb = verb; }
    void setSave(int s) { _sr = s; }
    int getVerbose() const { return _verb; }
    void set(cmd* command) { _cmd = command; }
    int ret() const { return _ret; }
    int addFunction(const string &def, const vector<string> &var, string name="");
    int addMesh(OFELI::Mesh* ms, const string& name);
    int getPar(int n, const string& msg, double& v);
    int getPar(int n, const string& msg, int& v);
    void print(const string &s);
    int getNbEq() const { return nb_ae+nb_ode+nb_pde; }

    int nb_fields, nb_fcts, nb_tabs, nb_meshes, nb_grids, nb_params, nb_vectors, nb_matrices;
    int nb_pde, nb_ode, nb_ae, nb_opt, nb_int, nb_eigen, nb_eq;
    vector<int> nb_dof;
    vector<OFELI::Vect<double> *> u;
    vector<double *> theParam;
    vector<string> Field, Param, Matr, AE, ODE, PDE;
    map<string,int> FieldName, ParamName, MatrixName, TabName, GridName, MeshName, FctName;
    map<string,int> AEName, ODEName, PDEName;
    map<string,Dat> dn;
    vector<eqType> FieldType;
    vector<DataSize> FieldSizeType;
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
    double obj, integral;
    bool ok;
    int iFct, iMesh, iGrid, iField, iMatrix, iTab, iParam, iAE, iODE, iPDE, iEq;
    vector<eqType> eq_type;
    vector<int> eqq;

 private:

    rita *_rita;
    int _size, _default_field;
    int _nb_args, _verb, _ret, _nb_dof, _nb, _sr;
    configure *_configure;
    cmd *_cmd;

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
