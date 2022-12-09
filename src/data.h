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
#include "OFELI.h"
#include "HVect.h"

namespace RITA {

class rita;
class configure;
class cmd;
class equa;

enum class DataType { NOTHING, PARAM, VECTOR, HVECTOR, MATRIX, GRID, MESH, TAB, FCT, AE, ODE, PDE };

struct odae {
   bool isSet, log, isFct;
   DataType type;
   vector<string> analytic, vars;
   vector<OFELI::Fct> theFct;
   OFELI::Vect<double> y, ph;
   OFELI::Vect<string> J;
   double init_time, final_time, time_step;
   OFELI::TimeScheme scheme;
   int size, vect, ind_fct, every;
   NonLinearIter nls;
   string fn, name, phase;
   odae();
   void setVars(int opt);
};

ostream& operator<<(ostream& s, const odae& e);


class data
{

 public:

    enum class eqType { AE, ODE, PDE, OPT, EIGEN, INTEGR };
    enum class DataSize { GIVEN_SIZE, GRID, NODES, ELEMENTS, SIDES, EDGES };
    enum class Storage { DENSE, SPARSE, SKYLINE, BAND, TRIDIAGONAL, DIAGONAL };
    struct Dat { int i; DataType dt; };

    data(rita *r, cmd *command, configure *config);
    ~data();
    int addVector(const string& name, double t=0., int n=1, string file="", bool opt=true);
    int remove(const string& name);
    int addMeshVector(const string& name, DataSize s, int nb_dof=1);
    int addGridVector(const string& name, int nb_dof=1);
    int addParam(const string& name, double value);
    int addAE(odae *ae, string name="");
    int addODE(odae *ode, string name="");
    int addPDE(equa *pde, string name="");
    int addTab(OFELI::Tabulation *tab, const string &name);
    void setTime(string s, double t);
    DataType getType(const string& s);
    int checkName(const string& name, const DataType& dt, int opt=0);
    int checkParam(const string& name, double& value);
    int checkParam(const string& name, int& value);
    int save();
    void setVerbose(int verb) { _verb = verb; }
    void setSave(int s) { _sr = s; }
    int getVerbose() const { return _verb; }
    void set(cmd* command) { _cmd = command; }
    int ret() const { return _ret; }
    int setDataExt(int key);
    int addFunction(const string &def, const vector<string> &var, string name="");
    int addMatrix(string name, int nr, int nc, string file="", string s="dense", bool opt=true);
    int addMesh(OFELI::Mesh* ms, const string& name);
    int getPar(int n, const string& msg, double& v);
    int getPar(int n, const string& msg, int& v);
    void print(const string &s);
    int getNbEq() const { return nb_ae+nb_ode+nb_pde; }
    void setTab2Grid(OFELI::Tabulation* tab);
    void setTab2Vector(OFELI::Tabulation* tab);
    int add2History(const string& s1, const string& s2, double t=0.0);
    int setHistory(const string& s1, const string& s2);
    void Summary();
  
  /*
   *  Each entity E:
   *  - iE is the index in the global array
   *  - theE[iE] is the pointer to entity number iE
   *  - NameE[iE] is the name of entity number iE
   *  - EName[name] is the index of the entity named name in the global array
   * 
   *  Entities are:
   *  Vector, Matrix, Mesh, Grid, Tab, Param, Fct, AE, ODE, PDE
   */
    int nb_vectors, nb_hvectors, nb_fcts, nb_tabs, nb_meshes, nb_grids, nb_params, nb_matrices;
    int nb_pde, nb_ode, nb_ae, nb_opt, nb_int, nb_eigen, nb_eq, save_freq;
    vector<int> nb_dof;
    vector<bool> aParam, aVector, aHVector, aMatrix, aGrid, aMesh, aTab, aFct, aAE, aODE, aPDE;
    vector<OFELI::Vect<double> *> theVector;
    vector<HVect *> theHVector;
    map<string,string> vect_hist;
    vector<double> VectorTime;
    vector<double *> theParam;
    vector<string> Vector, HVector, Param, Matr, AE, ODE, PDE;
    map<string,int> VectorName, HVectorName, ParamName, MatrixName, TabName, GridName, MeshName, FctName;
    map<string,string> Desc;
    map<string,int> AEName, ODEName, PDEName;
    map<string,Dat> dn;
    vector<eqType> VectorType;
    vector<DataSize> VectorSizeType;
    void setNodeBC(int code, string exp, double t, OFELI::Vect<double>& v);
    vector<int> VectorEquation;
    vector<OFELI::Grid *> theGrid;
    vector<OFELI::Mesh *> theMesh;
    vector<OFELI::Tabulation *> theTab;
    vector<OFELI::Fct *> theFct;
    vector<OFELI::Matrix<double> *> theMatrix;
    vector<odae *> theAE, theODE;
    vector<equa *> thePDE;
    vector<string> NameGrid, NameTab, NameMesh, NameFct, NameVector, NameParam, NameMatrix;
    double obj, integral;
    bool ok;
    int iFct, iMesh, iGrid, iVector, iHVector, iMatrix, iTab, iParam, iAE, iODE, iPDE, iEq;
    vector<eqType> eq_type;
    vector<int> eqq;
    int setDesc(const string& name, const string& desc);
    int setVector();
    int setMatrix();
    int setGrid();
    int setTab();
    int setFunction();
    int setDerivative();
    //    void Clear();
    void ListParams(int opt);
    void ListGrids(int opt);
    void ListMeshes(int opt);
    void ListVectors(int opt);
    void ListHVectors(int opt);
    void ListFunctions(int opt);
    void ListTabs(int opt);
    void ListMatrices(int opt);
    void ListAE(int opt);
    void ListODE(int opt);
    void ListPDE(int opt);

   int saveGmsh(const string &file, const Vect<double>& v);
   int saveVTK(const string &file, const Vect<double>& v);
   int saveTecplot(const string &file, const Vect<double>& v);


 private:

    rita *_rita;
    int _size, _default_vector;
    int _nb_args, _verb, _ret, _nb_dof, _nb, _sr;
    configure *_configure;
    cmd *_cmd;

    OFELI::Mesh *_theMesh;
    OFELI::Tabulation *_theTab;
    OFELI::Grid *_theGrid;
    OFELI::Fct *_theFct;
    OFELI::Vect<double> *_theVector;
    OFELI::Matrix<double> *_theMatrix;
    Vect<double> *_u;
    HVect *_hu;
    double *_theParam;
    void getHelp();
    int setNbDOF();
    map<string,OFELI::MatrixType> st_map = {{"dense",OFELI::DENSE},
                                            {"band",OFELI::BAND},
                                            {"skyline",OFELI::SKYLINE},
                                            {"sparse",OFELI::SPARSE}};

    map<string,int> ff = {{"gmsh",GMSH},
                          {"gnuplot",GNUPLOT},
                          {"vtk",VTK},
                          {"tecplot",TECPLOT},
                          {"matlab",MATLAB},
                          {"ofeli",OFELI_FF}};
    int setConfigure();
    vector<double> _xv;
};

} /* namespace RITA */
