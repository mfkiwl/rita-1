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

                         Implementation of class 'data'

  ==============================================================================*/

#include "data.h"
#include "configure.h"
#include "rita.h"
#include "equa.h"
#include "linear_algebra/Matrix.h"
#include "io/IOField.h"
#include "calc.h"

using std::cout;
using std::endl;
using std::map;


namespace RITA {

data::data(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _size(0), _default_vector(0), _verb(1), _configure(config), _cmd(command),
       _theMesh(nullptr), _theTab(nullptr), _theGrid(nullptr), _theFct(nullptr), 
       _theVector(nullptr), _theMatrix(nullptr), _u(nullptr), _hu(nullptr), _theParam(nullptr)
{
   nb_vectors = nb_params = nb_fcts = nb_meshes = nb_tabs = nb_grids = nb_matrices = 0;
   nb_pde = nb_ode = nb_ae = nb_int = nb_eigen = nb_eq = 0;
   iMesh = iVector = iHVector = iMatrix = iGrid = iParam = iFct = iTab = iAE = iODE = iPDE = iEq = 0;
   theParam.push_back(nullptr), theVector.push_back(nullptr), VectorTime.push_back(0.), theFct.push_back(nullptr);
   theMatrix.push_back(nullptr), theMesh.push_back(nullptr), theGrid.push_back(nullptr), theHVector.push_back(nullptr);
   theAE.push_back(nullptr), theODE.push_back(nullptr), thePDE.push_back(nullptr), theTab.push_back(nullptr);
   aParam.push_back(true), aVector.push_back(true), aHVector.push_back(true), aMatrix.push_back(true), aGrid.push_back(true);
   aMesh.push_back(true), aTab.push_back(true), aFct.push_back(true), aAE.push_back(true), aODE.push_back(true), 
   aPDE.push_back(true), AE.push_back(" "), ODE.push_back(" "), PDE.push_back(" "), NameFct.push_back(" ");
   NameMesh.push_back(" "), NameGrid.push_back(" "), NameTab.push_back(" "), Param.push_back(" "), Vector.push_back(" ");
   Matr.push_back(" "), HVector.push_back(" ");
   VectorType.push_back(eqType::AE);
   VectorSizeType.push_back(DataSize::GIVEN_SIZE);
   VectorEquation.push_back(0);
   eq_type.push_back(eqType::AE);
   eqq.push_back(0);
   _ret = 0;
   ok = false;
}


data::~data()
{
}


int data::setDataExt(int key)
{
   int ret = 0;
   string td = "", fn="", name="", str1="", str2="";
   switch (key) {

      case 200:
         ret = setGrid();
         break;

      case 201:
         ret = _rita->setMesh();
         break;

      case 202:
         ret = setVector();
         break;

      case 203:
         ret = setTab();
         break;

      case 204:
         ret = setFunction();
         break;

      case 205:
         ret = setMatrix();
         break;

      case 206:
         ret = save();
         break;

      case 207:
         if (_cmd->setNbArg(1,"Data name to remove.")) {
            _rita->msg("","Data name to remove.","",1);
            break;
         }
         _cmd->get(td);
         ret = remove(td);
         break;

      case 208:
         _cmd->get(name);
         _cmd->get(td);
         setDesc(name,td);
         break;

      case 209:
         _cmd->get(str1);
         _cmd->get(str2);
         setHistory(str1,str2);
         break;

      case 210:
         Summary();
         break;

      case 211:
         if (_cmd->setNbArg(1,"Type of data to list.")) {
            _rita->msg("","Missing data type to list.\nAvailable data: parameters, grids, meshes, vectors, hvectors"
                          "functions, tabulations, matrices","",1);
            break;
         }
         _cmd->get(td);
         if (td.substr(0,5)=="param")
            ListParams(1);
         else if (td.substr(0,4)=="grid")
            ListGrids(1);
         else if (td.substr(0,4)=="mesh")
            ListMeshes(1);
         else if (td.substr(0,4)=="vect")
            ListVectors(1);
         else if (td.substr(0,5)=="hvect")
            ListHVectors(1);
         else if (td.substr(0,5)=="funct")
            ListFunctions(1);
         else if (td.substr(0,3)=="tab")
            ListTabs(1);
         else if (td.substr(0,5)=="matri")
            ListMatrices(1);
         else
            _rita->msg("list>","Unknown data type "+td);
         ret = 0;
         break;

      case 212:
      case 213:
         if (_cmd->setNbArg(1,"Data name to be given.",1)) {
            _rita->msg("","Missing data to display.","",1);
            break;
         }
         if (!_cmd->get(fn))
            print(fn);
         break;
   }
   return ret;
}


int data::checkName(const string&   name,
                    const DataType& dt,
                    int             opt)
{
   int ret = 0;
   if (dn.size()==0)
      return 0;
   DataType d = dn[name].dt;
   if (d==dt) {
      switch (d) {

         case DataType::PARAM:
            ret = ParamName[name];
            if (aParam[ret]==false)
               ret = 0;
            break;

         case DataType::VECTOR:
            ret = VectorName[name];
            if (aVector[ret]==false)
               ret = 0;
            break;

         case DataType::HVECTOR:
            ret = HVectorName[name];
            if (aHVector[ret]==false)
               ret = 0;
            break;

         case DataType::MATRIX:
            ret = MatrixName[name];
            if (aMatrix[ret]==false)
               ret = 0;
            break;

         case DataType::GRID:
            ret = GridName[name];
            if (aGrid[ret]==false)
               ret = 0;
            break;

         case DataType::MESH:
            ret = MeshName[name];
            if (aMesh[ret]==false)
               ret = 0;
            break;

         case DataType::TAB:
            ret = TabName[name];
            if (aTab[ret]==false)
               ret = 0;
            break;

         case DataType::FCT:
            ret = FctName[name];
            if (aFct[ret]==false)
               ret = 0;
            break;

         case DataType::AE:
            ret = AEName[name];
            if (aAE[ret]==false)
               ret = 0;
            break;

         case DataType::ODE:
            ret = ODEName[name];
            if (aODE[ret]==false)
               ret = 0;
            break;

         case DataType::PDE:
            ret = PDEName[name];
            if (aPDE[ret]==false)
               ret = 0;
            break;

         case DataType::NOTHING:
            ret = -1;
            break;
      }
      return ret;
   }

   if (opt==0)
      return 0;

   switch (d) {

      case DataType::PARAM:
         _rita->msg("data>","Name "+name+" already used for a parameter.");
         return -2;

      case DataType::VECTOR:
         _rita->msg("data>","Name "+name+" already used for a vector.");
         return -2;

      case DataType::HVECTOR:
         _rita->msg("data>","Name "+name+" already used for a history vector.");
         return -2;

      case DataType::MATRIX:
         _rita->msg("data>","Name "+name+" already used for a matrix.");
         return -2;

      case DataType::GRID:
         _rita->msg("data>","Name "+name+" already used for a grid.");
         return -2;

      case DataType::MESH:
         _rita->msg("data>","Name "+name+" already used for a mesh.");
         return -2;

      case DataType::TAB:
         _rita->msg("data>","Name "+name+" already used for a tabulation.");
         return -2;

      case DataType::FCT:
         _rita->msg("data>","Name "+name+" already used for a function.");
         return -2;

      case DataType::AE:
         _rita->msg("data>","Name "+name+" already used for an algebraic equation.");
         return -2;

      case DataType::ODE:
         _rita->msg("data>","Name "+name+" already used for an ordinary differential equation.");
         return -2;

      case DataType::PDE:
         _rita->msg("data>","Name "+name+" already used for a partial differential equation.");
         return -2;

      case DataType::NOTHING:
         return 0;
   }
   return 0;
}


int data::add2History(const string& s1, const string& s2, double t)
{
   int k = checkName(s1,DataType::HVECTOR);
   if (k==0) {
      _rita->msg("history>","History Vector "+s1+" not found.");
      return 1;
   }
   _hu = theHVector[k];
   if (s2!=_hu->getVectorName()) {
      _rita->msg("history>","History Vector "+s2+" not found in history vector "+s1+".");
      return 1;
   }
   int j = checkName(s2,DataType::VECTOR);
   theHVector[k]->set(*theVector[j],t);
   return 0;
}


int data::setHistory(const string& s1, const string& s2)
{
   int j = checkName(s1,DataType::VECTOR);
   if (j==0) {
      _rita->msg("history>","Vector "+s1+" undefined.");
      return 1;
   }
   int k = checkName(s2,DataType::HVECTOR);
   if (k==0) {
      _hu = new HVect;
      HVector.push_back(s2);
      iHVector = ++nb_hvectors;
      theHVector.push_back(_hu);
      aHVector.push_back(true);
      Desc[s2] = " ";
   }
   else {
      _hu = theHVector[k];
      iHVector = k;
   }
   HVectorName[s2] = nb_hvectors;
   dn[s2].i = nb_hvectors;
   dn[s2].dt = DataType::HVECTOR;
   theHVector[iHVector]->setVectorName(s1);
   vect_hist[s1] = s2;
   return 0;
}


void data::setTab2Grid(OFELI::Tabulation* tab)
{
   double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;
   int nb = tab->getNbVar(1);
   OFELI::fct &f = tab->Funct[0];
   int nx = f.Np[0], ny = f.Np[1], nz = f.Np[2];
   switch (nb)
   {
      case 1:
         xmin = tab->getMinVar(1,1);
         xmax = tab->getMaxVar(1,1);
         _theGrid = new OFELI::Grid(xmin,xmax,nx);
         break;   

      case 2:
         xmin = tab->getMinVar(1,1);
         xmax = tab->getMaxVar(1,1);
         ymin = tab->getMinVar(1,2);
         ymax = tab->getMaxVar(1,2);
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
         break;   

      case 3:
         xmin = tab->getMinVar(1,1);
         xmax = tab->getMaxVar(1,1);
         ymin = tab->getMinVar(1,2);
         ymax = tab->getMaxVar(1,2);
         zmin = tab->getMinVar(1,3);
         zmax = tab->getMaxVar(1,3);
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
         break;   
   }
   iGrid = ++nb_grids;
   string name = "G-"+to_string(iGrid);
   NameGrid.push_back(name);
   theGrid.push_back(_theGrid);
   aGrid.push_back(true);
   GridName[name] = dn[name].i = nb_grids;
   dn[name].dt = DataType::GRID;
}


void data::setTab2Vector(OFELI::Tabulation* tab)
{
   int nb = tab->getNbVar(1);
   OFELI::fct &f = tab->Funct[0];
   int nx = f.Np[0], ny = f.Np[1], nz = f.Np[2];
   addGridVector(tab->getFunctName(1),1);
   OFELI::Vect<double> *v = theVector[iVector];

   switch (nb) {
      case 1:
         v->setSize(nx);
         for (int i=1; i<=nx; ++i)
            (*v)(i) = f.Val(i);
         break;   

      case 2:
         v->setSize(nx,ny);
         for (int i=1; i<=nx; ++i)
            for (int j=1; j<=ny; ++j)
               (*v)(i,j) = f.Val(i,j);
         break;   

      case 3:
         v->setSize(nx,ny,nz);
         for (int i=1; i<=nx; ++i)
            for (int j=1; j<=ny; ++j)
               for (int k=1; k<=nz; ++k)
                  (*v)(i,j,k) = f.Val(i,j,k);
         break;   
   }
}


int data::checkParam(const string& name,
                     double&       value)
{
   int r = checkName(name,DataType::PARAM,1);
   if (r==-1)
      return r;
   if (r>=0)
      value = *theParam[r];
   else
      r = -1;
   return r;
}


int data::checkParam(const string& name,
                     int&          value)
{
   int r = checkName(name,DataType::PARAM,1);
   if (r==-1)
      return r;
   if (r>=0)
      value = *theParam[r];
   else
      r = -1;
   return r;
}


int data::addAE(odae*  ae,
                string name)
{
   if (name=="")
      name = "AE-" + to_string(nb_pde+1);
   int i = AEName[name];
   if (i==0 || aAE[i]==false) {
      iAE = ++nb_ae;
      iEq = ++nb_eq;
      eq_type.push_back(eqType::AE);
      eqq.push_back(iEq);
      AE.push_back(name);
      theAE.push_back(ae);
      aAE.push_back(true);
      Desc[name] = " ";
   }
   else {
      iAE = i;
      theAE[iAE] = ae;
   }
   AE[iAE] = name;
   AEName[name] = dn[name].i = nb_ae;
   dn[name].dt = DataType::AE;
   ae->name = name;
   return iAE;
}


int data::addODE(odae*  ode,
                 string name)
{
   if (name=="")
      name = "ODE-" + to_string(nb_ode+1);
   int i = ODEName[name];
   if (i==0 || aODE[i]==false) {
      iODE = ++nb_ode;
      iEq = ++nb_eq;
      eq_type.push_back(eqType::ODE);
      eqq.push_back(iEq);
      ODE.push_back(name);
      theODE.push_back(ode);
      aODE.push_back(true);
      Desc[name] = " ";
   }
   else {
      iODE = i;
      theODE[iODE] = ode;
   }
   ODE[iODE] = name;
   ODEName[name] = dn[name].i = nb_ode;
   dn[name].dt = DataType::ODE;
   ode->name = name;
   return iODE;
}


int data::addPDE(equa*  pde,
                 string name)
{
   if (name=="")
      name = "PDE-" + to_string(nb_pde+1);
   int i = PDEName[name];
   if (PDEName[name]==0 || aPDE[i]==false) {
      iPDE = ++nb_pde;
      iEq = ++nb_eq;
      eq_type.push_back(eqType::PDE);
      eqq.push_back(iEq);
      PDE.push_back(name);
      thePDE.push_back(pde);
      aPDE.push_back(true);
      Desc[name] = " ";
   }
   else {
      iPDE = i;
      thePDE[iPDE] = pde;
   }
   PDEName[name] = dn[name].i = nb_pde;
   dn[name].dt = DataType::PDE;
   pde->name = name;
   return iPDE;
}


int data::addFunction(const string&         def,
                      const vector<string>& var,
                      string                name)
{
   if (name=="")
      name = "F-" + to_string(nb_fcts+1);
   int i = FctName[name];
   if (i==0 || aFct[i]==false) {
      _theFct = new OFELI::Fct(name,def,var);
      iFct = ++nb_fcts;
      NameFct.push_back(name);
      theFct.push_back(_theFct);
      aFct.push_back(true);
      Desc[name] = " ";
   }
   else {
      _rita->msg("data>","Function "+name+" redefined.");
      iFct = i;
      theFct[iFct] = _theFct;
   }
   NameFct[iFct] = name;
   FctName[name] = dn[name].i = nb_fcts;
   dn[name].dt = DataType::FCT;
   return iFct;
}


int data::addMesh(OFELI::Mesh*  ms,
                  const string& name)
{
   int i = MeshName[name];
   if (i==0 || aMesh[i]==false) {
      iMesh = ++nb_meshes;
      theMesh.push_back(ms);
      NameMesh.push_back(name);
      aMesh.push_back(true);
   }
   else {
      iMesh = i;
      theMesh[iMesh] = ms;
   }
   MeshName[name] = dn[name].i = nb_meshes;
   dn[name].dt = DataType::MESH;
   return iMesh;
}


int data::addParam(const string& name,
                   double        value)
{
   int i = ParamName[name];
   if (i==0 || aParam[i]==false) {
      _theParam = new double;
      theParam.push_back(_theParam);
      aParam.push_back(true);
      Param.push_back(name);
      iParam = ++nb_params;
      dn[name].dt = DataType::PARAM;
      ParamName[name] = dn[name].i = nb_params;
      Desc[name] = " ";
   }
   else {
      iParam = i;
      _theParam = theParam[iParam];
   }
   *_theParam = value;
   return iParam;
}


int data::addVector(const string& name,
                    double        t,
                    int           n,
                    string        file,
                    bool          opt)
{
   int i = VectorName[name];
   if (i==0 || aVector[i]==false) {
      if (n==0)
         _u = new OFELI::Vect<double>;
      else
         _u = new OFELI::Vect<double>(n);
      if (file!="") {
         OFELI::XMLParser xml(file,OFELI::XMLParser::FIELD);
         xml.get(*_u);
      }
      Vector.push_back(name);
      iVector = ++nb_vectors;
      theVector.push_back(_u);
      VectorTime.push_back(t);
      aVector.push_back(true);
      VectorSizeType.push_back(DataSize::GIVEN_SIZE);
      Desc[name] = " ";
      vect_hist[name] = "%$§&";
   }
   else {
      if (file!="") {
         _u = theVector[i];
         OFELI::XMLParser xml(file,OFELI::XMLParser::FIELD);
         xml.get(*_u);
      }
      iVector = i;
      *theVector[iVector] = *_u;
   }
   VectorName[name] = nb_vectors;
   dn[name].i = nb_vectors;
   dn[name].dt = DataType::VECTOR;
   if (!opt)
      _rita->_calc->setVector(name,_u);
   return iVector;
}


void data::setTime(string s,
                   double t)
{
   int it = 0;
   if ((it=checkName(s,DataType::VECTOR))) {
      VectorTime[it] = t;
   }
}


DataType data::getType(const string& s)
{
   if (checkName(s,DataType::PARAM))
      return DataType::PARAM;
   else if (checkName(s,DataType::VECTOR))
      return DataType::VECTOR;
   else if (checkName(s,DataType::MATRIX))
      return DataType::MATRIX;
   else if (checkName(s,DataType::GRID))
      return DataType::GRID;
   else if (checkName(s,DataType::MESH))
      return DataType::MESH;
   else if (checkName(s,DataType::FCT))
      return DataType::FCT;
   else if (checkName(s,DataType::TAB))
      return DataType::TAB;
   else
      return DataType::NOTHING;
}


int data::remove(const string& name)
{
   DataType d = getType(name);
   int i = 0;
   if (d==DataType::PARAM) {
      i = ParamName[name];
      if (i>0) {
         aParam[i] = false;
         nb_params--;
      }
   }
   else if (d==DataType::VECTOR) {
      i = VectorName[name];
      if (i>0) {
         aVector[i] = false;
         theVector[i]->clear();
         nb_vectors--;
      }
   }
   else if (d==DataType::MATRIX) {
      i = MatrixName[name];
      if (i>0) {
         theMatrix[i]->setSize(0);
         aMatrix[i] = false;
         nb_matrices--;
      }
   }
   else if (d==DataType::GRID) {
      i = GridName[name];
      if (i>0) {
         aGrid[i] = false;
         nb_grids--;
      }
   }
   else if (d==DataType::MESH) {
      i = MeshName[name];
      if (i>0) {
         aMesh[i] = false;
         nb_meshes--;
      }
   }
   else if (d==DataType::TAB) {
      i = TabName[name];
      if (i>0) {
         aTab[i] = false;
         nb_tabs--;
      }
   }
   else if (d==DataType::FCT) {
      i = FctName[name];
      if (i>0) {
         aFct[i] = false;
         nb_fcts--;
      }
   }
   else if (d==DataType::AE) {
      i = AEName[name];
      if (i>0) {
         aAE[i] = false;
         nb_ae--;
      }
   }
   else if (d==DataType::ODE) {
      i = ODEName[name];
      if (i>0) {
         aODE[i] = false;
         nb_ode--;
      }
   }
   else if (d==DataType::PDE) {
      i = PDEName[name];
      if (i>0) {
         aPDE[i] = false;
         nb_pde--;
      }
   }
   else {
      _rita->msg("","No data named "+name+" found");
      return 1;
   }
   return 0;
}


int data::addMatrix(string name,
                    int    nr,
                    int    nc,
                    string file,
                    string s,
                    bool   opt)
{
   if (checkName(name,DataType::MATRIX,1)<0)
      return 1;
   if (file!="")
      nr = 1, nc = 1;
   _theMatrix = new DMatrix<double>(nr,nc);
   if (file!="") {
      ifstream f(file.c_str());
      if (!f.good()) {
         _rita->msg("","File "+file+" not found");
         return 1;
      }
      OFELI::XMLParser xml(file,OFELI::XMLParser::MATRIX);
      xml.get(_theMatrix);
   }
   NameMatrix.push_back(name);
   if (MatrixName[name]==0) {
      Matr.push_back(name);
      iMatrix = ++nb_matrices;
      theMatrix.push_back(_theMatrix);
      aMatrix.push_back(true);
      Desc[name] = " ";
   }
   else {
      iMatrix = MatrixName[name];
      theMatrix[iMatrix] = _theMatrix;
   }
   MatrixName[name] = dn[name].i = nb_matrices;
   dn[name].dt = DataType::MATRIX;
   if (!opt)
      _rita->_calc->setMatrix(name,_theMatrix);
   return iMatrix;
}


int data::addTab(OFELI::Tabulation* tab,
                 const string&      name)
{
   if (checkName(name,DataType::TAB,1)<0)
      return 1;
   if (TabName[name]==0) {
      NameTab.push_back(name);
      iTab = ++nb_tabs;
      theTab.push_back(tab);
      aTab.push_back(true);
      Desc[name] = " ";
   }
   else {
      iMatrix = TabName[name];
      theTab[iTab] = tab;
   }
   TabName[name] = dn[name].i = nb_tabs;
   dn[name].dt = DataType::TAB;
   return iTab;
}


int data::addMeshVector(const string& name,
                        DataSize      s,
                        int           nb_dof)
{
   _theMesh = theMesh[iMesh];
   if (_theMesh==nullptr) {
      _rita->msg("data>","No mesh data available.");
      _ret = -1;
      return _ret;
   }
   _u = new OFELI::Vect<double>;
   _nb_dof = nb_dof;
   if (VectorName[name]==0) {
      iVector = ++nb_vectors;
      Vector.push_back(name);
      theVector.push_back(_u);
      aVector.push_back(true);
      VectorType.push_back(eqType::PDE);
      Desc[name] = " ";
   }
   else {
      iVector = VectorName[name];
      theVector[iVector] = _u;
      VectorType[iVector] = eqType::PDE;
   }
   VectorName[name] = dn[name].i = nb_vectors;
   dn[name].dt = DataType::VECTOR;

   if (s==DataSize::NODES) {
      if (_theMesh->getNbNodes()==0) {
         _rita->msg("","Mesh has no nodes");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==DataSize::ELEMENTS) {
      if (_theMesh->getNbElements()==0) {
         _rita->msg("","Mesh has no elements.");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==DataSize::SIDES) {
      if (_theMesh->getNbSides()==0) {
         _rita->msg("","Mesh has no sides");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==DataSize::EDGES) {
      if (_theMesh->getNbEdges()==0) {
         _rita->msg("","Mesh has no edges");
         _ret = 1;
         return _ret;
      }
   }
   _u->setName(name);
   if (s==DataSize::NODES)
      _u->setMesh(*_theMesh,NODE_DOF,_nb_dof);
   else if (s==DataSize::ELEMENTS)
      _u->setMesh(*_theMesh,ELEMENT_DOF,_nb_dof);
   else if (s==DataSize::SIDES)
      _u->setMesh(*_theMesh,SIDE_DOF,_nb_dof);
   else if (s==DataSize::EDGES)
      _u->setMesh(*_theMesh,EDGE_DOF,_nb_dof);
   return iVector;
}


int data::addGridVector(const string& name,
                        int           nb_dof)
{
   if (checkName(name,DataType::GRID,1)<0)
      return 1;
   _theGrid = theGrid[iGrid];
   if (_theGrid==nullptr) {
      _rita->msg("","No grid data available.");
      _ret = -1;
      return _ret;
   }
   _u = new OFELI::Vect<double>(*_theGrid);
   _nb_dof = nb_dof;
   _theGrid->setNbDOF(_nb_dof);
   if (VectorName[name]==0) {
      iVector = ++nb_vectors;
      Vector.push_back(name);
      theVector.push_back(_u);
      aVector.push_back(true);
      VectorType.push_back(eqType::PDE);
      Desc[name] = " ";
   }
   else {
      iVector = VectorName[name];
      theVector[iVector] = _u;
      VectorType[iVector] = eqType::PDE;
   }
   VectorName[name] = dn[name].i = nb_vectors;
   dn[name].dt = DataType::VECTOR;
   return iVector;
}


void data::getHelp()
{
   cout << "\nAvailable Commands in mode 'data':\n";
   cout << "grid:       Define a grid\n";
   cout << "mesh:       Define a finite element mesh (can be accessed from the top level)\n";
   cout << "vector:     Define a vector\n";
   cout << "matrix:     Define a matrix\n";
   cout << "tabulation: Define a tabulation\n";
   cout << "function:   Define a function\n";
   cout << "list:       List a specific entity\n";
   cout << "Global commands: \n";
   cout << "help or ?:  Display this help\n";
   cout << "set:        Set configuration data\n";
   cout << "print:      print a specific entity\n";
   cout << "save:       save a specific entity in file\n";
   cout << "data:       Summary of defined entities\n\n";
   cout << "end or <:   Back to higher level\n";
   cout << "exit:       Terminate execution\n" << endl;
}


void data::print(const string& s)
{
   int k = checkName(s,DataType::PARAM);
   if (k>0) {
      cout << "Parameter: " << s << " = " << *theParam[k] << endl;
      return;
   }
   k = checkName(s,DataType::VECTOR);
   if (k>0) {
      cout << *theVector[k];
      return;
   }
   k = checkName(s,DataType::HVECTOR);
   if (k>0) {
      cout << *theHVector[k];
      return;
   }
   k = checkName(s,DataType::MATRIX);
   if (k>0) {
      cout << "Matrix " << s << endl;
      for (int i=1; i<=int(theMatrix[k]->getNbRows()); ++i) {
         cout << "Row " << i << ": ";
         for (int j=1; j<=int(theMatrix[k]->getNbColumns()); ++j)
            cout << (*theMatrix[k])(i,j) << "  ";
         cout << endl;
      }
      return;
   }
   k = checkName(s,DataType::MESH);
   if (k>0) {
      cout << *theMesh[k];
      return;
   }
   k = checkName(s,DataType::GRID);
   if (k>0) {
      cout << *theGrid[k];
      return;
   }
   k = checkName(s,DataType::FCT);
   if (k>0) {
      cout << *theFct[k];
      return;
   }
   k =  checkName(s,DataType::TAB);
   if (k>0) {
      cout << *theTab[k];
      return;
   }
   k =  checkName(s,DataType::AE);
   if (k>0) {
      cout << *theAE[k];
      return;
   }
   k =  checkName(s,DataType::ODE);
   if (k>0) {
      cout << *theODE[k];
      return;
   }
   k =  checkName(s,DataType::PDE);
   if (k>0) {
      cout << *thePDE[k];
      return;
   }
   _rita->msg("","Entity/Data '"+s+"' not found");
}


int data::setConfigure()
{
   _configure->setVerbose(_verb);
   int ret = _configure->run();
   _verb = _configure->getVerbose();
   return ret;
}


int data::setGrid()
{
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, d1=-1, d2=-1, d3=-1, dim=0, nx=10, ny=10, nz=10;
   string name="G-"+to_string(nb_grids+1);
   static const vector<string> kw {"name","min","max","ne"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>grid>","Error in command.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
            name = _cmd->string_token(0);
            break;

         case   1:
            if (nb==0 || nb>3) {
               _rita->msg("data>grid>","Illegal number of arguments");
               return 1;
            }
            dim = d1 = nb;
            xmin = _cmd->double_token(0);
            if (d1>1)
               ymin = _cmd->double_token(1);
            if (d1>2)
               zmin = _cmd->double_token(2);
            break;

         case   2:
            if (nb==0 || nb>3) {
               _rita->msg("data>grid>","Illegal number of arguments");
               return 1;
            }
            dim = d2 = nb;
            xmax = _cmd->double_token(0);
            if (d2>1)
               ymax = _cmd->double_token(1);
            if (d2>2)
               zmax = _cmd->double_token(2);
            break;

         case   3:
            if (nb==0 || nb>3) {
               _rita->msg("data>grid>","Illegal number of arguments");
               return 1;
            }
            dim = d3 = nb;
            nx = _cmd->int_token(0);
            if (d3>1)
               ny = _cmd->int_token(1);
            if (d3>2)
               nz = _cmd->int_token(2);
            break;

         default:
            _rita->msg("data>grid>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if ((d1!=d2 || d1!=d3 || d2!=d3) && d1>=0 && d2>=0 && d3>=0) {
      _rita->msg("data>grid>","Dimensions do not match.");
      return 1;
   }
   if (xmin>=xmax || ymin>=ymax || zmin>=zmax) {
      _rita->msg("data>grid>","Domain definition is incorrect.");
      return 1;
   }
   if (nx<1 || ny<1 || nz<1) {
      _rita->msg("data>grid>","Number of subdivisions is incorrect.");
      return 1;
   }
   *_rita->ofh << "  grid  name=" << name;
   if (dim==1) {
      _theGrid = new OFELI::Grid(xmin,xmax,nx);
      *_rita->ofh << "  min=" << xmin << "  max=" << xmax << "  ne=" << nx;
   }
   else if (dim==2) {
      _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
      *_rita->ofh << "  min=" << xmin << "," << ymin << "  max=" << xmax << "," << ymax
                  << "  ne=" << nx << "," << ny;
   }
   else if (dim==3) {
      _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
      *_rita->ofh << "  min=" << xmin << "," << ymin << "," << zmin << "  max=" << xmax
                  << "," << ymax << "," << zmax << "  ne=" << nx << "," << ny << "," << nz;
   }
   *_rita->ofh << endl;
   theGrid.push_back(_theGrid);
   aGrid.push_back(true);
   NameGrid.push_back(name);
   iGrid++;
   nb_grids++;
   return 0;
}


int data::setVector()
{
   int nb=0, size=0;
   string name="vect-"+to_string(nb_vectors+1);
   static const vector<string> kw {"name","size","def$ine","set"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>vector>","No argument given\nAvailable arguments: name, size, define, set.","");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            size = _cmd->int_token(0);
            break;

         case 2:
            break;

         case 3:
            break;

         default:
            _rita->msg("data>vector>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (nb_args>0) {
      addVector(name,0.,size);
      *_rita->ofh << "  vector  size=" << size << "  name=" << name << endl;
   }
   return 0;
}


int data::setMatrix()
{
   string name="mat-"+to_string(nb_matrices+1);
   string file="", storage="dense";
   int nr=0, nc=0;
   static const vector<string> kw {"name","file","storage","nr","nc","def$ine","set"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>matrix>","No agument given.\nAvailable arguments:"
                                " name, file, storage, nr, nc, define, set.","");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg("=");
      switch (n) {

         case   0:
            name = _cmd->string_token();
            break;

         case   1:
            file = _cmd->string_token();
            break;

         case   2:
            storage = _cmd->string_token();
            break;

         case   3:
            nr = _cmd->int_token();
            if (nc==0)
               nc = nr;
            break;

         case   4:
            nc = _cmd->int_token();
            if (nr==0)
               nr = nc;
            break;

         default:
            _rita->msg("data>matrix>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (file=="" && nr==0) {
      _rita->msg("data>matrix>","Matrix must be given either by its size or read from a file.");
      return 1;
   }
   if (storage!="dense") {
      _rita->msg("data>matrix>","Only dense storage is allowed in the current release.");
      return 1;
   }
   addMatrix(name,nr,nc,file,storage,false);
   return 0;
}


int data::setFunction()
{
   _ret = 0;
   bool var_ok=false, def_ok=false, name_ok=false;
   string vv="", def="", name="f";
   int nb=1;
   vector<string> var;
   static const vector<string> kw {"name","var$iable","vect$or","nb","def","definition"};

   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<0) {
      _rita->msg("data>function>:","No argument given.\nAvailable arguments: name, variable, vector, nb, definition.","");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {

      switch(_cmd->getArg("=")) {

         case 0:
            name = _cmd->string_token();
            name_ok = true;
            break;

         case 1:
         case 2:
            vv = _cmd->string_token();
            var_ok = true;
            break;

         case 3:
            nb = _cmd->int_token();
            break;

         case 4:
         case 5:
            def = _cmd->string_token();
            def_ok = true;
            break;

         default:
            _rita->msg("data>function>","Unknown argument: "+_cmd->Arg());
            return 1;
      }
      if (!name_ok) {
         _rita->msg("data>function>","Missing function name.");
         return 1;
      }
   }

   if (nb_args==0) {
      _rita->msg("data>function>","No command argument given.");
      return 1;
   }
   else {
      if (!var_ok) {
         _rita->msg("data>function>","No variable defined.");
         return 1;
      }
      if (!def_ok) {
         _rita->msg("data>function>","No function definition given.");
         return 1;
      }
      var.clear();
      if (nb==1)
         var.push_back(vv);
      else {
         for (int i=0; i<nb; ++i)
            var.push_back(vv+to_string(i+1));
      }
      addFunction(def,var,name);
      if (name_ok)
         *_rita->ofh << "  function  name=" << name;
      for (const auto& v: var)
         *_rita->ofh << " var=" << v;
      *_rita->ofh << "  def=" << def << endl;
   }
   return 0;
}


int data::setNbDOF()
{
   if (_cmd->setNbArg(1,"Give number of degrees of freedom.")) {
      _rita->msg("data>>nbdof>","Missing value of nbdof.");
      return 1;
   }
   _ret = _cmd->get(_nb_dof);
   if (!_ret)
      *_rita->ofh << "  nbdof " << _nb_dof;
   return _ret;
}


int data::setTab()
{
   int dim1=0, dim2=0, dim3=0;
   string file="";
   string name = "T" + to_string(theTab.size()+1);
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, nx=10, ny=10, nz=10, grid_ok=0, file_ok=0, vector_ok=0;
   static const vector<string> kw {"name","file","min","max","ne","vector"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<0) {
      _rita->msg("data>tabulation>","Error in command.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            file = _cmd->string_token(0);
            file_ok = 1;
            break;

         case 2:
            dim1 = nb;
            xmin = _cmd->double_token(0);
            if (nb>1)
               ymin = _cmd->double_token(1);
            if (nb>2)
               zmin = _cmd->double_token(2);
            grid_ok += 1;
            break;

         case 3:
            dim2 = nb;
            xmax = _cmd->double_token(0);
            if (nb>1)
               ymax = _cmd->double_token(1);
            if (nb>2)
               zmax = _cmd->double_token(2);
            grid_ok += 10;
            break;

         case 4:
            dim3 = nb;
            nx = _cmd->int_token(0);
            if (nb>1)
               ny = _cmd->int_token(1);
            if (nb>2)
               nz = _cmd->int_token(2);
            grid_ok += 100;
            break;

         case 5:
	   //            fd = _cmd->string_token();
            vector_ok = true;
            break;

         default:
            _rita->msg("data>tabulation>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (!file_ok && grid_ok<111) {
      _rita->msg("data>tabulation>","No grid data given.");
      return 1;
   }
   if (!vector_ok && !file_ok) {
      _rita->msg("data>tabulation>","No associated vector given.");
      return 1;
   }
   if (!file_ok && (dim1!=dim2 || dim1!=dim3 || dim2!=dim3)) {
      _rita->msg("data>tabulation>","Incompatible space dimensions as given by grid data.");
      return 1;
   }
   *_rita->ofh << "  tabulation  name=" << name;
   if (file_ok) {
      ifstream ip(file);
      if (ip.is_open())
         ip.close();
      else {
         _rita->msg("data>tabulation>","Unable to open file: "+file);
         return 1;
      }
      _theTab = new OFELI::Tabulation(file);
      *_rita->ofh << "  file=" << file;
   }
   else {
      if (dim1==1)
         _theGrid = new OFELI::Grid(xmin,xmax,nx);
      else if (dim1==2)
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,nx,ny);
      if (dim1==3)
         _theGrid = new OFELI::Grid(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
   }
   NameTab.push_back(name);
   theTab.push_back(_theTab);
   aTab.push_back(true);
   nb_tabs++;
   *_rita->ofh << endl;
   return 0;
}


int data::getPar(int n, const string& msg, int& v)
{  
   string s="";
   if (n<0)
      _cmd->get(s);
   else
      s = _cmd->string_token(n);
   if (checkParam(s,v)<0) {
      _rita->msg(msg,"Undefined parameter: "+s);
      return 1;
   }
   return 0;
}


int data::getPar(int n, const string& msg, double& v)
{
   string s="";
   if (n<0)
      _cmd->get(s);
   else
      s = _cmd->string_token(n);
   if (_cmd->isNumeric(s)) {
      v = stof(s);
      return 0;
   }
   if (checkParam(s,v)<0) {
      _rita->msg(msg,"Undefined parameter: "+s);
      return 1;
   }
   return 0;
}


int data::setDesc(const string& name,
                  const string& desc)
{
   if (getType(name)==DataType::NOTHING)
      return 1;
   Desc[name] = desc;
   return 0;
}


int data::save()
{
   string name="", file="", format="ofeli", t="";
   double every=1, e=1;
   bool name_ok=false, file_ok=false;
   static const vector<string> kw {"name","file","format","every"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("save>","No argument given.\nAvailable arguments: name, file, format, every.","");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 0:
            if ((t=_cmd->string_token())!="")
               name = t;
            name_ok = true;
            break;

         case 1:
            if ((t=_cmd->string_token())!="")
               file = t;
            file_ok = true;
            break;

         case 2:
            if ((t=_cmd->string_token())!="")
               format = t;
            break;

         case 3:
            if ((e=_cmd->int_token())!=0)
               every = e;
            break;

         default:
            _rita->msg("data>save>","Unknown argument.");
            cout << "Usage: save name=na file=fi [format=ofeli] [frequency=1]\n";
            cout << "Available Arguments\n";
            cout << "name: Name of entity to save\n";
            cout << "file: File where to save data\n";
            cout << "format: File format\n";
            cout << "every: Frequency of saving (Default: 1)" << endl;
            _ret = 1;
            break;
      }
   }
   if (name_ok==false) {
      _rita->msg("save>","No data to save.");
      return 1;
   }
   if (file_ok==false) {
      _rita->msg("save>","No file given.");
      return 1;
   }

// Case of a parameter
   int k = checkName(name,DataType::PARAM);
   if (k>0) {
      _rita->msg("save>","Parameters cannot be saved.");
      return 1;
   }

// Save a vector
   k = checkName(name,DataType::VECTOR);
   if (k>0) {
      if (format=="ofeli") {
         OFELI::IOField ffo(file,OFELI::IOField::OUT);
         ffo.put(*theVector[k]);
         return 0;
      }
      else if (format=="gmsh") {
         int ret = saveGmsh(file,*theVector[k]);
         if (ret==1) {
            _rita->msg("data>save>","Cannot save in gmsh format a vector without mesh association.");
            return 1;
         }
      }
      else if (format=="vtk") { 
         int ret = saveVTK(file,*theVector[k]);
         if (ret==1) {
            _rita->msg("data>save>","Cannot save in gmsh format a vector without mesh association.");
            return 1;
         }
      }
      else if (format=="tecplot") { 
         int ret = saveTecplot(file,*theVector[k]);
         if (ret==1) {
            _rita->msg("data>save>","Cannot save in gmsh format a vector without mesh association.");
            return 1;
         }
      }
      else if (format=="gnuplot") {
         ofstream ffo(file);
         ffo << endl;
         return 0;
      }
   }

// Save a history vector
   k = checkName(name,DataType::HVECTOR);
   if (k>0) {
      if (format=="ofeli") {
         int ret = theHVector[k]->saveOFELI(file,every);
         return ret;
      }
      else if (format=="gmsh") {
         int ret = theHVector[k]->saveGmsh(file,every);
         if (ret==1) {
            _rita->msg("data>save>","Cannot save in gmsh format a vector without mesh association.");
            return 1;
         }
      }
      else if (format=="vtk") {
         int ret = theHVector[k]->saveVTK(file,every);
         if (ret==1) {
            _rita->msg("data>save>","Cannot save in vtk format a vector without mesh association.");
            return 1;
         }
      }
      else if (format=="tecplot") {
         int ret = theHVector[k]->saveTecplot(file,every);
         if (ret==1) {
            _rita->msg("data>save>","Cannot save in vtk format a vector without mesh association.");
            return 1;
         }
      }
      else if (format=="gnuplot") {
         int ret = theHVector[k]->saveGnuplot(file,every);
         return ret;
      }
   }

// Save a matrix
   k = checkName(name,DataType::MATRIX);
   if (k>0) {
      OFELI::saveMatrix(theMatrix[k],file);
      return 0;
   }

// Save a mesh
   k = checkName(name,DataType::MESH);
   if (k>0) {
      if (format=="ofeli") {
         saveMesh(file,*theMesh[k],OFELI_FF);
         return 0;
      }
      else if (format=="gmsh") {
         saveMesh(file,*theMesh[k],GMSH);
         return 0;
      }
      else if (format=="vtk") {
         saveMesh(file,*theMesh[k],VTK);
         return 0;
      }
      else if (format=="gnuplot") {
         saveMesh(file,*theMesh[k],GNUPLOT);
         return 0;
      }
   }

// Save a grid
   k = checkName(name,DataType::GRID);
   if (k>0) {
      _rita->msg("data>save>","This option is not yet implemented.");
      return 1;
   }

// Case of a function
   k = checkName(name,DataType::FCT);
   if (k>0) {
      _rita->msg("save>","Functions cannot be saved.");
      return 1;
   }

// Save a tabulation
   k = checkName(name,DataType::TAB);
   if (k>0) {
      _rita->msg("save>","This type of data cannot be saved in file.");
      return 0;
   }
   return 1;
}


void data::ListParams(int opt)
{
   if (opt && !nb_params) {
      cout << "No defined parameters." << endl;
      return;
   }
   cout << "Number of defined parameters: " << nb_params << endl;
   for (int i=1; i<theParam.size(); ++i)
      if (aParam[i]) {
         cout << Param[i] << " = " << *theParam[i] << endl;
         if (Desc[Param[i]]!=" ")
            cout << "Description: " << Desc[Param[i]] << endl;
      }
}


void data::ListVectors(int opt)
{
   if (opt && !nb_vectors) {
      cout << "No defined vectors." << endl;
      return;
   }
   cout << "Number of vectors: " << nb_vectors << endl;
   for (int k=1; k<theVector.size(); ++k) {
      if (aVector[k]) {
         if (VectorSizeType[k]==DataSize::GIVEN_SIZE)
            cout << "Vector: " << Vector[k] << ", Size: " << theVector[k]->size() << endl;
         else
            cout << "Vector: " << Vector[k] << ", Number of degrees of freedom: " << nb_dof[k] << endl;
      }
   }
}


void data::ListHVectors(int opt)
{
   if (opt && !nb_hvectors) {
      cout << "No defined history vectors." << endl;
      return;
   }
   cout << "Number of history vectors: " << nb_hvectors << endl;
   for (int k=1; k<theHVector.size(); ++k)
      cout << "History Vector: " << HVector[k] << ", Size: " << theHVector[k]->size() << endl;
}


void data::ListFunctions(int opt)
{
   if (opt && !nb_fcts) {
      cout << "No defined functions." << endl;
      return;
   }
   cout << "Number of functions: " << nb_fcts << endl;
   for (int k=1; k<theFct.size(); ++k) {
      if (aFct[k]) {
         Fct *f = theFct[k];
         cout << "Function: " << f->name << ", Variable(s): ";
         for (int j=0; j<int(f->nb_var)-1; ++j)
            cout << f->var[j] << ",";
         cout << f->var[theFct[k]->nb_var-1] << ", Definition: " << f->expr << endl;
      }
   }
}


void data::ListTabs(int opt)
{
   if (opt && !nb_tabs) {
      cout << "No defined tabulations." << endl;
      return;
   }
   cout << "Number of tabulations: " << nb_tabs << endl;
   for (int i=1; i<theTab.size(); ++i) {
      if (aTab[i])
         cout << "Tabulation: " << NameTab[i] << ", Nb. of variables: " 
              << theTab[i]->getNbVar(1) << ", Size: " << theTab[i]->getSize(1,1) << endl;
   }
}


void data::ListMatrices(int opt)
{
   if (opt && !nb_matrices) {
      cout << "No defined matrices." << endl;
      return;
   }
   cout << "Number of matrices: " << nb_matrices << endl;
   for (int i=1; i<theMatrix.size(); ++i) {
      if (aMatrix[i])
         cout << "Matrix: " << Matr[i] << ", Size: " << theMatrix[i]->getNbRows()
              << " x " << theMatrix[i]->getNbColumns() << endl;
   }
}


void data::ListGrids(int opt)
{
   if (opt && !nb_grids) {
      cout << "No defined grids." << endl;
      return;
   }
   cout << "Number of grids: " << nb_grids << endl;
   for (int k=1; k<theGrid.size(); ++k) {
      if (aGrid[k]) {
         OFELI::Grid *g = theGrid[k];
         cout << "Grid No.            " << k << endl;
         cout << "Grid name:          " << NameGrid[k] << endl;
         cout << "Space dimension:    " << g->getDim() << endl;
         if (g->getDim()==1) {
            cout << "Domain:                    (" << g->getX(1) << ","
                 << g->getX(theGrid[k]->getNx()+1) << ")" << endl;
            cout << "Number of intervals:    " << g->getNx() << endl;
         }
         else if (g->getDim()==2) {
            cout << "Domain:                    (" << g->getX(1) << ","
                 << g->getX(g->getNx()+1) << ")x(" << g->getY(1) << ","
                 << g->getY(g->getNy()+1) << ")" << endl;
            cout << "Number of intervals:    " << g->getNx() << " x " << g->getNy() << endl;
         }
         else if (g->getDim()==3) {
            cout << "Domain:                    (" << g->getX(1) << ","
                 << g->getX(g->getNx()+1) << ")x(" << g->getY(1) << ","
                 << g->getY(g->getNy()+1) << ")x(" << g->getZ(1) << ","
                 << g->getZ(g->getNz()+1) << ")" << endl;
            cout << "Number of intervals:    " << g->getNx() << " x " << g->getNy() << " x " << g->getNz() << endl;
         }
      }
   }
}


void data::ListMeshes(int opt)
{
   if (opt && !nb_meshes) {
      cout << "No defined meshes." << endl;
      return;
   }
   cout << "Number of meshes: " << nb_meshes << endl;
   for (int k=1; k<theMesh.size(); ++k) {
      if (aMesh[k]) {
         OFELI::Mesh *m = theMesh[k];
         cout << "Mesh No.            " << k << endl;
         cout << "Mesh name:          " << NameMesh[k] << endl;
         cout << "Number of nodes:    " << m->getNbNodes() << endl;
         cout << "Number of elements: " << m->getNbElements() << endl;
         cout << "Number of sides:    " << m->getNbSides() << endl;
      }
   }
}


void data::ListODE(int opt)
{
   if (opt && !nb_ode) {
      cout << "No defined ordinary differential equations." << endl;
      return;
   }
   cout << "Number of ordinary differential equations: " << nb_ode << endl;
   for (int k=1; k<theODE.size(); ++k) {
      if (aODE[k]) {
         odae *ode = theODE[k];
         cout << "ODE No.            " << k << endl;
         cout << "ODE name:          " << ODE[k] << endl;
         if (ode->size>1)
            cout << "Size:           " << ode->size << endl;
      }
   }
}


void data::ListPDE(int opt)
{
   if (opt && !nb_pde) {
      cout << "No defined partial differential equations." << endl;
      return;
   }
   cout << "Number of partial differential equations: " << nb_pde << endl;
   for (int k=1; k<thePDE.size(); ++k) {
      if (aPDE[k]) {
         equa *e = thePDE[k];
         cout << "  PDE No.            " << k << endl;
         cout << "  PDE name: " << e->name << endl;
         cout << "  PDE id: " << e->eq << endl;
         cout << "  PDE unknown vector(s): ";
         for (int i=0; i<e->nb_vectors-1; ++i)
            cout << e->fd[i].fn << ", ";
         cout << e->fd[nb_vectors-1].fn << endl;
      }
   }
}


void data::ListAE(int opt)
{
   if (opt && !nb_ae) {
      cout << "No defined algebraic equations." << endl;
      return;
   }
   cout << "Number of algebraic equations: " << nb_ae << endl;
   for (int k=1; k<=theAE.size(); ++k) {
      if (aAE[k]) {
         odae *ae = theAE[k];
         cout << "Algebraic System No.   " << k << endl;
         cout << "Algebraic System name: " << AE[k] << endl;
         if (ae->size>1)
            cout << "Size:           " << ae->size << endl;
         if (ae->isFct) {
            if (ae->size==1)
               cout << "Equation defined by function: " << ae->theFct[0].name << endl;
            else {
               for (int i=0; i<ae->size; ++i)
               cout << "Equation: " << i+1 << ", defined by function: " << ae->theFct[i].name << endl;
            }
         }
         else {
            if (ae->size==1) {
               cout << "Equation defined by: " << ae->theFct[0].expr << endl;
               cout << "Variable is          " << ae->theFct[0].var[0] << endl;
            }
            else {
               for (int i=0; i<ae->size; ++i) 
                  cout << "Equation: " << i+1 << ", defined by: " << ae->theFct[i].expr << endl;
               for (int i=0; i<ae->size; ++i)
                  cout << "Variable " << i+1 << " is " << ae->theFct[0].var[i] << endl;
            }
         }
      }
   }
}


void data::Summary()
{
   if (nb_params+nb_vectors+nb_matrices+nb_fcts+nb_tabs+nb_grids+nb_meshes+
      nb_ae+nb_ode+nb_pde==0) {
      cout << "No defined data." << endl;
      return;
   }
   cout << "SUMMARY OF DATA:" << endl;
   cout << "---------------------------------------------------------------" << endl;
   if (nb_params) {
      ListParams(0);
      cout << endl;
   }
   if (nb_vectors) {
      ListVectors(0);
      cout << endl;
   }
   if (nb_matrices) {
      ListMatrices(0);
      cout << endl;
   }
   if (nb_hvectors) {
      ListHVectors(0);
      cout << endl;
   }
   if (nb_fcts) {
      ListFunctions(0);
      cout << endl;
   }
   if (nb_tabs) {
      ListTabs(0);
      cout << endl;
   }
   if (nb_grids) {
      ListGrids(0);
      cout << endl;
   }
   if (nb_meshes) {
      ListMeshes(0);
      cout << endl;
   }
   if (nb_ae) {
      ListAE(0);
      cout << endl;
   }
   if (nb_ode) {
      ListODE(0);
      cout << endl;
   }
   if (nb_pde)
      ListPDE(0);
   cout << "---------------------------------------------------------------" << endl;
}


int data::saveGmsh(const string &file, const Vect<double>& v)
{
   using namespace OFELI;
   int nb_en=0;
   ofstream ff(file);

   if (v.WithMesh()==false)
      return 1;
   Mesh &ms = v.getMesh();
   int nb_dof = int(v.getNbDOF());
   ff << "View \"" << v.getName() << "\" {" << endl;
   switch (ms.getDim()) {

      case 1:
         element_loop(&ms) {
            ff << "SL(";
            ff << The_element(1)->getX() <<  ", 0., 0., "
               << The_element(2)->getX() <<  ", 0., 0. ) {" << endl;
            ff << v(The_element(1)->n(),1) << "," << v(The_element(2)->n(),1);
            ff << endl;
         }
         ff << "};" << endl;
         break;

      case 2:
         element_loop(&ms) {
            if (nb_dof==1)
               ff << 'S';
            else
               ff << 'V';
            int nb_en=int(The_element.getNbNodes());
            if (nb_en==3)
               ff << "T(";
            else if (nb_en==4)
               ff << "Q(";
            for (int k=1; k<nb_en; ++k)
               ff << The_element(k)->getX() << "," << The_element(k)->getY() << ",0.,";
            ff << The_element(nb_en)->getX() << "," << The_element(nb_en)->getY() << ",0.) {" << endl;
            for (int k=1; k<=nb_en; ++k) {
               ff << v(The_element(k)->n(),1);
               if (nb_dof > 1)
                  ff << "," << v(The_element(k)->n(),2) << ",0.0";
               if (k<nb_en)
                  ff << ",";
            }
            ff << "};" << endl;
         }
         ff << "};" << endl;
            break;

         case 3:
            element_loop(&ms) {
               if (nb_dof==1)
                  ff << 'S';
               else
                  ff << 'V';
               if ((nb_en=The_element.getNbNodes())==4)
                  ff << "S(";
               else if (nb_en==8)
                  ff << "H(";
               else if (nb_en==6)
                  ff << "I(";
               else if (nb_en==5)
                  ff << "Y(";
               for (int k=1; k<=nb_en-1; ++k)
                  ff << The_element(k)->getX() << ","
                     << The_element(k)->getY() << ","
                     << The_element(k)->getZ() << ",";
               ff << The_element(nb_en)->getX() << ","
                  << The_element(nb_en)->getY() << ","
                  << The_element(nb_en)->getZ() << ") {" << endl;
               for (int k=1; k<=nb_en; ++k) {
                  ff << v(The_element(k)->n(),1);
                  if (nb_dof > 1)
                     ff << "," << v(The_element(k)->n(),1) << ","
                        << v(The_element(k)->n(),3) << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;
      }
      return 0;
}


int data::saveVTK(const string &file, const Vect<double>& v)
{
   using namespace OFELI;
   map<int,int> shCode = {{LINE,2},{TRIANGLE,3},{QUADRILATERAL,4},{TETRAHEDRON,4},
                          {HEXAHEDRON,8},{PENTAHEDRON,6}};
   map<int,int> ShCode = {{LINE,3},{TRIANGLE,5},{QUADRILATERAL,9},{TETRAHEDRON,10},
                          {HEXAHEDRON,12},{PENTAHEDRON,13}};

   if (v.WithMesh()==false)
      return 1;
   int nb_dof = int(v.getNbDOF());
   Mesh &ms = v.getMesh();
   int sz=0;
   element_loop(&ms)
      sz += shCode[The_element.getShape()] + 1;
   ofstream ff(file.c_str());
   ff << setprecision(16) << std::scientific;
   ff << "# vtk DataFile Version 2.0\n# Imported from OFELI files\nASCII" << endl;
   ff << "DATASET UNSTRUCTURED_GRID\nPOINTS " << ms.getNbNodes() << " double" << endl;
   node_loop(&ms)
      ff << The_node.getX() << "  " << The_node.getY() << "  " << The_node.getZ() << endl;
   ff << "\nCELLS " << ms.getNbElements() << setw(10) << sz << endl;
   element_loop(&ms) {
      ff << setw(8) << shCode[The_element.getShape()];
      for (int i=1; i<=shCode[The_element.getShape()]; ++i)
         ff << setw(10) << The_element(i)->n()-1;
      ff << endl;
   }
   ff << "\nCELL_TYPES  " << ms.getNbElements() << endl;
   int k=0;
   element_loop(&ms) {
      ff << setw(4) << ShCode[The_element.getShape()];
      if (++k%30 == 0)
         ff << endl;
   }
   ff << "\nPOINT_DATA  " << ms.getNbNodes() << endl;

   if (nb_dof==1)
      ff << "SCALARS  " << v.getName()<< "  double  1\nLOOKUP_TABLE  default" << endl;
   else
      ff << "VECTORS  " << v.getName() << "  double" << endl;

   node_loop(&ms) {
      ff << v(node_label,1) << " ";
      if (nb_dof>1) {
         ff << v(node_label,2) << " ";
         if (nb_dof > 2)
            ff << v(node_label,3) << " ";
         else
            ff << 0. << " ";
      }
      ff << endl;
   }
   return 0;
}


int data::saveTecplot(const string& file, const Vect<double>& v)
{
   using namespace OFELI;
   map<int,string> shape = {{LINE,"LINESEG"},{QUADRILATERAL,"QUADRILATERAL"},{TRIANGLE,"TRIANGLE"},
                            {TETRAHEDRON,"TETRAHEDRON"},{HEXAHEDRON,"HEXAHEDRON"},{PENTAHEDRON,"HEXAHEDRON"}};
   if (v.WithMesh()==false)
      return 1;
   int nb_dof = int(v.getNbDOF());
   Mesh &ms = v.getMesh();
   ofstream ff(file.c_str());
   ff.setf(ios::right|ios::scientific);
   ff << "TITLE = \" \"\n\nVARIABLES = \"X\", \"Y\"";
   if (ms.getDim()==3)
      ff << ", \"Z\"";
   if (nb_dof==1)
      ff << ", \"T\"";
   else if (nb_dof==2)
      ff << ", \"UX\", \"UY\"";
   else if (nb_dof==3)
      ff << ", \"UX\", \"UY\", \"UZ\"";
   ff << "\n\nZONE T=\"" << "step-1" << "\", N=" << ms.getNbNodes() << ", E="
      << ms.getNbElements() << ", F=FEPOINT, ET=" << shape[ms.getShape()]
      << ", SOLUTIONTIME=" << v.getTime();
   ff << ", D=(1,";
   if (ms.getDim()>1)
      ff << "2,";
   if (ms.getDim()==3)
      ff << "3,";
   ff << "FECONNECT)";
   ff << endl;
   node_loop(&ms) {
     for (int i=1; i<=int(ms.getDim()); i++)
         ff << "  " << The_node.getCoord(i);
         for (int j=0; j<nb_dof; j++)
             ff << "  " << v[nb_dof*(node_label-1)+j];
         ff << endl;
   }
   element_loop(&ms) {
     for (int i=1; i<=int(The_element.getNbNodes()); ++i)
         ff << setw(10) << The_element(i)->n();
      ff << endl;
   }
   ff.close();
   return 0;
}


odae::odae()
     : isSet(false), log(false), isFct(false), vect(-1)
{
   every = 1;
}


void odae::setVars(int opt)
{
   if (opt)
      vars.push_back("t");
   if (size==1)
      vars.push_back(fn);
   else {
      for (int i=0; i<size; ++i)
         vars.push_back(fn+to_string(i+1));
   }
}


ostream& operator<<(ostream& s, const odae& e)
{
   map<NonLinearIter,string> Nls = {{BISECTION,"bisection"},
                                    {REGULA_FALSI,"regula-falsi"},
                                    {PICARD,"picard"},
                                    {SECANT,"secant"},
                                    {NEWTON,"newton"}};
   if (e.type==DataType::AE) {
      s << "Algebraic Equation Name: " << e.name << endl;
      s << "Size: " << e.size << endl;
      if (e.size==1)
         s << "Variable: " << e.fn << endl;
      else {
         s << "Variable(s): ";
         for (int i=0; i<e.size-1; ++i)
            s << e.vars[i] << ", ";
         s << e.vars[e.size-1] << endl;
      }
      s << "Nonlinear iteration procedure: " << Nls[e.nls] << endl;

   }
   else if (e.type==DataType::ODE) {
      s << "Ordinary Differential Equation Name: " << e.name << endl;
      s << "Size: " << e.size << endl;
      s << "Name of variable(s): " << e.fn << endl;
      s << "Time integration scheme: " << e.scheme << endl;
      s << "Initial time: " << e.init_time << endl;
      s << "Final time: " << e.final_time << endl;
      s << "Time step: " << e.time_step << endl;
      s << "Nonlinear iteration procedure: " << Nls[e.nls] << endl;
   }
   return s;
}

} /* namespace RITA */
