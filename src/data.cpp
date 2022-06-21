/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2022 Rachid Touzani

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

using std::cout;
using std::endl;
using std::map;

namespace RITA {

data::data(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _size(0), _default_field(0), _verb(1), _configure(config), _cmd(command),
       _theMesh(nullptr), _theTab(nullptr), _theGrid(nullptr), _theFct(nullptr), 
       _theVector(nullptr), _theMatrix(nullptr), _u(nullptr), _theParam(nullptr)
{
   nb_fields = nb_params = nb_fcts = nb_meshes = nb_tabs = nb_grids = nb_vectors = nb_matrices = 0;
   nb_pde = nb_ode = nb_ae = nb_int = nb_eigen = nb_eq = 0;
   iMesh = iField = iMatrix = iGrid = iParam = iFct = iTab = iAE = iODE = iPDE = iEq = 0;
   theParam.push_back(nullptr);
   u.push_back(nullptr);
   theFct.push_back(nullptr);
   theMatrix.push_back(nullptr);
   theMesh.push_back(nullptr);
   theGrid.push_back(nullptr);
   theAE.push_back(nullptr);
   theODE.push_back(nullptr);
   thePDE.push_back(nullptr);
   theTab.push_back(nullptr);
   AE.push_back(" ");
   ODE.push_back(" ");
   PDE.push_back(" ");
   fct_name.push_back(" ");
   mesh_name.push_back(" ");
   grid_name.push_back(" ");
   tab_name.push_back(" ");
   Param.push_back(" ");
   Field.push_back(" ");
   Matr.push_back(" ");
   FieldType.push_back(eqType::AE);
   FieldSizeType.push_back(DataSize::GIVEN_SIZE);
   FieldEquation.push_back(0);
   eq_type.push_back(eqType::AE);
   eqq.push_back(0);
   _u = nullptr;
   _ret = 0;
   ok = false;
}


data::~data()
{
/*   if (_theGrid!=nullptr)
      delete _theGrid;
   if (_theMesh!=nullptr)
      delete _theMesh;
   if (_theFct!=nullptr)
      delete _theFct;
   if (_theTab!=nullptr)
      delete _theTab;
   if (_theVector!=nullptr)
      delete _theVector;
//   if (_theMatrix!=nullptr)
//      delete _theMatrix;
//   if (_u!=nullptr)
//      delete _u;
   if (_theParam!=nullptr)
      delete _theParam;*/
}


int data::checkField(const string& name)
{
   int k = FieldName[name];
   if (k==0)
      return -1;
   return k;
}


int data::checkMatrix(const string& name)
{
   return checkName(name,DataType::MATRIX);
}


int data::checkMesh(const string& name)
{
   auto it = find(mesh_name.begin(),mesh_name.end(),name);
   if (it != mesh_name.end())
      return distance(mesh_name.begin(),it);
   else
      return -1;
}


int data::checkGrid(const string& name)
{
   auto it = find(grid_name.begin(),grid_name.end(),name);
   if (it != grid_name.end())
      return distance(grid_name.begin(),it);
   else
      return -1;
}


int data::checkFct(const string& name)
{
   int k = FctName[name];
   if (k==0)
      return -1;
   return k;
}


int data::checkName(const string&   name,
                    const DataType& dt,
                    int             opt)
{
   if (dn.size()==0)
      return 0;
   DataType d = dn[name].dt;
   if (d==dt) {
      switch (d) {

         case DataType::PARAM:
            return ParamName[name];

         case DataType::FIELD:
            return FieldName[name];

         case DataType::MATRIX:
            return MatrixName[name];

         case DataType::GRID:
            return GridName[name];

         case DataType::MESH:
            return MeshName[name];

         case DataType::TAB:
            return TabName[name];

         case DataType::FCT:
            return FctName[name];

         case DataType::AE:
            return AEName[name];

         case DataType::ODE:
            return ODEName[name];

         case DataType::PDE:
            return PDEName[name];
      }
   }
   if (opt==0)
      return 0;

   switch (d) {

      case DataType::PARAM:
         _rita->msg("data>","Name "+name+" already used for a parameter.");
         return -2;

      case DataType::FIELD:
         _rita->msg("data>","Name "+name+" already used for a field.");
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
   }
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
   grid_name.push_back(name);
   theGrid.push_back(_theGrid);
   GridName[name] = dn[name].i = nb_grids;
   dn[name].dt = DataType::GRID;
}
    
    
void data::setTab2Field(OFELI::Tabulation* tab)
{
   int nb = tab->getNbVar(1);
   OFELI::fct &f = tab->Funct[0];
   int nx = f.Np[0], ny = f.Np[1], nz = f.Np[2];
   addGridField(tab->getFunctName(1),1);
   OFELI::Vect<double> *v = u[iField];

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
   if (AEName[name]==0) {
      iAE = ++nb_ae;
      iEq = ++nb_eq;
      eq_type.push_back(eqType::AE);
      eqq.push_back(iEq);
      AE.push_back(name);
      theAE.push_back(ae);
   }
   else {
      iAE = AEName[name];
      theAE[iAE] = ae;
   }
   if (name=="")
      name = "AE-"+to_string(iAE);
   AE[iAE] = name;
   AEName[name] = dn[name].i = nb_ae;
   dn[name].dt = DataType::AE;
   return iAE;
}


int data::addODE(odae*  ode,
                 string name)
{
   if (ODEName[name]==0) {
      iODE = ++nb_ode;
      iEq = ++nb_eq;
      eq_type.push_back(eqType::ODE);
      eqq.push_back(iEq);
      ODE.push_back(name);
      theODE.push_back(ode);
   }
   else {
      iODE = ODEName[name];
      theODE[iODE] = ode;
   }
   if (name=="")
      name = "ODE-"+to_string(iODE);
   ODE[iODE] = name;
   ODEName[name] = dn[name].i = nb_ode;
   dn[name].dt = DataType::ODE;
   return iODE;
}


int data::addPDE(equa*  pde,
                 string name)
{
   if (PDEName[name]==0) {
      iPDE = ++nb_pde;
      iEq = ++nb_eq;
      eq_type.push_back(eqType::PDE);
      eqq.push_back(iEq);
      PDE.push_back(name);
      thePDE.push_back(pde);
   }
   else {
      iPDE = PDEName[name];
      thePDE[iPDE] = pde;
   }
   PDEName[name] = dn[name].i = nb_pde;
   dn[name].dt = DataType::PDE;
   return iPDE;
}


int data::addFunction(const string&         def,
                      const vector<string>& var,
                      string                name)
{
   if (name=="")
      name = "F-" + to_string(nb_fcts+1);
   _theFct = new OFELI::Fct(name,def,var);
   if (FctName[name]==0) {
      iFct = ++nb_fcts;
      fct_name.push_back(name);
      theFct.push_back(_theFct);
   }
   else {
      cout << "Function "+name+" is redefined." << endl;
      _rita->msg("data>","Function "+name+" redefined.");
      iFct = FctName[name];
      theFct[iFct] = _theFct;
   }
   fct_name[iFct] = name;
   FctName[name] = dn[name].i = nb_fcts;
   dn[name].dt = DataType::FCT;
   return iFct;
}


int data::addMesh(OFELI::Mesh*  ms,
                  const string& name)
{
   if (MeshName[name]==0) {
      iMesh = ++nb_meshes;
      theMesh.push_back(ms);
      mesh_name.push_back(name);
   }
   else {
      iMesh = MeshName[name];
      theMesh[iMesh] = ms;
   }
   MeshName[name] = dn[name].i = nb_meshes;
   dn[name].dt = DataType::MESH;
   return iMesh;
}


int data::addParam(const string& name,
                   double        value)
{
   _theParam = new double;
   *_theParam = value;
   if (ParamName[name]==0) {
      iParam = ++nb_params;
      Param.push_back(name);
      theParam.push_back(_theParam);
   }
   else {
      iParam = ParamName[name];
      theParam[iParam] = _theParam;
   }
   ParamName[name] = dn[name].i = nb_params;
   dn[name].dt = DataType::PARAM;
   return iParam;
}


int data::addField(const string& name,
                   int           n,
                   string        file)
{
   if (file!="")
      n = 1;
   _u = new OFELI::Vect<double>(n);
   if (file!="") {
      OFELI::XMLParser xml(file,OFELI::XMLParser::MATRIX);
      xml.get(*_u);
   }
   if (FieldName[name]==0) {
      Field.push_back(name);
      iField = ++nb_fields;
      u.push_back(_u);
      FieldSizeType.push_back(DataSize::GIVEN_SIZE);
   }
   else {
      iField = FieldName[name];
      u[iField] = _u;
   }
   FieldName[name] = nb_fields;
   dn[name].i = nb_fields;
   dn[name].dt = DataType::FIELD;
   return iField;
}


int data::addMatrix(string name,
                    int    nr,
                    int    nc,
                    string file,
                    string s)
{
   if (file!="")
      nr = 1, nc = 1;
   _theMatrix = new DMatrix<double>(nr,nc);
   if (file!="") {
      OFELI::XMLParser xml(file,OFELI::XMLParser::MATRIX);
      xml.get(_theMatrix);
   }
   matrix_name.push_back(name);
   if (MatrixName[name]==0) {
      Matr.push_back(name);
      iMatrix = ++nb_matrices;
      theMatrix.push_back(_theMatrix);
   }
   else {
      iMatrix = MatrixName[name];
      theMatrix[iMatrix] = _theMatrix;
   }
   MatrixName[name] = dn[name].i = nb_matrices;
   dn[name].dt = DataType::MATRIX;
   return iMatrix;
}


int data::addTab(OFELI::Tabulation* tab,
                 const string&      name)
{
   if (TabName[name]==0) {
      tab_name.push_back(name);
      iTab = ++nb_tabs;
      theTab.push_back(tab);
   }
   else {
      iMatrix = TabName[name];
      theTab[iTab] = tab;
   }
   TabName[name] = dn[name].i = nb_tabs;
   dn[name].dt = DataType::TAB;
   return iTab;
}


int data::addMeshField(const string& name,
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
   if (FieldName[name]==0) {
      iField = ++nb_fields;
      Field.push_back(name);
      u.push_back(_u);
      FieldType.push_back(eqType::PDE);
   }
   else {
      iField = FieldName[name];
      u[iField] = _u;
      FieldType[iField] = eqType::PDE;
   }
   FieldName[name] = dn[name].i = nb_fields;
   dn[name].dt = DataType::FIELD;

   if (s==DataSize::NODES) {
      if (_theMesh->getNbNodes()==0) {
         _rita->msg("data>","Mesh has no nodes");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==DataSize::ELEMENTS) {
      if (_theMesh->getNbElements()==0) {
         _rita->msg("data>","Mesh has no elements.");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==DataSize::SIDES) {
      if (_theMesh->getNbSides()==0) {
         _rita->msg("data>","Mesh has no sides");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==DataSize::EDGES) {
      if (_theMesh->getNbEdges()==0) {
         _rita->msg("data>","Mesh has no edges");
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
   return iField;
}


int data::addGridField(const string& name,
                       int           nb_dof)
{
   _theGrid = theGrid[iGrid];
   if (_theGrid==nullptr) {
      _rita->msg("data>","No grid data available.");
      _ret = -1;
      return _ret;
   }
   _u = new OFELI::Vect<double>(*_theGrid);
   _nb_dof = nb_dof;
   _theGrid->setNbDOF(_nb_dof);
   if (FieldName[name]==0) {
      iField = ++nb_fields;
      Field.push_back(name);
      u.push_back(_u);
      FieldType.push_back(eqType::PDE);
   }
   else {
      iField = FieldName[name];
      u[iField] = _u;
      FieldType[iField] = eqType::PDE;
   }
   FieldName[name] = dn[name].i = nb_fields;
   dn[name].dt = DataType::FIELD;
   return iField;
}


int data::run()
{
   int key = 0;
   string td = "", fn="";
   static const vector<string> kw {"grid","mesh","field","tab$ulation","func$tion",
                                   "vect$or","matr$ix","save","sum$mary","list"};
   *_rita->ofh << "data" << endl;
   while (1) {
      if (_cmd->readline("rita>data> ")<0)
         continue;

      switch (key=_cmd->getKW(kw,_rita->_gkw)) {

         case   0:
            _ret = setGrid();
            break;

         case   1:
            _ret = 10;
            return _ret;

         case   2:
            _ret = setField();
            break;

         case   3:
            _ret = setTab();
            break;

         case   4:
            _ret = setFunction();
            break;

         case   5:
            _ret = setVector();
            break;

         case   6:
            _ret = setMatrix();
            break;

         case   7:
            _ret = save();
            break;

         case   8:
            Summary();
            break;

         case   9:
            if (_cmd->setNbArg(1,"Type of data to list.")) {
               _rita->msg("data>list>","Missing data type to list.","",1);
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
            else if (td.substr(0,5)=="field")
               ListFields(1);
            else if (td.substr(0,5)=="funct")
               ListFunctions(1);
            else if (td.substr(0,3)=="tab")
               ListTabs(1);
            else if (td.substr(0,5)=="matri")
               ListMatrices(1);
            _ret = 0;
            break;

         case 100:
         case 101:
            _cmd->setNbArg(0);
            getHelp();
            break;

         case 102:
            _rita->getLicense();
            break;

         case 103:
            _ret = _configure->run();
            _verb = _configure->getVerbose();
            break;

         case 104:
         case 105:
            _rita->setParam();
            break;

         case 106:
            if (_cmd->setNbArg(1,"Data name to be given.",1)) {
               _rita->msg("print>","Missing data name.","",1);
               break;
            }
            if (!_cmd->get(fn))
               print(fn);
            break;

         case 107:
            Summary();
            break;

         case 108:
         case 109:
            _ret = 0;
            ok = true;
            *_rita->ofh << "  end" << endl;
            return _ret;

         default:
            _rita->msg("data>","Unknown Command "+_cmd->token(),
                       "Available commands: grid, mesh, field, tabulation, function, vector, matrix,"
                       " save, summary, list");
            break;
      }
   }
   return 0;
}


void data::getHelp()
{
   cout << "\nAvailable Commands in mode 'data':\n";
   cout << "grid:       Define a grid\n";
   cout << "mesh:       Define a finite element mesh (can be accessed from the top level)\n";
   cout << "field:      Define a field (or vector)\n";
   cout << "matrix:     Define a matrix\n";
   cout << "tabulation: Define a tabulation\n";
   cout << "function:   Define a function\n";
   cout << "matrix:     Define a matrix\n";
   cout << "list:       List a specific entity\n";
   cout << "Global commands: \n";
   cout << "help or ?:  Display this help\n";
   cout << "set:        Set configuration data\n";
   cout << "par or @:   Defined a parameter (constant)\n";
   cout << "print       print a specific entity\n";
   cout << "save        save a specific entity in file\n";
   cout << "summary:    Summary of defined entities\n\n";
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

   k = checkName(s,DataType::FIELD);
   if (k>0) {
      cout << *u[k];
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
   k = checkName(s,DataType::TAB);
   if (k>0) {
      cout << *theTab[k];
      return;
   }
   cout << "Data " << s << " not found." << endl;
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
   grid_name.push_back(name);
   iGrid++;
   nb_grids++;
   return 0;
}


int data::setVector()
{
   int nb=0;
   string name="vect-"+to_string(nb_vectors+1);
   static const vector<string> kw {"name","size","def$ine","set"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>vector>","Error in command.","Available arguments: name, size, define, set.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
	   //           size = _cmd->int_token(0);
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
      *_rita->ofh << endl;
   }
   theVector.push_back(_theVector);
   vector_name.push_back(name);
   nb_vectors++;
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
      _rita->msg("data>matrix>","Error in command.");
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
   if (file=="") {
      _rita->msg("data>matrix>","No file name given for matrix input.");
      return 1;
   }
   if (storage!="dense") {
      _rita->msg("data>matrix>","Only dense storage is allowed in the current release.");
      return 1;
   }

   addMatrix(name,nr,nc,file,storage);
   return 0;
}


int data::setField()
{
   int nb=0;
   double umin=0., umax=0.;
   string name="u", type="size", arg=" ", mn="";
   static const vector<string> kw {"name","size","mesh","grid","nbdof","type","uniform"};
   static const vector<string> types {"size","nodes","elements","sides","edges"};
   map<string,DataSize> tt {{"size",DataSize::GIVEN_SIZE}, {"nodes",DataSize::NODES}, {"elements",DataSize::ELEMENTS},
                            {"sides",DataSize::SIDES}, {"edges",DataSize::EDGES}};
   if (theMesh[0]!=nullptr)
      _theMesh = theMesh[0];
   if (_default_field==1) {
      iField = 0;
      FieldName.clear();
   }
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>field>","Error in command.",
                 "Available arguments: name, size, mesh, grid, nbdof, type, uniform.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 0:
            name = _cmd->string_token();
            break;

         case 1:
            _size = _cmd->int_token();
            type = "size";
            break;

         case 2:
            mn = _cmd->string_token();
	    //            mesh_ok = 1;
            break;

         case 3:
            _nb_dof = _cmd->int_token();
            break;

         case 4:
            _nb_dof = _cmd->int_token();
            break;

         case 5:
            type = _cmd->string_token();
            if (type!="size" && type!="nodes" && type!="elements" &&
                type!="sides" && type!="edges") {
               _rita->msg("data>field>","Unknown type: "+type);
               return 0;
            }
            break;

         case 6:
            umin = _cmd->double_token();
            if (nb>1)
               umax = _cmd->double_token();
	    //            uniform_ok = nb;
            break;

         default:
            _rita->msg("data>field>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (nb_args>0) {
      int k = 0;
      if (type=="size")
         addField(name,_size);
      else
         k = addMeshField(name,tt[type],_size);
      if (k==0) {
         if (_verb>1) {
            if (type=="size") 
               cout << "New field: " << name << ", size: " << _size << endl;
            else
               cout << "New field: " << name << ", Nb. of DOF: " << _nb_dof << endl;
            cout << "Total number of fields: " << nb_fields << endl;
         }
         *_rita->ofh << "  field  name=" << name;
         if (_size)
            *_rita->ofh << " size=" << _size;
         *_rita->ofh << " type=" << type << "  nbdof=" << _nb_dof << " min=" << umin << " max=" << umax;
         nb_dof[iField] = _nb_dof;
      /*      if (uniform==1) {
         _rita->msg("data>field>","Minimal and maximal values must be given for field.");
	 return 1;
	 }*/
         *_rita->ofh << endl;
      }
   }
   return 0;
}


int data::setFunction()
{
   _ret = 0;
   bool var_ok=false, def_ok=false, name_ok=false;
   string vv="", def="", name="f";
   int nb=1;
   vector<string> var;
   static const vector<string> kw {"name","var","field","nb","def","definition"};

   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<0) {
      _rita->msg("data>function>:","Error in command.");
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

/*
void data::Clear()
{
   *_rita->ofh << "  clear" << endl;
   nb_fields = 0;
   Field[0] = "u";
   FieldName["u"] = 0;
   FieldEquation[0] = 0;
   if (_verb)
      cout << "Field data cleared." << endl;
   if (u[_ifield].size()>0) {
      u[_ifield] = 0.;
   }
   if (bc[_ifield].size()>0) {
      bc[_ifield] = 0.;
      if (_verb)
         cout << "Boundary condition vector cleared." << endl;
   }
   if (sf[_ifield].size()>0) {
      sf[_ifield] = 0.;
      if (_verb)
         cout << "Surface force condition vector cleared." << endl;
   }
   if (bf[_ifield].size()>0) {
      bf[_ifield] = 0.;
      if (_verb)
         cout << "Body force vector cleared." << endl;
   }
   _ifield = 0;
   _ret = 90;
}*/


int data::setTab()
{
   int dim1=0, dim2=0, dim3=0;
   string file="";
   string name = "T" + to_string(theTab.size()+1);
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, nx=10, ny=10, nz=10, grid_ok=0, file_ok=0, field_ok=0;
   static const vector<string> kw {"name","file","min","max","ne","field"};
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
            field_ok = true;
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
   if (!field_ok && !file_ok) {
      _rita->msg("data>tabulation>","No associated field given.");
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
   tab_name.push_back(name);
   theTab.push_back(_theTab);
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


int data::save()
{
   string name="", file="", format="", t="";
   bool name_ok=false, file_ok=false, format_ok=true;
   static const vector<string> kw {"name","file","format"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("data>save>","No argument given.");
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
            format_ok = true;
            break;

         default:
            _rita->msg("data>save>","Unknown argument.");
            cout << "\nAvailable Arguments\n";
            cout << "name: Name of entity to save\n";
            cout << "file: File where to save datum" << endl;
            _ret = 1;
            break;
      }
   }
   if (name_ok==false) {
      _rita->msg("data>save>","No data to save.");
      return 1;
   }

   int k = checkName(name,DataType::FIELD);
   if (k>0) {
      OFELI::IOField ffo(file,OFELI::IOField::OUT);
      ffo.put(*u[k]);
      return 0;
   }
   k = checkName(name,DataType::MATRIX);
   if (k>0) {
      OFELI::saveMatrix(theMatrix[k],file);
      return 0;
   }
   k = checkName(name,DataType::MESH);
   if (k>0) {
      saveMesh(file,*theMesh[k],GMSH);
      return 0;
   }
   k = checkName(name,DataType::GRID);
   if (k>0) {
      _rita->msg("data>save>","This option is not yet implemented.");
      return 0;
   }
   k = checkName(name,DataType::FCT);
   if (k>0) {
      _rita->msg("data>save>","This type of data cannot be saved in file.");
      return 0;
   }
   k = checkName(name,DataType::TAB);
   if (k>0) {
      _rita->msg("data>save>","This type of data cannot be saved in file.");
      return 0;
   }
   _rita->msg("data>save>","No data entity "+name+" found.");
   return 1;
}


void data::ListParams(int opt)
{
   if (opt && !nb_params) {
      cout << "No defined parameters." << endl;
      return;
   }
   cout << "Number of defined parameters: " << nb_params << endl;
   for (int i=1; i<=nb_params; ++i)
      cout << Param[i] << " = " << *theParam[i] << endl;
}


void data::ListFields(int opt)
{
   if (opt && !nb_fields) {
      cout << "No defined fields." << endl;
      return;
   }
   cout << "Number of fields: " << nb_fields << endl;
   for (int k=1; k<=nb_fields; ++k) {
      if (FieldSizeType[k]==DataSize::GIVEN_SIZE)
         cout << "Field: " << Field[k] << ", Size: " << u[k]->size() << endl;
      else
         cout << "Field: " << Field[k] << ", Number of degrees of freedom: " << nb_dof[k] << endl;
   }
}


void data::ListFunctions(int opt)
{
   if (opt && !nb_fcts) {
      cout << "No defined functions." << endl;
      return;
   }
   cout << "Number of functions: " << nb_fcts << endl;
   for (int k=1; k<=nb_fcts; ++k) {
      Fct *f = theFct[k];
      cout << "Function: " << f->name << ", Variable(s): ";
      for (int j=0; j<int(f->nb_var)-1; ++j)
         cout << f->var[j] << ",";
      cout << f->var[theFct[k]->nb_var-1] << ", Definition: " << f->expr << endl;
   }
}


void data::ListTabs(int opt)
{
   if (opt && !nb_tabs) {
      cout << "No defined tabulations." << endl;
      return;
   }
   cout << "Number of tabulations: " << nb_tabs << endl;
   for (int i=1; i<=nb_tabs; ++i) {
      cout << "Tabulation: " << tab_name[i] << ", Nb. of variables: " 
           << theTab[i]->getNbVar(1) << ", Size: " << theTab[i]->getSize(1,1) << endl;
   }
}


void data::ListVectors(int opt)
{
   if (opt && !nb_tabs) {
      cout << "No defined vectors." << endl;
      return;
   }
   cout << "Number of vectors: " << nb_vectors << endl;
}


void data::ListMatrices(int opt)
{
   if (opt && !nb_matrices) {
      cout << "No defined matrices." << endl;
      return;
   }
   cout << "Number of matrices: " << nb_matrices << endl;
   for (int i=1; i<=nb_matrices; ++i) {
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
   for (int k=1; k<=nb_grids; ++k) {
      OFELI::Grid *g = theGrid[k];
      cout << "Grid No.            " << k << endl;
      cout << "Grid name:          " << grid_name[k] << endl;
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


void data::ListMeshes(int opt)
{
   if (opt && !nb_meshes) {
      cout << "No defined meshes." << endl;
      return;
   }
   cout << "Number of meshes: " << nb_meshes << endl;
   for (int k=1; k<=nb_meshes; ++k) {
      OFELI::Mesh *m = theMesh[k];
      cout << "Mesh No.            " << k << endl;
      cout << "Mesh name:          " << mesh_name[k] << endl;
      cout << "Number of nodes:    " << m->getNbNodes() << endl;
      cout << "Number of elements: " << m->getNbElements() << endl;
      cout << "Number of sides:    " << m->getNbSides() << endl;
   }
}


void data::ListODE(int opt)
{
   if (opt && !nb_ode) {
      cout << "No defined ordinary differential equations." << endl;
      return;
   }
   cout << "Number of ordinary differential equations: " << nb_ode << endl;
   for (int k=1; k<=nb_ode; ++k) {
      odae *ode = theODE[k];
      cout << "ODE No.            " << k << endl;
      cout << "ODE name:          " << ODE[k] << endl;
      if (ode->size>1)
         cout << "Size:           " << ode->size << endl;
   }
}


void data::ListPDE(int opt)
{
   if (opt && !nb_pde) {
      cout << "No defined partial differential equations." << endl;
      return;
   }
   cout << "Number of partial differential equations: " << nb_pde << endl;
   for (int k=1; k<=nb_pde; ++k) {
      equa *e = thePDE[k];
      cout << "PDE No.            " << k << endl;
      cout << "  PDE name: " << e->eq << endl;
      cout << "  PDE unknown field(s): ";
      for (int i=0; i<e->nb_fields-1; ++i)
         cout << e->fd[i].fn << ", ";
      cout << e->fd[nb_fields-1].fn << endl;
   }
}


void data::ListAE(int opt)
{
   if (opt && !nb_ae) {
      cout << "No defined algebraic equations." << endl;
      return;
   }
   cout << "Number of algebraic equations: " << nb_ae << endl;
   for (int k=1; k<=nb_ae; ++k) {
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


void data::Summary()
{
   if (nb_fields+nb_params+nb_fields+nb_matrices+nb_fcts+nb_tabs+nb_grids+nb_meshes+
       nb_ae+nb_ode+nb_pde==0) {
      cout << "No defined data." << endl;
      return;
   }
   cout << "SUMMARY OF DATA:" << endl;
   cout << "---------------------------------------------------------------" << endl;
   ListParams(0);
   cout << endl;
   ListFields(0);
   cout << endl;
   ListMatrices(0);
   cout << endl;
   ListFunctions(0);
   cout << endl;
   ListTabs(0);
   cout << endl;
   ListGrids(0);
   cout << endl;
   ListMeshes(0);
   cout << endl;
   ListAE(0);
   cout << endl;
   ListODE(0);
   cout << endl;
   ListPDE(0);
}

} /* namespace RITA */
