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
using std::pair;

namespace RITA {

data::data(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _size(0), _default_field(0), _verb(1), _configure(config), _cmd(command),
       _theMesh_alloc(0), _theTab_alloc(0), _theGrid_alloc(0), _theFct_alloc(0),
       _theVector_alloc(0), _theMatrix_alloc(0), _u_alloc(0), _theParam_alloc(0),
       _theMesh(nullptr), _theTab(nullptr), _theGrid(nullptr), _theFct(nullptr), 
       _theVector(nullptr), _theMatrix(nullptr)
{
   nb_fields = nb_params = nb_fcts = nb_meshes = nb_tabs = nb_grids = nb_vectors = nb_matrices = 0;
   nb_pde = nb_ode = nb_ae = nb_int = nb_eigen = nb_eq = 0;
   iMesh = iField = iGrid = iParam = iFct = iTab = iAE = iODE = iPDE = iEq = 0;
   _u = nullptr;
   _ret = 0;
   iMesh = iGrid = iField = -1;
   ok = false;
}


data::~data()
{
   if (_theGrid_alloc)
      delete _theGrid;
   if (_theMesh_alloc)
      delete _theMesh;
   if (_theFct_alloc)
      delete _theFct;
   if (_theTab_alloc)
      delete _theTab;
   if (_theVector_alloc)
      delete _theVector;
   if (_theMatrix_alloc)
      delete _theMatrix;
   if (_u_alloc)
      delete _u;
   if (_theParam_alloc)
      delete _theParam;
}


int data::checkField(const string& name) const
{
   auto it = find(Field.begin(),Field.end(),name);
   if (it != Field.end())
      return distance(Field.begin(),it);
   else
      return -1;
}


int data::checkMesh(const string& name) const
{
   auto it = find(mesh_name.begin(),mesh_name.end(),name);
   if (it != mesh_name.end())
      return distance(mesh_name.begin(),it);
   else
      return -1;
}


int data::checkGrid(const string& name) const
{
   auto it = find(grid_name.begin(),grid_name.end(),name);
   if (it != grid_name.end())
      return distance(grid_name.begin(),it);
   else
      return -1;
}


int data::checkFct(const string& name) const
{
   for (auto it=theFct.begin(); it!=theFct.end(); ++it) {
      if ((*it)->name==name)
         return distance(theFct.begin(),it);
   }
   return -1;
}


int data::checkParam(const string& name,
                     double&       value)
{
   if (_cmd->isNumeric(name,value))
      return 0;
   int i = ParamName.find(name)->second - 1;
   if (i>=0)
      value = *theParam[i];
   return i;
}


int data::checkParam(const string& name,
                     int&          value)
{
   if (_cmd->isNumeric(name,value))
      return 0;
   int i = ParamName.find(name)->second - 1;
   if (i>=0)
      value = *theParam[i];
   return i;
}


int data::addAE(const string& name,
                odae*         ae)
{
   theAE.push_back(ae);
   ae_name.push_back(name);
   iAE = nb_ae++, nb_eq++;
   iEq += iAE;
   theODE.push_back(nullptr);
   thePDE.push_back(nullptr);
   eq_type.push_back(ALGEBRAIC_EQ);
   return 0;
}


int data::addODE(const string& name,
                 odae*         ode)
{
   theODE.push_back(ode);
   ode_name.push_back(name);
   iODE = nb_ode++, nb_eq++;
   iEq += iODE;
   theAE.push_back(nullptr);
   thePDE.push_back(nullptr);
   eq_type.push_back(ODE_EQ);
   return 0;
}


int data::addPDE(const string& name,
                 equa*         pde)
{
   thePDE.push_back(pde);
   pde_name.push_back(name);
   iPDE = nb_pde++, nb_eq++;
   iEq += iPDE;
   theAE.push_back(nullptr);
   theODE.push_back(nullptr);
   eq_type.push_back(PDE_EQ);
   return 0;
}


int data::addFunction(const string&         name,
                      const string&         def,
                      const vector<string>& var)
{
   string nm = name;
   if (nm!="") {
      int k = checkFct(name);
      if (k>=0) {
         _rita->msg("rita>data>","Function "+nm+" exists already.");
         return 1;
      }
   }
   if (nm=="")
      nm = "F-" + to_string(nb_fcts+1);
   _theFct = new OFELI::Fct(nm,def,var);
   if (_theFct->check()) {
      _rita->msg("data>",_theFct->getErrorMessage());
      delete _theFct;
      return 1;
   }
   _theFct_alloc = 1;
   theFct.push_back(_theFct);
   nb_fcts++;
   return 0;
}


int data::addMesh(OFELI::Mesh*  ms,
                  const string& name)
{
   theMesh.push_back(ms);
   mesh_name.push_back(name);
   iMesh++;
   return ++nb_meshes;
}


int data::addParam(const string& name,
                   double        value)
{
   iParam = nb_params;
   int k = ParamName.find(name)->second - 1;
//cout<<"par: "<<name<<" ? "<<k<<endl;
   if (k>=0) {
      iParam = k;
      return -1;
   }
   Param.push_back(name);
   ParamName[name] = iParam + 1;
   nb_params++;
//cout<<"nb: "<<nb_params<<endl;
   _theParam = new double;
   *_theParam = value;
   _theParam_alloc = 1;
   theParam.push_back(_theParam);
   return iParam;
}


void data::addField(const string& name)
{
   bool new_field = true;
   int k = checkField(name);
   if (k>=0) {
      iField = k;
      new_field = false;
   }
   _size = 1;
   _nb_dof = 1;
   if (new_field) {
      Field.push_back(name);
      FieldName[name] = ++iField;
      FieldEquation.push_back(0);
   }
   _u = new OFELI::Vect<double>(1);
   _u_alloc = 1;
   _u->setName(name);
   u.push_back(_u);
   if (new_field)
      nb_fields++;
}


void data::addField(const string& name,
                    int           n)
{
   bool new_field = true;
   iField = nb_fields;
   int k = checkField(name);
   if (k>=0) {
      iField = k;
      new_field = false;
   }
   _size = n;
   _nb_dof = 1;
   if (new_field) {
      Field.push_back(name);
      FieldName[name] = iField;
      FieldSizeType.push_back(GIVEN_SIZE);
      FieldEquation.push_back(0);
   }
   _u = new OFELI::Vect<double>(n);
   _u_alloc = 1;
   _u->setName(name);
   u.push_back(_u);
   if (new_field)
      nb_fields++;
}


int data::addMeshField(const string& name,
                       dataSize      s,
                       int           nb_dof)
{
   bool new_field = true;
   _nb_dof = nb_dof;
   _theMesh = theMesh[iMesh];
   if (_theMesh==nullptr) {
      _rita->msg("data>","No mesh data available.");
      _ret = -1;
      return _ret;
   }
   iField = nb_fields;
   int k = checkField(name);
   if (k>=0) {
      iField = k;
      new_field = false;
   }
   if (new_field) {
      Field.push_back(name);
      FieldName[name] = iField;
      FieldType.push_back(PDE_EQ);
   }
   if (s==NODES) {
      if (_theMesh->getNbNodes()==0) {
         _rita->msg("data>","Mesh has no nodes");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==ELEMENTS) {
      if (_theMesh->getNbElements()==0) {
         _rita->msg("data>","Mesh has no elements.");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==SIDES) {
      if (_theMesh->getNbSides()==0) {
         _rita->msg("data>","Mesh has no sides");
         _ret = 1;
         return _ret;
      }
   }
   else if (s==EDGES) {
      if (_theMesh->getNbEdges()==0) {
         _rita->msg("data>","Mesh has no edges");
         _ret = 1;
         return _ret;
      }
   }
   if (new_field) {
      _u = new OFELI::Vect<double>;
      _u_alloc = 1;
      u.push_back(_u);
      nb_fields++;
   }
   _u->setName(name);
   if (s==NODES)
      _u->setMesh(*_theMesh,NODE_DOF,_nb_dof);
   else if (s==ELEMENTS)
      _u->setMesh(*_theMesh,ELEMENT_DOF,_nb_dof);
   else if (s==SIDES)
      _u->setMesh(*_theMesh,SIDE_DOF,_nb_dof);
   else if (s==EDGES)
      _u->setMesh(*_theMesh,EDGE_DOF,_nb_dof);
   return 0;
}


int data::addGridField(const string& name,
                       int           nb_dof)
{
   bool new_field = true;
   _theGrid = theGrid[iGrid];
   if (_theGrid==nullptr) {
      _rita->msg("data>","No grid data available.");
      _ret = -1;
      return _ret;
   }
   iField = nb_fields;
   int k = checkField(name);
   if (k>=0) {
      iField = k;
      new_field = false;
   }
   if (new_field) {
      Field.push_back(name);
      FieldName[name] = iField;
      FieldType.push_back(PDE_EQ);
      _u = new OFELI::Vect<double>(*_theGrid);
      _u_alloc = 1;
      nb_fields++;
   }
   _nb_dof = nb_dof;
   _theGrid->setNbDOF(_nb_dof);
   _u->setName(name);
   u.push_back(_u);
   return 0;
}


int data::run()
{
   int key = 0;
   string td = "", fn="";
   static const vector<string> kw {"grid","mesh","field","tab$ulation","func$tion",
                                   "vect$or","matr$ix","clear","print","save",
                                   "summary","list"};
   *_rita->ofh << "data" << endl;
   while (1) {
      _cmd->readline("rita>data> ");
      if (_nb < 0)
         continue;
      switch (key=_cmd->getKW(kw,_rita->_gkw)) {

         case 100:
         case 101:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands:\n";
            cout << "grid:       Define a grid\n";
            cout << "mesh:       Define a finite element mesh\n";
            cout << "field:      Define a field\n";
            cout << "tabulation: Define a tabulated function\n";
            cout << "function:   Define a function\n";
            cout << "vector:     Define a vector\n";
            cout << "matrix:     Define a matrix\n";
            cout << "print:      Print given field\n";
            cout << "save:       Save given field\n";
            cout << "list:       List a specific data entity\n";
            cout << "summary:    Summary of prescribed data" << endl;
            break;

         case 102:
            _ret = _rita->_configure->run();
            break;

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
	   //            Clear();
            break;

         case   8:
            if (_cmd->setNbArg(1,"Field name to be supplied.",1)) {
               _rita->msg("data>print>","Missing field name.","",1);
               break;
            }
            _ret = _cmd->get(fn);
            if (!_ret)
               print(fn);
            break;

         case   9:
            break;

         case  10:
            Summary();
            break;

         case  11:
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

         case 104:
         case 105:
            _rita->setParam();
            break;

         case 106:
         case 107:
            _ret = 0;
            ok = true;
            *_rita->ofh << "  end" << endl;
            return _ret;

         case -4:
            return 1;

         default:
            _rita->msg("data>","Unknown Command "+_cmd->token(),
                       "Available commands: grid, mesh, field, tabulation, function, vector, matrix,"
                       " summary, list");
            break;
       }
   }
   return 0;
}


void data::getHelp()
{
   cout << "In command data, the following data can be defined:\n";
   cout << "field: Define a field (unknown)\n";
   cout << "function: Define a function\n";
   cout << "mesh: Clear all defined fields\n";
   cout << "summary: Display this summary\n\n";
   cout << "Global commands: \n";
   cout << "help or ?: Display this help\n";
   cout << "set: Set configuration data\n";
   cout << "end or <: Back to higher level\n";
   cout << "exit: End the program\n" << endl;
}


void data::print(const string& s)
{
   int k = checkField(s);
   if (k>=0) {
      cout << u[k];
      return;
   }
   double d=0.;
   k = checkParam(s,d);
   if (k>=0) {
      cout << "Parameter: " << s << " = " << d << endl;
      return;
   }
   int i=0;
   k = checkParam(s,i);
   if (k>=0) {
      cout << "Parameter: " << s << " = " << i << endl;
      return;
   }
/*   k = checkVect(s);
   if (k>=0) {
      cout << theVector[k];
      return;
   }*/
/*   k = checkMatrix(s);
   if (k>=0) {
      cout << theMatrix[k];
      return;
   }*/
   k = checkMesh(s);
   if (k>=0) {
      cout << theMesh[k];
      return;
   }
   k = checkGrid(s);
   if (k>=0) {
      cout << theGrid[k];
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
   _theGrid_alloc = 1;
   theGrid.push_back(_theGrid);
   grid_name.push_back(name);
   iGrid++;
   nb_grids++;
   return 0;
}


int data::setParam()
{
   int nb=0;
   string name="par-"+to_string(nb_params+1);
   static const vector<string> kw {"name","def$ine","set"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<=0) {
      _rita->msg("data>parameter>","Error in command.","Available arguments: name, define, set.");
      return 1;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 0:
            name = _cmd->string_token(0);
            break;

         case 1:
            break;

         case 2:
            break;

         default:
            _rita->msg("data>parameter>","Unknown argument: "+kw[n]);
            return 1;
      }
   }
   if (nb_args>0)
      *_rita->ofh << endl;
   theParam.push_back(_theParam);
   ParamName[name] = iParam;
   nb_params++;
   return 0;
}


int data::setVector()
{
   int size=0, nb=0;
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
   static const vector<string> kw {"name","storage","size","def$ine","set"};
   theMatrix.push_back(_theMatrix);
   matrix_name.push_back(name);
   nb_matrices++;
   return 0;
}


int data::setField()
{
   int nb=0;
   double umin=0., umax=0.;
   string name="u", type="size", arg=" ", mn="";
   static const vector<string> kw {"name","size","mesh","grid","nbdof","type","uniform"};
   static const vector<string> types {"size","nodes","elements","sides","edges"};
   map<string,dataSize> tt {{"size",GIVEN_SIZE}, {"nodes",NODES}, {"elements",ELEMENTS},
                            {"sides",SIDES}, {"edges",EDGES}};
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
   int nb=1, ret=0;
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
      for (const auto& v: theFct) {
         if (v->name==name) {
            _rita->msg("data>function>","Function "+name+" exists already.");
            ret = 1;
         }
      }
      if (ret)
         return 1;
      else {
         if (nb==1)
            var.push_back(vv);
         else {
            for (int i=0; i<nb; ++i)
               var.push_back(vv+to_string(i+1));
         }
         addFunction(name,def,var);
         if (name_ok)
            *_rita->ofh << "  function  name=" << name;
         for (const auto& v: var)
            *_rita->ofh << " var=" << v;
         *_rita->ofh << "  def=" << def << endl;
      }
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
      _theTab_alloc = 1;
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
   if (checkParam(s,v)<0) {
      _rita->msg(msg,"Undefined parameter: "+s);
      return 1;
   }
   return 0;
}


void data::ListParams(int opt)
{
   if (opt && !nb_params) {
      cout << "No defined parameters." << endl;
      return;
   }
   cout << "Number of defined parameters: " << nb_params << endl;
   for (int i=0; i<nb_params; ++i)
      cout << Param[i] << " = " << *theParam[i] << endl;
}


void data::ListFields(int opt)
{
   if (opt && !nb_fields) {
      cout << "No defined fields." << endl;
      return;
   }
   cout << "Number of fields: " << nb_fields << endl;
   for (int i=0; i<nb_fields; ++i) {
      if (FieldSizeType[i]==GIVEN_SIZE)
         cout << "Field: " << Field[i] << ", Size: " << u[i]->size() << endl;
      else
         cout << "Field: " << Field[i] << ", Number of degrees of freedom: " << nb_dof[i] << endl;
   }
}


void data::ListFunctions(int opt)
{
   if (opt && !nb_fcts) {
      cout << "No defined functions." << endl;
      return;
   }
   cout << "Number of functions: " << nb_fcts << endl;
   for (const auto& v: theFct) {
      cout << "Function: " << v->name << ", Variable(s): ";
      for (size_t j=0; j<v->nb_var-1; ++j)
         cout << v->var[j] << ",";
      cout << v->var[v->nb_var-1] << ", Definition: " << v->expr << endl;
   }
}


void data::ListTabs(int opt)
{
   if (opt && !nb_tabs) {
      cout << "No defined tabulations." << endl;
      return;
   }
   cout << "Number of tabulations: " << nb_tabs << endl;
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
}


void data::ListGrids(int opt)
{
   if (opt && !nb_grids) {
      cout << "No defined grids." << endl;
      return;
   }
   cout << "Number of grids: " << nb_grids << endl;
   int nb = 0;
   for (const auto& v: theGrid) {
      cout << "Grid No.            " << ++nb << endl;
      cout << "Grid name:          " << grid_name[nb-1] << endl;
      cout << "Space dimension:    " << v->getDim() << endl;
      if (v->getDim()==1) {
         cout << "Domain:                    (" << v->getX(1) << ","
              << v->getX(v->getNx()+1) << ")" << endl;
         cout << "Number of x-intervals:    " << v->getNx() << endl;
      }
      else if (v->getDim()==2) {
         cout << "Domain:                    (" << v->getX(1) << ","
              << v->getX(v->getNx()+1) << ")x(" << v->getY(1) << ","
              << v->getY(v->getNy()+1) << ")" << endl;
         cout << "Number of x-intervals:    " << v->getNx() << endl;
         cout << "Number of y-intervals:    " << v->getNy() << endl;
      }
      else if (v->getDim()==3) {
         cout << "Domain:                    (" << v->getX(1) << ","
              << v->getX(v->getNx()+1) << ")x(" << v->getY(1) << ","
              << v->getY(v->getNy()+1) << ")x(" << v->getZ(1) << ","
              << v->getZ(v->getNz()+1) << ")" << endl;
         cout << "Number of x-intervals:    " << v->getNx() << endl;
         cout << "Number of y-intervals:    " << v->getNy() << endl;
         cout << "Number of z-intervals:    " << v->getNz() << endl;
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
   int nb = 0;
   for (const auto& v: theMesh) {
      cout << "Mesh No.            " << nb++ << endl;
      cout << "Mesh name:          " << mesh_name[nb-1] << endl;
      cout << "Space dimension:    " << v->getDim() << endl;
      cout << "Number of nodes:    " << v->getNbNodes() << endl;
      cout << "Number of elements: " << v->getNbElements() << endl;
      cout << "Number of sides:    " << v->getNbSides() << endl;
   }
}


void data::ListODE(int opt)
{
   if (opt && !nb_ode) {
      cout << "No defined ordinary differential equations." << endl;
      return;
   }
   cout << "Number of ordinary differential equations: " << nb_ode << endl;
   for (int i=0; i<nb_ode; ++i) {
   }
}


void data::ListPDE(int opt)
{
   if (opt && !nb_pde) {
      cout << "No defined partial differential equations." << endl;
      return;
   }
   cout << "Number of partial differential equations: " << nb_pde << endl;
   for (int i=0; i<nb_pde; ++i) {
   }
}


void data::ListAE(int opt)
{
   if (opt && !nb_ae) {
      cout << "No defined algebraic equations." << endl;
      return;
   }
   cout << "Number of algebraic equations: " << nb_ae << endl;
   for (int i=0; i<nb_ae; ++i) {
   }
}


void data::Summary()
{
   if (nb_fields+nb_params+nb_fields+nb_fcts+nb_tabs+nb_grids+nb_meshes+
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
