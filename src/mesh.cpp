/*==============================================================================

                                 r  i  t  a

            An environment for Modelling and Numerical Simulation

  ==============================================================================

    Copyright (C) 2021 - 2023 Rachid Touzani

    This file is part of rita.

    rita is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or

    rita is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

  ==============================================================================
 
                      Implementation of class 'mesh'

  ==============================================================================*/

#include <iostream>
#include <stdlib.h>
#include "mesh.h"
#include "mesh/saveMesh.h"
#include "rita.h"
#include "data.h"
#include "calc.h"

//#ifdef USE_GMSH
#include <gmsh.h>
//#endif

namespace RITA {

mesh::mesh(rita*      r,
           cmd*       command,
           configure* config)
     : _rita(r), _saved(false), _generated(false), _geo(false), _ret(0), _generator(0),
       _theMesh(nullptr), _theDomain(nullptr), _nb_Ccontour(0), _nb_Scontour(0), _nb_Vcontour(0),
       _nb_sub_domain(0), _nb_point(0), _nb_curve(0), _nb_surface(0), _nb_volume(0), 
       _configure(config), _cmd(command)
{
}


mesh::~mesh()
{
   if (_theMesh!=nullptr)
      delete _theMesh, _theMesh = nullptr;
   if (_theDomain!=nullptr)
      delete _theDomain, _theDomain = nullptr;
}


int mesh::run()
{
   string fn="";
   _nb_dof = 1;
   _data = _rita->_data;
   static const vector<string> kw {"1d","rect$angle","cube","point","curve","surface","volume","contour",
                                   "code","gen$erate","nbdof","list","plot","clear","save","read"};
#ifndef USE_GMSH
   _theDomain = new OFELI::Domain;
#endif
   while (1) {
      int nb_args = _cmd->readline(sPrompt+"mesh> ");
      if (nb_args < 0)
         continue;
      _key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
      if (_key>=200) {
         _data->setDataExt(_key);
         continue;
      }
      switch (_key) {

         case   0:
            set1D();
            break;

         case   1:
            setRectangle();
            break;

         case   2:
            setCube();
            break;

         case   3:
            setPoint();
            break;

         case   4:
            setCurve();
            break;

         case   5:
            setSurface();
            break;

         case   6:
            setVolume();
            break;

         case   7:
            setContour();
            break;

         case   8:
            setCode();
            break;

         case   9:
            Generate();
            break;

         case  10:
            setNbDOF();
            break;

         case  11:
            List();
            break;

         case  12:
            Plot();
            break;

         case  13:
            Clear();
            break;

         case  14:
            Save();
            break;

         case  15:
            Read();
            break;

         case 100:
         case 101:
            _cmd->setNbArg(0);
            cout << "\nAvailable commands:" << endl;
            cout << "1d        : Data to generate a 1-D mesh" << endl;
            cout << "rectangle : Data to mesh a rectangle" << endl;
            cout << "cube      : Data to mesh a cube (parallelepiped)" << endl;
            cout << "point     : Define a point" << endl;
            cout << "curve     : Define a curve" << endl;
            cout << "surface   : Define a surface" << endl;
            cout << "volume    : Define a volume" << endl;
            cout << "contour   : Define a contour as a sequence of curves or surfaces" << endl;
//   cout << "subdomain : Define a subdomain" << endl;
            cout << "code      : Set code for points, lines, surfaces, volumes" << endl;
            cout << "generate  : Generate mesh of a polygon" << endl;
            cout << "list      : List mesh data" << endl;
            cout << "plot      : Plot mesh" << endl;
            cout << "clear     : Clear mesh" << endl;
            cout << "read      : Read mesh from file" << endl;
            cout << "save      : Save mesh in file" << endl;
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
            *_rita->ofh << "  end" << endl;
             _ret = 0;
            return _ret;

         default:
            _ret = _rita->_calc->run();
            break;
      }
   }
   return 1;
}


void mesh::List()
{
   if (_generated==0) {
      cout << "Geometry data\n" << endl;
      cout << "Space dimension:    " << _theMesh->getDim() << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl;
   }
   else {
      cout << "Mesh data\n" << endl;
      cout << "Space dimension:    " << _theMesh->getDim() << endl;
      cout << "Number of nodes:    " << _theMesh->getNbNodes() << endl;
      cout << "Number of elements: " << _theMesh->getNbElements() << endl;
      cout << "Number of sides:    " << _theMesh->getNbSides() << endl;
   }
}


void mesh::setConfigure()
{
   if (_verb)
      cout << "Setting configuration parameters ..." << endl;
   _configure->setVerbose(_verb);
   _configure->run();
   _verb = _configure->getVerbose();
}


void mesh::set1D()
{
   string fn="";
   _dim = 1;
   _nb_dof = 1;
   double xmin=0., xmax=1.;
   int nb=0, ret=0, ne=10, cmin=0, cmax=0;
   _mesh_file = "rita-1d.m";
   _saved = false;
   _ret = 0;
   _generated = false;
   static const vector<string> kw {"domain","ne","codes","nbdof","save"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 100:
         case 101:
            cout << "1d [domain=m,M] [ne=n] [codes=c1,c2] [nbdof=d]\n"
                    "m, M: Extremal coordinates of interval ends. The default values are 0., 1.\n"
                    "n: Number of elements in the interval. Its default value is 10.\n"
                    "c1, c2: Codes associated to the first and the last node. These integer values\n"
                    "        are necessary to enforce boundary conditions. A code 0 (Default value) means\n"
                    "        no condition to prescribe.\n"
                    "d: Number of degrees of freedom associated to any generated node. Default value is 1." << endl << endl;
            _ret = 0;
            return;

         case   0:
            if (nb==1) {
               xmin = 0.;
               _ret = _data->getPar(0,"mesh>1d>",xmax);
               if (_ret)
                  return;
            }
            else if (nb==2) {
               _ret  = _data->getPar(0,"mesh>1d>",xmin);
               _ret += _data->getPar(1,"mesh>1d>",xmax);
               if (_ret)
                  return;
            }
            else {
               _rita->msg("mesh>1d>","This argument requires 1 or 2 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   1:
            _ret  = _data->getPar(0,"mesh>1d>",ne);
            if (_ret)
               return;
            break;

         case   2:
            if (nb==1)
               cmin = cmax = _cmd->int_token(0);
            else if (nb==2) {
               cmin = _cmd->int_token(0);
               cmax = _cmd->int_token(1);
            }
            else {
               _rita->msg("mesh>1d>","This argument requires 1 or 2 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   3:
            _nb_dof = _cmd->int_token(0);
            break;

         case   4:
            _mesh_file = _cmd->string_token(0);
            break;

         default:
            _rita->msg("mesh>1d>","Unknown argument: "+kw[n]);
            return;
      }
   }
   if (nb_args>0) {
      if (xmax<=xmin) {
         _rita->msg("mesh>1d>","Error in values of xmin and xmax: "+to_string(xmin)+", "+to_string(xmax));
         _ret = 1;
         return;
      }
      if (ne<2) {
         _rita->msg("mesh>1d>","Number of elements must be > 1");
         _ret = 1;
         return;
      }
      _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
      _theMesh->removeImposedDOF();
      _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
      _saved = false;
      _generator = 1;
      _generated = true;
      _data->NameMesh.push_back("M"+to_string(_data->theMesh.size()));
      _data->theMesh.push_back(_theMesh);
      *_rita->ofh << "  1d  domain=" << xmin << "," << xmax << "  codes=" << cmin
                  << "," << cmax << "  ne=" << ne << "  nbdof=" << _nb_dof << endl;
      if (_verb)
         cout << "1-D mesh complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
   }

   else {
      if (_verb) {
         cout << "Default interval: (0,1)\n";
         cout << "Default number of elements: 10\n";
         cout << "Default codes for end nodes: 0 0\n";
         cout << "Default nb of dof: 1" << endl;
      }
      *_rita->ofh << "  1d" << endl;
      while (1) {
         if (_cmd->readline(sPrompt+"mesh>1d> ")<0)
            continue;
         _key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
         if (_key>=200) {
            _data->setDataExt(_key);
            continue;
         }
         switch (_key) {

            case   0:
               if (_verb>1)
                  cout << "Setting interval bounds ..." << endl;
               if (_cmd->setNbArg(2,"Give xmin and xmax.")) {
                  _rita->msg("mesh>1d>domain>","Missing values of xmin and xmax.");
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>1d>domain>",xmin);
               _ret += _data->getPar(-1,"mesh>1d>domain>",xmin);
               if (_ret)
                  return;
               if (xmax<=xmin) {
                  _rita->msg("mesh>1d>domain>","Value of xmax: "+to_string(xmax)+" must be > xmin: "+to_string(xmin));
                  break;
               }
               if (!_ret)
                  *_rita->ofh << "    domain " << xmin << " " << xmax << endl;
               break;

            case   1:
               if (_verb>1)
                  cout << "Setting interval subdivision ..." << endl;
               if (_cmd->setNbArg(1,"Give number of elements.")) {
                  _rita->msg("mesh>1d>ne>","Missing number of elements.","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>1d>ne>",ne);
               if (_ret)
                  return;
               *_rita->ofh << "    ne " << ne << endl;
               break;

            case   2:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               if (_cmd->setNbArg(2,"Give cmin and cmax.")) {
                  _rita->msg("mesh>1d>codes>","Missing codes for end nodes.","",1);
                  break;
               }
               ret  = _cmd->get(cmin);
               ret += _cmd->get(cmax);
               if (!ret)
                  *_rita->ofh << "    codes " << cmin << " " << cmax << endl;
               break;

            case   3:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               if (_cmd->setNbArg(1,"Give number of dof per node.")) {
                  _rita->msg("mesh>1d>nbdof>","Missing number of dof.","",1);
                  break;
               }
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_rita->ofh << "    nbdof " << _nb_dof << endl;
               break;

            case   4:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  _rita->msg("mesh>1d>save>","You are trying to delete an existing mesh.",
                             "retype command 'save' to confirm.");
                  _saved = false;
                  break;
               }
               _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
               _theMesh->removeImposedDOF();
               _saved = false;
               _generator = 1;
               _generated = true;
               if (_cmd->getNbArgs()>0)
                  ret = _cmd->get(_mesh_file);
               if (!ret)
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
               break;

            case 100:
            case 101:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands:\n";
               cout << "domain   : Enter xmin and xmax\n";
               cout << "ne       : Enter number of elements to generate\n";
               cout << "codes    : Codes to associate to end nodes\n";
               cout << "nbdof    : Number of dof per node" << endl;
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
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               if (!_saved) {
                  if (_theMesh!=nullptr)
                     delete _theMesh, _theMesh = nullptr;
                  _theMesh = new OFELI::Mesh(xmin,xmax,ne,cmin,cmax,1,size_t(_nb_dof));
                  _theMesh->removeImposedDOF();
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
               }
               *_rita->ofh << "  end" << endl;
               if (_verb)
                  cout << "1-D mesh complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
               _ret = 90;
               return;

            case -2:
               break;

            default:
               _ret = _rita->_calc->run();
               break;
         }
      }
   }
}


void mesh::setRectangle()
{
   string fn="";
   double xmin=0., xmax=1., ymin=0., ymax=1.;
   int c[4] = {0,0,0,0}, cv[4] = {0,0,0,0};
   int nb=0, nx=10, ny=10, ret=0, ret1=0, ret2=0, ret3=0, ret4=0;
   _nb_dof = 1;
   _dim = 2;
   _ret = 0;
   _saved = false;
   const static string H = "rectangle [min=mx,my] [max=Mx,My] [ne=nx,ny]  [codes=c1,c2,c3,c4]\n"
                           "          [nbdof=d]\n"
                           "mx, my: coordinates of the lower left corner of the rectangle.\n"
                           "        The default values are 0., 0.\n"
                           "Mx, My: Coordinates of the upper right corner of the rectangle.\n"
                           "        The default values are 1., 1.\n"
                           "nx, ny: Number of elements in x and y direction respectively.\n"
                           "        Their default value is 10.\n"
                           "c1, c2, c3, c4: Codes associated to the nodes generated on the lines y=my,\n"
                           "                x=Mx, y=My, x=mx respectively. These integer values are necessary\n"
                           "                to enforce boundary conditions. A code 0 (Default value) means no\n"
                           "                condition to prescribe.\n"
                           "d: Number of degrees of freedom associated to any generated node. Default value is 1.\n";
   _mesh_file = "rita-rectangle.m";
   static const vector<string> kw {"min","max","ne","codes","nbdof"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case    0:
            if (nb==1) {
               _ret = _data->getPar(0,"mesh>rectangle>",xmin);
               if (_ret)
                  return;
               ymin = xmin;
            }
            if (nb==2) {
               _ret  = _data->getPar(0,"mesh>rectangle>",xmin);
               _ret += _data->getPar(1,"mesh>rectangle>",ymin);
               if (_ret)
                  return;
            }
            else {
               _rita->msg("mesh>rectangle>","This argument requires 1 or 2 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   1:
            if (nb==1) {
               _ret  = _data->getPar(0,"mesh>rectangle>",xmax);
               if (_ret)
                  return;
               ymax = xmax;
            }
            if (nb==2) {
               _ret  = _data->getPar(0,"mesh>rectangle>",xmax);
               _ret += _data->getPar(1,"mesh>rectangle>",ymax);
               if (_ret)
                  return;
            }
            else {
               _rita->msg("mesh>rectangle>","This argument requires 1 or 2 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   2:
            if (nb==1) {
               _ret  = _data->getPar(0,"mesh>rectangle>",nx);
               ny = nx;
            }
            else {
               _ret  = _data->getPar(0,"mesh>rectangle>",nx);
               _ret += _data->getPar(1,"mesh>rectangle>",ny);
               if (_ret)
                  return;
            }
            break;

         case   3:
            if (nb==1)
               c[0] = c[1] = c[2] = c[3] = _cmd->int_token(0);
            else if (nb==4) {
               c[0] = _cmd->int_token(0);
               c[1] = _cmd->int_token(1);
               c[2] = _cmd->int_token(2);
               c[3] = _cmd->int_token(3);
            }
            else {
               _rita->msg("mesh>rectangle>","This argument requires 1 or 4 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   4:
            _nb_dof = _cmd->int_token(0);
            break;

         case  100:
         case  101:
            cout << H << endl;
            _ret = 1;
            return;

         default:
            _rita->msg("mesh>rectangle>","Unknown argument: "+kw[n]);
            return;
      }
   }
   if (nb_args>0) {
      if (xmax<=xmin) {
         _rita->msg("mesh>rectangle>","xmax: "+to_string(xmax)+" must be > than xmin: "+ to_string(xmin));
         return;
      }
      if (ymax<=ymin) {
         _rita->msg("mesh>rectangle>","ymax: "+to_string(ymax)+" must be > ymin: "+to_string(ymin));
         return;
      }
      _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
      _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
      if (_verb)
         cout << "2-D mesh complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
      _generated = true;
      _generator = 1;
      *_rita->ofh << "  rectangle min=" << xmin << "," << ymin << " max=" << xmax << "," << ymax;
      *_rita->ofh << "  ne=" << nx << "," << ny << "  codes=" << c[0] << "," << c[1] << "," << c[2];
      *_rita->ofh << "," << c[3] << "  nbdof=" << _nb_dof << endl;
      _ret = 0;
      return;
   }

   else {
      if (_verb) {
         cout << "Default rectangle: (0,1)*(0,1)\n";
         cout << "Default node codes: 0\n";
         cout << "Default subdivision: 10*10\n";
         cout << "Default nb of dof: 1" << endl;
      }
      *_rita->ofh << "  rectangle" << endl;

      while (1) {
         if (_cmd->readline(sPrompt+"mesh>rectangle> ")<0)
            continue;
         _key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
         if (_key>=200) {
            _data->setDataExt(_key);
            continue;
         }
         switch (_key) {

            case   0:
               if (_verb>1)
                  cout << "Setting xmin and ymin ..." << endl;
               if (_cmd->setNbArg(2,"Give values of xmin and ymin.")) {
                  _rita->msg("mesh>rectangle>min>","Missing xmin and ymin values","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>rectangle>min>",xmin);
               _ret  = _data->getPar(-1,"mesh>rectangle>min>",ymin);
               if (_ret)
                  break;
               *_rita->ofh << "    min " << xmin << " " << ymin << endl;
               break;

            case   1:
               if (_verb>1)
                  cout << "Setting xmax and ymax ..." << endl;
               if (_cmd->setNbArg(2,"Give values of xmax and ymax.")) {
                  _rita->msg("mesh>rectangle>max>","Missing xmax and ymax values","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>rectangle>max>",xmax);
               _ret += _data->getPar(-1,"mesh>rectangle>max>",ymax);
               if (_ret)
                  break;
               *_rita->ofh << "    max " << xmax << " " << ymax << endl;
               break;

            case   2:
               if (_verb>1)
                  cout << "Setting mesh subdivisions ..." << endl;
               if (_cmd->setNbArg(2,"Give subdivision in the x and y-directions.")) {
                  _rita->msg("mesh>rectangle>ne>","Missing subdivisions in x and y directions","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>rectangle>ne>",nx);
               _ret += _data->getPar(-1,"mesh>rectangle>ne>",ny);
               if (_ret)
                  return;
               *_rita->ofh << "    ne " << nx << " " << ny << endl;
               break;

            case   3:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               if (_cmd->setNbArg(4,"Code to assign on the line y=ymin.")) {
                  _rita->msg("mesh>rectangle>codes>","Missing code to assign on the line y=ymin.","",1);
                  break;
               }
               ret1 = _cmd->get(c[0]);
               ret2 = _cmd->get(c[1]);
               ret3 = _cmd->get(c[2]);
               ret4 = _cmd->get(c[3]);
               if (!ret1 && !ret2 && !ret3 && !ret4) {
                  *_rita->ofh << "    codes " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << endl;
                  if (cv[0]==0)
                     cv[0] = c[0];
                  if (cv[1]==0)
                     cv[1] = c[0];
                  if (cv[1]==0)
                     cv[1] = c[1];
                  if (cv[2]==0)
                     cv[2] = c[1];
                  if (cv[2]==0)
                     cv[2] = c[2];
                  if (cv[3]==0)
                     cv[3] = c[2];
                  if (cv[3]==0)
                     cv[3] = c[3];
                  if (cv[0]==0)
                     cv[0] = c[3];
               }
               break;

            case   4:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               if (_cmd->setNbArg(1,"Give number of dof per node.")) {
                  _rita->msg("mesh>rectangle>nbdof>","Missing number of dof.","",1);
                  break;
               }
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_rita->ofh << "    nbdof " << _nb_dof << endl;
               break;

            case   5:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  _rita->msg("mesh>rectangle>save>","Trying to delete a created mesh.","retype command 'save' to confirm.");
                  _saved = false;
                  break;
               }
               _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
               _saved = true;
               _generator = 2;
               if (_cmd->getNbArgs()>0)
                  ret = _cmd->get(_mesh_file);
               if (!ret) {
                  _theMesh->put(_mesh_file);
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
                  _saved = true;
                  cout << "2-D mesh complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
               }
               break;

            case 100:
            case 101:
               _cmd->setNbArg(0);
               cout << H << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               _ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 106:
            case 107:
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               *_rita->ofh << "    end" << endl;
               if (!_saved) {
                  _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,nx,ny,c[3],c[1],c[0],c[2],TRIANGLE,size_t(_nb_dof));
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
                  _generator = 2;
                  _saved = true;
               }
               if (_verb && !_saved)
                  cout << "Mesh 'rectangle' complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
               _ret = 90;
               return;

            case -2:
               break;

            default:
               _ret = _rita->_calc->run();
               break;
         }
      }
   }
   _ret = 0;
   return;
}


void mesh::setCube()
{
   string fn="";
   double xmin=0., xmax=1., ymin=0., ymax=1., zmin=0., zmax=1.;
   int nb=0, ret=0, nx=10, ny=10, nz=10, cxmin=0, cxmax=0, cymin=0, cymax=0, czmin=0, czmax=0;
   _nb_dof = 1;
   _dim = 3;
   const static string H = "cube [min=mx,my,mz] [max=Mx,My,Mz] [ne=nx,ny,nz] [codes=cxm,cxM,cym,cyM,czm,czMs]\n"
                           "     [nbdof=d]\n"
                           "mx, my, mz: Minimal coordinates in each direction.\n"
                           "            The default values are 0., 0., 0.\n"
                           "Mx, My, Mz: Maximal coordinates in each direction.\n"
                           "            The default values are 1., 1.,1.\n"
                           "nx, ny, nz: Number of elements in x, y and z direction respectively.\n"
                           "            Their default value is 10.\n"
                           "cxm, cxM, cym, cyM, czm, czM: Codes associated to the nodes generated on the face x=mx\n" 
                           "x=Mx, y=my, y=My, z=mz, z=Mz respectively.\n"
                           "These integer values are necessary to enforce boundary\n"
                           "conditions. A code 0 (Default value) means no condition to prescribe.\n"
                           "d: Number of degrees of freedom associated to any generated node. Default value is 1.\n";
   *_rita->ofh << "  cube" << endl;
   _saved = false;
   _ret = 0;
   if (_verb) {
      cout << "Default cube: (0,1)*(0,1)*(0,1)\n";
      cout << "Default node codes: 0\n";
      cout << "Default subdivision: 10*10*10\n";
      cout << "Default nb of dof: 1" << endl;
   }
   static const vector<string> kw {"min","max","ne","codes","nbdof"};

   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case 100:
         case 101:
            cout << H << endl;
            _ret = 0;
           return;   

         case   0:
            if (nb==1) {
               _ret  = _data->getPar(0,"mesh>cube>",xmin);
               if (_ret)
                  return;
               ymin = zmin = xmin;
            }
            if (nb==3) {
               _ret  = _data->getPar(0,"mesh>cube>",xmin);
               _ret += _data->getPar(1,"mesh>cube>",ymin);
               _ret += _data->getPar(2,"mesh>cube>",zmin);
               if (_ret)
                  return;
            }
            else {
               _rita->msg("mesh>cube>","This argument requires 1 or 3 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   1:
            if (nb==1) {
               _ret  = _data->getPar(0,"mesh>cube>",xmax);
               if (_ret)
                  return;
               ymax = zmax = xmax;
            }
            if (nb==3) {
               _ret  = _data->getPar(0,"mesh>cube>",xmax);
               _ret += _data->getPar(1,"mesh>cube>",ymax);
               _ret += _data->getPar(2,"mesh>cube>",zmax);
               if (_ret)
                  return;
            }
            else {
               _rita->msg("mesh>cube>","This argument requires 1 or 3 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   2:
            if (nb==1) {
               _ret  = _data->getPar(0,"mesh>cube>",nx);
               if (_ret)
                  return;
               ny = nz = nx;
            } 
            else {
               _ret  = _data->getPar(0,"mesh>cube>",nx);
               _ret += _data->getPar(1,"mesh>cube>",ny);
               _ret += _data->getPar(2,"mesh>cube>",nz);
               if (_ret)
                  return;
            }
            break;

         case   3:
            if (nb==1)
               cxmin = cxmax = cymin = cymax = czmin = czmax = _cmd->int_token(0);
            else if (nb==6) {
               cxmin = _cmd->int_token(0);
               cxmax = _cmd->int_token(1);
               cymin = _cmd->int_token(2);
               cymax = _cmd->int_token(3);
               czmin = _cmd->int_token(4);
               czmax = _cmd->int_token(5);
            }
            else {
               _rita->msg("mesh>cube>","This argument requires 1 or 6 parameters.");
               _ret = 1;
               return;
            }
            break;

         case   4:
            _nb_dof = _cmd->int_token(0);
            break;

         case   5:
            _mesh_file = _cmd->string_token(0);
            break;

         default:
            _rita->msg("mesh>cube>","Unknown argument: "+kw[n]);
            return;
      }
   }
   if (nb_args>0) {
      if (xmax<=xmin) {
         _rita->msg("mesh>cube>","Value of xmax: "+to_string(xmax)+" must be > xmin: "+to_string(xmin));
         return;
      }
      if (ymax<=ymin) {
         _rita->msg("mesh>cube>","Value of ymax: "+to_string(ymax)+" must be > ymin: "+to_string(ymin));
         return;
      }
      if (zmax<=zmin) {
         _rita->msg("mesh>cube>","Value of zmax: "+to_string(zmax)+" must be > zmin: "+to_string(zmin));
         return;
      }
      if (!_saved) {
         _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,cxmax,cymin,
                                    cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
         _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
         if (_verb)
            cout << "3-D mesh complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
         _saved = _generated = true;
      }
      _generator = 1;
      _generated = true;
      *_rita->ofh << "  cube min=" << xmin << "," << ymin << "," << zmin << " max=" << xmax << "," << ymax << "," << zmax;
      *_rita->ofh << " ne=" << nx << "," << ny << "," << nz << " codes=" << cxmin << "," << cxmax << "," << cymin << ",";
      *_rita->ofh << cymax << "," << czmin << "," << czmax << " nbdof=" << _nb_dof << endl;
      _ret = 0;
      return;
   }

   else {
      while (1) {
         if (_cmd->readline(sPrompt+"mesh>cube> ")<0)
            continue;
         _key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
         if (_key>=200) {
            _data->setDataExt(_key);
            continue;
         }
         switch (_key) {

            case   0:
               if (_verb>1)
                  cout << "Setting xmin, ymin and zmin ..." << endl;
               if (_cmd->setNbArg(3,"Give values of xmin, ymin and zmin.")) {
                  _rita->msg("mesh>cube>min>","Missing values of xmin, ymin and zmin.","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>cube>min>",xmin);
               _ret += _data->getPar(-1,"mesh>cube>min>",ymin);
               _ret += _data->getPar(-1,"mesh>cube>min>",zmin);
               if (!ret)
                  *_rita->ofh << "    min " << xmin << " " << ymin << " " << zmin << endl;
               break;

            case   1:
               if (_verb>1)
               cout << "Setting xmax, ymax and zmax ..." << endl;
               if (_cmd->setNbArg(3,"Give values of xmax, ymax and zmax.")) {
                  _rita->msg("mesh>cube>min>","Missing values of xmax, ymax and zmax.","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>cube>max>",xmax);
               _ret += _data->getPar(-1,"mesh>cube>max>",ymax);
               _ret += _data->getPar(-1,"mesh>cube>max>",zmax);
               if (!_ret)  
                  *_rita->ofh << "    max " << xmax << " " << ymax << " " << zmax << endl;
               break;

            case   2:
               if (_verb>1)
                  cout << "Setting mesh subdivisions ..." << endl;
               if (_cmd->setNbArg(3,"Give subdivision in the x-, y- and z-directions.")) {
                  _rita->msg("mesh>cube>ne>","Missing subvdivisions in x, y and z directions","",1);
                  break;
               }
               _ret  = _data->getPar(-1,"mesh>cube>ne>",nx);
               _ret += _data->getPar(-1,"mesh>cube>ne>",ny);
               _ret += _data->getPar(-1,"mesh>cube>ne>",nz);
               if (!ret)
                  *_rita->ofh << "    ne " << nx << " " << ny << " " << nz << endl;
               break;
   
            case   3:
               if (_verb>1)
                  cout << "Setting boundary codes ..." << endl;
               if (_cmd->setNbArg(6,"Codes to assign to faces.")) {
                  _rita->msg("mesh>cube>codes>","Missing codes to assign to faces.","",1);
                  break;
               }
               ret  = _cmd->get(cxmin);
               ret += _cmd->get(cxmax);
               ret += _cmd->get(cymin);
               ret += _cmd->get(cymax);
               ret += _cmd->get(czmin);
               ret += _cmd->get(czmax);
               if (!ret)
                  *_rita->ofh << "    codes " << cxmin << " " << cxmax << " " << cymin << " " << cymax
                             << " " << czmin << " " << czmax << endl;
               break;

            case   4:
               if (_verb>1)
                  cout << "Setting number of degrees of freedom ..." << endl;
               if (_cmd->setNbArg(1,"Give number of dof per node.")) {
                  _rita->msg("mesh>cube>nbdof>","Missing number of dof.");
                  break;
               }
               ret = _cmd->get(_nb_dof);
               if (!ret)
                  *_rita->ofh << "    nbdof " << _nb_dof << endl;
               break;

            case   5:
               if (_verb>1)
                  cout << "Saving mesh in file ..." << endl;
               if (_saved) {
                  _rita->msg("mesh>cube>save>","Trying to delete a created mesh.",
                             "You are trying to delete an existing mesh.\n"
                             "retype command 'save' to confirm.");
                  _saved = false;
                  break;
               }
               _mesh_file = "rita-cube.m";
               if (_cmd->getNbArgs()>0)
                  _cmd->get(_mesh_file);
               _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,cxmax,cymin,
                                          cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
               _theMesh->put(_mesh_file);
               _data->NameMesh.push_back("M"+to_string(_data->theMesh.size()));
               _data->theMesh.push_back(_theMesh);
               _generator = 3;
               _generated = true;
               cout << "Mesh of cube complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
               _saved = true;
               break;

            case 100:
            case 101:
               _cmd->setNbArg(0);
               cout << H << endl;
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
               if (_verb>1)
                  cout << "Getting back to higher level ..." << endl;
               if (!_saved) {
                  _theMesh = new OFELI::Mesh(xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz,cxmin,
                                             cxmax,cymin,cymax,czmin,czmax,HEXAHEDRON,size_t(_nb_dof));
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
                  _saved = true;
               }
               *_rita->ofh << "    end" << endl;
               if (_verb && !_saved)
                  cout << "Mesh 'cube' complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
               _ret = 90;
               return;

            case -2:
               break;

            default:
               _ret = _rita->_calc->run();
               break;
         }
      }
   }
   _ret = 0;
   return;
}


void mesh::setCode()
{
   int nb=0, c=0, np=0, nc=0, ns=0, nv=0;
   int c_ok=0, points_ok=0, curves_ok=0, surfaces_ok=0, volumes_ok=0;
   vector<int> points, curves, surfaces, volumes;
   const static string H = "code value=v [points=p1,p2,...] [curves=c1,c2,...]  [surfaces=s1,s2,...] [volumes=v1,v2,...]\n\n"
                           "v: Code value (integer) to assign to an entity\n"
                           "p1, p2, ...: Points to which the code v is assigned\n"
                           "c1, c2, ...: Curves to which the code v is assigned\n"
                           "s1, s2, ...: Surfaces to which the code v is assigned\n"
                           "v1, v2, ...: Volumes to which the code v is assigned.";
   _ret = 0;
   if (_generator>0 && _generator<4) {
      _rita->msg("mesh>code>","Keyword not allowed for generated mesh");
      return;
   }
   _generator = 10;
   if (_verb>1) {
      cout << "Default codes for generated nodes:    0" << endl;
      cout << "Default codes for generated elements: 1" << endl;
   }

   static const vector<string> kw {"value","points","curves","surfaces","volumes"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<1) {
      _rita->msg("mesh>codes>"," No argument for command.",H);
      _ret = 1;
      return;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {
  
         case   0:
            c = _cmd->int_token(0);
            if (c<0)
               c = MY_RANDOM - c;
            c_ok++;
            break;

         case   1:
            np = nb;
            for (int j=0; j<np; ++j)
               points.push_back(_cmd->int_token(j));
            points_ok++;
            break;

         case   2:
            nc = nb;
            for (int j=0; j<nc; ++j)
               curves.push_back(_cmd->int_token(j));
            curves_ok++;
            break;

         case   3:
            ns = nb;
            for (int j=0; j<ns; ++j)
               surfaces.push_back(_cmd->int_token(j));
            surfaces_ok++;
            break;

         case   4:
            nv = nb;
            for (int j=0; j<nv; ++j)
               volumes.push_back(_cmd->int_token(j));
            volumes_ok++;
            break;

         case 100:
         case 101:
            cout << H << endl;
            _ret = 1;
            return;

         default:
            _rita->msg("mesh>code>","Unknown argument: "+_kw[n]);
            _ret = 1;
	    return;
      }
   }
   if (!c_ok) {
      _rita->msg("mesh>code>","Missing code value.");
      _ret = 1;
      return;
   }
   if (points_ok>1 || curves_ok>1 || surfaces_ok>1 || volumes_ok>1) {
      _rita->msg("mesh>code>","Each entity must be given only once.");
      _ret = 1;
      return;
   }
   if (np+nc+ns+nv==0) {
      _rita->msg("mesh>code>","At least one entity must be given.");
      _ret = 1;
      return;
   }
   *_rita->ofh << "  code  value=" << c;
   if (points_ok) {
      *_rita->ofh << "  points=";
      for (int j=0; j<np; ++j)
         _Pcode[c].l.push_back(points[j]), _Pcode[c].nb++;
      for (int j=0; j<np-1; ++j)
         *_rita->ofh << points[j] << ",";
      *_rita->ofh << points[np-1] << endl;
   }
   if (curves_ok) {
      *_rita->ofh << "  curves=";
      for (int j=0; j<nc; ++j)
         _Ccode[c].l.push_back(curves[j]), _Ccode[c].nb++;
      for (int j=0; j<nc-1; ++j)
         *_rita->ofh << curves[j] << ",";
      *_rita->ofh << curves[nc-1] << endl;
   }
   if (surfaces_ok) {
      *_rita->ofh << "  surfaces=";
      for (int j=0; j<ns; ++j)
         _Scode[c].l.push_back(surfaces[j]), _Scode[c].nb++;
      for (int j=0; j<ns-1; ++j)
         *_rita->ofh << surfaces[j] << ",";
      *_rita->ofh << surfaces[ns-1] << endl;
   }
   if (volumes_ok) {
      *_rita->ofh << "  volumes=";
      for (int j=0; j<nv; ++j)
         _Vcode[c].l.push_back(volumes[j]), _Vcode[c].nb++;
      for (int j=0; j<nv-1; j++)
         *_rita->ofh << volumes[j] << ",";
      *_rita->ofh << volumes[nv-1] << endl;
   }
   _ret = 0;
   return;
}


void mesh::setPoint()
{
   _ret = 0;
   int nb=0;
   bool n_ok=false, x_ok=false, y_ok=false, z_ok=false, h_ok=false, del_ok=false;
   static int it = 0;
   if (it==0) {
      _point.n = 1;
      _point.x = _point.y = _point.z = 0.;
      _point.h = 0.1;
   }
   it++;
   static const vector<string> kw {"label","n","coord","size"};
   static const string H = "point label=n coord=x,y,z size=h\n"
                           "n: Point's label\n"
                           "x, y, z: Point coordinates. If y and/or z are not given, their value is set to 0?\n"
                           "This can be the case for 1-D and 2-D.\n"
                           "h: Mesh size around the point.";
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
         case   1:
            _point.n = _cmd->int_token(0);
            n_ok = true;
            break;

         case   2:
            _point.x = _point.y = _point.z = 0.;
            if (nb==1) {
               _ret  = _data->getPar(0,"mesh>point>",_point.x);
               if (_ret)
                  return;
               x_ok = true;
            }
            else if (nb==2) {
               _ret  = _data->getPar(0,"mesh>point>",_point.x);
               _ret += _data->getPar(1,"mesh>point>",_point.y);
               if (_ret)
                  return;
               x_ok = y_ok = true;
            }
            else {
               _ret  = _data->getPar(0,"mesh>point>",_point.x);
               _ret += _data->getPar(1,"mesh>point>",_point.y);
               _ret += _data->getPar(2,"mesh>point>",_point.z);
               if (_ret)
                  return;
               x_ok = y_ok = z_ok = true;
            }
            break;

         case   3:
            _ret  = _data->getPar(0,"mesh>point>",_point.h);
            if (_ret)
               return;
            h_ok = true;
            break;

         case 100:
         case 101:
            cout << H << endl;
            _ret = 0;
            return;

         default:
            _rita->msg("mesh>point>","Unknown argument: "+_kw[n]);
            _ret = 1;
            return;
      }
   }
   if (del_ok) {
      if (!n_ok) {
         _rita->msg("mesh>point>","No point label given.");
         _ret = 1;
         return;
      }
      _points.erase(_point.n);
      *_rita->ofh << "  point  label="<< _point.n << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (!n_ok) {
      _point.n++;
      if (_verb>1)
         cout << "Assumed point label: " << _point.n << endl;
   }
   if (_point.n<=0) {
      _rita->msg("mesh>point>","Label must be positive.");
      _ret = 1;
      return;
   }
   if (!x_ok && _verb>1)
      cout << "Warning: No x-coordinate given. Assumed x-coordinate: " << _point.x << endl;
   if (!y_ok && _verb>1)
      cout << "Warning: No y-coordinate given. Assumed y-coordinate: " << _point.y << endl;
   if (!z_ok && _verb>1)
      cout << "Warning: No z-coordinate given. Assumed z-coordinate: " << _point.z << endl;
   if (!h_ok && _verb)
      cout << "Warning: No mesh size given. Assumed mesh size: " << _point.h << endl;
   if (_generator>0 && _generator<=3) {
      _rita->msg("mesh>point>:","Mesh already generated. Needs retyping command to confirm.");
      _generator = 10;
   }
   _points[_point.n] = _point;
   *_rita->ofh << "  point  label=" << _point.n << "  x=" << _point.x << "  y=" << _point.y
               << "  z=" << _point.z << "  size=" << _point.h << endl;
#ifndef USE_GMSH
   _theDomain->insertVertex(_point.x,_point.y,_point.z,_point.h,0);
#endif
   _ret = 0;
   return;
}


void mesh::setCurve()
{
   static int nn = 0;
   int nb=0, n1=0, n2=0, n3=0;
   bool n_ok=false, line_ok=false, circle_ok=false, del_ok=false;
   _ret = 0;
   const static string H = "curve label=n [line=n1,n2] [circle=n1,n2,n3]\n\n"
                           "n: Curve's label\n"
                           "n1, n2: Labels of points defining a line\n"
                           "n1, n2, n3: Labels of points defining a circular arc, n1 and n2 are the two \n"
                           "            points defining the extremities of the arc and n3 is the point defining\n"
                           "            the center of the circular arc.";
   if (_generator>0 && _generator<=3) {
      _rita->msg("mesh>curve>","Mesh already generated. Needs retyping command to confirm.");
      _generator = 0;
      return;
   }
   _generator = 10;
   static const vector<string> kw {"label","n","line","circle","del$ete"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("mesh>curve>","No argument for command.",H);
      _ret = 1;
      return;
   }
   if (nb_args<1) {
      _rita->msg("mesh>curve>","No argument for command.");
      _ret = 1;
      return;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
         case   1:
            nn = _cmd->int_token(0);
            n_ok = true;
            break;

         case   2:
            if (nb!=2) {
               _rita->msg("mesh>curve>","Illegal number of arguments for a line");
               _ret = 1;
               return;
            }
            n1 = _cmd->int_token(0);
            n2 = _cmd->int_token(1);
            line_ok = true;
            break;

         case   3:
            if (nb!=3) {
               _rita->msg("mesh>curve>","Illegal number of arguments for a circle");
               _ret = 1;
               return;
            }
            n1 = _cmd->int_token(0);
            n2 = _cmd->int_token(1);
            n3 = _cmd->int_token(2);
            circle_ok = true;
            break;

         case   4:
            del_ok = true;
            break;

         case 100:
         case 101:
            cout << H << endl;
            _ret = 0;
            return;

         default:
            _rita->msg("mesh>curve>","Unknown argument: "+_kw[n]);
            _ret = 1;
            return;
      }
   }

   if (del_ok) {
      if (!n_ok) {
         _rita->msg("mesh>curve>","No curve label given.");
         _ret = 1;
         return;
      }
      _curve.erase(nn);
      *_rita->ofh << "  curve  label="<< nn << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (nn<=0) {
      _rita->msg("mesh>curve>","Label must be positive.");
      _ret = 1;
      return;
   }
   if (!n_ok) {
      nn++;
      if (_verb>1)
         cout << "Assumed curve label: " << nn << endl;
   }
   if (!line_ok && !circle_ok) {
      _rita->msg("mesh>curve>","A line or a circle must be defined.");
      _ret = 1;
      return;
   }
   if (line_ok && circle_ok) {
      _rita->msg("mesh>curve>","Curve cannot be defined as a line and a circle simultaneously.");
      _ret = 1;
      return;
   }
   _curve[nn].nb = nn;
   *_rita->ofh << "  curve label=" << nn;
   if (line_ok) {
      if (_points.find(n1)==_points.end()) {
         _rita->msg("mesh>curve>","Undefined end point: "+to_string(n1));
         _ret = 1;
         return;
      }
      if (_points.find(n2)==_points.end()) {
         _rita->msg("mesh>curve>","Undefined end point: "+to_string(n2));
         _ret = 1;
         return;
      }
      _curve[nn].l.clear(); _curve[nn].nb = 0;
      _curve[nn].l.push_back(n1); _curve[nn].nb++;
      _curve[nn].l.push_back(n2); _curve[nn].nb++;
      _curve[nn].type = 1;
      *_rita->ofh << "  line=" << n1 << "," << n2 << endl;
   }
   if (circle_ok) {
      if (_points.find(n1)==_points.end()) {
         _rita->msg("mesh>curve>","Undefined point: "+to_string(n1));
         _ret = 1;
         return;
      }
      if (_points.find(n2)==_points.end()) {
         _rita->msg("mesh>curve>","Undefined point: "+to_string(n2));
         _ret = 1;
         return;
      }
      if (_points.find(n3)==_points.end()) {
         _rita->msg("mesh>curve>","Undefined point: "+to_string(n3));
         _ret = 1;
         return;
      }
      _curve[nn].l.clear();
      _curve[nn].l.push_back(n1);
      _curve[nn].l.push_back(n2);
      _curve[nn].l.push_back(n3);
      _curve[nn].nb = 3;
      _curve[nn].type = 2;
      *_rita->ofh << "  circle=" << n1 << "," << n2 << "," << n3 << endl;
   }
#ifndef USE_GMSH
   _theDomain->insertLine(_curve[nn].l[0],_curve[nn].l[1],0);
#endif
   _ret = 0;
   return;
}


void mesh::setContour()
{
   int nb=0, nn=0, s=1, n1=0, n2=0;
   vector<int> curv, surf;
   int n_ok=0, curves_ok=0, surfaces_ok=0, del_ok=0;
   const static string H = "contour label=n curve=c1,c2,...\n\n"
                           "n: Contour's label\n"
                           "c1, c2, ...: Labels of curves defining the contour";
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      _rita->msg("mesh>contour>","Mesh already generated.");
      _generator = 0;
      return;
   }
   _generator = 10;
   static const vector<string> kw {"label","n","curv$es","surf$aces"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args==0) {
      _rita->msg("mesh>contour>","No argument for command.",H);
      _ret = 1;
      return;
   }
   if (nb_args<1) {
      _rita->msg("mesh>contour>","No argument for command.");
      _ret = 1;
      return;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
         case   1:
            nn = _cmd->int_token(0);
            n_ok++;
            break;

         case   2:
            for (int i=0; i<nb; i++)
               curv.push_back(_cmd->int_token(i));
            curves_ok++;
            break;

         case   3:
            for (int j=0; j<nb; j++)
               surf.push_back(_cmd->int_token(j));
            surfaces_ok++;
            break;

         case   4:
            del_ok++;
            break;

         case 100:
         case 101:
            cout << H << endl;
            _ret = 0;
            return;

         default:
            _rita->msg("mesh>contour>","Unknown argument: "+_kw[n]);
            _ret = 1;
	    return;
      }
   }
   if (!n_ok) {
      _rita->msg("mesh>contour>","Missing contour label.");
      _ret = 1;
      return;
   }
   if (n_ok>1) {
      _rita->msg("mesh>contour>","Contour label cannot be given more than once.");
      _ret = 1;
      return;
   }
   if (del_ok) {
      if (del_ok>1) {
         _rita->msg("mesh>contour>","Keyword delete cannot be given more than once.");
         _ret = 1;
         return;
      }
      _Ccontour.erase(nn);
      _Scontour.erase(nn);
      *_rita->ofh << "  contour  label="<< nn << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (nn==0) {
      cout << "Contour label cannot be zero." << endl;
      _rita->msg("mesh>contour>","Contour label cannot be zero.");
      _ret = 1;
      return;
   }
   if (!curves_ok && !surfaces_ok) {
      _rita->msg("mesh>contour>","Contour must be defined either by curves or surfaces.");
      _ret = 1;
      return;
   }
   if (curves_ok && surfaces_ok) {
      _rita->msg("mesh>contour>","Contour cannot be defined by curves or surfaces simultaneously.");
      _ret = 1;
      return;
   }
   s = 1;
   if (nn<0)
      s = -1, nn = -nn;
   if (curves_ok) {
      if (curves_ok>1) {
         _rita->msg("mesh>contour>","A curve contour already defined.");
         _ret = 1;
         return;
      }
      _Ccontour[nn].nb = 0;
      for (int i=0; i<nb; ++i) {
         int l = curv[i];
         if (_curve.find(l)==_curve.end()) {
            _rita->msg("mesh>contour>","Undefined curve: "+to_string(l));
            _ret = 1;
            return;
         }
         _Ccontour[nn].l.push_back(l), _Ccontour[nn].nb++;
      }
      n1 = _curve[_Ccontour[nn].l[0]].l[0];
      if (s==-1)
         n1 = _curve[_Ccontour[nn].l[0]].l[1];
      n2 = _curve[_Ccontour[nn].l[_Ccontour[nn].nb-1]].l[1];
      if (s==-1)
         n2 = _curve[_Ccontour[nn].l[_Ccontour[nn].nb-1]].l[0];
      if (n1!=n2) {
         _rita->msg("mesh>contour>:","Non closed contour.");
         _ret = 1;
         return;
      }	 
      *_rita->ofh << "  contour  label=" << nn << "  curves=";
      for (int i=0; i<nb-1; ++i)
         *_rita->ofh << curv[i] << ",";
      *_rita->ofh << curv[nb-1] << endl;
      _nb_Ccontour = _Ccontour.size();
   }
   if (surfaces_ok) {
      if (surfaces_ok>1)
         _rita->msg("mesh>contour>","A surface contour already defined.");
      _Scontour[nn].nb = 0;
      for (int i=0; i<nb; ++i) {
         int l = surf[i];
         if (_surface.find(l)==_surface.end()) {
            _rita->msg("mesh>contour>","Undefined surface: "+to_string(l));
            _ret = 1;
            return;
         }
         _Scontour[nn].l.push_back(l); _Scontour[nn].nb++;
      }
      n1 = _curve[_Scontour[nn].l[0]].l[0];
      if (s==-1)
         n1 = _surface[_Scontour[nn].l[0]].l[1];
      n2 = _surface[_Scontour[nn].l[_Scontour[nn].nb-1]].l[1];
      if (s==-1)
         n2 = _surface[_Scontour[nn].l[_Scontour[nn].nb-1]].l[0];
      if (n1 != n2) {
         _rita->msg("mesh>contour>","Non closed contour.");
         _ret = 1;
         return;
      }
      *_rita->ofh << "  contour  label=" << nn << "  surfaces=";
      for (int i=0; i<nb-1; ++i)
         *_rita->ofh << surf[i] << ",";
      *_rita->ofh << curv[nb-1] << endl;
      _nb_Scontour = _Scontour.size();
   }
   _ret = 0;
   return;
}


void mesh::setSurface()
{
   static int nn = 1;
   int nb=0, n=_surface.size()+1;
   bool n_ok=false, cont_ok=false, del_ok=false;
   vector<int> cont;
   const static string H = "surface label=n contours=c1,c2,...\n"
                           "n: Surface's label\n"
                           "c1, c2, ...: Labels of contours defining the surface.";
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      _rita->msg("mesh>surface>","Mesh already generated. Needs retyping command to confirm.");
      _generator = 0;
      return;
   }
   _generator = 10;
   static const vector<string> kw {"label","n","contours"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   if (nb_args<1) {
      _rita->msg("mesh>surface>","No argument for command.",H);
      _ret = 1;
      return;
   }
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArgs(nb);
      switch (n) {

         case   0:
         case   1:
            nn = _cmd->int_token(0);
            n_ok = true;
            break;

         case   2:
            for (int j=0; j<nb; ++j)
               cont.push_back(_cmd->int_token(j));
            cont_ok = true;
            break;

         case   3:
            del_ok = true;
            break;

         case 100:
         case 101:
            cout << H << endl;
            _ret = 0;
            return;

         default:
            _rita->msg("mesh>surface>","Unknown argument: "+_kw[n]);
            _ret = 1;
	    return;
      }
   }
   if (del_ok) {
      if (!n_ok) {
         _rita->msg("mesh>surface>","No surface label given.");
         _ret = 1;
         return;
      }
      _surface.erase(nn);
      *_rita->ofh << "  surface  label="<< nn << "  delete" << endl;
      _ret = 0;
      return;
   }
   if (!n_ok) {
      nn++;
      if (_verb)
         cout << "Assumed surface label: " << nn << endl;
   }
   if (nn==0) {
      _rita->msg("mesh>surface>","Label 0 is forbidden");
      _ret = 1;
      return;
   }
   if (!cont_ok) {
      _rita->msg("mesh>surface>","Surface contours must be given.");
      _ret = 1;
      return;
   }
   _surface[nn].nb = 0;
   for (int i=0; i<nb; ++i) {
      int c = cont[i];
      if (_Ccontour.find(c)==_Ccontour.end()) {
         _rita->msg("mesh>surface>","Undefined curve contour: "+to_string(n));
         _ret = 1;
         return;
      }
      _surface[nn].l.push_back(c); _surface[nn].nb++;
   }
   *_rita->ofh << "  surface  label=" << nn << "  contours=";
   for (int i=0; i<nb-1; ++i)
      *_rita->ofh << cont[i] << ",";
   *_rita->ofh << cont[nb-1] << endl;
   _nb_surface = _surface.size();
   _ret = 0;
   return;
}


void mesh::setVolume()
{
   _ret = 0;
   _generator = 10;
   _nb_volume = 0;
   _rita->msg("mesh>volume>","Volume generation not implemented yet.");
}


void mesh::setSubDomain()
{
   _ret = 0;
   string fn="";
   if (_generator>0 && _generator<=3) {
      _rita->msg("mesh>vertex>","Mesh already generated. Needs retyping command to confirm.");
      _generator = 0;
      return;
   }
   _generator = 10;
   int ret=0;
   *_rita->ofh << "   subdomain " << endl;   
   _sd.ln = 1, _sd.orientation = 1, _sd.code = 1;
   if (_verb) {
      cout << "Default line in subdomain: 1\n";
      cout << "Default orientation: 1\n";
      cout << "Default code: 10" << endl;
   }

   static const vector<string> kw {"line","orient$ation"};
   while (1) {
      if (_cmd->readline(sPrompt+"mesh>subdomain> ")<0)
         continue;
      _key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
      if (_key>=200) {
         _data->setDataExt(_key);
         continue;
      }
      switch (_key) {

         case   0:
            if (_verb>1)
               cout << "Setting line ..." << endl;
            if (_cmd->setNbArg(1,"Give label of a line in subdomain.")) {
               _rita->msg("mesh>subdomain>line>","Missing label of vertex in subdomain.","",1);
               break;
            }
            ret = _cmd->get(_sd.ln);
            if (!ret)
               *_rita->ofh << "    line " << _sd.ln << endl;
            break;

         case   1:
            if (_verb>1)
               cout << "Setting orientation ..." << endl;
            if (_cmd->setNbArg(1,"Give orientation of subdomain (1/-1).")) {
               _rita->msg("mesh>subdomain>orientation>","Missing orientation.","",1);
               break;
            }
            ret = _cmd->get(_sd.orientation);
            if (!ret)
               *_rita->ofh << "    orientation " << _sd.orientation << endl;   
            break;

         case   2:
            if (_verb>1)
               cout << "Setting code ..." << endl;
            if (_cmd->setNbArg(1,"Give code to associate to subdomain")) {
               _rita->msg("mesh>subdomain>code>","Missing code to associate to subdomain.","",1);
               break;
            }
            ret = _cmd->get(_sd.code);
            if (!ret)
               *_rita->ofh << "    code " << _sd.code << endl;   
            break;

         case   3:
            if (_verb>1)
               cout << "Saving subdomain ..." << endl;
            _cmd->setNbArg(0);
    //        _theDomain->insertSubDomain(_sd.ln,_sd.orientation,_sd.code);
            _nb_sub_domain++;
            _subdomains.push_back(_sd);
            cout << "Subdomain Added." << endl;
            cout << "Line " << _sd.ln << ", Orientation " << _sd.orientation
                 << ", Code " << _sd.code << endl;
            cout << "Total number of subdomains: " << _nb_sub_domain << endl;
            *_rita->ofh << "    save" << endl;
            _saved = true;
            return;

         case   4:
         case   5:
            if (_verb>1)
               cout << "Getting back to mesh menu ..." << endl;
            _cmd->setNbArg(0);
            if (!_saved)
               _rita->msg("mesh>subdomain>","Subdomain not saved.","Type again 'subdomain' and save generated one");
            *_rita->ofh << "    end" << endl;
            _ret = 0;
            return;

         case 100:
         case 101:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands\n";
            cout << "line        : Label for a line in subdomain\n";
            cout << "orientation : Orientation\n";
            cout << "code        : Code to associate to subdomain" << endl;
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
            _ret = 100;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            _rita->msg("mesh>subdomain>","Unknown command "+_cmd->token(),
                       "Available commands: lines");
            break;
       }
   }
   _ret = 0;
   return;
}


void mesh::saveGeo(const string& file)
{
   _generator = 10;
   ofstream s(file.c_str());
   for (auto const& v: _points)
      s << "Point(" << v.first << ") = {" << v.second.x << ", " << v.second.y
        << ", " << v.second.z << ", " << v.second.h << "};" << endl;

   for (auto const& v: _Pcode) {
      s << "Physical Point(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }
                   
   for (auto const& v: _curve) {
      Entity curve = v.second;
      if (curve.type==1)
         s << "Line(" << v.first << ") = {" << curve.l[0] << ", " << curve.l[1] << "};" << endl;
      else if (curve.type==2)
         s << "Circle(" << v.first << ") = {" << curve.l[0] << ", " << curve.l[2] << ", " 
           << curve.l[1] << "};" << endl;
   }

   for (auto const& v: _Ccode) {
      Entity code = v.second;
      s << "Physical Curve(" << v.first << ") = {";
      for (int i=0; i<code.nb-1; ++i)
         s << code.l[i] << ", ";
      s << code.l[code.nb-1] << "};" << endl;
   }

   for (auto const& v: _Ccontour) {
      s << "Curve Loop(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _surface) {
      s << "Plane Surface(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _Scode) {
      s << "Physical Surface(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }

   for (auto const& v: _Vcode) {
      s << "Physical Volume(" << v.first << ") = {";
      for (int i=0; i<v.second.nb-1; ++i)
         s << v.second.l[i] << ", ";
      s << v.second.l[v.second.nb-1] << "};" << endl;
   }
   _geo = true;
}


void mesh::saveDomain(const string& file)
{
   ofstream id(file.c_str());
   id << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>\n<OFELI_File>\n<info>" << endl;
   id << "<title>Domain file created by rita</title>" << endl;
   id << "<date></date>\n<author></author>\n</info>" << endl;
   id << "<Project name=\"\">" << endl;
   id << "  <Domain dim=\"2\">" << endl;
   for (int i=0; i<_nb_point; ++i)
      id << "    <vertex>" << _points[i].x << " " << _points[i].y << " "
         << 0 << " " << _points[i].h << "</vertex>" << endl;
   for (int i=1; i<=_nb_curve; ++i)
      id << "    <line>" << _curve[i].l[0] << " " << _curve[i].l[1] << " "
         << 0 << "</line>" << endl;
/*   for (int i=0; i<_nb_circle; ++i)
      id << "    <circle>" << _circles[i].pt1 << " " << _circles[i].c << " " << _circles[i].pt2
         << " " << 0 << "</circle>" << endl;*/
   for (int i=0; i<_nb_sub_domain; ++i)
      id << "    <Subdomain>" << _subdomains[i].ln << " " << _subdomains[i].orientation << " "
         << _subdomains[i].code << "</Subdomain>" << endl;
   id << "  </Domain>\n</Project>\n</OFELI_File>" << endl;
}


void mesh::Generate()
{
   _ret = 0;
   if (_generator>0 && _generator<=3) {
      if (_generator==1) {
         _rita->msg("mesh>generate>","A 1-D mesh has already been generated.");
      }
      else if (_generator==2) {
         _rita->msg("mesh>generate>","A 2-D mesh has already been generated.");
      }
      else if (_generator==3) {
         _rita->msg("mesh>generate>","A 3-D mesh has already been generated.");
      }
      return;
   }
   _generator = 10;
   if (_theMesh!=nullptr)
      delete _theMesh, _theMesh = nullptr;
   _mesh_file = "rita.m";
   
// Save geo gmsh file and generate gmsh mesh
   if (!_geo)
      saveGeo("rita.geo");

#ifdef USE_GMSH
   if (_verb)
      cout << "Starting mesh generation using Gmsh ..." << endl;
   gmsh::initialize();
   gmsh::option::setNumber("General.Terminal",1);
   gmsh::model::add("rita");
   for (auto const& v: _points)
      gmsh::model::geo::addPoint(v.second.x,v.second.y,v.second.z,v.second.h,v.first);
   for (auto const& v: _curve) {
      if (v.second.type==1)
         gmsh::model::geo::addLine(v.second.l[0],v.second.l[1],v.first);
      else if (v.second.type==2)
         gmsh::model::geo::addCircleArc(v.second.l[0],v.second.l[2],v.second.l[1],v.first);
   }
   for (auto const& v: _Ccontour)
      gmsh::model::geo::addCurveLoop(v.second.l,v.first);
   for (auto const& v: _surface)
      gmsh::model::geo::addPlaneSurface(v.second.l,v.first);
   for (auto const& v: _Pcode)
      gmsh::model::addPhysicalGroup(0,v.second.l,v.first);
   for (auto const& v: _Ccode)
      gmsh::model::addPhysicalGroup(1,v.second.l,v.first);
   for (auto const& v: _Scode)
      gmsh::model::addPhysicalGroup(2,v.second.l,v.first);
   for (auto const& v: _Vcode)
      gmsh::model::addPhysicalGroup(3,v.second.l,v.first);
   gmsh::model::setPhysicalName(2,4,"Domain");
   gmsh::model::geo::synchronize();
   gmsh::model::mesh::generate(2);
   gmsh::write("rita.msh");
   gmsh::finalize();
   if (_verb)
      cout << "Gmsh mesh generation complete." << endl;
   _theMesh = new OFELI::Mesh("rita.msh",false,NODE_DOF,_nb_dof);
   _generated = true;
#else
   _theDomain->setDim(_dim);
   _theDomain->setNbDOF(size_t(_nb_dof));
   saveDomain("rita.dom");
   _theDomain->genMesh(_mesh_file);
   _theMesh = new OFELI::Mesh(_mesh_file,false,NODE_DOF,2);
   _generated = true;
   _generator = 4;
#endif
   *_rita->ofh << "  generate" << endl;
   _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
   if (_verb)
      cout << "Mesh complete. Mesh Name: M-"+to_string(_data->theMesh.size()-1) << endl;
   _ret = 0;
   return;
}


void mesh::setNbDOF()
{
   int n = 0;
   if (_cmd->setNbArg(1,"Give number of degrees of freedom per node.")) {
      _rita->msg("mesh>nbdof>","Missing number of degrees of freedom.","",1);
      _ret = 1;
      return;
   }
   _ret = _cmd->get(n);
   if (!_ret) {
      *_rita->ofh << "    nbdof " << n << endl;
      _nb_dof = n;
   }
}


void mesh::Plot()
{
   string fn="";
   if (_theMesh==nullptr) {
      _rita->msg("mesh>plot>","No mesh to plot.");
      return;
   }
   _cmd->setNbArg(0);
   string file = "rita.msh";
   _theMesh->save(file);
   string com = "gmsh " + file;
   if (system(com.c_str())) {
      _rita->msg("mesh>plot>","Unrecognizable system command.");
      _ret = 1;
      return;
   }
   remove(file.c_str());
   _ret = 0;
}


void mesh::Read()
{
   int ret=0;
   ifstream ip;
   string dom_file, file, bamg_file, geo_file, out_file, Cmd, msh_file, fn="";
   int mesh_ok=0, geo_ok=0, gmsh_ok=0;
   static const vector<string> kw {"mesh","geo","gmsh"};
   _cmd->set(kw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case   0:
            file = _cmd->string_token();
            mesh_ok++;
            break;

         case   1:
            file = _cmd->string_token();
            geo_ok++;
            break;

         case   2:
            file = _cmd->string_token();
            gmsh_ok++;
            break;

         default:
            _rita->msg("mesh>read>","Unknown argument: "+_kw[n]);
            _ret = 1;
            return;
      }
   }
   if (nb_args>0) {
      if (mesh_ok+geo_ok+gmsh_ok==0) {
         _rita->msg("mesh>read>","No input file provided.");
         _ret = 1;
         return;
      }
      if (mesh_ok+geo_ok+gmsh_ok>1) {
         _rita->msg("mesh>read>","Only one input file must be provided.");
         _ret = 1;
         return;
      }
      ip.open(file);
      if (ip.is_open())
         ip.close();
      else {
         _rita->msg("mesh>read>","Unable to open file: "+file);
         _ret = 1;
         return;
      }
      *_rita->ofh << "  read";
      if (mesh_ok) {
         if (file.substr(file.find_last_of(".")+1)!="m") {
            _rita->msg("mesh>read>","File extension must be \".m\"");
            _ret = 1;
            return;
         }
         _theMesh = new OFELI::Mesh(file);
         if (_theMesh->getNbNodes()==0) {
            _rita->msg("mesh>read>","Empty mesh");
            _ret = 1;
            return;
         }
         _mesh_file = file;
         _nb_dof = _theMesh->getNbDOF() / _theMesh->getNbNodes();
         *_rita->ofh << "  mesh=" << file;
         _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
         _generator = 1;
         _generated = true;
      }
      else if (geo_ok) {
         if (file.substr(file.find_last_of(".")+1)!="geo") {
            _rita->msg("mesh>read>","File extension must be \".geo\"");
            _ret = 1;
            return;
         }
         msh_file = file.substr(0,file.rfind(".")) + "msh";
         Cmd = "gmsh -2 " + file + " -o " + msh_file;
         if (system(Cmd.c_str())) {
            _rita->msg("mesh>read>:","Unrecognizable system command.");
            _ret = 1;
            return;
         }
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file,GMSH);
         _data->NameMesh.push_back("M"+to_string(_data->theMesh.size()));
         _data->theMesh.push_back(_theMesh);
         _generated = true;
         *_rita->ofh << "  geo=" << file;
      }
      else if (gmsh_ok) {
         if (file.substr(file.find_last_of(".")+1)!="msh") {
            _rita->msg("mesh>read>","File extension must be \".msh\"");
            _ret = 1;
            return;
         }
         _theMesh = new OFELI::Mesh;
         _theMesh->get(file,GMSH);
         _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
         _generator = 1;
         _generated = true;
         *_rita->ofh << "  gmsh=" << file;
      }
      *_rita->ofh << endl;
   }
   else {
      _cmd->setNbArg(0);
      while (1) {
         if (_cmd->readline(sPrompt+"mesh>read> ")<0)
            continue;

         _key = _cmd->getKW(kw,_rita->_gkw,_rita->_data_kw);
         if (_key>=200) {
            _data->setDataExt(_key);
            continue;
         }
         switch (_key) {

            case   0:
               if (_cmd->setNbArg(1,"Give OFELI mesh file name.")) {
                  _rita->msg("mesh>read>mesh>","Missing Mesh file name.","",1);
                  break;
               }
               ret = _cmd->get(file);
               ip.open(file);
               if (ip.is_open()) {
                  ip.close();
                  _theMesh = new OFELI::Mesh(file);
                  *_rita->ofh << "  read mesh " << file << endl;
                  if (_theMesh->getNbNodes()==0) {
                     _rita->msg("mesh>read>mesh>","Empty mesh");
                     _ret = 1;
                     return;
                  }
                  _nb_dof = _theMesh->getNbDOF() / _theMesh->getNbNodes();
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
                  _generator = 1;
                  _generated = true;
                  _ret = 90;
               }
               else {
                  _rita->msg("mesh>read>mesh>","Unable to open file: "+file);
                  _generated = false;
                  _ret = 0;
               }
               return;

            case   1:
               if (_cmd->setNbArg(1,"Give geo gmsh file name.")) {
                  _rita->msg("mesh>read>geo>","Missing geo file name.","",1);
                  break;
               }
               ret = _cmd->get(file);
               if (!ret) {
                  if (file.substr(file.find_last_of(".")+1)!=".geo") {
                     _rita->msg("mesh>read>geo>","File extension must be \".geo\"");
                     _ret = 1;
                     break;
                  }
                  msh_file = file.substr(0,file.rfind(".")) + "msh";
                  Cmd = "gmsh -2 " + file + " -o " + msh_file;
                  if (system(Cmd.c_str())) {
                     _rita->msg("mesh>read>geo>","Unrecognizable system command.");
                     _ret = 1;
                     break;
                  }
                  _theMesh = new OFELI::Mesh;
                  _theMesh->get(msh_file,GMSH);
                  *_rita->ofh << "  read geo " << file << endl;
                  _data->NameMesh.push_back("M"+to_string(_data->theMesh.size()));
                  _data->theMesh.push_back(_theMesh);
                  _geo = true;
                  _generated = false;
               }
               _ret = 90;
               return;

            case 2:
               if (_cmd->setNbArg(1,"Give gmsh file name where to save mesh.")) {
                  _rita->msg("mesh>read>gmsh>","Missing gmsh file name.","",1);
                  break;
               }
               ret = _cmd->get(file);
               if (!ret) {
                  if (file.substr(file.find_last_of(".")+1)!=".msh") {
                     _rita->msg("mesh>read>gmsh>","File extension must be \".msh\"");
                     _ret = 1;
                     break;
                  }
                  _theMesh = new OFELI::Mesh;
                  _theMesh->get(file,GMSH);
                  _data->addMesh(_theMesh,"M-"+to_string(_data->nb_meshes+1));
                  *_rita->ofh << "  read gmsh " << file << endl;
               }
               _generated = true;
               _generator = 1;
               _ret = 90;
               return;

            case 100:
            case 101:
               _cmd->setNbArg(0);
               cout << "\nAvailable Commands\n";
               cout << "mesh     : Read mesh in OFELI mesh file\n";
               cout << "geo      : Read mesh in OFELI mesh file\n";
               cout << "gmsh     : Read mesh in gmsh file" << endl;
               break;

            case 102:
               _rita->getLicense();
               break;

            case 103:
               _ret = _configure->run();
               _verb = _configure->getVerbose();
               break;

            case 104:
               if (_cmd->setNbArg(1,"Data name to be given.",1)) {
                  _rita->msg("print>","Missing data name.","",1);
                  break;
               }
               if (!_cmd->get(fn))
                  _data->print(fn);
               break;

            case 105:
               _data->Summary();
               break;

            case 106:
            case 107:
               _ret = 0;
               return;

            case -2:
            case -3:
            case -4:
               break;

            default:
               _rita->msg("mesh>read>","Unknown command "+_cmd->token(),
                          "Available commands: mesh, geo, gmsh");
               break;
         }
      }
   }
   _ret = 0;
   return;
}


void mesh::Save()
{
   string domain_f="rita.dom", geo_f="rita.geo", mesh_f="rita.m", gmsh_f="rita.msh";
   string t="", vtk_f="rita.vtk", gnuplot_f="rita-gpl.dat", matlab_f="rita-matlab.m";
   string tecplot_f="rita-tecplot.dat";
   int domain_ok=0, geo_ok=0, mesh_ok=0, gmsh_ok=0, vtk_ok=0, gnuplot_ok=0, matlab_ok=0, tecplot_ok=0;
   _ret = 0;
   static const vector<string> kw {"domain","geo","mesh","gmsh","vtk","gnuplot","matlab","tecplot"};
   _cmd->set(kw,_rita->_gkw);
   int nb_args = _cmd->getNbArgs();
   for (int i=0; i<nb_args; ++i) {
      int n = _cmd->getArg();
      switch (n) {

         case 0:
#ifndef USE_GMSH
            if ((t=_cmd->string_token())!="")
               domain_f = t;
            domain_ok++;
#endif
            break;

         case 1:
            if ((t=_cmd->string_token())!="")
               geo_f = t;
            geo_ok++;
            break;

         case 2:
            if ((t=_cmd->string_token())!="")
               mesh_f = t;
            mesh_ok++;
            break;

         case 3:
            if ((t=_cmd->string_token())!="")
               gmsh_f = t;
            gmsh_ok++;
            break;

         case 4:
            if ((t=_cmd->string_token())!="")
               vtk_f = t;
            vtk_ok++;
            break;

         case 5:
            if ((t=_cmd->string_token())!="")
               gnuplot_f = t;
            gnuplot_ok++;
            break;

         case 6:
            if ((t=_cmd->string_token())!="")
               matlab_f = t;
            matlab_ok++;
            break;

         case 7:
            if ((t=_cmd->string_token())!="")
               tecplot_f = t;
            tecplot_ok++;
            break;

         case 100:
         case 101:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands\n";
            cout << "domain:  Save domain file\n";
            cout << "geo:     Save geo file\n";
            cout << "mesh:    Save mesh in OFELI format\n";
            cout << "gmsh:    Save mesh in gmsh format\n";
            cout << "vtk:     Save mesh in vtk format\n";
            cout << "gnuplot: Save mesh in gnuplot format\n";
            cout << "matlab:  Save mesh in matlab format\n";
            cout << "tecplot: Save mesh in tecplot format\n" << endl;
            break;

         default:
            _rita->msg("mesh>save>","Unknown argument.");
            _ret = 1;
	    break;
      }
   }
   if (domain_ok+geo_ok+mesh_ok+gmsh_ok+vtk_ok+gnuplot_ok+matlab_ok+tecplot_ok==0) {
      _rita->msg("mesh>save>","Nothing to save.");
      _ret = 1;
      return;
   }
   if (mesh_ok+gmsh_ok+vtk_ok+gnuplot_ok+matlab_ok+tecplot_ok>0 && !_generated) {
      _rita->msg("mesh>save>","No generated mesh to be saved");
      _ret = 1;
      return;
   }
   *_rita->ofh << "  save";
   if (domain_ok) {
      saveDomain(domain_f);
      *_rita->ofh << "  domain=" << domain_f;
   }
   if (geo_ok) {
      saveGeo(geo_f);
      *_rita->ofh << "  geo=" << geo_f;
   }
   if (mesh_ok) {
      _mesh_file = mesh_f;
      _theMesh->put(mesh_f);
      *_rita->ofh << "  mesh=" << mesh_f;
   }
   if (gmsh_ok) {
      saveMesh(gmsh_f,*_theMesh,GMSH);
      *_rita->ofh << "  gmsh=" << gmsh_f;
   }
   if (vtk_ok) {
      saveMesh(vtk_f,*_theMesh,VTK);
      *_rita->ofh << "  vtk=" << vtk_f;
   }
   if (gnuplot_ok) {
      saveMesh(gnuplot_f,*_theMesh,GNUPLOT);
      *_rita->ofh << "  gnuplot=" << gnuplot_f;
   }
   if (matlab_ok) {
      saveMesh(matlab_f,*_theMesh,MATLAB);
      *_rita->ofh << "  matlab=" << matlab_f;
   }
   if (tecplot_ok) {
      saveMesh(tecplot_f,*_theMesh,TECPLOT);
      *_rita->ofh << "  tecplot=" << tecplot_f;
   }
   *_rita->ofh << endl;
   return;
}

 /*
void mesh::Save()
{
   int ret=0;
   string file;
   *_rita->ofh << "  save" << endl;
   while (1) {
      if (_cmd->readline(sPrompt+"mesh>save> ")<0)
         continue;

      switch (_key=_cmd->getKW(_kw_save)) {

         case 0:
         case 1:
            _cmd->setNbArg(0);
            cout << "\nAvailable Commands\n";
#ifndef USE_GMSH
            cout << "domain  : Save Domain in a Domain format file\n";
#endif
            cout << "geo     : Save geometry in gmsh geo file\n";
            cout << "mesh    : Save mesh in OFELI mesh file\n";
            cout << "gmsh    : Save mesh in gmsh file\n";
            cout << "vtk     : Save mesh in vtk file\n";
            cout << "gnuplot : Save mesh in gnuplot file\n";
            cout << "matlab  : Save mesh in Matlab(R) file\n";
            cout << "tecplot : Save mesh in Tecplot(R) graphics file" << endl;
            break;

         case 2:
            setConfigure();
            break;

         case 3:
#ifndef USE_GMSH
            file = "rita.dom";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            if (!ret)
               saveDomain(file);
            *_rita->ofh << "    domain " << file << endl;
            _ret = 0;
#endif
            return;

         case 4:
            file = "rita.geo";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            if (!ret)
               *_rita->ofh << "    geo " << file << endl;
            saveGeo(file);
            _ret = 0;
            return;

         case 5:
            _mesh_file = "rita.m";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(_mesh_file);
            if (_generated==false) {
               cout << "Error: No mesh generated to be saved." << endl;
               *_ofl << "In rita>mesh>save>mesh>: No mesh generated to be saved" << endl;
               return;
            }
            _theMesh->put(_mesh_file);
            if (!ret)
               *_rita->ofh << "    mesh " << _mesh_file << endl;
            if (_verb)
               cout << "Mesh saved in file: " << _mesh_file << endl;
            _ret = 90;
            return;

         case 6:
            file = "rita.msh";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            if (_generated==false) {
               cout << "Error: No mesh generated to be saved." << endl;
               *_ofl << "In rita>mesh>save>msh>: No mesh generated to be saved" << endl;
               return;
            }
            saveMesh(file,*_theMesh,VTK);
            if (!ret)
               *_rita->ofh << "    gmsh " << file << endl;
            if (_verb)
               cout << "Mesh saved in file: " << file << endl;
            _ret = 90;
            return;

         case 7:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>vtk>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita.vtk";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,VTK);
            if (!ret)
               *_rita->ofh << "    vtk " << file << endl;
            _ret = 0;
            return;

         case 8:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>gnuplot>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita-gpl.dat";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,GNUPLOT);
            if (!ret)
               *_rita->ofh << "    gnuplot " << file << endl;
            _ret = 0;
            return;

         case 9:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>matlab>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita-matlab.m";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,MATLAB);
            if (!ret)
               *_rita->ofh << "    matlab " << file << endl;
            _ret = 0;
            return;

         case 10:
            if (_generated==false) {
               cout << "Error: No mesh generated." << endl;
               *_ofl << "In rita>mesh>save>tecplot>: No mesh generated." << endl;
               _ret = 1;
               return;
            }
            file = "rita-tecplot.dat";
            if (_cmd->getNbArgs()==1)
               ret = _cmd->get(file);
            saveMesh(file,*_theMesh,TECPLOT);
            if (!ret)
               *_rita->ofh << "    tecplot " << file << endl;
            _ret = 0;
            return;

         case 11:
         case 12:
            _cmd->setNbArg(0);
            if (!_saved) {
               cout << "Warning: Nothing saved." << endl;
               *_ofl << "In rita>mesh>save>: Nothing saved." << endl;
            }
            *_rita->ofh << "    end" << endl;
            _ret = 0;
            return;

         case 13:
         case 14:
            _ret = 100;
            return;

         case 15:
            _ret = 200;
            return;

         case -2:
         case -3:
         case -4:
            break;

         default:
            cout << "Unknown Command " << _cmd->token() << endl;
#ifdef USE_GMSH
            cout << "Available commands: geo, mesh, gmsh, vtk, gnuplot, matlab, tecplot, end, <" << endl;
#else
            cout << "Available commands: domain, geo, mesh, gmsh, vtk, gnuplot, matlab, tecplot, end, <" << endl;
#endif
            cout << "Global commands:    help, ?, set, quit, exit" << endl;
            *_ofl << "In rita>mesh>save>: Unknown command " << _cmd->token() << endl;
            break;
       }
   }

   if (_generator < 4) {
      cout << "Mesh already saved in file." << endl;
      *_ofl << "In rita>mesh>save>: Mesh already saved in file." << endl;
      _ret = 1;
      return;
   }
   if (_cmd->setNbArg(1,"Give name to file where to save."))
      return;
   if (_theMesh==nullptr) {
      cout << "Mesh was not properly created." << endl;
      *_ofl << "In rita>mesh>save>: mesh was not properly created" << endl;
      _ret = 1;
      return;
   }
   _cmd->get(file);
   _theMesh->save(file);
   _ret = 0;
   return;
}*/


void mesh::Clear()
{
   _cmd->setNbArg(0);
   if (_theMesh!=nullptr)
      delete _theMesh;
   _theMesh = nullptr;
   _saved = false;
   _dim = 1;
   _ret = 0;
/*   if (_theDomain)
      delete _theDomain;
   _theDomain = nullptr;*/
   _nb_point = _nb_curve = _nb_sub_domain = 0;
   _points.clear();
   _curve.clear();
   _surface.clear();
   _volume.clear();
   _subdomains.clear();
   *_rita->ofh << "  clear" << endl;
   if (_verb)
      cout << "Mesh cleared." << endl;
}

} /* namespace RITA */
