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

                  Definition and implementation of class 'HVect'

  ==============================================================================*/

#pragma once


#include "linear_algebra/Vect_impl.h"
using std::vector;
using OFELI::Vect;

namespace RITA {

class HVect
{

 public:

    HVect() : nt(0) { }

    HVect(int n) { nt = n; }

    void set(Vect<double>& v, double t)
    {
       nt++;
       ts.push_back(t);
       v.setTime(t);
       vs.push_back(v);
    }

    Vect<double> *get(int n)
    {
       if (n<=nt)
          return &vs[n-1];
       else
          return nullptr; 
    }

    Vect<double> *get(double t)
    {
       for (int i=0; i<nt; ++i) {
          if (ts[i]==t)
             return &vs[i];
       }
       return nullptr;
    }

    int size() const { return nt; }

    string getVectorName() const { return name; }

    void setVectorName(string s) { name = s; }

    double getTime(int i) const { return ts[i-1]; } 

/// \brief Destructor
    ~HVect() {}

   int saveOFELI(const string &file, int e=1)
   {
      OFELI::IOField ff(file,OFELI::IOField::OUT);
      for (int n=0; n<nt; n+=e)
         ff.put(vs[n]);
      return 0;
   }

   int saveGnuplot(const string &file, int e=1)
   {
      ofstream ff(file);
      for (int n=0; n<nt; n+=e) {
         Vect<double> &v = vs[n];
         ff << ts[n];
         for (int i=0; i<v.size(); ++i)
            ff << "  " << v[i];
         ff << endl;
      }
      return 0;
   }

   int saveGmsh(const string &file, int e=1)
   {
      using namespace OFELI;
      int nb_en=0;
      ofstream ff(file);

      Vect<double> &v = *get(1);
      if (v.WithMesh()==false)
         return 1;
      Mesh &ms = v.getMesh();
      int nb_dof = v.getNbDOF();

      ff << "View \"" << v.getName() << "\" {" << endl;
      switch (ms.getDim()) {

         case 1:
            element_loop(&ms) {
               ff << "SL(";
               ff << The_element(1)->getX() <<  ", 0., 0., "
                  << The_element(2)->getX() <<  ", 0., 0. ) {" << endl;
               for (int n=0; n<nt; n+=e) {
                  Vect<double> &v = *get(n+1);
                  ff << v(The_element(1)->n(),1) << "," << v(The_element(2)->n(),1);
                  if (n<nt-1 && n+e<nt)
                     ff << ",";
                  ff << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;

         case 2:
            element_loop(&ms) {
               if (nb_dof==1)
                  ff << 'S';
               else
                  ff << 'V';
               if ((nb_en=The_element.getNbNodes())==3)
                  ff << "T(";
               else if (nb_en==4)
                  ff << "Q(";
               for (int k=1; k<nb_en; ++k)
                  ff << The_element(k)->getX() << "," << The_element(k)->getY() << ",0.,";
               ff << The_element(nb_en)->getX() << "," << The_element(nb_en)->getY() << ",0.) {" << endl;
               for (int n=0; n<nt; n+=e) {
                  Vect<double> &v = *get(n+1);
                  for (int k=1; k<nb_en; ++k) {
                     ff << v(The_element(k)->n(),1);
                     if (nb_dof > 1)
                        ff << "," << v(The_element(k)->n(),2) << ",0.0";
                     ff << "," << endl;
                  }
                  ff << v(The_element(nb_en)->n(),1);
                  if (nb_dof > 1)
                     ff << "," << v(The_element(nb_en)->n(),2) << ",0.0";
                  if (n<nt-1 && n+e<nt)
                     ff << ",";
                  ff << endl;
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
               for (int k=1; k<nb_en; ++k)
                  ff << The_element(k)->getX() << "," << The_element(k)->getY() << ","
                     << The_element(k)->getZ() << ",";
               ff << The_element(nb_en)->getX() << "," << The_element(nb_en)->getY() << "," 
                  << The_element(nb_en)->getZ() << ") {" << endl;
               for (int n=0; n<nt; n+=e) {
                  Vect<double> &v = *get(n+1);
                  for (int k=1; k<nb_en; ++k) {
                     ff << v(The_element(k)->n(),1);
                     if (nb_dof>1)
                        ff << "," << v(The_element(k)->n(),2) << "," << v(The_element(k)->n(),3);
                     ff << "," << endl;
                  }
                  ff << v(The_element(nb_en)->n(),1);
                  if (nb_dof>1)
                     ff << "," << v(The_element(nb_en)->n(),2) << "," << v(The_element(nb_en)->n(),3);
                  if (n<nt-1 && n+e<nt)
                     ff << ",";
                  ff << endl;
               }
               ff << "};" << endl;
            }
            ff << "};" << endl;
            break;
      }
      return 0;
   }

   int saveVTK(const string &file, int e=1)
   {
      using namespace OFELI;
      string of="", proj=file.substr(0,file.rfind("."));
      map<int,int> nen = {{LINE,2},{TRIANGLE,3},{QUADRILATERAL,4},{TETRAHEDRON,4},
                          {HEXAHEDRON,8},{PENTAHEDRON,6}};
      map<int,int> ShCode = {{LINE,3},{TRIANGLE,5},{QUADRILATERAL,9},{TETRAHEDRON,10},
                             {HEXAHEDRON,12},{PENTAHEDRON,13}};

      Vect<double> &v = *get(1);
      if (v.WithMesh()==false)
         return 1;
      int nb_dof = v.getNbDOF();
      Mesh &ms = v.getMesh();

      int sz=0;
      element_loop(&ms)
         sz += nen[The_element.getShape()] + 1;
      int nn = 0;
      for (int n=0; n<nt; n+=e) {
         Vect<double> &v = *get(n+1);
         of = proj + ".vtk";
         if (nt>1)
            of = proj + "-" + zeros(nn++) + ".vtk";
         cout << "   Storing time step " << n << " in file " << of << endl;
         ofstream pf(of.c_str());
         pf << setprecision(16) << std::scientific;
         pf << "# vtk DataFile Version 2.0\n# Imported from OFELI files\nASCII" << endl;
         pf << "DATASET UNSTRUCTURED_GRID\nPOINTS " << ms.getNbNodes() << " double" << endl;
         node_loop(&ms)
            pf << The_node.getX() << "  " << The_node.getY() << "  " << The_node.getZ() << endl;
         pf << "\nCELLS " << ms.getNbElements() << setw(10) << sz << endl;
         element_loop(&ms) {
            pf << setw(8) << nen[The_element.getShape()];
            for (int i=1; i<=nen[The_element.getShape()]; i++)
               pf << setw(10) << The_element(i)->n()-1;
            pf << endl;
         }
         pf << "\nCELL_TYPES  " << ms.getNbElements() << endl;
         int k=0;
         element_loop(&ms) {
            pf << setw(4) << ShCode[The_element.getShape()];
            if (++k%30 == 0)
               pf << endl;
         }
         pf << "\nPOINT_DATA  " << ms.getNbNodes() << endl;
         if (nb_dof==1)
            pf << "SCALARS  " << v.getName() << "  double  1\nLOOKUP_TABLE  default" << endl;
         else
            pf << "VECTORS  " << v.getName() << "  double" << endl;

         node_loop(&ms) {
            pf << v(node_label,1) << " ";
            if (nb_dof>1) {
               pf << v(node_label,2) << " ";
               if (nb_dof>2)
                  pf << v(node_label,3) << " ";
               else
                  pf << 0. << " ";
            }
            pf << endl;
         }
      }
      return 0;
   }

   int saveTecplot(const string& file, int e=1)
   {
      using namespace OFELI;
      map<int,string> shape = {{LINE,"LINESEG"},{QUADRILATERAL,"QUADRILATERAL"},{TRIANGLE,"TRIANGLE"},
                               {TETRAHEDRON,"TETRAHEDRON"},{HEXAHEDRON,"HEXAHEDRON"},{PENTAHEDRON,"HEXAHEDRON"}};
      Vect<double> &v = *get(1);
      if (v.WithMesh()==false)
         return 1;
      int nb_dof = v.getNbDOF();
      Mesh &ms = v.getMesh();

      ofstream ff(file.c_str());
      ff.setf(ios::right|ios::scientific);
      int count = 0;
      for (int n=0; n<nt; n+=e) {
         Vect<double> &v = *get(n+1);
         if (count==0) {
            ff << "TITLE = \" \"\n" << endl;
            ff << "VARIABLES = \"X\", \"Y\"";
            if (ms.getDim()==3)
               ff << ", \"Z\"";
            if (nb_dof==1)
               ff << ", \"T\"";
            else if (nb_dof==2)
               ff << ", \"UX\", \"UY\"";
            else if (nb_dof==3)
               ff << ", \"UX\", \"UY\", \"UZ\"";
            ff << endl;
         }
         ff << "\nZONE T=\"" << "step-" << n << "\", N=" << ms.getNbNodes() << ", E="
            << ms.getNbElements() << ", F=FEPOINT, ET=" << shape[ms.getShape()]
            << ", SOLUTIONTIME=" << v.getTime();
         if (count) {
            ff << ", D=(1,";
            if (ms.getDim()>1)
               ff << "2,";
            if (ms.getDim()==3)
               ff << "3,";
            ff << "FECONNECT)";
         }
         ff << endl;
         node_loop(&ms) {
            if (count==0)
               for (int i=1; i<=ms.getDim(); i++)
                  ff << "  " << The_node.getCoord(i);
            for (int j=1; j<=nb_dof; j++)
               ff << "  " << v(node_label,j);
            ff << endl;
         }
         if (count==0) {
            element_loop(&ms) {
               for (int i=1; i<=The_element.getNbNodes(); ++i)
                  ff << setw(10) << The_element(i)->n();
               ff << endl;
            }
         }
         count++;
      }
      ff.close();
      return 0;
   }

 private:

    int nt;
    string name;
    vector<Vect<double> > vs;
    vector<double> ts;

};


inline ostream &operator<<(ostream&     s,
                           const HVect& v)
{
   s << "Contents of history vector:" << endl;
   s << "Name:                " << v.getVectorName() << endl;
   s << "Number of time steps: " << v.size() << endl;
   if (v.size())
      s << "Time interval: [" << v.getTime(1) << "," << v.getTime(v.size()) << "]" << endl;
   return s;
}

} /* namespace RITA */
