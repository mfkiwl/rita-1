#!/bin/sh

# test_rita.sh 
# A script to test rita

echo "==============================================="
echo "Testing rita Tutorial ..."
echo "==============================================="

cd tutorial/ae

echo "-----------------------------------------------"
echo "Test solution of algebraic equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then

rita example1.rita
rita example2.rita
rita example3.rita
fi

echo "-----------------------------------------------------------"
echo "Test solution of ordinary differential equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../ode
rita example1.rita
rita example2.rita
fi

echo "----------------------------------------------------------"
echo "Test solution of partial differential equations (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../pde
rita example1.rita
rita example2.rita
rita example3.rita
rita example4.rita
rita example5.rita
fi

echo "-------------------------------------------------"
echo "Test solution of optimization problems (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../optim
rita example1.rita
rita example2.rita
rita example3.rita
rita example4.rita
fi

echo "-------------------------------------"
echo "Test numerical integration (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../integration
rita example1.rita
rita example2.rita
fi

echo "-------------------------------------"
echo "Test calculator (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
cd ../calc
rita example1.rita
rita example2.rita
fi

echo "-----------------------------------"
echo "Remove all created files (y/n) ? \c"
read ans
if test "$ans" = "y" ; then
echo "Cleaning ..."
cd ../ae
rm -f *.dat *.sol .rita.log .rita.his
cd ../ode
rm -f *.dat *.sol .rita.log .rita.his
cd ../pde
rm -f *.dat *.bamg *.dom rita.m *.pos *.geo *.vtk rita-1d.m example*.m *.msh *.sol .rita.log .rita.his
cd ../optim
rm -f *.dat *.sol .rita.log .rita.his
cd ../integration
rm -f *.dat *.sol .rita.log .rita.his
cd ../calc
rm -f *.dat *.sol .rita.log .rita.his
fi

cd ..
echo "========================================================="
echo "rita Testing complete"
