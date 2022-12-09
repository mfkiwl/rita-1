#! /bin/sh
# install.sh

PREFIX="/usr/local"

echo "Execution of this script my require administrator privileges"
echo "In this case, type: sudo install.sh"

mkdir -p $PREFIX/bin
cp bin/rita $PREFIX/bin/rita
cp bin/gmsh $PREFIX/bin/.
mkdir -p $PREFIX/lib
cp -a lib/* $PREFIX/lib/.
mkdir -p $PREFIX/share/rita
cp -r doc/ $PREFIX/share/rita/
cp -r tutorial/ $PREFIX/share/rita/tutorial/

echo "Installation complete"
