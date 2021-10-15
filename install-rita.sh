#! /bin/sh
# install.sh

PREFIX="/usr/local"

mkdir -p $PREFIX/bin
cp bin/rita $PREFIX/bin/rita
cp test_rita.sh $PREFIX/bin/test_rita.sh
cp bin/gmsh $PREFIX/bin/.
mkdir -p $PREFIX/lib
cp -a lib/* $PREFIX/lib/.
mkdir -p $PREFIX/share/rita
cp -r doc/ $PREFIX/share/rita/doc/
cp -r tutorial/ $PREFIX/share/rita/tutorial/

echo "Installation complete"
echo "You can execute tutorial by typing test_rita.sh"
