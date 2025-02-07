#! /bin/sh
# pack.sh 
# A script to build rita distribution

if [ $# -eq 0 ]; then
   echo 1>&2 "Usage: pack.sh <release> [MacOSX|Linux64|Win64]"
   exit 2
fi

RELEASE=$1
PREFIX="/usr/local"
SYSTEM=$2
RELEASE=rita-${RELEASE}

if [${SYSTEM} -eq ""]; then
   tar --exclude='./build' --exclude='.DS_Store' --exclude='.git' --exclude='.gitattributes' --exclude='.vscode' -czf rita.tar.gz rita/
   mv rita ${RELEASE}
   tar --exclude='./build' --exclude='.DS_Store' --exclude='.git' --exclude='.gitattributes' --exclude='.vscode' --exclude='./pack.sh' --exclude='install.sh' -czf ${RELEASE}-src.tar.gz ${RELEASE}
   mv ${RELEASE} rita
fi

if [${SYSTEM} -ne "MacOSX"] &&
   [${SYSTEM} -ne "Linux64"] &&
   [${SYSTEM} -ne "Win64"]; then
     echo "Error: Unavailable Platform $MACHINE"
   exit
fi

echo "Preparing release $RELEASE..."
echo "-----------------------------"

case "$SYSTEM" in

    MacOSX)
        echo "Creating Compiled Package ..."
        mkdir -p ${RELEASE}
        cp rita/install.sh ${RELEASE}/.
        cd ${RELEASE}
        mkdir -p bin
        cp ${PREFIX}/bin/rita bin/rita
        cp ${PREFIX}/bin/gmsh bin/gmsh
        mkdir -p doc
        cp -r ${PREFIX}/share/rita/doc/ doc/
        mkdir -p tutorial
        cp -r ${PREFIX}/share/rita/tutorial/ tutorial/
        find tutorial/ -name "CMakeLists.txt" -exec rm {} \;
        mkdir -p lib
        cp -a ${PREFIX}/lib/libgmsh* lib/
	cd ../
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}
	/bin/rm -rf ${RELEASE}
	;;

    Linux64)
        echo "Creating Compiled Package ..."
        mkdir -p ${RELEASE}
        cp rita/install.sh ${RELEASE}/.
        cd ${RELEASE}
        mkdir -p bin
        cp ${PREFIX}/bin/rita bin/rita
        cp ${PREFIX}/bin/gmsh bin/gmsh
        mkdir -p doc
        cp -r ${PREFIX}/share/rita/doc/ doc/
        mkdir -p tutorial
        cp -r ${PREFIX}/share/rita/tutorial/ tutorial/
        find tutorial/ -name "CMakeLists.txt" -exec rm {} \;
        mkdir -p lib
        cp -a ${PREFIX}/lib/libgmsh* lib/
	cd ../
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}
	/bin/rm -rf ${RELEASE}
        ;;

    Win64)
        PREFIX="/usr/local"
        ;;

esac

echo "Packing complete"
