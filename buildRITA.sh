#! /bin/sh
# buildRITA.sh

if [ $# -ne 2 ]; then
   echo 1>&2 "Usage: buildRITA.sh <release> <MacOSX|Linux64|Linux32|Win32|Win64>"
   exit 2
fi

RELEASE=$1
PREFIX="/usr/local"
SYSTEM=$2

if [${SYSTEM} != "MacOSX"] &&
   [${SYSTEM} != "Linux64"] &&
   [${SYSTEM} != "Linux32"] &&
   [${SYSTEM} != "Win32"] &&
   [${SYSTEM} != "Win64"]; then
     echo "Error: Unavailable Platform $MACHINE"
   exit
fi

echo "Preparing release $RELEASE..."
echo "--------------------------"

RELEASE=rita-${RELEASE}

case "$SYSTEM" in

    MacOSX)
        echo "Creating Source Package ..."
        tar --exclude='./build' --exclude='.git' --exclude='.gitattributes' --exclude='.vscode' -czf rita.tar.gz rita/
        mv rita ${RELEASE}
        tar --exclude='./build' --exclude='.git' --exclude='.gitattributes' --exclude='.vscode' --exclude='./buildOFELI.sh' --exclude='install.sh' -czf ${RELEASE}-src.tar.gz ${RELEASE}
        mv ${RELEASE} rita

        echo "Creating Compiled Package ..."
        mkdir -p ${RELEASE}-${SYSTEM}
        cp rita/install.sh ${RELEASE}-${SYSTEM}/.
        cd ${RELEASE}-${SYSTEM}
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
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}-${SYSTEM}
	/bin/rm -rf ${RELEASE}-${SYSTEM}
       ;;

    Linux64)
        echo "Creating Compiled Package ..."
        mkdir -p ${RELEASE}-${SYSTEM}
        cp rita/install.sh ${RELEASE}-${SYSTEM}/.
        cd ${RELEASE}-${SYSTEM}
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
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}-${SYSTEM}
	/bin/rm -rf ${RELEASE}-${SYSTEM}
        ;;

    Linux32)
        echo "Creating Compiled Package ..."
        mkdir -p ${RELEASE}-${SYSTEM}
        cp rita/install.sh ${RELEASE}-${SYSTEM}/.
        cd ${RELEASE}-${SYSTEM}
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
        tar czf ${RELEASE}-${SYSTEM}.tar.gz ${RELEASE}-${SYSTEM}
	/bin/rm -rf ${RELEASE}-${SYSTEM}
        ;;

    Win32)
        PREFIX="/usr/local"
        ;;

esac

echo "Packing complete"
