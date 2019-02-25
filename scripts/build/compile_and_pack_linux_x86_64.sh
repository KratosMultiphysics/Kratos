#!/bin/bash

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
#to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 "

clear

# export AMDAPPSSDKROOT=/home/compile-kratos/Libraries/AMD-APP-SDK-v2.8-RC-lnx64
# export AMDAPPSDKSAMPLESROOT=/home/compile-kratos/Libraries/AMD-APP-SDK-v2.8-RC-lnx64
# export LD_LIBRARY_PATH=$AMDAPPSDKSAMPLESROOT/lib/x86_64:$LD_LIBRARY_PATH

#you may want to decomment this the first time you compile
rm CMakeCache.txt

################################################################################
#### EXTENDED CONFIGURATION
################################################################################

# Here we are goin to define all version numbers, svn revisions and dates
export SVN_ADDRESS="https://svn.cimne.upc.edu/p/kratos/kratos"
export SVN_WORKDIR="../"

export KRATOS_MAJOR="4"															      # Kratos major
export KRATOS_MINOR="4"															      # Kratos minor
export KRATOS_ARCHY="64"														      # Adress model
export KRATOS_NEWSD="Fixes in StructuralApplication"      # Changelog of the version
export KRATOS_TVERS="1.0.7"                               # Version of the interface
export KRATOS_TDATE="$(date '+%Y-%m-%d')"									# Today, or some other day

export GIDINT_SRC=$HOME/gid_interfaces/standard_14_07_14	# Superproblemtype path ( path to the problemtype that contains Kratos and DEM)
export KRATOS_INSTALL_DIR=$HOME/KratosRelease							# Kratos install dir
export KRATOS_PACK=$HOME/KratosPack												# Kratos pack dir ( where the tar would be )

export GIDINT_KTS=kratos.gid/Kratos/kratos.gid						# Path to Kratos Problemtype
export GIDINT_DEM=kratos.gid/DEM/kratos.gid								# Path to DEM Problemtype

# The sub-problemtypes In this case, Kratos and DEM
GIDINT_ARR=( $GIDINT_KTS $GIDINT_DEM )

# Extract the subversion number
export SVN_REV=`svn info ${SVN_WORKDIR} | grep '^Revision:' | sed -e 's/^Revision: //'`

## Set Promblemtype version and date
export PROBLEMTYPE_VERSION="${KRATOS_MAJOR}.${KRATOS_MINOR}.${SVN_REV}"

################################################################################
#### COMPILE PARAMETERS AND DIRECTORIES
################################################################################

export BLA_VENDOR="ATLAS"                                 # Blas vendor, we use atlas. Leave empty to use system default
export BOOST_ROOT=$HOME/boost_1_55_0                      # Boost directory file, Can be replaced by -DBOOT_ROOT in cmake
export DEPLOY_PATH=$HOME/masterdisc/www/data              # Where the tar.gz is going to be uploaded at the end. Can be omited

rm -rf $PROBLEMTYPE_PATH

# Set this so the system knows where atlas is installed. It is not needed inless you set BLA_VENDOR="ATLAS"
export LD_LIBRARY_PATH=/usr/lib/altas-base/atlas/:$LD_LIBRARY_PATH

################################################################################
#### CONFIGURE KRATOS
################################################################################

cmake ..  																		                                 \
-DCMAKE_C_COMPILER=/usr/bin/gcc 												                       \
-DCMAKE_CXX_COMPILER=/usr/bin/g++												                       \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3   " 								               \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -msse3   " 									                 \
-DBOOST_LIBRARYDIR="${BOOST_ROOT}/stage/lib"									                 \
-DPYTHON_EXECUTABLE="/usr/bin/python${PYTHON_VERSION_S}.${PYTHON_VERSION_M}"                    \
-DCMAKE_BUILD_TYPE=Release  													                         \
-DINCOMPRESSIBLE_FLUID_APPLICATION=ON  											                   \
-DMESHING_APPLICATION=ON 														                           \
-DEXTERNAL_SOLVERS_APPLICATION=ON												                       \
-DSTRUCTURAL_APPLICATION=ON 													                         \
-DCONVECTION_DIFFUSION_APPLICATION=ON 											                   \
-DFLUID_DYNAMICS_APPLICATION=ON 												                       \
-DMESH_MOVING_APPLICATION=ON 															                             \
-DDEM_APPLICATION=ON 															                             \
-DFSI_APPLICATION=ON 															                             \
-DULF_APPLICATION=ON															                             \
-DMIXED_ELEMENT_APPLICATION=ON													                       \
-DTHERMO_MECHANICAL_APPLICATION=ON												                     \
-DBLOOD_FLOW_APPLICATION=OFF													                         \
-DUSE_INTEL_GREATER_THAN_12=TRUE 												                       \
-DMETIS_APPLICATION=ON															                           \
-DPARMETIS_ROOT_DIR="/home/odin/compiled_libraries/ParMetis-3.2.0" 				     \
-DTRILINOS_APPLICATION=ON														                           \
-DTRILINOS_ROOT="/usr/local/trilinos-11.2.3"									                 \
-DINSTALL_EMBEDDED_PYTHON=ON 													                         \
-DINCLUDE_PASTIX=OFF 															                             \
-DPASTIX_INSTALL_DIR="/home/odin/Downloads/pastix_release_3923/install"			   \
-DSCOTCH_INSTALL_DIR="/home/odin/Downloads/scotch_6.0.0rc17/lib/" 				     \
-DEXCLUDE_ITSOL=ON 																                             \
-DKRATOS_INSTALL_PREFIX="$KRATOS_INSTALL_DIR" 									               \
-DINSTALL_EMBEDDED_PYTHON=ON 													                         \
-DINSTALL_PYTHON_FILES=ON														                           \
-DSOLID_MECHANICS_APPLICATION=ON 												                       \
-DSWIMMING_DEM_APPLICATION=ON													                         \
-DINCLUDE_FEAST=ON                                                                                                                                                     \

################################################################################
#### COMPILE THE CODE
################################################################################

make
make install

# If your machine has multiple cores, you may want to compile in parallel
# consider ~3GB of ram memory per core
# make -j4
# make install

################################################################################
#### CHECK THE RESULTS
################################################################################
#TODO

################################################################################
#### BUILD THE PROBLEMTYPE
################################################################################
cp -r ${GIDINT_SRC}/kratos.gid ${KRATOS_PACK}/

for GIDINT_DIR in $GIDINT_ARR
do
	cp -r $KRATOS_INSTALL_DIR ${KRATOS_PACK}/${GIDINT_DIR}/kratos
done

################################################################################
#### GENERATE THE TARBALL
################################################################################
for GIDINT_DIR in $GIDINT_ARR
do
	# Copying additional files to the install directory
	cp /usr/lib/atlas-base/atlas/*.so.3gf           ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/
	cp /usr/lib/atlas-base/atlas/*.so               ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/
	cp /usr/lib/atlas-base/*.so.3gf 			          ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/
	cp /usr/lib/atlas-base/*.so 					          ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/

	cp /lib/x86_64-linux-gnu/libcrypto.so.1.0.0 	  ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/libcrypto.so.1.0.0
	cp /lib/x86_64-linux-gnu/libssl.so.1.0.0 		    ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/libssl.so.1.0.0
	cp /usr/lib/x86_64-linux-gnu/libgfortran.so.3 	${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/libgfortran.so.3

	cp /lib/x86_64-linux-gnu/libz.so.1 				      ${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/libz.so.1

	# Coping aditional trillinos missing libraries
	cp /usr/local/trilinos-11.2.3/lib/libgaleri.so 	${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs/libgaleri.so

	# Redistibutable MPI support (you may probably want to comment this line unless you have a Redistibutable openmpi)
	cp -r /home/odin/openmpi_libs/ 					        ${KRATOS_PACK}/${GIDINT_DIR}/openmpi

	# Gratinc execution permission to script files
	chmod 777 ${KRATOS_PACK}/${GIDINT_DIR}/*.bat

	# Creating a link to the libpython
	export old_cd=`pwd`
	cd "${KRATOS_PACK}/${GIDINT_DIR}/kratos/libs"
	ln -s ../libpython3.3m.so.1.0 libpython3.3m.so.1.0
	ln -s libatlas.so.3gf.0 libatlas.so.3gf
	cd "$old_cd"

	# Copy Python zip file
	mkdir ${KRATOS_PACK}/${GIDINT_DIR}/kratos/lib
	cp -r /home/odin/Python-3.3.4/Lib ${KRATOS_PACK}/${GIDINT_DIR}/kratos/lib/python3.3

	# find ${PROBLEMTYPE_PATH} -type d -name .svn -exec rm -rf {} \;
	find ${KRATOS_PACK}/${GIDINT_DIR} -type d -name .svn 	-exec rm -rf {} \;
	find ${KRATOS_PACK}/${GIDINT_DIR} -name "*~" 			-exec rm -f  {} \;
	find ${KRATOS_PACK}/${GIDINT_DIR} -name "lib*.a" 		-exec rm -f  {} \;
done

# Setup kratos.xml, since it seems that no one is keeping it up to date...
sed -i 's/<Version>[^<]*<\/Version>/<Version>'"${PROBLEMTYPE_VERSION}"'<\/Version>/g' 	${KRATOS_PACK}/kratos.gid/kratos.xml

export TAR_NAME="kratos-${KRATOS_MAJOR}.${KRATOS_MINOR}.${SVN_REV}-linux-64.tgz"
export TGZ_NAME="${KRATOS_PACK}/kratos-${KRATOS_MAJOR}.${KRATOS_MINOR}.${SVN_REV}-linux-64.tgz"

# Creating the file
export old_cd=`pwd`
cd "${KRATOS_PACK}"
tar -cvzf ${TGZ_NAME} kratos.gid > ${KRATOS_PACK}/tar.log
cd "$old_cd"
echo ${TGZ_NAME} created.

# Replacing the older release
sudo rm ${DEPLOY_PATH}/*-linux-64.tgz
sudo cp ${KRATOS_PACK}/${TAR_NAME} ${DEPLOY_PATH}/${TAR_NAME}
