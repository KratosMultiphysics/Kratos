*******************************************************************************
PREREQUISITES:
*******************************************************************************
  the full compilation of the Kratos depends on a number of external libraries that should be present (and findable) on the system.

INSTALL PREREQUISITES THE EASY WAY: (install all what you may need)

  - on "deb based" Linux systems (ex ubuntu) --> run the following command to install all what you may need
              .............
  - on "rpm based" Linux systems (ex Fedora, RedHat, Suse)  --> run the following command on the command line
	      .............
  - Manual install of external libraries ... see detailed instructions at the end of the current file

Please note that the prerequisites described here also contain the mpi lib. If mpi version
is not needed nothing happens if installing mpi fails

*******************************************************************************
COMPILING AND INSTALLING
*******************************************************************************
This "cmake_build" directory is designed for building Kratos using the cmake tool.
After downloading the Kratos, the "cmake build" directory contains an example configuration file
named
    example_configure.sh.do_not_touch

to compile first of all open a terminal windows, go to the dir "cmake_build" and copy the example file 
to your own version, say "local_compile.sh"

    cp example_configure.sh.do_not_touch local_compile.sh

and edit local_compile.sh ... to customarize it following your taste.

REMARK: it is highly advisable that the "local_compile.sh" file is located 
in the directory cmake_build. If another directory is chosen one should change
"cmake .." with "cmake /path/to/kratos/root"

After modifying the configuration file with the local paths run (assuming the file you copied to is local_compile.sh)

    sh local_compile.sh

hopefully you should enjoy the kratos compilation with CMAKE


*******************************************************************************
MANUALLY INSTALLING LIBRARIES
*******************************************************************************
This section is needed for "power users" which want to customarize their version of the
Kratos or ... for users which do not have superuser rights on a given system.

The Kratos code is divided in applications so that an user can switch on some features.
Depending on the set of applications selected some applications may not be needed.
In the following we detail the libraries to be installed in order of importance

REQUIRED LIBRARIES:
  Boost Library: the boost library and the boost-python libraries should be available on the system.
  KRATOS CAN NOT COMPILE without boost. Installation instructions can be found at www.boost.org.
  After installation the user should inform the kratos compilation of the position of the boost library
  by adding at the beginning of your configuration script (or in your .bashrc)

           export BOOST_ROOT=/path/to/your/boost
	
  Python library: this library is normally included in the python-dev package and needs to be present as a shared-library (.so or .dll).
  If more than one library is present in the system, the user can control which one to link to by (for example) setting
  
		-DPYTHON_LIBRARY="/usr/lib64/libpython2.6.so"
		-DPYTHON_INCLUDE_DIR="/usr/include/python2.6"
   

OPTIONAL LIBRARIES 
   BLAS and LAPACK: 	--> Needed for "External Solvers Application" and for "Trilinos App"
		    this libraries are the most troublesome as they need linking a fortran lib :-(
                    we STRONGLY advise the user to avoid using (.a) version of the libs in the compilation of 
		    the Kratos (and in compiling trilinos)since linking statically brings countless problems 
		    in linking due to the need of explicitely linking in the fortran libs.
                    This problem is automatically solved by using (.so) version of the libs.

		    Please note that on the Kratos webpage it is possible to download a modified version of the 
		    ref BLAS which compiles to ".so". The modified version simply replaces the Netlib's Makefile 
		    with a CMake file.

		    In order to control which library to link to, the user should add to the config script something like:

		      -DBLAS_LIBRARIES="/home/rrossi/mn/SharedBlas/libmy_shared_blas.so" 	\
		      -DLAPACK_LIBRARIES="/usr/lib/lapack/liblapack.so"	\

		    Thus telling explicitely which lib to link too.

                    in order to use the Superlu solver, the user should control explicitely the fortran 
		    mangling by adding something like 
		     
		      -DKRATOS_SUPERLU_FORTRAN_MANGLING="-DAdd_"

		    Available options are (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
		    if the var is not prescribed it will be autommatically added "Add_" which works for gfortran

   MKL: 		-->Needed for MKL solver application
		    Installation path is controlled by

		      -DMKLSOLVER_INCLUDE_DIR="/opt/intel/Compiler/11.1/072/mkl/include"		\
		      -DMKLSOLVER_LIB_DIR="/opt/intel/Compiler/11.1/072/mkl/lib/em64t"		\


   ParMETIS LIBRARY: 	-->Needed for MPI
		    metis library (< 4.0, the latest one is not yet supported) is needed for the mpi version of Kratos.
		    The user should set the following vars:

		      -DPARMETIS_ROOT_DIR="/home/rrossi/Libraries/ParMetis-3.2.0" 		\

			if you want to prescribe exactly the library to which you wish to link the user should set the following vars:
			
		      -DPARMETIS_LIBRARY="/home/rrossi/Libraries/ParMetis-3.2.0" 		\
			
			
   TRILINOS LIBRARY: 	-->Needed for MPI
		    Trilinos provides the basis for kratos mpi-parallel capabilities. 
		    Installation path is controlled by 
		      -DTRILINOS_ROOT="/home/rrossi/Libraries/trilinos-10.8.3/"

		    please note that Kratos uses a feature of trilinos which is present only since version 10.8.3, so please
		    do not use older versions. If this is a problem set to off the MeshingApplication and in this case any 
		    version can be used


