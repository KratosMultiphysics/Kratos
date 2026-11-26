GiDPost 2.13
============

gidpost is a set of functions (library) for writing postprocess results for GiD in ASCII, binary compressed or HDF5 format.

This software is copyrighted by CIMNE www.gidsimulation.com. The software can be used freely under the terms described in license.terms, all terms described there apply to all files associated with the software unless explicitly disclaimed in individual files. Particular terms apply to the third party code, "cfortran.h", which has its own distribution policy (please read the "cfortran.doc" for this code).
This description assumes that the reader is familiar with the postprocess terminology.
For further details please check the online help available in GiD ( Postprocess data files chapter).

The library was implemented taking into account two of the must widely used development environments: C/C++ and FORTRAN.

**Look into the *doc* folder for the library manual**

## How to build the gidpost library (Linux - gcc):

Prerequisites: if you want to compile gidpost with HDF5 support, make sure the libraries are installed, for instance with

```shell
$ sudo apt install hdf5-helpers libhdf5-dev libhdf5-cpp-11
$ sudo apt install hdf5-tools hdfview    ;# recommended tools to manage and look into hdf5 files
```

Now, to build gidpost:

```shell
$ cd .../gidpost
$ mkdir build-linux
$ cd build-linux
$ cmake -DENABLE_HDF5=ON -DENABLE_FORTRAN_EXAMPLES=ON ..     ;# gfortran is needed to ENABLE_FORTRAN_EXAMPLES (by default it's off)
# more options are: -DENABLE_SHARED_LIBS=ON -DENABLE_EXAMPLES=ON -DENABLE_PARALLEL_EXAMPLE=ON
# extra flag -DUSE_PKGCONFIG=false   to use pkg-config and pkg_check_modules() instead of cmake's find_package()
# by default, if environemt variable VCPKG_ROOT is defined, then USE_PKGCONFIG = true
$ make
$ cd examples
$ ./testc -help      ;# to view format options
$ ./testc -f hdf5    ;# to write a gid post file in hdf5 format
$ ./testf90          ;# fortran 90 example writing hdf5 gid post file
```

## How to build gidpost as python module (Linux - gcc):

Requirements:
* Python
* [SWIG](https://www.swig.org/)

which can be installed with your favorite package manager.

Once you have configured and build the gidpost library then:

```shell
$ cd .../gidpost/gidpost-swig/
$ make
$ make test
```

For further information please read the inside [README](gidpost-swig/README.md)

## How to build the gidpost library (Linux - nvidia hpc sdk ):

Using nVidia HPC sdk https://developer.nvidia.com/hpc-sdk:

First, remember to add the environment variables to your $HOME/.bashrc: 
https://docs.nvidia.com/hpc-sdk/hpc-sdk-install-guide/index.html#install-linux-end-usr-env-settings

```shell
$ NVARCH=`uname -s`_`uname -m`; export NVARCH
$ NVCOMPILERS=/opt/nvidia/hpc_sdk; export NVCOMPILERS
$ MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/23.5/compilers/man; export MANPATH
$ PATH=$NVCOMPILERS/$NVARCH/23.5/compilers/bin:$PATH; export PATH
```

Once done, log in again and:

```shell
$ cd .../gidpost
$ mkdir build-linux
$ cd build-linux
$ cmake -DENABLE_HDF5=ON -DENABLE_FORTRAN_EXAMPLES=ON \
-DMAKE_CXX_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/23.5/compilers/bin/nvc++ \
-DMAKE_C_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/23.5/compilers/bin/nvc \
-DMAKE_Fortran_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/23.5/compilers/bin/nvfortran \
..          ;# you can also use 'ccmake/cmake-gui ..' with 'Advanced options' enabled.
$ make
$ cd examples
$ ./testc -help      ;# to view format options
$ ./testc -f hdf5    ;# to write a gid post file in hdf5 format
$ ./testf90          ;# fortran 90 example writing hdf5 gid post file
```


ChangeLog:
==========

*From version 2.12 to 2.13*

* CMake changes: 
    * Added support for pkg-config, extra flag `-DUSE_PKGCONFIG=false` to use pkg-config and pkg_check_modules() instead of cmake's find_package().
    * By default, if environment variable `VCPKG_ROOT` is defined, then `USE_PKGCONFIG = true`.
    * Added `shlwapi` dependency for Windows+vcpkg".
* Added prefix to hash functions to avoid collision with POSIX ones.
* Allow writing `*MESH + Elements* without *Coords* block, i.e. these elements use already written global coordinates, for instance:
```c++
  GiD_fBeginMeshColor( fdm, "Spheres", GiD_3D, GiD_Sphere, 1, 1.0, 0.7, 0.3 );
  /* empty or no coordinates = uses global coordinates */
  // GiD_fBeginCoordinates( fdm);
  // GiD_fEndCoordinates( fdm);
  /* elements */
  GiD_fBeginElements( fdm);
  ...
  GiD_fEndElements( fdm);
  GiD_fEndMesh( fdm);
  ```

*From version 2.11 to 2.12*

* CMake changes: HDF5 to be optional, Intel fortran compiler detection,
* some refactoring about GiD_fOpenPostResultFile() in Fortran API, utf-8 encoding.
* GiD_fWriteResultBlock() added unit_name option, and optional NULL parameters.
* corrected some warnings and errors, for instance, when using hdf5 1.8.x.

*From version 2.10 to 2.11*

* HDF5 several corrections when using SWMR mode.
* removed deprecated examples
* some refactoring
* debug build mode provides more error messages.

*From version 2.9 to 2.10*

* HDF5 writes using *Single Write Multiple Read*, when using hdf5 version 1.10 or higher, i.e. GiD and other programs (if they have SWMR support) can read the hdf5 file when gidpost is writing it. https://docs.hdfgroup.org/archive/support/HDF5/Tutor/swmr.html

*From version 2.8 to 2.9*

* GiD_fWriteCoordinates*Block() Added coordinates block function to write array of coordinates with a single call.
* GiD_fWriteElements*Block() Added elements block function to write array of element's connectivities with a single call.
* GiD_fWriteResultBlock() Added result block function to write array of result values with a single call.
* GiD_PostGetVersion() added to return library version, and if hdf5 is enabled, hdf5 library version too.
* GiD_PostIsThreadSafe( format) : function which returns if format is thread-safe, i.e. several files with this format can be opened and written at the same time. For instance, some HDF5 compilations are not thread-safe, and can not be used to write several hdf5 files at the same time from the same program.
* gidpost-swig/ : to use gidpost library as python module.
* several corrected bugs in library.
* updated some C and Fortran examples to use the new GiD_fWrite*Block functions.
* Updated `examples/CMakeLists.txt` to include and build more examples.
* updated documentation.

*From version 2.7 to 2.8*

* CMake updated and corrected to compile fortran90 example.
* Added functions to allow the user to add mesh and result properties / attributes:
  * User defined properties defined inside Mesh or Result blocks
    * HDF5: stored as properties/attributes (Name, value) of the current Mesh/N or Result/N folder
    * ASCII / raw binary: stored as comments

    `# Name: value`

    * Define the macro COMPASSIS_USER_ATTRIBUTES_FORMAT to have it like compassis wants:

    `# ResultUserDefined \"%s\" \"%s\"`      or
    `# ResultUserDefined \"%s\" %s`

* int GiD_WriteMeshUserAttribute(GP_CONST char * Name,GP_CONST char * Value);
* int GiD_fWriteMeshUserAttribute( GiD_FILE fd, GP_CONST char *Name, GP_CONST char *Value );
* int GiD_WriteResultUserAttribute( GP_CONST char *Name, GP_CONST char *Value );
* int GiD_fWriteResultUserAttribute( GiD_FILE fd, GP_CONST char *Name, GP_CONST char *Value );

*From version 2.6 to 2.7*

- Added support for MeshGroups ( BeginMeshGroup(), End..., BeginOnMeshGroup(), End...) in HDF5 as used in GiD, useful for dynamic meshes, refinement, multi-stage simulation, etc.
- All non implemented functions for HDF5 now return -1, i.e. all GiD_fXXX() functions.

*From version 2.5 to 2.6*

- Added support for ComplexMatrix Result types.
- More exhaustive example and some corrections.

*From version 2.4 to 2.5*

- Default format to print real numbers changed from "%g" to "%.9g" increasing the amount of significant digits to be printed in ASCII.
- New function GiD_PostSetFormatReal to allow change the default format of real numbers.

*From version 2.3 to 2.4*

- Allow write OnNurbsSurface results for iso-geometrical analysis

*From version 2.1 to 2.3*

- Allow use HD5 from FORTRAN interface
- Fixed bug with mesh group

*From version 2.0 to 2.1*

- Allow write complex scalar and complex vector results
- Write units of mesh or results
- Enhanced FORTRAN 77 and 90 interface to be more compatible.

*From version 1.7.2 to 2.0*

- New set of functions "GiD_fxxx" that allow specify the file to write, in case of have multiple files of mesh or results.
- Library rewritten in C avoiding C++ to be more compatible linking with other languages.
- Optional define of HDF5 to enable write postprocess files using the HDF5 library

*Old changes*

2008/09/30
* Updated cfortran to release 4.4

2008/02/25
* changing to release 1.8
* M build/Jamroot, source/gidpost.cpp source/gidpostInt.h: removed	C++ dependencies.

2008/01/29
* source/gidpost.cpp:
  - GiD_OpenPostMeshFile should not write a header. Thanks to	Janosch Stascheit <janosch.stascheit@rub.de>
  - invalid access in GiD_BeginMeshGroup: should ensure that a valid	mesh file object is available so etMeshFile() should be called before writting. Thanks to Janosch Stascheit	<janosch.stascheit@rub.de> for providing the fix.

2008/01/18
* source/gidpost.cpp: GiD_BeginOnMeshGroup should be "OnGroup	\"%s\"" instead of "Group \"%s\"". Thanks to Janosch Stascheit <janosch.stascheit@rub.de>
2007/06/13
* source/gidpost.cpp: GiD_WritePlainDefMatrix can be written in	binary files.

2007/06/11
* GiD_Sphere and GiD_Circle

2007/06/11
* source/gidpost.cc: Write2DMatrix is now available on binary	format, Z component is set to 0.

2007/01/23:
* source/gidpostfor.c: bug reported by a user:	GiD_Begin3DMatResult must be ncomp = SetupComponents(6, not 4	Factoring fortran declaration in gidpostfor.h
	Defining symbols for both gcc and icc.

2006/06/25
* doc/*:
* sources/gidpost.cpp,gidpost.h: included pyramid element.

2006/06/25
* tag 1.7: rel_1_7
* examples/testpost.c: fixed a bug, AsciiZipped must has the mesh file in ascii format.

2006/05/25
* gidpost.h: dllexport.

2006/05/18 
* all: GiD_BeginMeshColor, GiD_BeginMeshGroup, GiD_EndMeshGroup, GiD_BeginOnMeshGroup, GiD_EndOnMeshGroup

2005/10/24
* doc:
* gidpost.h : documented GiD_ResultDescriptionDim
2005/10/21
* gidpost* : added support to new type specification for result groups. The type could be Type:dim, where Type is one valid result type and dim is a number giving a valid dimension. Updated fortran interface for recently added functions.

2005/10/07
* all: new version 1.60. "Result Groups" can now be written. The macro USE_CONST can be used to enable (if defined) or desable (if not defined) the use of const char* arguments.

2005/09/22
 * all: revert change from 'char *' to 'const char *'

2005/09/22
* gidpost.cpp: flag_begin_values = 0; must be set in GiD_BeginResultGroup. 

2005/09/21
* gidpost.h:
  * gidpost.cpp: values_location_{min,max} not needed, CBufferValue
  * gidpostInt.cpp: now controls the types of results written within 
  * gidpostInt.h: a ResultGroup.

2005/09/16
* gidpost.cpp: in GiD_WriteVectorModule, when writing if we are	inside a ResultGroup block the module value is ignored.
* gidpostInt.cpp: in CBufferValues::WriteValues, buffer overflow is checked after checking for change in location.

2005/09/15
* gidpost.cpp:
* gidpostInt.cpp:
* gidpostInt.h: fixed a bug when writing results in a "Result Group", vectors can be 3 or 4 component long.
* examples/testpost.c : bug when passing location id in result group sample code.

2005/06/27
* source/gidpost.{h,cpp} including Prism element
* doc/gidpost.{html,subst} including Prism element documentation.
* all: change to version 1.52

2005/05/09
* all: constification, version change from 1.5 to 1.51

2005/01/04
* source/gidpostInt.cpp: fixed a bug when writing 3D vectors.

2005/01/07
* source/gidpost.cpp added const char and a filter to remove double quotes from names which can cause problems within gid. 

2003/07/29
* doc/gidpost.html   : commented the GiD_ResultUnit
* source/gidpost.h   : interface because it is not 
* source/gidpostfor.c: supported yet inside GiD.

2003/07/28
* doc/gidpost.html: updated documentation for release 1.5
* binary/gidpost.lib: binary release for windows
2003/07/15
* gidpost.h:
* gidpost.cpp:
* gidpostfor.c: new function 'GiD_WriteCoordinates2D' (to write	coordinates in 2D but only in ASCII format, the function can be	used in binary but the library provide the 'z' coordinate a zero.

2003/07/15
* gidpost.h:
* gidpost.cpp: removed 'char * UnitName' argument , new functions:

	 GiD_BeginResultHeader, GiD_ResultRange, GiD_ResultComponents,
	 GiD_ResultUnit, GiD_BeginResultGroup, GiD_ResultDescription,
	 GiD_ResultValues, GiD_FlushPostFile. Validation in debug mode.

* gidpostInt.h:
* gidpostInt.cpp: new member

	CPostFile::WriteValues( int id, int n, double * buffer), new class CBufferValues to write the values in a result group, validation in debug mode.

* gidpostfor.c: updated fortran interface for the new functionality.
* testpost.c:
* testpostfor.f: added test code for the new functionality.
