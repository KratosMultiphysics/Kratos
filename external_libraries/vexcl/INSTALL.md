# Installation

Since VexCL is header-only library, installation is straightforward: you just
need to copy vexcl folder somewhere and tell your compiler to scan it for
include files.

* [Gentoo Linux ebuild](https://github.com/ddemidov/ebuilds/blob/master/dev-util/vexcl)
* [Arch Linux PKGBUILD](https://aur.archlinux.org/packages/vexcl-git)


## Building boost on OSX

_Here are instructions that worked for at least one of VexCL users on a MacOSX
system (thanks, @d-meiser!):_

Because VexCL requires c++11 features I had to tweak the boost libraries on my
system a bit. The following instructions worked for me. I build boost version
1.53. Download boost sources from sourceforge
(http://sourceforge.net/projects/boost/files/boost/1.53.0/). In the boost
source directory do:

    ./bootstrap.sh --with-toolset=clang --prefix=${BOOST_DIR}
    ./b2 link=shared threading=multi variant=release -s NO_COMPRESSION=1 \
        --layout=system --without-mpi --prefix=${BOOST_DIR} \
        toolset=clang \
        cxxflags="-std=c++ -stdlib=libc++" \
        linkflags="-stdlib=libc++" \
        install

`${BOOST_DIR}` is the boost installation directory. I install in my home
directory just so it's easier for me to uninstall and so there are no
unintended interactions with the system boost libraries. The `cxxflags` and
`linkflags` were necessary to run the tests. Without them the boost unit test
framework segfaults


## Making cmake find boost headers and libraries

It took me a while to figure out how to make cmake find the right boost
libraries instead of the system libraries. The following command line worked
for me

     cmake \
        -DBoost_NO_BOOST_CMAKE:BOOL=TRUE \
        -DBoost_NO_SYSTEM_PATHS:BOOL=TRUE \
        -DBOOST_ROOT=${BOOST_DIR} \
        -DUSE_LIBCPP:BOOL=TRUE \
        ${VEXCL_SOURCE_DIR}

`${BOOST_DIR}` is the location of your boost library, e.g. the one installed
here. `${VEXCL_SOURCE_DIR}` is the location of the VexCL sources. Note that
without the `-DBoost_NO_BOOST_CMAKE:BOOL=TRUE` flag the `-DBOOST_ROOT` option
was being ignored. Without `-DUSE_LIBCPP:BOOL=TRUE` the `tuple` header was not
found.
