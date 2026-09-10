//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// The Kratos dense/sparse linear-algebra interface (Matrix, Vector,
// DenseMatrix, BoundedMatrix, array_1d-compatible operations, CompressedMatrix,
// prod/trans/noalias/..., the proxies and the lazy factories), resolved at
// configure time through the KRATOS_LINEAR_ALGEBRA_BACKEND CMake option:
// "eigen" (which defines KRATOS_USE_EIGEN_BACKEND) selects the Eigen-backed
// types and the pure-Eigen implementation of the uBLAS idiom, "ublas" (the
// default) keeps the boost::numeric::ublas types. This is the header to
// include; it is the counterpart of spaces/default_spaces.h for the containers.
//
// NOTE: the backend define changes the meaning of the aliases and hence the
// mangled names of everything instantiated with them. The define is set
// globally by the root CMakeLists.txt; never mix binaries compiled with
// different KRATOS_LINEAR_ALGEBRA_BACKEND values.

// WIP: TODO

//#ifdef KRATOS_USE_EIGEN_BACKEND
//#include "includes/eigen_interface.h"
//#else
#include "includes/ublas_interface.h"
//#endif
