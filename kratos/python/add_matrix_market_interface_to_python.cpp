//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/matrix_market_interface.h"
#include "includes/ublas_interface.h"
#include "includes/ublas_complex_interface.h"
#ifdef KRATOS_USE_EIGEN_BACKEND
#include "includes/kratos_eigen_interface.h"
#endif
#include "python/add_matrix_market_interface_to_python.h"

namespace Kratos::Python
{
void  AddMatrixMarketInterfaceToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix <Kratos::CompressedMatrix>);
    m.def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix <Kratos::ComplexCompressedMatrix>);
    m.def("WriteMatrixMarketMatrix", WriteMatrixMarketMatrix <Kratos::CompressedMatrix>);
    m.def("WriteMatrixMarketMatrix", WriteMatrixMarketMatrix <Kratos::ComplexCompressedMatrix>);

#ifdef KRATOS_USE_EIGEN_BACKEND
    // System-matrix variant of the Eigen sparse backend, so it can be
    // exercised from python without a uBLAS<->Eigen conversion. There is no
    // Eigen WriteMatrixMarketMatrix: it relies on the uBLAS-only nested
    // iterator1/iterator2 concept, which the Eigen wrapper does not mirror.
    m.def("ReadMatrixMarketMatrix", ReadMatrixMarketMatrix <Kratos::EigenCompressedMatrix<double>>);
#endif

    m.def("ReadMatrixMarketVector", ReadMatrixMarketVector <Kratos::Vector>);
    m.def("ReadMatrixMarketVector", ReadMatrixMarketVector <Kratos::ComplexVector>);
    m.def("WriteMatrixMarketVector", WriteMatrixMarketVector <Kratos::Vector>);
    m.def("WriteMatrixMarketVector", WriteMatrixMarketVector <Kratos::ComplexVector>);

#ifdef KRATOS_USE_EIGEN_BACKEND
    m.def("ReadMatrixMarketVector", ReadMatrixMarketVector <Kratos::EigenVector<double>>);
    m.def("WriteMatrixMarketVector", WriteMatrixMarketVector <Kratos::EigenVector<double>>);
#endif

}

}  // namespace Kratos::Python.

