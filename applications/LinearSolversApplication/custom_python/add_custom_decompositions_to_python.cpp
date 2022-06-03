//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "spaces/ublas_space.h"
#include "utilities/dense_svd_decomposition.h"

// Application includes
#include "custom_python/add_custom_decompositions_to_python.h"
#include "custom_decompositions/eigen_dense_bdc_svd_decomposition.h"
#include "custom_decompositions/eigen_dense_jacobi_svd_decomposition.h"

namespace Kratos {
namespace Python {

void AddCustomDecompositionsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef typename DenseSpaceType::MatrixType MatrixType;
    typedef typename DenseSpaceType::VectorType VectorType;
    typedef DenseSingularValueDecomposition<DenseSpaceType> BaseSVDType;

    typedef EigenDenseBDCSVD<DenseSpaceType> BDCSVDType;
    py::class_<BDCSVDType, typename BDCSVDType::Pointer, BaseSVDType>(m,"EigenDenseBDCSVD")
        .def(py::init<>())
        .def("Compute", [](BDCSVDType& rBDCSVD, MatrixType& rInputMatrix, Parameters Settings){rBDCSVD.Compute(rInputMatrix, Settings);})
        .def("Compute", [](BDCSVDType& rBDCSVD, MatrixType& rInputMatrix, VectorType& rVectorS, MatrixType& rMatrixU, MatrixType& rMatrixV, Parameters Settings){rBDCSVD.Compute(rInputMatrix, rVectorS, rMatrixU, rMatrixV, Settings);})
        .def("SingularValues", &BDCSVDType::SingularValues)
        .def("MatrixU", &BDCSVDType::MatrixU)
        .def("MatrixV", &BDCSVDType::MatrixV)
        .def("Rank", &BDCSVDType::Rank)
        .def("NonZeroSingularValues", &BDCSVDType::NonZeroSingularValues)
        ;

    typedef EigenDenseJacobiSVD<DenseSpaceType> JacobiSVDType;
    py::class_<JacobiSVDType, typename JacobiSVDType::Pointer, BaseSVDType>(m,"EigenDenseJacobiSVD")
        .def(py::init<>())
        .def("Compute", [](JacobiSVDType& rJacobiSVD, MatrixType& rInputMatrix, Parameters Settings){rJacobiSVD.Compute(rInputMatrix, Settings);})
        .def("Compute", [](JacobiSVDType& rJacobiSVD, MatrixType& rInputMatrix, VectorType& rVectorS, MatrixType& rMatrixU, MatrixType& rMatrixV, Parameters Settings){rJacobiSVD.Compute(rInputMatrix, rVectorS, rMatrixU, rMatrixV, Settings);})
        .def("SingularValues", &JacobiSVDType::SingularValues)
        .def("MatrixU", &JacobiSVDType::MatrixU)
        .def("MatrixV", &JacobiSVDType::MatrixV)
        .def("Rank", &JacobiSVDType::Rank)
        .def("NonZeroSingularValues", &JacobiSVDType::NonZeroSingularValues)
        ;

}

} // namespace Python

} // namespace Kratos
