//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "solving_strategies/convergence_accelerators/convergence_accelerator.h"
#include "spaces/ublas_space.h"
#include "utilities/dense_qr_decomposition.h"
#include "utilities/dense_svd_decomposition.h"

// Application includes
#include "custom_python/add_convergence_accelerators_to_python.h"
#include "custom_utilities/constant_relaxation_convergence_accelerator.h"
#include "custom_utilities/mvqn_convergence_accelerator.hpp"
#include "custom_utilities/mvqn_randomized_svd_convergence_accelerator.h"
#include "custom_utilities/mvqn_recursive_convergence_accelerator.hpp"
#include "custom_utilities/aitken_convergence_accelerator.hpp"
#include "custom_utilities/ibqn_mvqn_convergence_accelerator.h"
#include "custom_utilities/ibqn_mvqn_randomized_svd_convergence_accelerator.h"

namespace Kratos
{

namespace Python
{

void AddConvergenceAcceleratorsToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef ConvergenceAccelerator<SparseSpaceType, DenseSpaceType> BaseConvergenceAcceleratorType;

    // Constant relaxation convergence accelerator
    typedef ConstantRelaxationConvergenceAccelerator<SparseSpaceType, DenseSpaceType> ConstantRelaxationConvergenceAcceleratorType;
    typedef typename ConstantRelaxationConvergenceAcceleratorType::Pointer ConstantRelaxationConvergenceAcceleratorPointerType;
    py::class_<ConstantRelaxationConvergenceAcceleratorType, ConstantRelaxationConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "ConstantRelaxationConvergenceAccelerator")
        .def(py::init<double>())
        .def(py::init<Parameters>())
    ;

    // Aitken convergence accelerator
    typedef AitkenConvergenceAccelerator<SparseSpaceType, DenseSpaceType> AitkenConvergenceAcceleratorType;
    typedef typename AitkenConvergenceAcceleratorType::Pointer AitkenConvergenceAcceleratorPointerType;
    py::class_<AitkenConvergenceAcceleratorType, AitkenConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "AitkenConvergenceAccelerator")
        .def(py::init<double>())
        .def(py::init<Parameters>())
    ;

    // MVQN convergence accelerator
    typedef MVQNFullJacobianConvergenceAccelerator<SparseSpaceType, DenseSpaceType> MVQNFullJacobianConvergenceAcceleratorType;
    typedef typename MVQNFullJacobianConvergenceAcceleratorType::Pointer MVQNFullJacobianConvergenceAcceleratorPointerType;
    py::class_<MVQNFullJacobianConvergenceAcceleratorType, MVQNFullJacobianConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "MVQNFullJacobianConvergenceAccelerator")
        .def(py::init<Parameters>())
    ;

    // MVQN randomized SVD convergence accelerator
    typedef typename DenseQRDecomposition<DenseSpaceType>::Pointer DenseQRPointerType;
    typedef typename DenseSingularValueDecomposition<DenseSpaceType>::Pointer DenseSVDPointerType;
    typedef MVQNRandomizedSVDConvergenceAccelerator<SparseSpaceType, DenseSpaceType> MVQNRandomizedSVDConvergenceAcceleratorType;
    typedef typename MVQNRandomizedSVDConvergenceAcceleratorType::Pointer MVQNRandomizedSVDConvergenceAcceleratorPointerType;
    py::class_<MVQNRandomizedSVDConvergenceAcceleratorType, MVQNRandomizedSVDConvergenceAcceleratorPointerType, MVQNFullJacobianConvergenceAcceleratorType>(m, "MVQNRandomizedSVDConvergenceAccelerator")
        .def(py::init<DenseQRPointerType, DenseSVDPointerType, Parameters>())
    ;

    // IBQN-MVQN convergence accelerator
    typedef IBQNMVQNConvergenceAccelerator<SparseSpaceType, DenseSpaceType> IBQNMVQNConvergenceAcceleratorType;
    typedef typename IBQNMVQNConvergenceAcceleratorType::Pointer IBQNMVQNConvergenceAcceleratorPointerType;
    py::class_<IBQNMVQNConvergenceAcceleratorType, IBQNMVQNConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "IBQNMVQNConvergenceAccelerator")
        .def(py::init<Parameters>())
        .def("UpdateSolutionLeft", &IBQNMVQNConvergenceAcceleratorType::UpdateSolutionLeft)
        .def("UpdateSolutionRight", &IBQNMVQNConvergenceAcceleratorType::UpdateSolutionRight)
    ;

    // IBQN-MVQN randomized SVD convergence accelerator
    typedef IBQNMVQNRandomizedSVDConvergenceAccelerator<SparseSpaceType, DenseSpaceType> IBQNMVQNRandomizedSVDConvergenceAcceleratorType;
    typedef typename IBQNMVQNRandomizedSVDConvergenceAcceleratorType::Pointer IBQNMVQNRandomizedSVDConvergenceAcceleratorPointerType;
    py::class_<IBQNMVQNRandomizedSVDConvergenceAcceleratorType, IBQNMVQNRandomizedSVDConvergenceAcceleratorPointerType, IBQNMVQNConvergenceAcceleratorType>(m, "IBQNMVQNRandomizedSVDConvergenceAccelerator")
        .def(py::init<DenseQRPointerType, DenseSVDPointerType, Parameters>())
    ;

    // MVQN recursive convergence accelerator
    typedef MVQNRecursiveJacobianConvergenceAccelerator<SparseSpaceType, DenseSpaceType> MVQNRecursiveJacobianConvergenceAcceleratorType;
    typedef typename MVQNRecursiveJacobianConvergenceAcceleratorType::Pointer MVQNRecursiveJacobianConvergenceAcceleratorPointerType;
    py::class_<MVQNRecursiveJacobianConvergenceAcceleratorType, MVQNRecursiveJacobianConvergenceAcceleratorPointerType, BaseConvergenceAcceleratorType>(m, "MVQNRecursiveJacobianConvergenceAccelerator")
        .def(py::init<Parameters>())
        .def(py::init<double, unsigned int, double>())
    ;
}

}  // namespace Python.

} // Namespace Kratos
