//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif                                    // KRATOS_USE_AMATRIX

// External includes
#include "pybind11/pybind11.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "processes/process.h"
#include "custom_utilities/solver_settings.h"

#include "spaces/ublas_space.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

//schemes
#include "custom_strategies/schemes/residualbased_predictorcorrector_velocity_bossak_scheme_turbulent_with_no_node_displacement.h"


namespace Kratos
{
namespace Python
{

void AddCustomStrategiesToPython(pybind11::module &m)
{
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    py::class_<
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement<SparseSpaceType, LocalSpaceType>,
        typename ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement<SparseSpaceType, LocalSpaceType>::Pointer,
        ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent<SparseSpaceType, LocalSpaceType>>(m, "ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulentWithNoNodeDisplacement")
    .def(py::init<double, double, unsigned int, Process::Pointer>())
    .def(py::init<double, double, unsigned int, double, Process::Pointer>())
    .def(py::init<double, double, unsigned int>())                        // constructor without a turbulence model
    .def(py::init<double, unsigned int, const Kratos::Variable<int> &>()) // constructor without a turbulence model for periodic boundary conditions
    ;

}

} // namespace Python.

} // Namespace Kratos
