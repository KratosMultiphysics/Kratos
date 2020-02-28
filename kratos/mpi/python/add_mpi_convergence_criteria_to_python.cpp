
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:     Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "mpi/python/add_mpi_convergence_criteria_to_python.h"

#include "mpi/spaces/amgcl_mpi_space.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/convergencecriterias/convergence_criteria.h"

// ConvergenceCriterias
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"


namespace Kratos {
namespace Python {

void AddMPIConvergenceCriteriasToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Amgcl type definitions
    typedef amgcl::backend::builtin<double> Backend;
    typedef amgcl::mpi::distributed_matrix<Backend> amgcl_mpi_matrix;
    typedef typename Backend::vector amgcl_mpi_vector;

    // Type definitions
    typedef AmgclMPISpace<amgcl_mpi_matrix, amgcl_mpi_vector> MPISparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> MPILocalSpaceType;

    // Base ConvergenceCriteria
    typedef ConvergenceCriteria< MPISparseSpaceType, MPILocalSpaceType > MPIConvergenceCriteria;

    typedef typename ConvergenceCriteria< MPISparseSpaceType, MPILocalSpaceType > ::Pointer MPIConvergenceCriteriaPointer;

    py::class_< MPIConvergenceCriteria, MPIConvergenceCriteriaPointer > (m,"MPIConvergenceCriteria")
        .def(py::init<>())
        .def("SetActualizeRHSFlag", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::SetActualizeRHSFlag)
        .def("GetActualizeRHSflag", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::GetActualizeRHSflag)
        .def("PreCriteria", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::PreCriteria)
        .def("PostCriteria", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::PostCriteria)
        .def("Initialize", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::Initialize)
        .def("InitializeSolutionStep", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::InitializeSolutionStep)
        .def("FinalizeSolutionStep", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::FinalizeSolutionStep)
        .def("Check", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::Check)
        .def("SetEchoLevel", &ConvergenceCriteria<MPISparseSpaceType, MPILocalSpaceType >::SetEchoLevel)
        ;

    // ConvergenceCriterias
    py::class_<And_Criteria<MPISparseSpaceType, MPILocalSpaceType >,
            typename And_Criteria<MPISparseSpaceType, MPILocalSpaceType >::Pointer,
            MPIConvergenceCriteria>(m,"MPIAndCriteria")
        .def(py::init<MPIConvergenceCriteriaPointer, MPIConvergenceCriteriaPointer > ());

    py::class_<Or_Criteria<MPISparseSpaceType, MPILocalSpaceType >,
            typename Or_Criteria<MPISparseSpaceType, MPILocalSpaceType >::Pointer,
            MPIConvergenceCriteria>(m,"MPIOrCriteria")
        .def(py::init<MPIConvergenceCriteriaPointer, MPIConvergenceCriteriaPointer > ());
}

} // namespace Python.
} // Namespace Kratos
