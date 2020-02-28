
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
#include "mpi/python/add_mpi_schemes_to_python.h"

#include "mpi/spaces/amgcl_mpi_space.h"
#include "spaces/ublas_space.h"

#include "solving_strategies/schemes/scheme.h"

// Schemes
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"


namespace Kratos {
namespace Python {

void AddMPISchemesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Amgcl type definitions
    typedef amgcl::backend::builtin<double> Backend;
    typedef amgcl::mpi::distributed_matrix<Backend> amgcl_mpi_matrix;
    typedef typename Backend::vector amgcl_mpi_vector;

    // Type definitions
    typedef AmgclMPISpace<amgcl_mpi_matrix, amgcl_mpi_vector> MPISparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> MPILocalSpaceType;
    typedef Scheme< MPISparseSpaceType, MPILocalSpaceType > MPISchemeType;

    // Base Scheme
    py::class_< MPISchemeType, typename MPISchemeType::Pointer >(m, "MPIScheme")
        .def(py::init< >() )
        .def("Initialize", &MPISchemeType::Initialize )
        .def("SchemeIsInitialized", &MPISchemeType::SchemeIsInitialized )
        .def("ElementsAreInitialized", &MPISchemeType::ElementsAreInitialized )
        .def("ConditionsAreInitialized", &MPISchemeType::ConditionsAreInitialized )
        .def("InitializeElements", &MPISchemeType::InitializeElements )
        .def("InitializeConditions", &MPISchemeType::InitializeConditions )
        .def("InitializeSolutionStep", &MPISchemeType::InitializeSolutionStep )
        .def("FinalizeSolutionStep", &MPISchemeType::FinalizeSolutionStep )
        .def("InitializeNonLinIteration", &MPISchemeType::InitializeNonLinIteration )
        .def("FinalizeNonLinIteration", &MPISchemeType::FinalizeNonLinIteration )
        .def("Predict", &MPISchemeType::Predict )
        .def("Update", &MPISchemeType::Update )
        .def("CalculateOutputData", &MPISchemeType::CalculateOutputData )
        .def("Clean", &MPISchemeType::Clean )
        .def("Clear", &MPISchemeType::Clear )
        .def("Check", &MPISchemeType::Check )
        ;

    // Schemes
    // there is an issue with the DofUpdater
    // py::class_<ResidualBasedIncrementalUpdateStaticScheme< MPISparseSpaceType, MPILocalSpaceType>,
    //     typename ResidualBasedIncrementalUpdateStaticScheme< MPISparseSpaceType, MPILocalSpaceType>::Pointer,
    //     MPISchemeType >(m,"MPIResidualBasedIncrementalUpdateStaticScheme")
    //        .def(py::init< >() );

    // py::class_<ResidualBasedBossakDisplacementScheme< MPISparseSpaceType, MPILocalSpaceType>,
    //     typename ResidualBasedBossakDisplacementScheme< MPISparseSpaceType, MPILocalSpaceType>::Pointer,
    //     MPISchemeType > (m,"MPIResidualBasedBossakDisplacementScheme")
    //     .def(py::init<double >())
    //     ;

    // py::class_<ResidualBasedBDFDisplacementScheme< MPISparseSpaceType, MPILocalSpaceType>,
    //     typename ResidualBasedBDFDisplacementScheme< MPISparseSpaceType, MPILocalSpaceType>::Pointer,
    //     MPISchemeType > (m,"MPIResidualBasedBDFDisplacementScheme")
    //     .def(py::init<  >())
    //     .def(py::init <const std::size_t>())
    //     ;

    // py::class_<ResidualBasedBDFCustomScheme< MPISparseSpaceType, MPILocalSpaceType>,
    //     typename ResidualBasedBDFCustomScheme< MPISparseSpaceType, MPILocalSpaceType>::Pointer,
    //     MPISchemeType > (m,"MPIResidualBasedBDFCustomScheme")
    //     .def(py::init<  >())
    //     .def(py::init <const std::size_t>())
    //     .def(py::init <const std::size_t, Parameters>())
    //     ;
}

} // namespace Python.
} // Namespace Kratos
