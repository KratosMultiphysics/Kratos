//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"
#include "spaces/default_spaces.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/potential_flow_residualbased_incrementalupdate_static_scheme.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{

void AddCustomStrategiesToPython(pybind11::module &m)
{
    namespace py = pybind11;
    typedef TDefaultSparseSpace<double> SparseSpaceType;
    typedef TDefaultDenseSpace<double> LocalSpaceType;
    typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedIncrementalUpdateStaticSchemeType;

    //********************************************************************
    typedef PotentialFlowResidualBasedIncrementalUpdateStaticScheme<SparseSpaceType, LocalSpaceType> PotentialFlowResidualBasedIncrementalUpdateStaticScheme;
    py::class_<PotentialFlowResidualBasedIncrementalUpdateStaticScheme, typename PotentialFlowResidualBasedIncrementalUpdateStaticScheme::Pointer, ResidualBasedIncrementalUpdateStaticSchemeType>(m, "PotentialFlowResidualBasedIncrementalUpdateStaticScheme")
        .def(py::init<Parameters>())
        .def(py::init< >() )
        ;
}

} // namespace Python.

} // Namespace Kratos
