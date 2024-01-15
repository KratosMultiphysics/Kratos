//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"

#include "processes/process.h"
#include "custom_utilities/solver_settings.h"

#include "spaces/ublas_space.h"

//schemes
#include "custom_strategies/residualbased_incrementalupdate_static_scheme_mod.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{

void AddCustomStrategiesToPython(pybind11::module &m)
{
    namespace py = pybind11;
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
   
    py::class_< ResidualBasedIncrementalUpdateStaticSchemeMod< SparseSpaceType, LocalSpaceType>,
        typename ResidualBasedIncrementalUpdateStaticSchemeMod< SparseSpaceType, LocalSpaceType>::Pointer,
        BaseSchemeType >
        (m, "ResidualBasedIncrementalUpdateStaticSchemeMod")
        .def(py::init<Parameters >() )
        .def(py::init< >()
        );

}

} // namespace Python.

} // Namespace Kratos
