//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Schemes
#include "custom_strategies/schemes/residualbased_incrementalupdate_wetting_scheme.h"


namespace Kratos
{

namespace Python
{

  void AddCustomStrategiesToPython(pybind11::module& m)
  {
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    // Schemes
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    typedef ResidualBasedIncrementalUpdateWettingScheme<SparseSpaceType, LocalSpaceType> ResidualBasedIncrementalUpdateWettingScheme;
    py::class_<
        ResidualBasedIncrementalUpdateWettingScheme,
        typename ResidualBasedIncrementalUpdateWettingScheme::Pointer,
        BaseSchemeType>
        (m, "ResidualBasedIncrementalUpdateWettingScheme")
        .def(py::init<>())
        .def(py::init<Process::Pointer>())
        ;

  }

}  // namespace Python.

} // Namespace Kratos
