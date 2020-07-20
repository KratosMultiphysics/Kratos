//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//                  Altug Emiroglu, http://github.com/emiroglu
//
//


// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "custom_python/add_custom_utilities_to_python.h"

// Utilities
#include "custom_utilities/rom_residuals_utility.h"

namespace Kratos {
namespace Python {


using namespace pybind11;

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme<SparseSpaceType, LocalSpaceType> BaseSchemeType;

    class_<RomResidualsUtility, typename RomResidualsUtility::Pointer>(m, "RomResidualsUtility")
    .def(init<ModelPart&, Parameters, BaseSchemeType::Pointer>()) // 
    .def("GetResiduals",&RomResidualsUtility::Calculate) //
    ;     

}

} // namespace Python.
} // Namespace Kratos
