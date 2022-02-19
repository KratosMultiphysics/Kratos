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
//
//  Contributors:   Altug Emiroglu, http://github.com/emiroglu
//
//


// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"

// Utilities

// Project includes

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rom_residuals_utility.h"
#include "custom_utilities/rom_auxiliary_utilities.h"

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

    class_<RomAuxiliaryUtilities>(m, "RomAuxiliaryUtilities")
        .def_static("SetHRomComputingModelPart", &RomAuxiliaryUtilities::SetHRomComputingModelPart)
        .def_static("SetHRomVolumetricVisualizationModelPart", &RomAuxiliaryUtilities::SetHRomVolumetricVisualizationModelPart)
        .def_static("GetHRomConditionParentsIds", &RomAuxiliaryUtilities::GetHRomConditionParentsIds)
        .def_static("GetHRomMinimumConditionsIds", &RomAuxiliaryUtilities::GetHRomMinimumConditionsIds)
        .def_static("ProjectRomSolutionIncrementToNodes", &RomAuxiliaryUtilities::ProjectRomSolutionIncrementToNodes)
        ;
}

} // namespace Python.
} // Namespace Kratos
