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
//


// System includes

// External includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


// Project includes
#include "includes/define.h"

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/rom_residuals_utility.h"
#include "custom_utilities/rom_auxiliary_utilities.h"
#include "custom_utilities/base_encoder_decoder.h"
//#include "custom_utilities/global_linear_encoder_decoder.h"

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
    .def("GetProjectedResidualsOntoPhi",&RomResidualsUtility::GetProjectedResidualsOntoPhi) //
    .def("GetProjectedResidualsOntoPsi",&RomResidualsUtility::GetProjectedResidualsOntoPsi) //
    ;

    class_<RomAuxiliaryUtilities>(m, "RomAuxiliaryUtilities")
        .def_static("SetHRomComputingModelPart", &RomAuxiliaryUtilities::SetHRomComputingModelPart)
        .def_static("SetHRomVolumetricVisualizationModelPart", &RomAuxiliaryUtilities::SetHRomVolumetricVisualizationModelPart)
        .def_static("GetHRomConditionParentsIds", &RomAuxiliaryUtilities::GetHRomConditionParentsIds)
        .def_static("GetNodalNeighbouringElementIdsNotInHRom", &RomAuxiliaryUtilities::GetNodalNeighbouringElementIdsNotInHRom)
        .def_static("GetNodalNeighbouringElementIds", &RomAuxiliaryUtilities::GetNodalNeighbouringElementIds)
        .def_static("GetConditionIdsNotInHRomModelPart", &RomAuxiliaryUtilities::GetConditionIdsNotInHRomModelPart)
        .def_static("GetElementIdsNotInHRomModelPart", &RomAuxiliaryUtilities::GetElementIdsNotInHRomModelPart)
        .def_static("GetHRomMinimumConditionsIds", &RomAuxiliaryUtilities::GetHRomMinimumConditionsIds)
        .def_static("ProjectRomSolutionIncrementToNodes", &RomAuxiliaryUtilities::ProjectRomSolutionIncrementToNodes)
        .def_static("GetElementIdsInModelPart", &RomAuxiliaryUtilities::GetElementIdsInModelPart)
        .def_static("GetConditionIdsInModelPart", &RomAuxiliaryUtilities::GetConditionIdsInModelPart)
        ;

    class_<BaseEncoderDecoder, typename BaseEncoderDecoder::Pointer>(m, "BaseEncoderDecoder")
    .def(init<>()) //
    //.def(init<Parameters>()) //
    ;

    class_<GlobalLinearEncoderDecoder, typename GlobalLinearEncoderDecoder::Pointer, BaseEncoderDecoder>(m, "GlobalLinearEncoderDecoder")
    .def(init<>()) //
    .def("SetNodalBasis",&GlobalLinearEncoderDecoder::SetNodalBasis)
    //.def(init<Parameters>()) //
    ;

}

} // namespace Python.
} // Namespace Kratos
