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
    .def("GetProjectedResidualsOntoJPhi",&RomResidualsUtility::GetProjectedResidualsOntoJPhi) //
    ;

    class_<RomAuxiliaryUtilities>(m, "RomAuxiliaryUtilities")
        .def_static("SetHRomComputingModelPart", &RomAuxiliaryUtilities::SetHRomComputingModelPart)
        .def_static("SetHRomComputingModelPartWithLists", &RomAuxiliaryUtilities::SetHRomComputingModelPartWithLists)
        .def_static("SetHRomComputingModelPartWithNeighbours", &RomAuxiliaryUtilities::SetHRomComputingModelPartWithNeighbours)
        .def_static("SetHRomVolumetricVisualizationModelPart", &RomAuxiliaryUtilities::SetHRomVolumetricVisualizationModelPart)
        .def_static("GetHRomConditionParentsIds", [](ModelPart& rModelPart, const std::vector<IndexType>& rConditionIds) {
                return RomAuxiliaryUtilities::GetHRomConditionParentsIds(rModelPart, rConditionIds);})
        .def_static("GetHRomConditionParentsIds", [](const ModelPart& rModelPart, const std::map<std::string, std::map<IndexType, double>>& rHRomWeights) {
                return RomAuxiliaryUtilities::GetHRomConditionParentsIds(rModelPart, rHRomWeights);})
        .def_static("GetNodalNeighbouringElementIdsNotInHRom", &RomAuxiliaryUtilities::GetNodalNeighbouringElementIdsNotInHRom)
        .def_static("GetNodalNeighbouringElementIds", [](Kratos::ModelPart& rModelPart, Kratos::ModelPart& rGivenModelPart) {
                return Kratos::RomAuxiliaryUtilities::GetNodalNeighbouringElementIds(rModelPart, rGivenModelPart);})
        .def_static("GetNodalNeighbouringElementIds", [](Kratos::ModelPart& rModelPart, const std::vector<Kratos::IndexType>& rNodeIds, bool retrieveSingleNeighbour) {
                return Kratos::RomAuxiliaryUtilities::GetNodalNeighbouringElementIds(rModelPart, rNodeIds, retrieveSingleNeighbour);})
        .def_static("GetNodalNeighbouringConditionIds", &RomAuxiliaryUtilities::GetNodalNeighbouringConditionIds)
        .def_static("GetConditionIdsNotInHRomModelPart", &RomAuxiliaryUtilities::GetConditionIdsNotInHRomModelPart)
        .def_static("GetElementIdsNotInHRomModelPart", &RomAuxiliaryUtilities::GetElementIdsNotInHRomModelPart)
        .def_static("GetHRomMinimumConditionsIds", &RomAuxiliaryUtilities::GetHRomMinimumConditionsIds)
        .def_static("ProjectRomSolutionIncrementToNodes", &RomAuxiliaryUtilities::ProjectRomSolutionIncrementToNodes)
        .def_static("ProjectToReducedBasis", &RomAuxiliaryUtilities::ProjectToReducedBasis)
        .def_static("GetElementIdsInModelPart", &RomAuxiliaryUtilities::GetElementIdsInModelPart)
        .def_static("GetConditionIdsInModelPart", &RomAuxiliaryUtilities::GetConditionIdsInModelPart)
        ;
}

} // namespace Python.
} // Namespace Kratos
