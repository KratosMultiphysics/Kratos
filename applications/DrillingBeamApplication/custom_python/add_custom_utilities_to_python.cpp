//
//  Main authors:    Miguel Angel Celigueta   maceli@cimne.upc.edu
//
//

// System includes

// External includes

// Project includes

#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/beam_solids_utility.h"
#include "custom_utilities/post_utilities.h"

#include "includes/define.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void AddCustomUtilitiesToPython(pybind11::module &m)
{

    py::class_<PostUtilities, PostUtilities::Pointer>(m, "PostUtilities", py::module_local())
        .def(py::init<>())
        .def("CreateSkinForBeam", &PostUtilities::CreateSkinForBeam)
        ;

    py::class_<BeamSolidsUtility, BeamSolidsUtility::Pointer>(m, "BeamSolidsUtility", py::module_local())
        .def(py::init<ModelPart&>())
        .def("ComputeEdgePointOfBeamSolids", &BeamSolidsUtility::ComputeEdgePointOfBeamSolids)        
        .def("ComputeDirectiveLineOfBeamSolids", &BeamSolidsUtility::ComputeDirectiveLineOfBeamSolids)
        .def("ComputeReactionsOfBeamSolids", &BeamSolidsUtility::ComputeReactionsOfBeamSolids)        
        ;
}

} // namespace Python.

} // Namespace Kratos