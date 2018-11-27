//
// Author: Salva Latorre latorre@cimne.upc.edu
//

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/flex_wrapper.h"
#include "custom_utilities/nvidia_flex_pre_utilities.h"

namespace Kratos {

    namespace Python {

        namespace py = pybind11;

        void AddCustomUtilitiesToPython(pybind11::module& m) {

        py::class_<FlexWrapper, FlexWrapper::Pointer>(m, "FlexWrapper")
            .def(py::init<ModelPart&, ModelPart&, ParticleCreatorDestructor&, Parameters>())
            .def("UpdateFlex", &FlexWrapper::UpdateFlex)
            .def("TransferDataFromFlexToKratos", &FlexWrapper::TransferDataFromFlexToKratos)
            .def("SolveTimeSteps", &FlexWrapper::SolveTimeSteps)
            ;

        py::class_<NvidiaFlexPreUtilities, NvidiaFlexPreUtilities::Pointer >(m, "NvidiaFlexPreUtilities")
        .def(py::init<>())
        .def(py::init<ModelPart&>())
        .def("RemoveSpheresInitiallyIndentedWithFEM", &NvidiaFlexPreUtilities::RemoveSpheresInitiallyIndentedWithFEM)
        ;

        }
    }  // namespace Python
} // Namespace Kratos
