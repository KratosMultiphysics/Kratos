// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "modeler/modeler.h"
#include "custom_python/add_custom_modelers_to_python.h"

#ifdef INCLUDE_MMG
#include "custom_modelers/mmg/mmg_modeler.h"
#endif

namespace Kratos::Python
{

namespace py = pybind11;

void AddCustomModelersToPython(py::module& m)
{
#ifdef INCLUDE_MMG
    py::class_<MmgModeler<MMGLibrary::MMG2D>,
               MmgModeler<MMGLibrary::MMG2D>::Pointer,
               Modeler>(m, "MmgModeler2D")
        .def(py::init<Model&, Parameters>())
        .def("SetupModelPart", &MmgModeler<MMGLibrary::MMG2D>::SetupModelPart)
        ;

    py::class_<MmgModeler<MMGLibrary::MMG3D>,
               MmgModeler<MMGLibrary::MMG3D>::Pointer,
               Modeler>(m, "MmgModeler3D")
        .def(py::init<Model&, Parameters>())
        .def("SetupModelPart", &MmgModeler<MMGLibrary::MMG3D>::SetupModelPart)
        ;

    py::class_<MmgModeler<MMGLibrary::MMGS>,
               MmgModeler<MMGLibrary::MMGS>::Pointer,
               Modeler>(m, "MmgModelerSurface")
        .def(py::init<Model&, Parameters>())
        .def("SetupModelPart", &MmgModeler<MMGLibrary::MMGS>::SetupModelPart)
        ;
#endif
}

} // namespace Kratos::Python
