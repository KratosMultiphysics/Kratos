/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "iga_application_variables.h"


namespace Kratos {
namespace Python {

void AddCustomUtilitiesToPython(
    pybind11::module& m)
{
    using namespace pybind11::literals;
    using namespace pybind11;

    class_< IgaFlags > iga_flags = class_< IgaFlags >(m, "IgaFlags")
        .def(init<>())
        ;

    iga_flags.attr("FIX_DISPLACEMENT_X") = IgaFlags::FIX_DISPLACEMENT_X;
    iga_flags.attr("FIX_DISPLACEMENT_Y") = IgaFlags::FIX_DISPLACEMENT_Y;
    iga_flags.attr("FIX_DISPLACEMENT_Z") = IgaFlags::FIX_DISPLACEMENT_Z;
    iga_flags.attr("FIX_ROTATION_X") = IgaFlags::FIX_ROTATION_X;
    iga_flags.attr("FIX_ROTATION_Y") = IgaFlags::FIX_ROTATION_Y;
    iga_flags.attr("FIX_ROTATION_Z") = IgaFlags::FIX_ROTATION_Z;

    pybind11::class_<BrepJsonIO, typename BrepJsonIO::Pointer>(m, "BrepJsonIO")
        .def(pybind11::init<>())
        ;

    pybind11::class_<NurbsBrepModeler, typename NurbsBrepModeler::Pointer>(m, "NurbsBrepModeler")
        .def(pybind11::init<ModelPart&>())
        .def("ImportGeometry", &NurbsBrepModeler::ImportGeometry)
        .def("ImportModelPart", &NurbsBrepModeler::ImportModelPart)
        .def("GetInterfaceConditionsDEM", &NurbsBrepModeler::GetInterfaceConditionsDEM)
        ;
}

} // namespace Python
} // Namespace Kratos
