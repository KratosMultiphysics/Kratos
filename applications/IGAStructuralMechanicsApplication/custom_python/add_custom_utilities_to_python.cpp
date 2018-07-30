//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/iga_flags.h"

namespace Kratos
{

namespace Python
{


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    // Adding the flags that can be used for application of enforcements on conditions
    class_< IGAFlags > iga_flags = class_< IGAFlags >(m, "IGAFlags")
        .def(init<>())
        ;

    iga_flags.attr("FIX_DISPLACEMENT_X") = IGAFlags::FIX_DISPLACEMENT_X;
    iga_flags.attr("FIX_DISPLACEMENT_Y") = IGAFlags::FIX_DISPLACEMENT_Y;
    iga_flags.attr("FIX_DISPLACEMENT_Z") = IGAFlags::FIX_DISPLACEMENT_Z;
    iga_flags.attr("FIX_ROTATION_X") = IGAFlags::FIX_ROTATION_X;
    iga_flags.attr("FIX_ROTATION_Y") = IGAFlags::FIX_ROTATION_Y;
    iga_flags.attr("FIX_ROTATION_Z") = IGAFlags::FIX_ROTATION_Z;

}

}  // namespace Python.

} // Namespace Kratos
