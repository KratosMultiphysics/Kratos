//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// System includes

// Project includes

//Application includes
#include "custom_python/add_custom_saturation_laws_to_python.h"
#include "custom_saturation/saturation_law.hpp"
#include "custom_saturation/saturation_law_wrapper.hpp"

//Saturation laws
#include "custom_saturation/brooksandcorey_law.hpp"
#include "custom_saturation/vangenuchten_law.hpp"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void  AddCustomSaturationLawsToPython(pybind11::module& m)
{
    // Module local to avoid conflicts with other apps
    py::class_< BrooksAndCoreyLaw, BrooksAndCoreyLaw::Pointer, SaturationLaw >
    (m, "BrooksAndCoreyLaw", py::module_local())
    .def( py::init<>() );
    py::class_< VanGenuchtenLaw, VanGenuchtenLaw::Pointer, SaturationLaw >
    (m, "VanGenuchtenLaw", py::module_local())
    .def( py::init<>() ) ;

}

}  // namespace Python.
}  // namespace Kratos.
