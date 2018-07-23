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


// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/poro_condition_utilities.hpp"
#include "custom_utilities/poro_element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "custom_utilities/fracture_propagation_3D_utilities.hpp"
#include "custom_utilities/fracture_propagation_2D_utilities.hpp"
#include "custom_utilities/nonlocal_damage_utilities.hpp"
#include "custom_utilities/nonlocal_damage_2D_utilities.hpp"
#include "custom_utilities/nonlocal_damage_3D_utilities.hpp"


namespace Kratos
{
	
namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m) 
{
    using namespace pybind11;

    class_< FracturePropagation3DUtilities > 
    (m, "FracturePropagation3DUtilities")
    .def(init<>())
    .def("CheckFracturePropagation",&FracturePropagation3DUtilities::CheckFracturePropagation)
    .def("MappingModelParts",&FracturePropagation3DUtilities::MappingModelParts);
    
    class_< FracturePropagation2DUtilities >
    (m, "FracturePropagation2DUtilities")
    .def(init<>())
    .def("CheckFracturePropagation",&FracturePropagation2DUtilities::CheckFracturePropagation)
    .def("MappingModelParts",&FracturePropagation2DUtilities::MappingModelParts);
}

}  // namespace Python.
} // Namespace Kratos
