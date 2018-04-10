//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
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
