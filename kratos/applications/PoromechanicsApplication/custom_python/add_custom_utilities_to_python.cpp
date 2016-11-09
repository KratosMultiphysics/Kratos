//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/condition_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/interface_element_utilities.hpp"
#include "custom_utilities/fracture_propagation_2D_utilities.hpp"
#include "custom_utilities/nonlocal_damage_2D_utilities.hpp"


namespace Kratos
{
	
namespace Python
{

void  AddCustomUtilitiesToPython() 
{
    using namespace boost::python;
    
    class_< FracturePropagation2DUtilities > ("FracturePropagation2DUtilities", init<>())
    .def("CheckFracturePropagation",&FracturePropagation2DUtilities::CheckFracturePropagation)
    .def("MappingModelParts",&FracturePropagation2DUtilities::MappingModelParts)
    ;
}

}  // namespace Python.
} // Namespace Kratos
