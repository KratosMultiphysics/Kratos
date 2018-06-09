// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"


//Utilities
#include "custom_utilities/formfinding_io_utility.h"


namespace Kratos
{
namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<FormfindingIOUtility>(m,"FormfindingIOUtility")
    .def(init<ModelPart&, const Parameters>())
    .def("PrintModelPart",&FormfindingIOUtility::PrintModelPart)
    .def("ReadPrestressData",&FormfindingIOUtility::ReadPrestressData )
    .def("PrintPrestressData",&FormfindingIOUtility::PrintPrestressData )
    ;

}

}  // namespace Python.  

} // Namespace Kratos

