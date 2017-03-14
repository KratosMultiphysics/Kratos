// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_mappers_to_python.h"
#include "custom_utilities/interface_preprocess.h"

namespace Kratos
{

namespace Python
{
void  AddCustomMappersToPython()
{

    using namespace boost::python;

    class_<InterfacePreprocessCondition>("InterfacePreprocessCondition", init<>())
    .def("GenerateInterfacePart",&InterfacePreprocessCondition::GenerateInterfacePart)
    .def("GenerateLine2NInterfacePart",&InterfacePreprocessCondition::GenerateLine2NInterfacePart)
    .def("GenerateLine3NInterfacePart",&InterfacePreprocessCondition::GenerateLine3NInterfacePart)
    .def("GenerateTriangle3NInterfacePart",&InterfacePreprocessCondition::GenerateTriangle3NInterfacePart)
    .def("GenerateTriangle6NInterfacePart",&InterfacePreprocessCondition::GenerateTriangle6NInterfacePart)
    .def("GenerateQuadrilateral4NInterfacePart",&InterfacePreprocessCondition::GenerateQuadrilateral4NInterfacePart)
    .def("GenerateQuadrilateral8NInterfacePart",&InterfacePreprocessCondition::GenerateQuadrilateral8NInterfacePart)
    .def("GenerateQuadrilateral9NInterfacePart",&InterfacePreprocessCondition::GenerateQuadrilateral9NInterfacePart)
    ;
}

}  // namespace Python.

} // Namespace Kratos

