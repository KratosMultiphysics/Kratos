//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "add_modeler_to_python.h"
#include "modeler/modeler.h"
#include "modeler/edge_swapping_2d_modeler.h"
#include "modeler/connectivity_preserve_modeler.h"


namespace Kratos
{

namespace Python
{

using namespace pybind11;

void GenerateModelPart(Modeler& GM, ModelPart& origin_model_part, ModelPart& destination_model_part, const char* ElementName, const char* ConditionName)
{
    GM.GenerateModelPart(origin_model_part, destination_model_part,
                         KratosComponents<Element>::Get(ElementName),
                         KratosComponents<Condition>::Get(ConditionName));

}

void GenerateMesh(Modeler& GM, ModelPart& model_part, const char* ElementName, const char* ConditionName)
{
    GM.GenerateMesh(model_part,
                    KratosComponents<Element>::Get(ElementName),
                    KratosComponents<Condition>::Get(ConditionName));

}


void  AddModelerToPython(pybind11::module& m)
{
    class_<Modeler, Modeler::Pointer>(m,"Modeler")
    .def(init<>())
    .def("GenerateModelPart",&GenerateModelPart)
    .def("GenerateMesh",&GenerateMesh)
    .def("GenerateNodes",&Modeler::GenerateNodes)
    .def("__str__", PrintObject<Modeler>)
    ;

    class_<ConnectivityPreserveModeler,ConnectivityPreserveModeler::Pointer,Modeler>(m,"ConnectivityPreserveModeler")
    .def(init< >())
            ;


    class_< EdgeSwapping2DModeler, EdgeSwapping2DModeler::Pointer, Modeler >(m,"EdgeSwapping2DModeler")
            .def(init< >())
            .def("ReGenerateMesh",&EdgeSwapping2DModeler::Remesh)
    ;
}

}  // namespace Python.

} // Namespace Kratos
