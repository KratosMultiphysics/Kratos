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

#include "modeler/modeler_factory.h"

namespace Kratos
{

namespace Python
{

namespace py = pybind11;

void GenerateMesh(Modeler& GM, ModelPart& model_part, const std::string& rElementName, const std::string& rConditionName)
{
    GM.GenerateMesh(model_part,
                    KratosComponents<Element>::Get(rElementName),
                    KratosComponents<Condition>::Get(rConditionName));

}

void GeneratePartialModelPart(ConnectivityPreserveModeler& GM, ModelPart& origin_model_part, ModelPart& destination_model_part, const std::string& rName)
{
    if (KratosComponents<Element>::Has(rName)) {
        GM.GenerateModelPart(origin_model_part, destination_model_part,
                             KratosComponents<Element>::Get(rName));
    }
    else if (KratosComponents<Condition>::Has(rName)) {
        GM.GenerateModelPart(origin_model_part, destination_model_part,
                             KratosComponents<Condition>::Get(rName));
    }
    else {
        KRATOS_ERROR << "Unknown Element/Condition name " << rName << "." << std::endl;
    }
}

void  AddModelerToPython(pybind11::module& m)
{
    m.def("CreateModeler", &ModelerFactory::Create);
    m.def("HasModeler", &ModelerFactory::Has);

    py::class_<Modeler, Modeler::Pointer>(m,"Modeler")
    .def(py::init<>())
    .def(py::init<Model&, Parameters>())
    // Modeler Stages Initialize
    .def("SetupGeometryModel", &Modeler::SetupGeometryModel)
    .def("PrepareGeometryModel", &Modeler::PrepareGeometryModel)
    .def("SetupModelPart", &Modeler::SetupModelPart)
    // Additional Old Functions
    .def("GenerateModelPart",
        [] (Modeler& rModeler, ModelPart& origin_model_part, ModelPart& destination_model_part, const std::string& rElementName, const std::string& rConditionName)
        {rModeler.GenerateModelPart(origin_model_part, destination_model_part,
            KratosComponents<Element>::Get(rElementName), KratosComponents<Condition>::Get(rConditionName));})
    .def("GenerateMesh",&GenerateMesh)
    .def("GenerateNodes",&Modeler::GenerateNodes)
    .def("__str__", PrintObject<Modeler>)
    ;

    py::class_<ConnectivityPreserveModeler,ConnectivityPreserveModeler::Pointer,Modeler>(m,"ConnectivityPreserveModeler")
    .def(py::init< >())
    .def("GenerateModelPart",
        [] (Modeler& rModeler, ModelPart& origin_model_part, ModelPart& destination_model_part, const std::string& rElementName, const std::string& rConditionName)
        {rModeler.GenerateModelPart(origin_model_part, destination_model_part,
            KratosComponents<Element>::Get(rElementName), KratosComponents<Condition>::Get(rConditionName));})
    .def("GenerateModelPart",&GeneratePartialModelPart)
    ;

    py::class_< EdgeSwapping2DModeler, EdgeSwapping2DModeler::Pointer, Modeler >(m,"EdgeSwapping2DModeler")
            .def(py::init< >())
            .def("ReGenerateMesh",&EdgeSwapping2DModeler::Remesh)
    ;
}

}  // namespace Python.

} // Namespace Kratos
