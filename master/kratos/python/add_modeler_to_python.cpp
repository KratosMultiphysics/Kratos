//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "add_modeler_to_python.h"
#include "modeler/modeler_factory.h"
#include "modeler/edge_swapping_2d_modeler.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "modeler/create_entities_from_geometries_modeler.h"
#include "modeler/serial_model_part_combinator_modeler.h"
#include "modeler/duplicate_mesh_modeler.h"
#include "modeler/copy_properties_modeler.h"
#include "modeler/combine_model_part_modeler.h"
#include "modeler/voxel_mesh_generator_modeler.h"
#include "modeler/clean_up_problematic_triangles_modeler.h"

namespace Kratos::Python
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
    } else if (KratosComponents<Condition>::Has(rName)) {
        GM.GenerateModelPart(origin_model_part, destination_model_part,
                             KratosComponents<Condition>::Get(rName));
    } else {
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
    .def("Create", &Modeler::Create)
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

    py::class_< SerialModelPartCombinatorModeler, SerialModelPartCombinatorModeler::Pointer, Modeler >(m,"SerialModelPartCombinatorModeler")
        .def(py::init< >())
        .def(py::init<Model&, Parameters>())
    ;

    py::class_< DuplicateMeshModeler, DuplicateMeshModeler::Pointer, Modeler >(m,"DuplicateMeshModeler")
        .def(py::init<ModelPart&>())
    ;

    py::class_< CopyPropertiesModeler, CopyPropertiesModeler::Pointer, Modeler >(m,"CopyPropertiesModeler")
        .def(py::init<Model&, Parameters>())
        .def(py::init<ModelPart&, ModelPart&>())
    ;

    py::class_< CombineModelPartModeler, CombineModelPartModeler::Pointer, Modeler >(m,"CombineModelPartModeler")
        .def(py::init<Model&, Parameters>())
    ;

    py::class_< CreateEntitiesFromGeometriesModeler, CreateEntitiesFromGeometriesModeler::Pointer, Modeler >(m, "CreateEntitiesFromGeometriesModeler")
        .def(py::init<Model&, Parameters>())
    ;

    py::class_<VoxelMeshGeneratorModeler, VoxelMeshGeneratorModeler::Pointer, Modeler>(m, "VoxelMeshGeneratorModeler")
        .def(py::init<Model &, Parameters>())
    ;

    py::class_<CleanUpProblematicTrianglesModeler, CleanUpProblematicTrianglesModeler::Pointer, Modeler>(m, "CleanUpProblematicTrianglesModeler")
        .def(py::init<Model&, Parameters>())
    ;
}

}  // namespace Kratos::Python.