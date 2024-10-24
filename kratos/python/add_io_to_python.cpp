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

// #define JSON_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/io.h"

#include "includes/model_part_io.h"
#include "includes/reorder_consecutive_model_part_io.h"
#include "input_output/stl_io.h"
#include "input_output/obj_io.h"
#include "includes/gid_io.h"
#include "python/add_io_to_python.h"
#include "containers/flags.h"

// Outputs
#include "input_output/vtk_output.h"
#include "input_output/unv_output.h"
#include "input_output/cad_json_input.h"
#include "input_output/vtu_output.h"

#ifdef JSON_INCLUDED
#include "includes/json_io.h"
#endif

namespace Kratos::Python
{

void (GidIO<>::*pointer_to_io_read_single_properties)(Properties& rThisProperties ) = &IO::ReadProperties;

void (GidIO<>::*pointer_to_io_read_properties)(IO::PropertiesContainerType& rThisProperties ) = &IO::ReadProperties;

void ReadInitialValues1(IO& IO, IO::NodesContainerType& rThisNodes, IO::ElementsContainerType& rThisElements, IO::ConditionsContainerType& rThisConditions){ IO.ReadInitialValues(rThisNodes, rThisElements, rThisConditions);}
void ReadInitialValues2(IO& IO, ModelPart& rThisModelPart){ IO.ReadInitialValues(rThisModelPart);}

void  AddIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<IO, IO::Pointer> io_python_interface = py::class_<IO, IO::Pointer>(m,"IO")
    .def("ReadNode",&IO::ReadNode)
    .def("ReadNodes",&IO::ReadNodes)
    .def("WriteNodes",&IO::WriteNodes)
    .def("ReadProperties",[](IO& self, Properties& rThisProperties){self.ReadProperties(rThisProperties);})
    .def("ReadProperties",[](IO& self,IO::PropertiesContainerType& rThisProperties){self.ReadProperties(rThisProperties);})
//     .def("ReadProperties",pointer_to_io_read_properties)
    .def("ReadGeometries",&IO::ReadGeometries)
    .def("WriteGeometries",&IO::WriteGeometries)
    .def("ReadElements",&IO::ReadElements)
    .def("WriteElements",&IO::WriteElements)
    .def("ReadConditions",&IO::ReadConditions)
    .def("WriteConditions",&IO::WriteConditions)
    .def("ReadInitialValues",&ReadInitialValues1)
    .def("ReadInitialValues",&ReadInitialValues2)
    .def("ReadMesh",&IO::ReadMesh)
    .def("ReadModelPart",&IO::ReadModelPart)
    .def("WriteModelPart", py::overload_cast<const ModelPart&>(&IO::WriteModelPart)) // overload_cast can be removed once the legacy version is removed
    ;

    io_python_interface.attr("READ") = IO::READ;
    io_python_interface.attr("WRITE") =IO::WRITE;
    io_python_interface.attr("APPEND") = IO::APPEND;
    io_python_interface.attr("IGNORE_VARIABLES_ERROR" ) = IO::IGNORE_VARIABLES_ERROR;
    io_python_interface.attr("SKIP_TIMER" ) = IO::SKIP_TIMER;
    io_python_interface.attr("MESH_ONLY" ) = IO::MESH_ONLY;
    io_python_interface.attr("SCIENTIFIC_PRECISION" ) = IO::SCIENTIFIC_PRECISION;

    py::class_<ModelPartIO, ModelPartIO::Pointer, IO>(
       m, "ModelPartIO")
        .def(py::init<const std::filesystem::path&>())
        .def(py::init<const std::filesystem::path&, const Flags>())
    ;


    py::class_<ReorderConsecutiveModelPartIO, ReorderConsecutiveModelPartIO::Pointer, ModelPartIO>(m,"ReorderConsecutiveModelPartIO")
        .def(py::init<const std::filesystem::path&>())
        .def(py::init<const std::filesystem::path&, const Flags>())
    ;

#ifdef JSON_INCLUDED
    py::class_<KratosJsonIO, KratosJsonIO::Pointer, IO>(m,
         "JsonIO",init<std::string const&>())
        .def(py::init<std::string const&, const Flags>())
    ;
#endif

    py::class_<GidIO<>, GidIO<>::Pointer, IO>(m, "GidIO")
    .def(py::init<std::string const&, const GiD_PostMode, const MultiFileFlag, const WriteDeformedMeshFlag, const WriteConditionsFlag>())
    .def(py::init<std::string const&, const GiD_PostMode, const MultiFileFlag, const WriteDeformedMeshFlag, const WriteConditionsFlag, const bool>())
    //.def(py::init<std::string const&>())
    .def("WriteMesh",[](GidIO<>& dummy, GidIO<>::MeshType& rThisMesh){dummy.WriteMesh( rThisMesh );})
    .def("WriteNodeMesh",[](GidIO<>& dummy, GidIO<>::MeshType& rThisMesh){dummy.WriteNodeMesh( rThisMesh );})
    .def("WriteSphereMesh",[](GidIO<>& dummy, const GidIO<>::MeshType& rThisMesh){dummy.WriteSphereMesh( rThisMesh );})
    .def("WriteCircleMesh",[](GidIO<>& dummy, const GidIO<>::MeshType& rThisMesh){dummy.WriteCircleMesh( rThisMesh );})
    .def("WriteClusterMesh",[](GidIO<>& dummy, const GidIO<>::MeshType& rThisMesh){dummy.WriteClusterMesh( rThisMesh );})
    .def("InitializeMesh",&GidIO<>::InitializeMesh)
    .def("FinalizeMesh",&GidIO<>::FinalizeMesh)
    .def("InitializeResults",&GidIO<>::InitializeResults)
    .def("FinalizeResults",&GidIO<>::FinalizeResults)
    .def("WriteNodalResults",[](GidIO<>& dummy, Variable<bool> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteNodalResults(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    .def("WriteNodalResults",[](GidIO<>& dummy, Variable<double> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteNodalResults(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    .def("WriteNodalResults",[](GidIO<>& dummy, Variable<int> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteNodalResults(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    .def("WriteNodalResults",[](GidIO<>& dummy, Variable<array_1d<double, 3>> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteNodalResults(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    .def("WriteNodalResults",[](GidIO<>& dummy, Variable<Vector> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteNodalResults(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    .def("WriteNodalResults",[](GidIO<>& dummy, Variable<Matrix> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteNodalResults(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    .def("WriteLocalAxesOnNodes",[](GidIO<>& dummy, Variable<array_1d<double, 3>> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag, const std::size_t SolutionStepNumber){dummy.WriteLocalAxesOnNodes(rVariable, rNodes, SolutionTag, SolutionStepNumber);})
    // NonHistorical
    .def("WriteNodalFlags",[](GidIO<>& dummy, const Kratos::Flags& rFlag, const std::string& rFlagName, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalFlags(rFlag, rFlagName, rNodes, SolutionTag);})
    .def("WriteNodalResultsNonHistorical",[](GidIO<>& dummy, Variable<bool> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalResultsNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("WriteNodalResultsNonHistorical",[](GidIO<>& dummy, Variable<double> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalResultsNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("WriteNodalResultsNonHistorical",[](GidIO<>& dummy, Variable<int> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalResultsNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("WriteNodalResultsNonHistorical",[](GidIO<>& dummy, Variable<array_1d<double, 3>> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalResultsNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("WriteNodalResultsNonHistorical",[](GidIO<>& dummy, Variable<Vector> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalResultsNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("WriteNodalResultsNonHistorical",[](GidIO<>& dummy, Variable<Matrix> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteNodalResultsNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("WriteLocalAxesOnNodesNonHistorical",[](GidIO<>& dummy, Variable<array_1d<double, 3>> const& rVariable, const GidIO<>::NodesContainerType& rNodes, const double SolutionTag){dummy.WriteLocalAxesOnNodesNonHistorical(rVariable, rNodes, SolutionTag);})
    .def("PrintFlagsOnGaussPoints", [](GidIO<>& dummy, const Kratos::Flags& rFlag, const std::string& rFlagName, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintFlagsOnGaussPoints( rFlag, rFlagName, rModelPart, SolutionTag);})
    .def("PrintOnGaussPoints", [](GidIO<>& dummy, const Variable<bool>& rVariable, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );})
    .def("PrintOnGaussPoints", [](GidIO<>& dummy, const Variable<int>& rVariable, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );})
    .def("PrintOnGaussPoints", [](GidIO<>& dummy, const Variable<double>& rVariable, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );})
    .def("PrintOnGaussPoints", [](GidIO<>& dummy, const Variable<array_1d<double,3>>& rVariable, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );})
    .def("PrintOnGaussPoints", [](GidIO<>& dummy, const Variable<Vector>& rVariable, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );})
    .def("PrintOnGaussPoints", [](GidIO<>& dummy, const Variable<Matrix>& rVariable, const ModelPart& rModelPart, const double SolutionTag){dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );})
    .def("Flush",&GidIO<>::Flush)
    .def("ChangeOutputName",&GidIO<>::ChangeOutputName)
    .def("CloseResultFile",&GidIO<>::CloseResultFile)
    ;

    py::enum_<GiD_PostMode>(m,"GiDPostMode")
    .value("GiD_PostAscii", GiD_PostAscii)
    .value("GiD_PostAsciiZipped", GiD_PostAsciiZipped)
    .value("GiD_PostBinary", GiD_PostBinary)
    .value("GiD_PostHDF5", GiD_PostHDF5)
    ;

    py::enum_<WriteDeformedMeshFlag>(m,"WriteDeformedMeshFlag")
    .value("WriteDeformed", WriteDeformed)
    .value("WriteUndeformed", WriteUndeformed)
    ;

    py::enum_<WriteConditionsFlag>(m,"WriteConditionsFlag")
    .value("WriteConditions",WriteConditions)
    .value("WriteElementsOnly",WriteElementsOnly)
    .value("WriteConditionsOnly",WriteConditionsOnly)
    ;

    py::enum_<MultiFileFlag>(m,"MultiFileFlag")
    .value("SingleFile",SingleFile)
    .value("MultipleFiles",MultipleFiles)
    ;


    py::class_<VtkOutput, VtkOutput::Pointer, IO>(m, "VtkOutput")
        .def(py::init< ModelPart&>())
        .def(py::init< ModelPart&, Parameters >())
        .def("PrintOutput", &VtkOutput::PrintOutput, py::arg("output_filename")="")
        .def_static("GetDefaultParameters", &VtkOutput::GetDefaultParameters)
        ;

    py::class_<UnvOutput, UnvOutput::Pointer>(m, "UnvOutput")
        .def(py::init<ModelPart&, const std::string &>())
        .def("InitializeMesh", &UnvOutput::InitializeOutputFile)
        .def("WriteMesh", &UnvOutput::WriteMesh)
        .def("PrintOutput", (void (UnvOutput::*)(const Variable<bool>&, const double)) &UnvOutput::WriteNodalResults)
        .def("PrintOutput", (void (UnvOutput::*)(const Variable<int>&, const double)) &UnvOutput::WriteNodalResults)
        .def("PrintOutput", (void (UnvOutput::*)(const Variable<double>&, const double)) &UnvOutput::WriteNodalResults)
        .def("PrintOutput", (void (UnvOutput::*)(const Variable<array_1d<double,3>>&, const double)) &UnvOutput::WriteNodalResults)
        ;


    py::class_<StlIO, StlIO::Pointer, IO>(m, "StlIO")
        .def(py::init<std::filesystem::path const& >())
        .def(py::init<std::filesystem::path const&, Parameters>())
        ;

    py::class_<ObjIO, ObjIO::Pointer, IO>(m, "ObjIO")
        .def(py::init<std::filesystem::path const& >())
        .def(py::init<std::filesystem::path const&, Parameters>())
        ;

    // Import of CAD models to the model part
    py::class_<CadJsonInput<>, CadJsonInput<>::Pointer>(m, "CadJsonInput")
        .def(py::init<const std::string &>())
        .def(py::init<const std::string&, std::size_t>())
        .def(py::init<Parameters>())
        .def(py::init<Parameters, std::size_t>())
        .def("ReadModelPart", &CadJsonInput<>::ReadModelPart)
        ;

    auto vtu_output = py::class_<VtuOutput, VtuOutput::Pointer>(m, "VtuOutput");

    py::enum_<VtuOutput::WriterFormat>(vtu_output, "WriterFormat")
        .value("ASCII", VtuOutput::WriterFormat::ASCII)
        .value("BINARY", VtuOutput::WriterFormat::BINARY)
        .export_values();

    vtu_output
        .def(py::init<ModelPart&, const bool, const VtuOutput::WriterFormat, const std::size_t>(), py::arg("model_part"), py::arg("is_initial_configuration_considered") = true, py::arg("binary_output") = VtuOutput::WriterFormat::BINARY, py::arg("precision") = 9)
        .def("AddHistoricalVariable", &VtuOutput::AddHistoricalVariable<int>, py::arg("int_variable"))
        .def("AddHistoricalVariable", &VtuOutput::AddHistoricalVariable<double>, py::arg("double_variable"))
        .def("AddHistoricalVariable", &VtuOutput::AddHistoricalVariable<array_1d<double, 3>>, py::arg("array3_variable"))
        .def("AddHistoricalVariable", &VtuOutput::AddHistoricalVariable<array_1d<double, 4>>, py::arg("array4_variable"))
        .def("AddHistoricalVariable", &VtuOutput::AddHistoricalVariable<array_1d<double, 6>>, py::arg("array6_variable"))
        .def("AddHistoricalVariable", &VtuOutput::AddHistoricalVariable<array_1d<double, 9>>, py::arg("array9_variable"))
        .def("AddNonHistoricalVariable", &VtuOutput::AddNonHistoricalVariable<int>, py::arg("int_variable"), py::arg("entity_flags"))
        .def("AddNonHistoricalVariable", &VtuOutput::AddNonHistoricalVariable<double>, py::arg("double_variable"), py::arg("entity_flags"))
        .def("AddNonHistoricalVariable", &VtuOutput::AddNonHistoricalVariable<array_1d<double, 3>>, py::arg("array3_variable"), py::arg("entity_flags"))
        .def("AddNonHistoricalVariable", &VtuOutput::AddNonHistoricalVariable<array_1d<double, 4>>, py::arg("array4_variable"), py::arg("entity_flags"))
        .def("AddNonHistoricalVariable", &VtuOutput::AddNonHistoricalVariable<array_1d<double, 6>>, py::arg("array6_variable"), py::arg("entity_flags"))
        .def("AddNonHistoricalVariable", &VtuOutput::AddNonHistoricalVariable<array_1d<double, 9>>, py::arg("array9_variable"), py::arg("entity_flags"))
        .def("AddFlagVariable", &VtuOutput::AddFlagVariable, py::arg("flag_variable_name"), py::arg("flag"), py::arg("entity_flags"))
        .def("AddContainerExpression", &VtuOutput::AddContainerExpression<ModelPart::NodesContainerType>, py::arg("container_expression_name"), py::arg("container_expression"))
        .def("AddContainerExpression", &VtuOutput::AddContainerExpression<ModelPart::ConditionsContainerType>, py::arg("container_expression_name"), py::arg("container_expression"))
        .def("AddContainerExpression", &VtuOutput::AddContainerExpression<ModelPart::ElementsContainerType>, py::arg("container_expression_name"), py::arg("container_expression"))
        .def("ClearHistoricalVariables", &VtuOutput::ClearHistoricalVariables)
        .def("ClearNodalNonHistoricalVariables", &VtuOutput::ClearNodalNonHistoricalVariables)
        .def("ClearNodalFlags", &VtuOutput::ClearNodalFlags)
        .def("ClearCellNonHistoricalVariables", &VtuOutput::ClearCellNonHistoricalVariables)
        .def("ClearCellFlags", &VtuOutput::ClearCellFlags)
        .def("ClearNodalContainerExpressions", &VtuOutput::ClearNodalContainerExpressions)
        .def("ClearCellContainerExpressions", &VtuOutput::ClearCellContainerExpressions)
        .def("GetModelPart", &VtuOutput::GetModelPart, py::return_value_policy::reference)
        .def("PrintOutput", &VtuOutput::PrintOutput, py::arg("output_file_name_prefix"))
        ;

    vtu_output.attr("NODES") = VtuOutput::NODES;
    vtu_output.attr("CONDITIONS") = VtuOutput::CONDITIONS;
    vtu_output.attr("ELEMENTS") = VtuOutput::ELEMENTS;
}
}  // namespace Kratos::Python.


