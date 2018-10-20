//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
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
#include "includes/gid_io.h"
#include "python/add_io_to_python.h"
#include "containers/flags.h"

#ifdef JSON_INCLUDED
#include "includes/json_io.h"
#endif

namespace Kratos
{

namespace Python
{

void (GidIO<>::*pointer_to_io_read_single_properties)(Properties& rThisProperties ) = &IO::ReadProperties;

void (GidIO<>::*pointer_to_io_read_properties)(IO::PropertiesContainerType& rThisProperties ) = &IO::ReadProperties;


void WriteNodeMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    dummy.WriteNodeMesh( rThisMesh );
}

void WriteSphereMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    dummy.WriteSphereMesh( rThisMesh );
}

void WriteCircleMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    dummy.WriteCircleMesh( rThisMesh );
}

void WriteClusterMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    dummy.WriteClusterMesh( rThisMesh );
}

void WriteMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
{
    dummy.WriteMesh( rThisMesh );
}


void FlagsPrintOnGaussPoints( GidIO<>& dummy, Kratos::Flags rFlag, std::string rFlagName,
                               ModelPart& rModelPart, double SolutionTag )
{
    dummy.PrintFlagsOnGaussPoints( rFlag, rFlagName, rModelPart, SolutionTag );
}

void DoublePrintOnGaussPoints( GidIO<>& dummy, const Variable<double>& rVariable,
                               ModelPart& rModelPart, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );
}

void IntPrintOnGaussPoints( GidIO<>& dummy, const Variable<int>& rVariable,
                               ModelPart& rModelPart, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );
}

void Array1DPrintOnGaussPoints( GidIO<>& dummy, const Variable<array_1d<double,3> >& rVariable,
                                ModelPart& rModelPart, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );
}

void VectorPrintOnGaussPoints( GidIO<>& dummy, const Variable<Vector>& rVariable,
                               ModelPart& rModelPart, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );
}

void MatrixPrintOnGaussPoints( GidIO<>& dummy, const Variable<Matrix>& rVariable,
                               ModelPart& rModelPart, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, rModelPart, SolutionTag );
}


void (GidIO<>::*pointer_to_bool_write_nodal_results)( Variable<bool> const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_double_write_nodal_results)( Variable<double> const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_int_write_nodal_results)( Variable<int> const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_array1d_write_nodal_results)(
    Variable<array_1d<double, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes,
    double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_vector_write_nodal_results)(
    Variable<Vector > const& rVariable, GidIO<>::NodesContainerType& rNodes, 
    double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*pointer_to_matrix_write_nodal_results)(
    Variable<Matrix > const& rVariable, GidIO<>::NodesContainerType& rNodes, 
    double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;
void (GidIO<>::*local_axes_write_nodal_results)( Variable<array_1d<double, 3> > const& rVariable,
        GidIO<>::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GidIO<>::WriteLocalAxesOnNodes;

/////////////////////////////////////////////////////////////
/// NON-HISTORICAL DATABASE                               ///
////////////////////////////////////////////////////////////
void (GidIO<>::*pointer_to_flags_write_nodal_results_NH)( Kratos::Flags rFlag, std::string rFlagName,  GidIO<>::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIO<>::WriteNodalFlags;
        
void (GidIO<>::*pointer_to_bool_write_nodal_results_NH)( Variable<bool> const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIO<>::WriteNodalResultsNonHistorical;
void (GidIO<>::*pointer_to_double_write_nodal_results_NH)( Variable<double> const& rVariable,  GidIO<>::NodesContainerType& rNodes, double SolutionTag) 
        = &GidIO<>::WriteNodalResultsNonHistorical;
void (GidIO<>::*pointer_to_array1d_write_nodal_results_NH)(Variable<array_1d<double, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag) = &GidIO<>::WriteNodalResultsNonHistorical;
void (GidIO<>::*pointer_to_vector_write_nodal_results_NH)(Variable<Vector > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag) 
        = &GidIO<>::WriteNodalResultsNonHistorical;
void (GidIO<>::*pointer_to_matrix_write_nodal_results_NH)(Variable<Matrix > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIO<>::WriteNodalResultsNonHistorical;
void (GidIO<>::*local_axes_write_nodal_results_NH)( Variable<array_1d<double, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag) 
    = &GidIO<>::WriteLocalAxesOnNodesNonHistorical;

//         void (GidIO::*pointer_to_double_cond_print_on_gauss_points)(const Variable<double>& rVariable,
//               ModelPart& rModelPart, double SolutionTag) = &GidIO::CondPrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_double_print_on_gauss_points)(const Variable<double >& rVariable,
//               ModelPart& rModelPart, double SolutionTag)
//                 = &GidIO<>::PrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_matrix_print_on_gauss_points)(const Variable<Matrix >& rVariable,
//               ModelPart& rModelPart, double SolutionTag)
//                 = &GidIO<>::PrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_vector_print_on_gauss_points)(const Variable<Vector >& rVariable,
//               ModelPart& rModelPart, double SolutionTag)
//                 = &GidIO<>::PrintOnGaussPoints;

void ReadInitialValues1(IO& IO, IO::NodesContainerType& rThisNodes, IO::ElementsContainerType& rThisElements, IO::ConditionsContainerType& rThisConditions){ IO.ReadInitialValues(rThisNodes, rThisElements, rThisConditions);}
void ReadInitialValues2(IO& IO, ModelPart& rThisModelPart){ IO.ReadInitialValues(rThisModelPart);}


        
        
void  AddIOToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<IO, IO::Pointer> io_python_interface = class_<IO, IO::Pointer>(m,"IO")
    .def("ReadNode",&IO::ReadNode)
    .def("ReadNodes",&IO::ReadNodes)
    .def("WriteNodes",&IO::WriteNodes)
    .def("ReadProperties",[](IO& self, Properties& rThisProperties){self.ReadProperties(rThisProperties);})
    .def("ReadProperties",[](IO& self,IO::PropertiesContainerType& rThisProperties){self.ReadProperties(rThisProperties);})
//     .def("ReadProperties",pointer_to_io_read_properties)
    .def("ReadElements",&IO::ReadElements)
    .def("WriteElements",&IO::WriteElements)
    .def("ReadConditions",&IO::ReadConditions)
    .def("WriteConditions",&IO::WriteConditions)
    .def("ReadInitialValues",&ReadInitialValues1)
    .def("ReadInitialValues",&ReadInitialValues2)
    .def("ReadMesh",&IO::ReadMesh)
    .def("ReadModelPart",&IO::ReadModelPart)
    .def("WriteModelPart",&IO::WriteModelPart)
    ;
    io_python_interface.attr("READ") = IO::READ;
    io_python_interface.attr("WRITE") =IO::WRITE;
    io_python_interface.attr("APPEND") = IO::APPEND;
    io_python_interface.attr("IGNORE_VARIABLES_ERROR" ) = IO::IGNORE_VARIABLES_ERROR;
    io_python_interface.attr("SKIP_TIMER" ) = IO::SKIP_TIMER;
    io_python_interface.attr("MESH_ONLY" ) = IO::MESH_ONLY;
 

    
    class_<ModelPartIO, ModelPartIO::Pointer, IO>(
       m, "ModelPartIO")
        .def(init<std::string const&>())
        .def(init<std::string const&, const Flags>())
    ;

    
    class_<ReorderConsecutiveModelPartIO, ReorderConsecutiveModelPartIO::Pointer, ModelPartIO>(m,"ReorderConsecutiveModelPartIO")
        .def(init<std::string const&>())
        .def(init<std::string const&, const Flags>())
    ;
#ifdef JSON_INCLUDED
    class_<KratosJsonIO, KratosJsonIO::Pointer, IO>(m,
         "JsonIO",init<std::string const&>())
        .def(init<std::string const&, const Flags>())
    ;
#endif
    
    /*
    class_<ModelPartIO, ModelPartIO::Pointer, bases<IO>,  boost::noncopyable>(
        "ModelPartIOWithSkipReadFlag",init<std::string const&, bool const>())
    ;*/

    class_<GidIO<>, GidIO<>::Pointer, IO>(m,
        "GidIO")
    .def(init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>())
    //.def(init<std::string const&>())
    .def("WriteMesh",WriteMesh)
    .def("WriteNodeMesh",WriteNodeMesh)
    .def("WriteSphereMesh",WriteSphereMesh)
    .def("WriteCircleMesh",WriteCircleMesh)
    .def("WriteClusterMesh",WriteClusterMesh)

    .def("InitializeMesh",&GidIO<>::InitializeMesh)
    .def("FinalizeMesh",&GidIO<>::FinalizeMesh)

    .def("InitializeResults",&GidIO<>::InitializeResults)
    .def("FinalizeResults",&GidIO<>::FinalizeResults)

    .def("WriteNodalResults",pointer_to_bool_write_nodal_results)
    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
    .def("WriteNodalResults",pointer_to_int_write_nodal_results)
    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)
    .def("WriteNodalResults",pointer_to_vector_write_nodal_results)
    .def("WriteNodalResults",pointer_to_matrix_write_nodal_results)

    .def("WriteLocalAxesOnNodes",local_axes_write_nodal_results)
    // NonHistorical
    .def("WriteNodalFlags",pointer_to_flags_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_bool_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_double_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_array1d_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_vector_write_nodal_results_NH)
    .def("WriteNodalResultsNonHistorical",pointer_to_matrix_write_nodal_results_NH)
    .def("WriteLocalAxesOnNodesNonHistorical",local_axes_write_nodal_results_NH)

//                     .def("PrintOnGaussPoints", pointer_to_double_print_on_gauss_points)
    .def("PrintFlagsOnGaussPoints", FlagsPrintOnGaussPoints)
    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", IntPrintOnGaussPoints)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)


//                     .def("PrintOnGaussPoints", pointer_to_vector_print_on_gauss_points)
//                     .def("PrintOnGaussPoints", pointer_to_matrix_print_on_gauss_points)
//                     .def("CondPrintOnGaussPoints", pointer_to_double_cond_print_on_gauss_points)

    .def("Flush",&GidIO<>::Flush)
    .def("ChangeOutputName",&GidIO<>::ChangeOutputName)
    .def("CloseResultFile",&GidIO<>::CloseResultFile)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

    enum_<GiD_PostMode>(m,"GiDPostMode")
    .value("GiD_PostAscii", GiD_PostAscii)
    .value("GiD_PostAsciiZipped", GiD_PostAsciiZipped)
    .value("GiD_PostBinary", GiD_PostBinary)
    .value("GiD_PostHDF5", GiD_PostHDF5)
    ;

    enum_<WriteDeformedMeshFlag>(m,"WriteDeformedMeshFlag")
    .value("WriteDeformed", WriteDeformed)
    .value("WriteUndeformed", WriteUndeformed)
    ;

    enum_<WriteConditionsFlag>(m,"WriteConditionsFlag")
    .value("WriteConditions",WriteConditions)
    .value("WriteElementsOnly",WriteElementsOnly)
    .value("WriteConditionsOnly",WriteConditionsOnly)
    ;

    enum_<MultiFileFlag>(m,"MultiFileFlag")
    .value("SingleFile",SingleFile)
    .value("MultipleFiles",MultipleFiles)
    ;
}
}  // namespace Python.

} // Namespace Kratos

    
