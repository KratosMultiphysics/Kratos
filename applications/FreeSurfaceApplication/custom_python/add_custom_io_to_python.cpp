//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_utilities/edgebased_gid_io.h"
#include "python/add_io_to_python.h"
#include "custom_python/add_custom_io_to_python.h"

namespace Kratos
{
namespace Python
{
typedef GidIO<EdgebasedGidGaussPointsContainer,EdgebasedGidMeshContainer> GidIOType;
typedef GidIO<> GidIOBaseType;

void WriteNodeMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
    KRATOS_WATCH("writing Node Mesh");
    dummy.WriteNodeMesh( rThisMesh );
}
void WriteMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
    dummy.WriteMesh( rThisMesh );
}

void DoublePrintOnGaussPoints( GidIOType& dummy, const Variable<double>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void VectorPrintOnGaussPoints( GidIOType& dummy, const Variable<Vector>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void MatrixPrintOnGaussPoints( GidIOType& dummy, const Variable<Matrix>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void (GidIOType::*pointer_to_double_write_nodal_results)(
    Variable<double> const& rVariable,
    GidIOType::NodesContainerType& rNodes, double SolutionTag,
    std::size_t SolutionStepNumber ) = &GidIOType::WriteNodalResults;
void (GidIOType::*pointer_to_array1d_write_nodal_results)(
    Variable<array_1d<double, 3> > const& rVariable,
    GidIOType::NodesContainerType& rNodes, double SolutionTag,
    std::size_t SolutionStepNumber) = &GidIOType::WriteNodalResults;


//         void (GidIOType::*pointer_to_double_write_nodal_results)( Variable<double> const& rVariable,
//               GidIOType::NodesContainerType& rNodes, double SolutionTag,
//               std::size_t SolutionStepNumber ) = &GidIO::WriteNodalResults;
//         void (GidIOType::*pointer_to_array1d_write_nodal_results)(
//               Variable<array_1d<double, 3> > const& rVariable, GidIOType::NodesContainerType& rNodes,
//               double SolutionTag, std::size_t SolutionStepNumber) = &GidIOType::WriteNodalResults;

//         void (GidIOType::*pointer_to_double_print_on_gauss_points)(const Variable<double >& rVariable,
//               ModelPart& r_model_part, double SolutionTag)
//                 = &GidIOType::PrintOnGaussPoints;
//         void (GidIOType::*pointer_to_matrix_print_on_gauss_points)(const Variable<Matrix >& rVariable,
//               ModelPart& r_model_part, double SolutionTag)
//                 = &GidIOType::PrintOnGaussPoints;
//         void (GidIOType::*pointer_to_vector_print_on_gauss_points)(const Variable<Vector >& rVariable,
//               ModelPart& r_model_part, double SolutionTag)
//                 = &GidIOType::PrintOnGaussPoints;

void  AddCustomIOToPython(pybind11::module& pymodule)
{


    namespace py = pybind11;

    py::class_<GidIOType, GidIOType::Pointer>(pymodule, "EdgebasedGidIO")
        .def(py::init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>())
    //.def(py::init<std::string const&>())
    .def("WriteMesh",WriteMesh)
    .def("WriteNodeMesh",WriteNodeMesh)

    .def("InitializeMesh",&GidIOType::InitializeMesh)
    .def("FinalizeMesh",&GidIOType::FinalizeMesh)

    .def("InitializeResults",&GidIOType::InitializeResults)
    .def("FinalizeResults",&GidIOType::FinalizeResults)

    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)
//                     .def("WriteNodalResults",pointer_to_vector_write_nodal_results)
//                     .def("WriteNodalResults",pointer_to_matrix_write_nodal_results)

//                     .def("PrintOnGaussPoints", pointer_to_double_print_on_gauss_points)
//                     .def("PrintOnGaussPoints", pointer_to_vector_print_on_gauss_points)
//                     .def("PrintOnGaussPoints", pointer_to_matrix_print_on_gauss_points)
    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)
//                     .def("CondPrintOnGaussPoints", pointer_to_double_cond_print_on_gauss_points)

    .def("Flush",&GidIOType::Flush)
    .def("ChangeOutputName",&GidIOType::ChangeOutputName)
    .def("CloseResultFile",&GidIOType::CloseResultFile)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

}
}  // namespace Python.

} // Namespace Kratos

