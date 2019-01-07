//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes

// #include "includes/gid_io.h"
#include "custom_utilities/custom_gid_io.h"
#include "custom_utilities/custom_gid_mesh_container.h"
#include "python/add_io_to_python.h"
#include "custom_python/add_custom_io_to_python.h"

namespace Kratos
{
namespace Python
{
typedef GidIO<TrilinosGidGaussPointsContainer,TrilinosGidMeshContainer> GidIOType;
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
void DoublePrintOnGaussPointsIndexed( GidIOType& dummy, const Variable<double>& rVariable,
                                      ModelPart& r_model_part, double SolutionTag, int value_index = 0 )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag, value_index );
}

void Array1DPrintOnGaussPointsIndexed( GidIOType& dummy, const Variable<array_1d<double,3> >rVariable, ModelPart& r_model_part, double SolutionTag, unsigned int value_index = 0 )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag, value_index );
}

void VectorPrintOnGaussPointsIndexed( GidIOType& dummy, const Variable<Vector>& rVariable,
                                      ModelPart& r_model_part, double SolutionTag, int value_index = 0 )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag, value_index );
}

void MatrixPrintOnGaussPointsIndexed( GidIOType& dummy, const Variable<Matrix>& rVariable,
                                      ModelPart& r_model_part, double SolutionTag, int value_index = 0 )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag, value_index );
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

void  AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<GidIOType, GidIOType::Pointer, IO>(m,"TrilinosGidIO")
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
    .def("PrintOnGaussPoints", DoublePrintOnGaussPointsIndexed)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPointsIndexed)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPointsIndexed)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPointsIndexed)
//                     .def("CondPrintOnGaussPoints", pointer_to_double_cond_print_on_gauss_points)

    .def("Flush",&GidIOType::Flush)
    .def("CloseResultFile",&GidIOType::CloseResultFile)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

}
}  // namespace Python.

} // Namespace Kratos

