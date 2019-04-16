// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors: Pooyan Dadvand
//                Josep Maria Carbonell Puigbo
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "custom_io/pfem_gid_io.h"
#include "python/add_io_to_python.h"
#include "custom_python/add_custom_io_to_python.h"
#include "includes/datafile_io.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;
    
typedef GidIO<PfemGidGaussPointsContainer,PfemGidMeshContainer> GidIOType;
typedef GidIO<> GidIOBaseType;

void WriteNodeMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
    KRATOS_INFO("WriteNodeMesh") << "Writing Node Mesh" << std::endl;
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

void  AddCustomIOToPython(pybind11::module& m)
{
    py::class_<GidIOType, GidIOType::Pointer,DatafileIO>(m,
        "PFEMGidIO",py::init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>());
    .def(py::init<std::string const&,
         std::string const&,
         std::string const&,
         std::string const&,
         std::string const&,
         std::string const&,
         GiD_PostMode,
         MultiFileFlag,
         WriteDeformedMeshFlag,
         WriteConditionsFlag>())
    .def("WriteMesh",WriteMesh)
    .def("WriteNodeMesh",WriteNodeMesh)
    .def("InitializeMesh",&GidIOType::InitializeMesh)
    .def("FinalizeMesh",&GidIOType::FinalizeMesh)
    .def("InitializeResults",&GidIOType::InitializeResults)
    .def("FinalizeResults",&GidIOType::FinalizeResults)
    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)
    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)
    .def("Flush",&GidIOType::Flush)
    .def("CloseResultFile",&GidIOType::CloseResultFile)
    ;

}
}  // namespace Python.

} // Namespace Kratos

