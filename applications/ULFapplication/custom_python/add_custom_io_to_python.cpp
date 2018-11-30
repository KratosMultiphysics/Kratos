/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: janosch $
//   Date:                $Date: 2008-04-29 15:05:55 $
//   Revision:            $Revision: 1.4 $
//
//


// System includes

// External includes
//#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
// #include "includes/datafile_io.h"
// #include "includes/gid_io.h"
#include "custom_utilities/pfem_gid_io.h"
#include "python/add_io_to_python.h"
#include "custom_python/add_custom_io_to_python.h"

namespace Kratos
{
namespace Python
{
typedef GidIO<PfemGidGaussPointsContainer,PfemGidMeshContainer> GidIOType;
typedef GidIO<> GidIOBaseType;

void WriteNodeMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
    KRATOS_WATCH("writing Node Mesh");
    dummy.WriteNodeMesh( rThisMesh );
}
void WriteMesh( GidIOType& dummy, GidIOType::MeshType& rThisMesh )
{
//            dummy.WriteMesh( rThisMesh );
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

void  AddCustomIOToPython()
{


    using namespace boost::python;

    py::class_<GidIOType, GidIOType::Pointer, bases<DatafileIO>, boost::noncopyable>(
        "PFEMGidIO",init<std::string const&, GiD_PostMode,
        MultiFileFlag,
        WriteDeformedMeshFlag,
        WriteConditionsFlag>())
    //.def(py::init<std::string const&>())
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
    .def("CloseResultFile",&GidIOType::CloseResultFile)
    //.def("",&DatafileIO::)
    //.def(self_ns::str(self))
    ;

}
}  // namespace Python.

} // Namespace Kratos

