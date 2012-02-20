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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-06-20 17:38:25 $
//   Revision:            $Revision: 1.24 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/datafile_io.h"
#include "includes/model_part_io.h"
#include "includes/gid_io.h"
#include "python/add_io_to_python.h"

namespace Kratos
{
	
    namespace Python
    {

      void (GidIO<>::*pointer_to_io_read_single_properties)(Properties& rThisProperties ) = &IO::ReadProperties;

      void (GidIO<>::*pointer_to_io_read_properties)(IO::PropertiesContainerType& rThisProperties ) = &IO::ReadProperties;

        std::string PrintDatafileIO(DatafileIO const& rDataFileIO)
        {
            std::stringstream buffer;
            rDataFileIO.PrintInfo(buffer);
            buffer << std::endl;
            rDataFileIO.PrintData(buffer);

            return buffer.str();
        }
        
        void WriteNodeMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
        {
            KRATOS_WATCH("writing Node Mesh");
            dummy.WriteNodeMesh( rThisMesh );
        }

        void WriteSphereMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
        {
            KRATOS_WATCH("writing Sphere Mesh");
            dummy.WriteSphereMesh( rThisMesh );
        }
        
        
        void WriteMesh( GidIO<>& dummy, GidIO<>::MeshType& rThisMesh )
        {
            dummy.WriteMesh( rThisMesh );
        }
        
        
        void DoublePrintOnGaussPoints( GidIO<>& dummy, const Variable<double>& rVariable,
                                       ModelPart& r_model_part, double SolutionTag )
        {
            dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
        }

        void Array1DPrintOnGaussPoints( GidIO<>& dummy, const Variable<array_1d<double,3> >& rVariable,
                                       ModelPart& r_model_part, double SolutionTag )
        {
            dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
        }

        void VectorPrintOnGaussPoints( GidIO<>& dummy, const Variable<Vector>& rVariable,
                                       ModelPart& r_model_part, double SolutionTag )
        {
            dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
        }
        
        void MatrixPrintOnGaussPoints( GidIO<>& dummy, const Variable<Matrix>& rVariable,
                                       ModelPart& r_model_part, double SolutionTag )
        {
            dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
        }
        

         void (GidIO<>::*pointer_to_bool_write_nodal_results)( Variable<bool> const& rVariable,
              GidIO<>::NodesContainerType& rNodes, double SolutionTag, 
              std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults; 
        void (GidIO<>::*pointer_to_double_write_nodal_results)( Variable<double> const& rVariable,
              GidIO<>::NodesContainerType& rNodes, double SolutionTag, 
              std::size_t SolutionStepNumber ) = &GidIO<>::WriteNodalResults;
        void (GidIO<>::*pointer_to_array1d_write_nodal_results)( 
              Variable<array_1d<double, 3> > const& rVariable, GidIO<>::NodesContainerType& rNodes, 
              double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;
        
//         void (GidIO::*pointer_to_vector_write_nodal_results)(Variable<Vector > const& rVariable, GidIO::NodesContainerType& rNodes, double SolutionTag, std::size_t SolutionStepNumber) = &GidIO::WriteNodalResults;
        void (GidIO<>::*pointer_to_matrix_write_nodal_results)(Variable<Matrix > const& rVariable, GidIO<>::NodesContainerType& rNodes, double SolutionTag, std::size_t SolutionStepNumber) = &GidIO<>::WriteNodalResults;

//         void (GidIO::*pointer_to_double_cond_print_on_gauss_points)(const Variable<double>& rVariable,
//               ModelPart& r_model_part, double SolutionTag) = &GidIO::CondPrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_double_print_on_gauss_points)(const Variable<double >& rVariable,
//               ModelPart& r_model_part, double SolutionTag) 
//                 = &GidIO<>::PrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_matrix_print_on_gauss_points)(const Variable<Matrix >& rVariable,
//               ModelPart& r_model_part, double SolutionTag) 
//                 = &GidIO<>::PrintOnGaussPoints;
//         void (GidIO<>::*pointer_to_vector_print_on_gauss_points)(const Variable<Vector >& rVariable,
//               ModelPart& r_model_part, double SolutionTag) 
//                 = &GidIO<>::PrintOnGaussPoints;

        void  AddIOToPython()
        {
            using namespace boost::python;
            
            class_<IO, IO::Pointer, boost::noncopyable>("IO")
                    .def("ReadNode",&IO::ReadNode)
                    .def("ReadNodes",&IO::ReadNodes)
                    .def("WriteNodes",&IO::WriteNodes)
                    .def("ReadProperties",pointer_to_io_read_single_properties)
                    .def("ReadProperties",pointer_to_io_read_properties)
                    .def("ReadElements",&IO::ReadElements)
                    .def("WriteElements",&IO::WriteElements)
                    .def("ReadConditions",&IO::ReadConditions)
                    .def("ReadInitialValues",&IO::ReadInitialValues)
                    .def("ReadMesh",&IO::ReadMesh)
                    .def("ReadModelPart",&IO::ReadModelPart)
                    .def("__str__",PrintDatafileIO)
	      ;

            class_<DatafileIO, DatafileIO::Pointer, bases<IO>,  boost::noncopyable>(
                    "DatafileIO",init<std::string const&>())
                    .def(init<std::string const&, 
                         std::string const&, 
                         std::string const&, 
                         std::string const&, 
                         std::string const&>())
                    //.def("",&DatafileIO::)
                    //.def(self_ns::str(self))
                    ;
            
            class_<ModelPartIO, ModelPartIO::Pointer, bases<IO>,  boost::noncopyable>(
                    "ModelPartIO",init<std::string const&>())
                    ;
            
            class_<GidIO<>, GidIO<>::Pointer, bases<DatafileIO>, boost::noncopyable>(
                    "GidIO",init<std::string const&, GiD_PostMode,
                    MultiFileFlag,
                    WriteDeformedMeshFlag,
                    WriteConditionsFlag>())
                    //.def(init<std::string const&>())
                    .def(init<std::string const&, 
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
                    .def("WriteSphereMesh",WriteSphereMesh)
                    
                    .def("InitializeMesh",&GidIO<>::InitializeMesh)
                    .def("FinalizeMesh",&GidIO<>::FinalizeMesh)
                    
                    .def("InitializeResults",&GidIO<>::InitializeResults)
                    .def("FinalizeResults",&GidIO<>::FinalizeResults)
                    
                     .def("WriteNodalResults",pointer_to_bool_write_nodal_results) 
                    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
                    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)

//                    .def("WriteNodalResults",pointer_to_vector_write_nodal_results)
                    .def("WriteNodalResults",pointer_to_matrix_write_nodal_results)
                    
//                     .def("PrintOnGaussPoints", pointer_to_double_print_on_gauss_points)
                    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
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

            enum_<GiD_PostMode>("GiDPostMode")
                    .value("GiD_PostAscii", GiD_PostAscii)
                    .value("GiD_PostAsciiZipped", GiD_PostAsciiZipped)
                    .value("GiD_PostBinary", GiD_PostBinary)
                    ;
            
            enum_<WriteDeformedMeshFlag>("WriteDeformedMeshFlag")
                    .value("WriteDeformed", WriteDeformed)
                    .value("WriteUndeformed", WriteUndeformed)
                    ;
            
            enum_<WriteConditionsFlag>("WriteConditionsFlag") 
                    .value("WriteConditions",WriteConditions)
                    .value("WriteElementsOnly",WriteElementsOnly)
                    .value("WriteConditionsOnly",WriteConditionsOnly)
                    ;
            
            enum_<MultiFileFlag>("MultiFileFlag")
                    .value("SingleFile",SingleFile)
                    .value("MultipleFiles",MultipleFiles)
                    ;
        }	
    }  // namespace Python.

} // Namespace Kratos

