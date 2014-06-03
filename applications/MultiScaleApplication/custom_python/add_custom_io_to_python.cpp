/*
==============================================================================
KratosMultiScaleApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-12-30 12:00:00 $
//   Revision:            $Revision: 1.00 $
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
#include "../custom_utilities/gid_io_extended.h"

#include "add_custom_io_to_python.h"
#include "custom_utilities/gid_extended_gauss_point_container.h"

namespace Kratos
{

namespace Python
{



typedef GidExtendedGaussPointsContainer GaussPointContainterType;

typedef GidIOExtended<GaussPointContainterType> GiDIOType;



void (GiDIOType::*pointer_to_io_read_single_properties)(Properties& rThisProperties ) = &IO::ReadProperties;

void (GiDIOType::*pointer_to_io_read_properties)(IO::PropertiesContainerType& rThisProperties ) = &IO::ReadProperties;

std::string PrintDatafileIO(DatafileIO const& rDataFileIO)
{
    std::stringstream buffer;
    rDataFileIO.PrintInfo(buffer);
    buffer << std::endl;
    rDataFileIO.PrintData(buffer);

    return buffer.str();
}

void WriteNodeMesh( GiDIOType& dummy, GiDIOType::MeshType& rThisMesh )
{
    std::cout<<"start printing nodes mesh "<<std::endl;
    dummy.WriteNodeMesh( rThisMesh );
    std::cout<<"end printing nodes mesh "<<std::endl;
}

void WriteSphereMesh( GiDIOType& dummy, GiDIOType::MeshType& rThisMesh )
{
    dummy.WriteSphereMesh( rThisMesh );
}


void WriteMesh( GiDIOType& dummy, GiDIOType::MeshType& rThisMesh )
{
    dummy.WriteMesh( rThisMesh );
}


void DoublePrintOnGaussPoints( GiDIOType& dummy, const Variable<double>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void Array1DPrintOnGaussPoints( GiDIOType& dummy, const Variable<array_1d<double,3> >& rVariable,
                                ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void VectorPrintOnGaussPoints( GiDIOType& dummy, const Variable<Vector>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}

void MatrixPrintOnGaussPoints( GiDIOType& dummy, const Variable<Matrix>& rVariable,
                               ModelPart& r_model_part, double SolutionTag )
{
    dummy.PrintOnGaussPoints( rVariable, r_model_part, SolutionTag );
}


void (GiDIOType::*pointer_to_bool_write_nodal_results)( Variable<bool> const& rVariable,
        GiDIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GiDIOType::WriteNodalResults;
		
void (GiDIOType::*pointer_to_double_write_nodal_results)( Variable<double> const& rVariable,
        GiDIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GiDIOType::WriteNodalResults;
		
void (GiDIOType::*pointer_to_array1d_write_nodal_results)(
    Variable<array_1d<double, 3> > const& rVariable, GiDIOType::NodesContainerType& rNodes,
    double SolutionTag, std::size_t SolutionStepNumber) = &GiDIOType::WriteNodalResults;

void (GiDIOType::*pointer_to_matrix_write_nodal_results)(
	Variable<Matrix > const& rVariable, GiDIOType::NodesContainerType& rNodes, 
	double SolutionTag, std::size_t SolutionStepNumber) = &GiDIOType::WriteNodalResults;

void (GiDIOType::*local_axes_write_nodal_results)( Variable<array_1d<double, 3> > const& rVariable,
        GiDIOType::NodesContainerType& rNodes, double SolutionTag,
        std::size_t SolutionStepNumber ) = &GiDIOType::WriteLocalAxesOnNodes;

void  AddCustomIOToPython()
{
    using namespace boost::python;

  //  class_<IO, IO::Pointer, boost::noncopyable>("IO_Extended")
  //  .def("ReadNode",&IO::ReadNode)
  //  .def("ReadNodes",&IO::ReadNodes)
  //  .def("WriteNodes",&IO::WriteNodes)
  //  .def("ReadProperties",pointer_to_io_read_single_properties)
  //  .def("ReadProperties",pointer_to_io_read_properties)
  //  .def("ReadElements",&IO::ReadElements)
  //  .def("WriteElements",&IO::WriteElements)
  //  .def("ReadConditions",&IO::ReadConditions)
  //  .def("ReadInitialValues",&IO::ReadInitialValues)
  //  .def("ReadMesh",&IO::ReadMesh)
  //  .def("ReadModelPart",&IO::ReadModelPart)
  //  .def("__str__",PrintDatafileIO)
  //  ;

  //  class_<DatafileIO, DatafileIO::Pointer, bases<IO>,  boost::noncopyable>(
  //      "DatafileIO_Extended",
		//init<std::string const&>())
  //  .def(init<std::string const&,
  //       std::string const&,
  //       std::string const&,
  //       std::string const&,
  //       std::string const&>())
  //  ;

  //  class_<ModelPartIO, ModelPartIO::Pointer, bases<IO>,  boost::noncopyable>(
  //      "ModelPartIO_Extended",
		//init<std::string const&>())
  //  ;

    class_<GiDIOType, GiDIOType::Pointer, bases<DatafileIO>, boost::noncopyable>(
        "GidIO_Extended",
		init<std::string const&, GiD_PostMode,
			MultiFileFlag,
			WriteDeformedMeshFlag,
			WriteConditionsFlag>())
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

    .def("InitializeMesh",&GiDIOType::InitializeMesh)
    .def("FinalizeMesh",&GiDIOType::FinalizeMesh)

    .def("InitializeResults",&GiDIOType::InitializeResults)
    .def("FinalizeResults",&GiDIOType::FinalizeResults)

    .def("WriteNodalResults",pointer_to_bool_write_nodal_results)
    .def("WriteNodalResults",pointer_to_double_write_nodal_results)
    .def("WriteNodalResults",pointer_to_array1d_write_nodal_results)

    .def("WriteNodalResults",pointer_to_matrix_write_nodal_results)

    .def("WriteLocalAxesOnNodes",local_axes_write_nodal_results)

    .def("PrintOnGaussPoints", DoublePrintOnGaussPoints)
    .def("PrintOnGaussPoints", Array1DPrintOnGaussPoints)
    .def("PrintOnGaussPoints", VectorPrintOnGaussPoints)
    .def("PrintOnGaussPoints", MatrixPrintOnGaussPoints)

    .def("Flush",&GiDIOType::Flush)
    .def("ChangeOutputName",&GiDIOType::ChangeOutputName)
    .def("CloseResultFile",&GiDIOType::CloseResultFile)
    ;

 //   enum_<GiD_PostMode>("GiDPostMode")
 //   .value("GiD_PostAscii", GiD_PostAscii)
 //   .value("GiD_PostAsciiZipped", GiD_PostAsciiZipped)
 //   .value("GiD_PostBinary", GiD_PostBinary)
	//.value("GiD_PostHDF5", GiD_PostHDF5)
 //   ;

 //   enum_<WriteDeformedMeshFlag>("WriteDeformedMeshFlag")
 //   .value("WriteDeformed", WriteDeformed)
 //   .value("WriteUndeformed", WriteUndeformed)
 //   ;

 //   enum_<WriteConditionsFlag>("WriteConditionsFlag")
 //   .value("WriteConditions",WriteConditions)
 //   .value("WriteElementsOnly",WriteElementsOnly)
 //   .value("WriteConditionsOnly",WriteConditionsOnly)
 //   ;

 //   enum_<MultiFileFlag>("MultiFileFlag")
 //   .value("SingleFile",SingleFile)
 //   .value("MultipleFiles",MultipleFiles)
    ;
}
}

}

