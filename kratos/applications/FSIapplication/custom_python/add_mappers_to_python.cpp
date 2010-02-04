//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:42 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_mappers_to_python.h"
#include "custom_utilities/NMPointsMapper.h"
#include "custom_utilities/AdvancedNMPointsMapper.h"
#include "custom_utilities/InterfacePreprocess.h"

#include "custom_utilities/shared_points_mapper.h" 
#include "includes/node.h"

namespace Kratos
{
	
namespace Python
{
  void  AddMappersToPython()
  {
	
	using namespace boost::python;

	  class_<SharedPointsMapper >("SharedPointsMapper",
		 init< const ModelPart::NodesContainerType&, const  ModelPart::NodesContainerType&, double>()) 
	  .def("ScalarMap",&SharedPointsMapper::ScalarMap)
	  .def("VectorMap",&SharedPointsMapper::VectorMap)
		;

	  class_<NMPointsMapper>("NMPointsMapper", init<ModelPart&, ModelPart&>())
		.def("FindNeighbours",&NMPointsMapper::FindNeighbours)
		.def("ScalarMap",&NMPointsMapper::ScalarMap< double >)
		.def("VectorMap",&NMPointsMapper::VectorMap< array_1d<double,3> > )
		;

	  class_<AdvancedNMPointsMapper>("AdvancedNMPointsMapper", init<const ModelPart&, ModelPart&>())
		.def("FindNeighbours",&AdvancedNMPointsMapper::FindNeighbours)
		.def("ScalarMap",&AdvancedNMPointsMapper::ScalarMap)
                .def("VectorMap",&AdvancedNMPointsMapper::VectorMap)
		;

          class_<InterfacePreprocess>("InterfacePreprocess", init<>())
                .def("GenerateInterfacePart",&InterfacePreprocess::GenerateInterfacePart)
                ;
  }
	
}  // namespace Python.

} // Namespace Kratos

