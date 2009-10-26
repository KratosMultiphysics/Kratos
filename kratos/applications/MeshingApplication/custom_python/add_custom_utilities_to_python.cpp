//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "custom_utilities/projection.h"
#include "custom_utilities/GenerateModelPartUtilities.h"


namespace Kratos
{
	
namespace Python
{

	void GenerateModelTemperaturePart(GenerateModelPartUtilities& GM,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
	{
		if(domain_size == 2)
		{
			GM.GenerateModelPart(origin_model_part, destination_model_part,
				KratosComponents<Element>::Get("ConvDiff2D"),
				KratosComponents<Condition>::Get("ThermalFace2D")	); 
		}
		else if(domain_size == 3)
		{
			GM.GenerateModelPart(origin_model_part, destination_model_part,
				KratosComponents<Element>::Get("ConvDiff3D"),
				KratosComponents<Condition>::Get("ThermalFace3D")	); 
		}
	}

  void  AddCustomUtilitiesToPython()
  {
	
	using namespace boost::python;
	
	  class_<MeshTransfer<2> >("MeshTransfer2D",init< >()) 
	  	.def("DirectModelPartInterpolation", &MeshTransfer<2>::DirectInterpolation)
		.def("DirectScalarVarInterpolation", &MeshTransfer<2>::DirectVariableInterpolation<double>)
		.def("DirectVectorialVarInterpolation", &MeshTransfer<2>::DirectVariableInterpolation< array_1d < double,3 > >);
	  
	 class_<MeshTransfer<3> >("MeshTransfer3D",init< >()) 
	  	.def("DirectModelPartInterpolation", &MeshTransfer<3>::DirectInterpolation)
		.def("DirectScalarVarInterpolation", &MeshTransfer<3>::DirectVariableInterpolation<double>)
		.def("DirectVectorialVarInterpolation", &MeshTransfer<3>::DirectVariableInterpolation< array_1d < double,3 > >);

	 class_<GenerateModelPartUtilities >("GenerateModelPartUtilities",init< >()) 
	   .def("GenerateModelTemperaturePart", GenerateModelTemperaturePart );
	



  }

	
}  // namespace Python.

} // Namespace Kratos

