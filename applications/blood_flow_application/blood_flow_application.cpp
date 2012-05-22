//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_3d_2.h"
#include "blood_flow_application.h"
#include "includes/variables.h"
#include "custom_elements/artery_element.h"

namespace Kratos
{
	//Example
 	KRATOS_CREATE_VARIABLE(double, FLOW)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WORK );
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosBloodFlowApplication::KratosBloodFlowApplication(): 
	   mArteryElement(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>()))))
 	{}
 	
 	void KratosBloodFlowApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosBloodFlowApplication... " << std::endl;
 
 		KRATOS_REGISTER_VARIABLE( FLOW )
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WORK )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

		KRATOS_REGISTER_ELEMENT("ArteryElement", mArteryElement);
 
 	}

}  // namespace Kratos.


