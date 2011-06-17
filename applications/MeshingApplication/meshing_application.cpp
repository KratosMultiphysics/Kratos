//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-22 16:37:10 $
//   Revision:            $Revision: 1.6 $
//
// 



// System includes
  

// External includes   


// Project includes
#include "includes/define.h"
 
#include "meshing_application.h"
#include "includes/variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"


namespace Kratos
{
	//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE)
	KRATOS_CREATE_VARIABLE(double, COUNTER);

	KratosMeshingApplication::KratosMeshingApplication():
	mTestElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
	mTestElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))
	{}


	void KratosMeshingApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing Kratos MeshingApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE(COUNTER);
		
		KRATOS_REGISTER_ELEMENT("TestElement2D", mTestElement2D);
		KRATOS_REGISTER_ELEMENT("TestElement3D", mTestElement3D);
	}


	
}  // namespace Kratos.


