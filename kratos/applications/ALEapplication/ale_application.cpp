//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-05-06 15:11:09 $
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
#include "geometries/line_2d.h"
#include "ale_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

	KratosALEApplication::KratosALEApplication():
		mLaplacianMeshMovingElem2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mLaplacianMeshMovingElem3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))
	{}
	
	void KratosALEApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosALEApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
//		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
//		KRATOS_REGISTER_VARIABLE(NODAL_AREA);
	
		KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElem2D", mLaplacianMeshMovingElem2D);
		KRATOS_REGISTER_ELEMENT("LaplacianMeshMovingElemt3D", mLaplacianMeshMovingElem3D);

	}

}  // namespace Kratos.


