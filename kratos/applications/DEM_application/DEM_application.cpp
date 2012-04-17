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


#include "includes/variables.h"
#include "DEM_application.h"
#include "geometries/point_3d.h"

namespace Kratos
{
	/*
	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);*/


	KratosDEMApplication::KratosDEMApplication():
	mSphericParticle(0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
	{}
        
	void KratosDEMApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosDEMApplication... " << std::endl;
		KRATOS_REGISTER_ELEMENT("SphericParticle", mSphericParticle)
	}

}  // namespace Kratos.


