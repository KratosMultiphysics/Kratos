//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jmora $ rrossi $
//   Date:                $Date: 2010-02-02 $
//   Revision:            $Revision: 1.4 $
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
#include "kMagnetostatic.h"
#include "includes/variables.h"


namespace Kratos
{

	KRATOS_CREATE_VARIABLE(double,  CONDUCTIVITY)
	KRATOS_CREATE_VARIABLE(double,  SPECIFIC_HEAT)
	KRATOS_CREATE_VARIABLE(double,  HEAT_FLUX)	
	KRATOS_CREATE_VARIABLE(double,  TEMP_CONV_PROJ)	

	KRATOS_CREATE_VARIABLE(double,  AMBIENT_TEMPERATURE)	
	KRATOS_CREATE_VARIABLE(double,  EMISSIVITY)	
	KRATOS_CREATE_VARIABLE(double,  CONVECTION_COEFFICIENT)	
	KRATOS_CREATE_VARIABLE(double,  FACE_HEAT_FLUX)	
		
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

	// for electromagnetic applications
	// for kMagnetostatic application
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_PERMEABILITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(COERCIVITY)
	KRATOS_CREATE_VARIABLE(double, MAGNETOSTATIC_POTENTIAL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAGNETOSTATIC_VECTOR_POTENTIAL)

	KRATOS_CREATE_VARIABLE(double, MAGNETOSTATIC_POINT_CURRENT)
	KRATOS_CREATE_VARIABLE(double, MAGNETOSTATIC_SURFACE_CURRENT)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAGNETOSTATIC_VOLUME_CURRENT)

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_FIELD_INTENSITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_FLUX_DENSITY)
	KRATOS_CREATE_VARIABLE(double, INFINIT_COEFFICIENT)



	KratosR1MagnetostaticApplication::KratosR1MagnetostaticApplication():
		mMagnetostatic2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mMagnetostatic3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mPointCurrent2D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
		mMfield2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
		mMfield3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
	{}

		
	void KratosR1MagnetostaticApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosR1MagnetostaticApplication... " << std::endl;

		KRATOS_REGISTER_VARIABLE( BDF_COEFFICIENTS );
		//KRATOS_REGISTER_VARIABLE( NODAL_AREA)
		KRATOS_REGISTER_VARIABLE( AUX_INDEX)
		KRATOS_REGISTER_VARIABLE( CONDUCTIVITY)
		KRATOS_REGISTER_VARIABLE( SPECIFIC_HEAT)
		KRATOS_REGISTER_VARIABLE( HEAT_FLUX)
		KRATOS_REGISTER_VARIABLE( TEMP_CONV_PROJ)

		KRATOS_REGISTER_VARIABLE(AMBIENT_TEMPERATURE)	
		KRATOS_REGISTER_VARIABLE(EMISSIVITY)	
		KRATOS_REGISTER_VARIABLE(CONVECTION_COEFFICIENT)	
		KRATOS_REGISTER_VARIABLE(FACE_HEAT_FLUX)
				
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

		// for electromagnetic applications
		// for kMagnetostatic application
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_PERMEABILITY)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(COERCIVITY)
		KRATOS_REGISTER_VARIABLE(MAGNETOSTATIC_POTENTIAL)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAGNETOSTATIC_VECTOR_POTENTIAL)
		KRATOS_REGISTER_VARIABLE(MAGNETOSTATIC_POINT_CURRENT)
		KRATOS_REGISTER_VARIABLE(MAGNETOSTATIC_SURFACE_CURRENT)
	    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAGNETOSTATIC_VOLUME_CURRENT)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_FIELD_INTENSITY)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_FLUX_DENSITY)
		KRATOS_REGISTER_VARIABLE(INFINIT_COEFFICIENT)


		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Magnetostatic2D", mMagnetostatic2D);
		KRATOS_REGISTER_ELEMENT("Magnetostatic3D", mMagnetostatic3D);
		KRATOS_REGISTER_CONDITION("PointCurrent2D", mPointCurrent2D);
		KRATOS_REGISTER_CONDITION("Mfield2D", mMfield2D);
		KRATOS_REGISTER_CONDITION("Mfield3D", mMfield3D);


	}
	
}  // namespace Kratos.


