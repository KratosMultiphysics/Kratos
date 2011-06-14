//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-15 15:41:36 $
//   Revision:            $Revision: 1.6 $
//
// 



// System includes


// External includes 
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "convection_diffusion_application.h"
#include "includes/variables.h"


namespace Kratos
{
	
	
//	KRATOS_CREATE_VARIABLE( Vector, BDF_COEFFICIENTS )
	//KRATOS_CREATE_VARIABLE(double, NODAL_AREA)
//	KRATOS_CREATE_VARIABLE(int, AUX_INDEX)
//	KRATOS_CREATE_VARIABLE(double,  CONDUCTIVITY)
//	KRATOS_CREATE_VARIABLE(double,  SPECIFIC_HEAT)
//	KRATOS_CREATE_VARIABLE(double,  HEAT_FLUX)
//	KRATOS_CREATE_VARIABLE(double,  TEMP_CONV_PROJ)

	//Added by Pavel and Annelie
//	KRATOS_CREATE_VARIABLE(double,  ENTHALPY)
	KRATOS_CREATE_VARIABLE(double,  LATENT_HEAT)	
	KRATOS_CREATE_VARIABLE(double,  MELT_TEMPERATURE_1)
	KRATOS_CREATE_VARIABLE(double,  MELT_TEMPERATURE_2)

	KRATOS_CREATE_VARIABLE(double,  AMBIENT_TEMPERATURE)	
//	KRATOS_CREATE_VARIABLE(double,  EMISSIVITY)
	KRATOS_CREATE_VARIABLE(double,  CONVECTION_COEFFICIENT)	
//	KRATOS_CREATE_VARIABLE(double,  FACE_HEAT_FLUX)
			
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

	KratosConvectionDiffusionApplication::KratosConvectionDiffusionApplication():
		mConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        	mConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mConvDiffChangeOfPhase2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mThermalFace2D(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
		mThermalFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
	{}


		
	void KratosConvectionDiffusionApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosConvectionDiffusionApplication... " << std::endl;
		

//		KRATOS_REGISTER_VARIABLE( BDF_COEFFICIENTS );
		//KRATOS_REGISTER_VARIABLE( NODAL_AREA)
//		KRATOS_REGISTER_VARIABLE( AUX_INDEX)
//		KRATOS_REGISTER_VARIABLE( CONDUCTIVITY)
//		KRATOS_REGISTER_VARIABLE( SPECIFIC_HEAT)
//		KRATOS_REGISTER_VARIABLE( HEAT_FLUX)
//		KRATOS_REGISTER_VARIABLE( TEMP_CONV_PROJ)

		//Added by Pavel and Annelie
//		KRATOS_REGISTER_VARIABLE(ENTHALPY)
		KRATOS_REGISTER_VARIABLE(LATENT_HEAT)	
		KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_1)
		KRATOS_REGISTER_VARIABLE(MELT_TEMPERATURE_2)

		KRATOS_REGISTER_VARIABLE(AMBIENT_TEMPERATURE)	
//		KRATOS_REGISTER_VARIABLE(EMISSIVITY)
		KRATOS_REGISTER_VARIABLE(CONVECTION_COEFFICIENT)	
//		KRATOS_REGISTER_VARIABLE(FACE_HEAT_FLUX)
					
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("ConvDiff2D", mConvDiff2D);
		KRATOS_REGISTER_ELEMENT("ConvDiff3D", mConvDiff3D);

		KRATOS_REGISTER_ELEMENT("ConvDiffChangeOfPhase2D", mConvDiffChangeOfPhase2D);
	
		KRATOS_REGISTER_CONDITION("ThermalFace2D", mThermalFace2D);
		KRATOS_REGISTER_CONDITION("ThermalFace3D", mThermalFace3D);
		
		// Registering elements and conditions here
//		KRATOS_REGISTER_ELEMENT("ConvDiff2D", msConvDiff2D);
//		KRATOS_REGISTER_ELEMENT("ConvDiff3D", msConvDiff3D);

//		KRATOS_REGISTER_ELEMENT("ConvDiff2DChangeOfPhase", mConvDiff2DChangeOfPhase);
		
//		KRATOS_REGISTER_CONDITION("ThermalFace2D", msThermalFace2D);
//		KRATOS_REGISTER_CONDITION("ThermalFace3D", msThermalFace3D);
	}




/*
	// Initializing static members
	const ConvDiff2D  KratosConvectionDiffusionApplication::msConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));

	const ConvDiff3D  KratosConvectionDiffusionApplication::msConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));
	
		const ThermalFace2D  KratosConvectionDiffusionApplication::msThermalFace2D(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>()))));

		const ThermalFace3D  KratosConvectionDiffusionApplication::msThermalFace3D(0, Element::GeometryType::Pointer(new Triangle3D<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));
*/		
}  // namespace Kratos.


