//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.19 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"
#include "incompressible_fluid_application.h"
#include "includes/variables.h"


namespace Kratos
{

//	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
//	//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
//	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)

	KRATOS_CREATE_VARIABLE( double, MACH_NUMBER )
	KRATOS_CREATE_VARIABLE( double, PRESSURE_COEFFICIENT )
//	KRATOS_CREATE_VARIABLE( double, PRESSURE_OLD_IT )
//	KRATOS_CREATE_VARIABLE( Vector, BDF_COEFFICIENTS )
//	KRATOS_CREATE_VARIABLE(double, NODAL_MASS)
//	KRATOS_CREATE_VARIABLE(int, AUX_INDEX)
//	KRATOS_CREATE_VARIABLE(double, EXTERNAL_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, DIAMETER)
	KRATOS_CREATE_VARIABLE(double, PERMEABILITY_INV)


	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)


	KratosIncompressibleFluidApplication::KratosIncompressibleFluidApplication():	
		mFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mFluid2DCoupled(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mFluid3DCoupled(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4,Node<3>())))),
                mFluid3DNeumann(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mNDFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> 
                >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mNDFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>
          	>(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mNDFluid2DCrankNicolson(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> 
                >(Element::GeometryType::PointsArrayType(3, Node<3>())))),


		mASGS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),

		mASGSCompressible2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mASGSCOMPPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),


		mMonolithic2DNeumann(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),

		mFluid2DGLS_expl(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mFluid2DGLS_expl_comp(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		
		mFluid2DSplit(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
	{}

	void KratosIncompressibleFluidApplication::Register()
	{
		std::cout << "Initializing KratosIncompressibleFluidApplication... " << std::endl;
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosIncompressibleFluidApplication...Register completed " << std::endl;

//		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
//		//KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
//		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)


		KRATOS_REGISTER_VARIABLE(  MACH_NUMBER )
		KRATOS_REGISTER_VARIABLE(  PRESSURE_COEFFICIENT )
//		KRATOS_REGISTER_VARIABLE( PRESSURE_OLD_IT )
//		KRATOS_REGISTER_VARIABLE( BDF_COEFFICIENTS )
//		KRATOS_REGISTER_VARIABLE( NODAL_MASS)
//		KRATOS_REGISTER_VARIABLE( AUX_INDEX)
//		KRATOS_REGISTER_VARIABLE( EXTERNAL_PRESSURE)
		KRATOS_REGISTER_VARIABLE( DIAMETER)
		KRATOS_REGISTER_VARIABLE( PERMEABILITY_INV)

		KRATOS_REGISTER_VARIABLE( DENSITY_AIR )

		KRATOS_REGISTER_VARIABLE( ARRHENIUS)


		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)

		std::cout << "Initializing KratosIncompressibleFluidApplication...variables succesfully registered " << std::endl;

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("Fluid3D", mFluid3D);
		KRATOS_REGISTER_ELEMENT("Fluid2D", mFluid2D);

		KRATOS_REGISTER_ELEMENT("Fluid2DCoupled", mFluid2DCoupled);
		KRATOS_REGISTER_ELEMENT("Fluid3DCoupled", mFluid3DCoupled);
		
		KRATOS_REGISTER_CONDITION("Fluid3DNeumann", mFluid3DNeumann);
		KRATOS_REGISTER_CONDITION("Monolithic2DNeumann", mMonolithic2DNeumann);
		
		KRATOS_REGISTER_ELEMENT("NDFluid2D", mNDFluid2D);
		KRATOS_REGISTER_ELEMENT("NDFluid3D", mNDFluid3D);

		KRATOS_REGISTER_ELEMENT("NDFluid2DCrankNicolson", mNDFluid2DCrankNicolson);


		KRATOS_REGISTER_ELEMENT("ASGS2D", mASGS2D);
		KRATOS_REGISTER_ELEMENT("ASGS3D", mASGS3D);	
		//KRATOS_REGISTER_ELEMENT("ASGSPRDC", mASGSPRDC);	
// 		KRATOS_REGISTER_ELEMENT("ASGSCompressible2D", mASGSCompressible2D);
// 		KRATOS_REGISTER_ELEMENT("ASGSCOMPPRDC2D", mASGSCOMPPRDC2D);

	


		KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl", mFluid2DGLS_expl);		
		KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl_comp", mFluid2DGLS_expl_comp);	

		KRATOS_REGISTER_ELEMENT("Fluid2DSplit", mFluid2DSplit);		

		std::cout << "Initializing KratosIncompressibleFluidApplication...elements succesfully registered " << std::endl;
	
	}





/*	// Initializing static members
	const Fluid3D  KratosIncompressibleFluidApplication::msFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));
	const Fluid2D  KratosIncompressibleFluidApplication::msFluid2D(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));

	const Fluid2DCoupled  KratosIncompressibleFluidApplication::msFluid2DCoupled(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));
	const Fluid3DCoupled  KratosIncompressibleFluidApplication::msFluid3DCoupled(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));
	*/
}  // namespace Kratos.


