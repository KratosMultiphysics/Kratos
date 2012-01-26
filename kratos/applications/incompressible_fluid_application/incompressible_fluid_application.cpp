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
#include "geometries/point_3d.h"


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
//	KRATOS_CREATE_VARIABLE(double, DIAMETER)
	KRATOS_CREATE_VARIABLE(double, PERMEABILITY_INV)
	//for disabling elements
	KRATOS_CREATE_VARIABLE(int, DISABLE)
	
	KRATOS_CREATE_VARIABLE(int, ACTIVATE_TAU2)


	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL1)


	KratosIncompressibleFluidApplication::KratosIncompressibleFluidApplication():
//                mFluid2Dlevelset(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
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
                mASGSPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		//mASGS2DPARTICLE(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),



		mASGSCompressible2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mASGSCompressible3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mASGSCOMPPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mASGSCOMPPRDC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),

		mMonolithic2DNeumann(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
		//mMonolithic3DNeumann(0, Element::GeometryType::Pointer(new Point3D <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),

		mFluid2DGLS_expl(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mFluid3DGLS_expl(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),

		mFluid2DGLS(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),

		mFluid2DGLS_expl_comp(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mFluid3DGLS_expl_comp(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),

		mProjDirichletCond(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
			
		mFluid2DSplit(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),

		mNoNewtonianASGS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mNoNewtonianASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mBinghamNonNewtonianASGS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mBinghamNonNewtonianASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),

		mNoSlipCondition2D(0, Element::GeometryType::Pointer(new Geometry <Node<3>  >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
// 		mNoSlipFractStep(0, Element::GeometryType::Pointer(new Geometry <Node<3>  >(Element::GeometryType::PointsArrayType(2, Node<3>()))))

		mExplicitASGSCompressible2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
		mExplicitASGSCOMPPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),


		mExplicitASGSCompressible3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mExplicitASGSCOMPPRDC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))

		//mABC2D(0, Element::GeometryType::Pointer(new Point3D <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>()))))




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
		//for disabling elements
		KRATOS_REGISTER_VARIABLE(DISABLE)
		
		KRATOS_REGISTER_VARIABLE( ACTIVATE_TAU2)


		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL1)

		std::cout << "Initializing KratosIncompressibleFluidApplication...variables succesfully registered " << std::endl;

		// Registering elements and conditions here
//		KRATOS_REGISTER_ELEMENT("Fluid2Dlevelset", mFluid2Dlevelset);
		KRATOS_REGISTER_ELEMENT("Fluid3D", mFluid3D);
		KRATOS_REGISTER_ELEMENT("Fluid2D", mFluid2D);

		KRATOS_REGISTER_ELEMENT("Fluid2DCoupled", mFluid2DCoupled);
		KRATOS_REGISTER_ELEMENT("Fluid3DCoupled", mFluid3DCoupled);
		
		KRATOS_REGISTER_CONDITION("Fluid3DNeumann", mFluid3DNeumann);
		KRATOS_REGISTER_CONDITION("Monolithic2DNeumann", mMonolithic2DNeumann);
		//KRATOS_REGISTER_CONDITION("Monolithic3DNeumann", mMonolithic3DNeumann);
		
		KRATOS_REGISTER_ELEMENT("NDFluid2D", mNDFluid2D);
		KRATOS_REGISTER_ELEMENT("NDFluid3D", mNDFluid3D);

		KRATOS_REGISTER_ELEMENT("NDFluid2DCrankNicolson", mNDFluid2DCrankNicolson);


		KRATOS_REGISTER_ELEMENT("ASGS2D", mASGS2D);
		KRATOS_REGISTER_ELEMENT("ASGS3D", mASGS3D);	
		KRATOS_REGISTER_ELEMENT("ASGSPRDC2D", mASGSPRDC2D);	
 		KRATOS_REGISTER_ELEMENT("ASGSCompressible2D", mASGSCompressible2D);
 		KRATOS_REGISTER_ELEMENT("ASGSCompressible3D", mASGSCompressible3D);
 		KRATOS_REGISTER_ELEMENT("ASGSCOMPPRDC2D", mASGSCOMPPRDC2D);
 		KRATOS_REGISTER_ELEMENT("ASGSCOMPPRDC3D", mASGSCOMPPRDC3D);

	
	//	KRATOS_REGISTER_ELEMENT("ASGS2DPARTICLE", mASGS2DPARTICLE);

		KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl", mFluid2DGLS_expl);
		KRATOS_REGISTER_ELEMENT("Fluid3DGLS_expl", mFluid3DGLS_expl);		

		KRATOS_REGISTER_ELEMENT("Fluid2DGLS", mFluid2DGLS);

		KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl_comp", mFluid2DGLS_expl_comp);	
		KRATOS_REGISTER_ELEMENT("Fluid3DGLS_expl_comp", mFluid3DGLS_expl_comp);	


		KRATOS_REGISTER_CONDITION("ProjDirichletCond", mProjDirichletCond);
		KRATOS_REGISTER_CONDITION("NoSlipCondition2D", mNoSlipCondition2D);
// 		KRATOS_REGISTER_CONDITION("NoSlipFractStep", mNoSlipFractStep);

		//KRATOS_REGISTER_CONDITION("ABC2D", mABC2D);


		KRATOS_REGISTER_ELEMENT("Fluid2DSplit", mFluid2DSplit);
		
		KRATOS_REGISTER_ELEMENT("NoNewtonianASGS2D", mNoNewtonianASGS2D);	
		KRATOS_REGISTER_ELEMENT("NoNewtonianASGS3D", mNoNewtonianASGS3D);		
		KRATOS_REGISTER_ELEMENT("BinghamNonNewtonianASGS2D", mBinghamNonNewtonianASGS2D);
		KRATOS_REGISTER_ELEMENT("BinghamNonNewtonianASGS3D", mBinghamNonNewtonianASGS3D);


		KRATOS_REGISTER_ELEMENT("ExplicitASGSCompressible2D", mExplicitASGSCompressible2D);
		KRATOS_REGISTER_ELEMENT("ExplicitASGSCOMPPRDC2D", mExplicitASGSCOMPPRDC2D);

		
		KRATOS_REGISTER_ELEMENT("ExplicitASGSCompressible3D", mExplicitASGSCompressible3D);
 		KRATOS_REGISTER_ELEMENT("ExplicitASGSCOMPPRDC3D", mExplicitASGSCOMPPRDC3D);		
		
		std::cout << "Initializing KratosIncompressibleFluidApplication...elements succesfully registered " << std::endl;
	
	}





/*	// Initializing static members
	const Fluid3D  KratosIncompressibleFluidApplication::msFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));
	const Fluid2D  KratosIncompressibleFluidApplication::msFluid2D(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));

	const Fluid2DCoupled  KratosIncompressibleFluidApplication::msFluid2DCoupled(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));
	const Fluid3DCoupled  KratosIncompressibleFluidApplication::msFluid3DCoupled(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));
	*/
}  // namespace Kratos.


