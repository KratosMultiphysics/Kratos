// Kratos Multi-Physics
//
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.




// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"
#include "incompressible_fluid_application.h"
#include "includes/c2c_variables.h"
#include "includes/variables.h"
#include "geometries/point_3d.h"


namespace Kratos
{

//	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
//	//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
//	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)

KRATOS_CREATE_VARIABLE( double, MACH_NUMBER )
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
KRATOS_CREATE_VARIABLE(double, ENRICHED_PRESSURE)
KRATOS_CREATE_VARIABLE(double, ENRICHED_PRESSURE_IT)
KRATOS_CREATE_VARIABLE(double, INLET_VELOCITY)
KRATOS_CREATE_VARIABLE(double, INLET_PRESSURE)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL1)


KratosIncompressibleFluidApplication::KratosIncompressibleFluidApplication():
    KratosApplication("IncompressibleFluidApplication"),
//                mFluid2Dlevelset(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
    mFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mFluid2DCoupled(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mFluid3DCoupled(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFluid3DNeumann(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mNDFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>
               >(Element::GeometryType::PointsArrayType(3 )))),
    mNDFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>
               >(Element::GeometryType::PointsArrayType(4 )))),
    mNDFluid2DCrankNicolson(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>
                            >(Element::GeometryType::PointsArrayType(3 )))),



    mASGS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
mVP_PRECOND2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mASGSPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),




    mASGSCompressible2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
mASGS3D_ENR(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
 		mASGS3D_COMP_ENR(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
    mASGSCompressible3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
    mASGSCOMPPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mASGSCOMPPRDC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),

    mMonolithic2DNeumann(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2 )))),
   mMonolithic3DNeumann(0, Element::GeometryType::Pointer(new Point3D <Node<3> >(Element::GeometryType::PointsArrayType(1 )))),

    mFluid2DGLS_expl(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mFluid3DGLS_expl(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),

    mFluid2DGLS(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),

    mFluid2DGLS_expl_comp(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mFluid3DGLS_expl_comp(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),

    mProjDirichletCond(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mProjDirichletCond3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),

    mNoNewtonianASGS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mNoNewtonianASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
    mBinghamNonNewtonianASGS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mBinghamNonNewtonianASGS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),

    mNoSlipCondition2D(0, Element::GeometryType::Pointer(new Geometry <Node<3>  >(Element::GeometryType::PointsArrayType(2 )))),
// 		mNoSlipFractStep(0, Element::GeometryType::Pointer(new Geometry <Node<3>  >(Element::GeometryType::PointsArrayType(2 ))))

    mExplicitASGSCompressible2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),
    mExplicitASGSCOMPPRDC2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 )))),


    mExplicitASGSCompressible3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 )))),
    mExplicitASGSCOMPPRDC3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 ))))


    //mABC2D(0, Element::GeometryType::Pointer(new Point3D <Node<3> >(Element::GeometryType::PointsArrayType(1 ))))
    //mExplicitHydro2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 ))))




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
    KRATOS_REGISTER_VARIABLE(FRONT_MEETING)
    KRATOS_REGISTER_VARIABLE(INLET_VELOCITY)
	KRATOS_REGISTER_VARIABLE(INLET_PRESSURE)

    KRATOS_REGISTER_VARIABLE(ENRICHED_PRESSURE)
    KRATOS_REGISTER_VARIABLE(ENRICHED_PRESSURE_IT)
    KRATOS_REGISTER_VARIABLE( ACTIVATE_TAU2)


		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
		//KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
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
	KRATOS_REGISTER_CONDITION("Monolithic3DNeumann", mMonolithic3DNeumann);

	KRATOS_REGISTER_ELEMENT("NDFluid2D", mNDFluid2D);
	KRATOS_REGISTER_ELEMENT("NDFluid3D", mNDFluid3D);


		KRATOS_REGISTER_ELEMENT("ASGS2D", mASGS2D);
		KRATOS_REGISTER_ELEMENT("VP_PRECOND2D", mVP_PRECOND2D);
		KRATOS_REGISTER_ELEMENT("ASGS3D", mASGS3D);
		KRATOS_REGISTER_ELEMENT("ASGS3D_ENR", mASGS3D_ENR);
 		KRATOS_REGISTER_ELEMENT("ASGS3D_COMP_ENR", mASGS3D_COMP_ENR);
		KRATOS_REGISTER_ELEMENT("ASGSPRDC2D", mASGSPRDC2D);
 		KRATOS_REGISTER_ELEMENT("ASGSCompressible2D", mASGSCompressible2D);
 		KRATOS_REGISTER_ELEMENT("ASGSCompressible3D", mASGSCompressible3D);
 		KRATOS_REGISTER_ELEMENT("ASGSCOMPPRDC2D", mASGSCOMPPRDC2D);
 		KRATOS_REGISTER_ELEMENT("ASGSCOMPPRDC3D", mASGSCOMPPRDC3D);

    KRATOS_REGISTER_ELEMENT("NDFluid2DCrankNicolson", mNDFluid2DCrankNicolson);


    //	KRATOS_REGISTER_ELEMENT("ASGS2DPARTICLE", mASGS2DPARTICLE);

    KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl", mFluid2DGLS_expl);
    KRATOS_REGISTER_ELEMENT("Fluid3DGLS_expl", mFluid3DGLS_expl);

    KRATOS_REGISTER_ELEMENT("Fluid2DGLS", mFluid2DGLS);

    KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl_comp", mFluid2DGLS_expl_comp);
    KRATOS_REGISTER_ELEMENT("Fluid3DGLS_expl_comp", mFluid3DGLS_expl_comp);


    KRATOS_REGISTER_CONDITION("ProjDirichletCond", mProjDirichletCond);
    KRATOS_REGISTER_CONDITION("ProjDirichletCond3D", mProjDirichletCond3D);
    KRATOS_REGISTER_CONDITION("NoSlipCondition2D", mNoSlipCondition2D);
// 		KRATOS_REGISTER_CONDITION("NoSlipFractStep", mNoSlipFractStep);

    //KRATOS_REGISTER_CONDITION("ABC2D", mABC2D);


    KRATOS_REGISTER_ELEMENT("NoNewtonianASGS2D", mNoNewtonianASGS2D);
    KRATOS_REGISTER_ELEMENT("NoNewtonianASGS3D", mNoNewtonianASGS3D);
    KRATOS_REGISTER_ELEMENT("BinghamNonNewtonianASGS2D", mBinghamNonNewtonianASGS2D);
    KRATOS_REGISTER_ELEMENT("BinghamNonNewtonianASGS3D", mBinghamNonNewtonianASGS3D);


    KRATOS_REGISTER_ELEMENT("ExplicitASGSCompressible2D", mExplicitASGSCompressible2D);
    KRATOS_REGISTER_ELEMENT("ExplicitASGSCOMPPRDC2D", mExplicitASGSCOMPPRDC2D);


    KRATOS_REGISTER_ELEMENT("ExplicitASGSCompressible3D", mExplicitASGSCompressible3D);
    KRATOS_REGISTER_ELEMENT("ExplicitASGSCOMPPRDC3D", mExplicitASGSCOMPPRDC3D);


    //KRATOS_REGISTER_ELEMENT("ExplicitHydro2D", mExplicitHydro2D);

    std::cout << "Initializing KratosIncompressibleFluidApplication...elements succesfully registered " << std::endl;

}





/*	// Initializing static members
	const Fluid3D  KratosIncompressibleFluidApplication::msFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 ))));
	const Fluid2D  KratosIncompressibleFluidApplication::msFluid2D(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3 ))));

	const Fluid2DCoupled  KratosIncompressibleFluidApplication::msFluid2DCoupled(0, Element::GeometryType::Pointer(new Triangle2D<Node<3> >(Element::GeometryType::PointsArrayType(3 ))));
	const Fluid3DCoupled  KratosIncompressibleFluidApplication::msFluid3DCoupled(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4 ))));
	*/
}  // namespace Kratos.
