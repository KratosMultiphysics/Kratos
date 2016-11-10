//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d.h"
#include "geometries/line_2d_3.h"
#include "fluid_dynamics_application.h"
#include "includes/variables.h"


namespace Kratos
{

KratosFluidDynamicsApplication::KratosFluidDynamicsApplication():
    mVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoFluidVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStationaryStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mStationaryStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFractionalStepDiscontinuous2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mFractionalStepDiscontinuous3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSpalartAllmaras2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3))),GeometryData::GI_GAUSS_2),
    mSpalartAllmaras3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))),GeometryData::GI_GAUSS_2),
    mWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mFSWernerWengleWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSWernerWengleWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mFSGeneralizedWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSGeneralizedWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mWallConditionDiscontinuous2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mWallConditionDiscontinuous3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMonolithicWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mMonolithicWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mStokesWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mFSPeriodicCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSPeriodicCondition3D(0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSPeriodicConditionEdge2D(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mFSPeriodicConditionEdge3D(0, Element::GeometryType::Pointer( new Quadrilateral3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mDPGVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDPGVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mBinghamVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mBinghamVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mBinghamFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mBinghamFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mBinghamFractionalStepDiscontinuous2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mBinghamFractionalStepDiscontinuous3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mHerschelBulkleyVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mHerschelBulkleyVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStokes3DTwoFluid(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    // Navier-Stokes symbolic elements
    mNavierStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mNavierStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    // Embedded Navier-Stokes symbolic elements
    mEmbeddedNavierStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedNavierStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}

void KratosFluidDynamicsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosFluidDynamicsApplication... " << std::endl;

    // Register Variables (defined in fluid_dynamics_application_variables.h)
    KRATOS_REGISTER_VARIABLE(PATCH_INDEX);
    KRATOS_REGISTER_VARIABLE(TAUONE);
    KRATOS_REGISTER_VARIABLE(TAUTWO);
    KRATOS_REGISTER_VARIABLE(PRESSURE_MASSMATRIX_COEFFICIENT)

//    KRATOS_REGISTER_VARIABLE(Y_WALL);
    KRATOS_REGISTER_VARIABLE(SUBSCALE_PRESSURE);
    KRATOS_REGISTER_VARIABLE(C_DES);
//    KRATOS_REGISTER_VARIABLE(C_SMAGORINSKY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VORTICITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY);
    
    // Non-Newtonian constitutive relations
    KRATOS_REGISTER_VARIABLE(REGULARIZATION_COEFFICIENT);

    KRATOS_REGISTER_VARIABLE(BINGHAM_SMOOTHER);
    KRATOS_REGISTER_VARIABLE(GEL_STRENGTH);

    KRATOS_REGISTER_VARIABLE(Q_VALUE);
    KRATOS_REGISTER_VARIABLE(VORTICITY_MAGNITUDE);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT);
    KRATOS_REGISTER_VARIABLE(NODAL_WEIGHTS);


    // Register Elements
    KRATOS_REGISTER_ELEMENT("VMS2D3N",mVMS2D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("VMS3D4N",mVMS3D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("VMS2D",mVMS2D);
    KRATOS_REGISTER_ELEMENT("VMS3D",mVMS3D);
    KRATOS_REGISTER_ELEMENT("TwoFluidVMS3D",mTwoFluidVMS3D);

    KRATOS_REGISTER_ELEMENT("StationaryStokes2D",mStationaryStokes2D);
    KRATOS_REGISTER_ELEMENT("StationaryStokes3D",mStationaryStokes3D);

    KRATOS_REGISTER_ELEMENT("FractionalStep2D3N",mFractionalStep2D); //same as just below but with a standardized name
    KRATOS_REGISTER_ELEMENT("FractionalStep3D4N",mFractionalStep3D); //same as just below but with a standardized name
    KRATOS_REGISTER_ELEMENT("FractionalStep2D",mFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("FractionalStep3D",mFractionalStep3D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuous2D",mFractionalStepDiscontinuous2D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuous3D",mFractionalStepDiscontinuous3D);
	
    KRATOS_REGISTER_ELEMENT("SpalartAllmaras2D",mSpalartAllmaras2D);
    KRATOS_REGISTER_ELEMENT("SpalartAllmaras3D",mSpalartAllmaras3D);
    
    KRATOS_REGISTER_ELEMENT("DPGVMS2D",mDPGVMS2D);
    KRATOS_REGISTER_ELEMENT("DPGVMS3D",mDPGVMS3D);

    KRATOS_REGISTER_ELEMENT("BinghamVMS2D",mBinghamVMS2D);
    KRATOS_REGISTER_ELEMENT("BinghamVMS3D",mBinghamVMS3D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStep2D",mBinghamFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStep3D",mBinghamFractionalStep3D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStepDiscontinuous2D",mBinghamFractionalStepDiscontinuous2D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStepDiscontinuous3D",mBinghamFractionalStepDiscontinuous3D);

    KRATOS_REGISTER_ELEMENT("HerschelBulkleyVMS2D",mHerschelBulkleyVMS2D);
    KRATOS_REGISTER_ELEMENT("HerschelBulkleyVMS3D",mHerschelBulkleyVMS3D);
    
    KRATOS_REGISTER_ELEMENT("Stokes3D4N",mStokes3D);
    KRATOS_REGISTER_ELEMENT("StokesTwoFluid3D4N",mStokes3DTwoFluid);
    
    // Navier-Stokes symbolic elements
    KRATOS_REGISTER_ELEMENT("NavierStokes2D3N",mNavierStokes2D);
    KRATOS_REGISTER_ELEMENT("NavierStokes3D4N",mNavierStokes3D);
    KRATOS_REGISTER_ELEMENT("EmbeddedNavierStokes2D3N",mEmbeddedNavierStokes2D);
    KRATOS_REGISTER_ELEMENT("EmbeddedNavierStokes3D4N",mEmbeddedNavierStokes3D);

    // Register Conditions
    KRATOS_REGISTER_CONDITION("WallCondition2D2N",mWallCondition2D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("WallCondition3D3N",mWallCondition3D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("WallCondition2D",mWallCondition2D);
    KRATOS_REGISTER_CONDITION("WallCondition3D",mWallCondition3D);
    KRATOS_REGISTER_CONDITION("FSWernerWengleWallCondition2D",mFSWernerWengleWallCondition2D);
    KRATOS_REGISTER_CONDITION("FSWernerWengleWallCondition3D",mFSWernerWengleWallCondition3D);
    KRATOS_REGISTER_CONDITION("FSGeneralizedWallCondition2D",mFSGeneralizedWallCondition2D);
    KRATOS_REGISTER_CONDITION("FSGeneralizedWallCondition3D",mFSGeneralizedWallCondition3D);
    KRATOS_REGISTER_CONDITION("WallConditionDiscontinuous2D",mWallConditionDiscontinuous2D);
    KRATOS_REGISTER_CONDITION("WallConditionDiscontinuous3D",mWallConditionDiscontinuous3D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition2D",mMonolithicWallCondition2D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition3D",mMonolithicWallCondition3D);
    KRATOS_REGISTER_CONDITION("StokesWallCondition3D",mStokesWallCondition3D);
    KRATOS_REGISTER_CONDITION("FSPeriodicCondition2D",mFSPeriodicCondition2D);
    KRATOS_REGISTER_CONDITION("FSPeriodicCondition3D",mFSPeriodicCondition3D);
    KRATOS_REGISTER_CONDITION("FSPeriodicConditionEdge2D",mFSPeriodicConditionEdge2D);
    KRATOS_REGISTER_CONDITION("FSPeriodicConditionEdge3D",mFSPeriodicConditionEdge3D);
}


}  // namespace Kratos.


