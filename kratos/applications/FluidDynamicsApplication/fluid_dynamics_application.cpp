//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-11-11 $
//   Revision:            $Revision: 1.0 $
//
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
    mVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mBinghamVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mBinghamVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mDynamicVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))),GeometryData::GI_GAUSS_2),
    mDynamicVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))),GeometryData::GI_GAUSS_2),
    mDynamicVMS2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))),GeometryData::GI_GAUSS_2),
    mDynamicVMS3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>()))),GeometryData::GI_GAUSS_2),
    mTwoFluidVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mStationaryStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mStationaryStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mFractionalStepDiscontinuous2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mFractionalStepDiscontinuous3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mSpalartAllmaras2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))),GeometryData::GI_GAUSS_2),
    mSpalartAllmaras3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))),GeometryData::GI_GAUSS_2),
    mWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mWernerWengleWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mWernerWengleWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mGeneralizedWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mGeneralizedWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mWallConditionDiscontinuous2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mWallConditionDiscontinuous3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mMonolithicWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mMonolithicWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mFSPeriodicCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mFSPeriodicCondition3D(0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mFSPeriodicConditionEdge2D(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mFSPeriodicConditionEdge3D(0, Element::GeometryType::Pointer( new Quadrilateral3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mLinearWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>())))),
    mDPGVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mDPGVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mFractionalStepBingham2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mFractionalStepBingham3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))), 
    mFractionalStepDiscontinuousBingham2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mFractionalStepDiscontinuousBingham3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))) 
{}

void KratosFluidDynamicsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosFluidDynamicsApplication... " << std::endl;

    // Register Variables (defined in fluid_dynamics_application_variables.h)
    KRATOS_REGISTER_VARIABLE(PATCH_INDEX);
    KRATOS_REGISTER_VARIABLE(TRACK_SUBSCALES);
    KRATOS_REGISTER_VARIABLE(TAUONE);
    KRATOS_REGISTER_VARIABLE(TAUTWO);
//    KRATOS_REGISTER_VARIABLE(Y_WALL);
    KRATOS_REGISTER_VARIABLE(SUBSCALE_PRESSURE);
    KRATOS_REGISTER_VARIABLE(C_DES);
//    KRATOS_REGISTER_VARIABLE(C_SMAGORINSKY);
    KRATOS_REGISTER_VARIABLE(SUBSCALE_VELOCITY);
    KRATOS_REGISTER_VARIABLE(VORTICITY);
    KRATOS_REGISTER_VARIABLE(COARSE_VELOCITY);
    
    KRATOS_REGISTER_VARIABLE(BINGHAM_SMOOTHER)


    // Register Elements
    KRATOS_REGISTER_ELEMENT("VMS2D",mVMS2D);
    KRATOS_REGISTER_ELEMENT("VMS3D",mVMS3D);
    KRATOS_REGISTER_ELEMENT("BinghamVMS2D",mBinghamVMS2D);
    KRATOS_REGISTER_ELEMENT("BinghamVMS3D",mBinghamVMS3D);
    KRATOS_REGISTER_ELEMENT("DynamicVMS2D",mDynamicVMS2D);
    KRATOS_REGISTER_ELEMENT("DynamicVMS3D",mDynamicVMS3D);
    KRATOS_REGISTER_ELEMENT("DynamicVMS2D4N",mDynamicVMS2D4N);
    KRATOS_REGISTER_ELEMENT("DynamicVMS3D8N",mDynamicVMS3D8N);
    KRATOS_REGISTER_ELEMENT("TwoFluidVMS3D",mTwoFluidVMS3D);

    KRATOS_REGISTER_ELEMENT("StationaryStokes2D",mStationaryStokes2D);
    KRATOS_REGISTER_ELEMENT("StationaryStokes3D",mStationaryStokes3D);

    KRATOS_REGISTER_ELEMENT("FractionalStep2D",mFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("FractionalStep3D",mFractionalStep3D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuous2D",mFractionalStepDiscontinuous2D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuous3D",mFractionalStepDiscontinuous3D);
	
    KRATOS_REGISTER_ELEMENT("SpalartAllmaras2D",mSpalartAllmaras2D);
    KRATOS_REGISTER_ELEMENT("SpalartAllmaras3D",mSpalartAllmaras3D);
    
    KRATOS_REGISTER_ELEMENT("DPGVMS2D",mDPGVMS2D);
    KRATOS_REGISTER_ELEMENT("DPGVMS3D",mDPGVMS3D);  
    
    KRATOS_REGISTER_ELEMENT("FractionalStepBingham2D",mFractionalStepBingham2D);
    KRATOS_REGISTER_ELEMENT("FractionalStepBingham3D",mFractionalStepBingham3D);  
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuousBingham2D",mFractionalStepDiscontinuousBingham2D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuousBingham3D",mFractionalStepDiscontinuousBingham3D); 
    

    // Register Conditions
    KRATOS_REGISTER_CONDITION("WallCondition2D",mWallCondition2D);
    KRATOS_REGISTER_CONDITION("WallCondition3D",mWallCondition3D);
    KRATOS_REGISTER_CONDITION("WernerWengleWallCondition2D",mWernerWengleWallCondition2D);
    KRATOS_REGISTER_CONDITION("WernerWengleWallCondition3D",mWernerWengleWallCondition3D);
    KRATOS_REGISTER_CONDITION("GeneralizedWallCondition2D",mGeneralizedWallCondition2D);
    KRATOS_REGISTER_CONDITION("GeneralizedWallCondition3D",mGeneralizedWallCondition3D);
    KRATOS_REGISTER_CONDITION("WallConditionDiscontinuous2D",mWallConditionDiscontinuous2D);
    KRATOS_REGISTER_CONDITION("WallConditionDiscontinuous3D",mWallConditionDiscontinuous3D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition2D",mMonolithicWallCondition2D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition3D",mMonolithicWallCondition3D);
    KRATOS_REGISTER_CONDITION("FSPeriodicCondition2D",mFSPeriodicCondition2D);
    KRATOS_REGISTER_CONDITION("FSPeriodicCondition3D",mFSPeriodicCondition3D);
    KRATOS_REGISTER_CONDITION("FSPeriodicConditionEdge2D",mFSPeriodicConditionEdge2D);
    KRATOS_REGISTER_CONDITION("FSPeriodicConditionEdge3D",mFSPeriodicConditionEdge3D);
    KRATOS_REGISTER_CONDITION("LinearWallCondition3D",mLinearWallCondition3D);    
}


}  // namespace Kratos.


