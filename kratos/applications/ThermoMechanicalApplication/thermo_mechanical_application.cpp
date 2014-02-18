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
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "thermo_mechanical_application.h"
#include "includes/variables.h"


namespace Kratos
{
//Example
// 	KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//
KRATOS_CREATE_VARIABLE(int, NODE_PROPERTY_ID)
KRATOS_CREATE_VARIABLE(double,  HTC)
KRATOS_CREATE_VARIABLE(int, REF_ID)
KRATOS_CREATE_VARIABLE(double, PARTICLE_RADIUS)
KRATOS_CREATE_VARIABLE(double, POSETIVE_DISTANCE)
KRATOS_CREATE_VARIABLE(double, NAGATIVE_DISTANCE)
KRATOS_CREATE_VARIABLE(bool, IS_ESCAPED)
KRATOS_CREATE_VARIABLE(int, IS_SOLIDIFIED)
Kratos::Variable<double> SOLIDFRACTION( "SOLID FRACTION" );
Kratos::Variable<double> SOLIDIF_TIME( "SOLIDIF TIME" );
Kratos::Variable<double> SOLIDIF_MODULUS( "SOLIDIF MODULUS (cm)" );
Kratos::Variable<double> FILLTIME( "FILLTIME (s)" );
//KRATOS_CREATE_VARIABLE(double, FILLTIME );
KRATOS_CREATE_VARIABLE(double, MACRO_POROSITY )    
Kratos::Variable<double> MAX_VEL( "MAX VEL (m/s)" );
KRATOS_CREATE_VARIABLE(int, IS_GRAVITY_FILLING)

KRATOS_CREATE_VARIABLE(double, VOLUME_FRACTION ) 

KRATOS_CREATE_VARIABLE(double, KAPPA ) 
KRATOS_CREATE_VARIABLE(double, EPSILON ) 

KratosThermoMechanicalApplication::KratosThermoMechanicalApplication():
// 		mElem2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
/*		mMonolithic2DNeumann(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),*/
    mHeatContact2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
    mHeatContact3D(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
    mThermalFace2D(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
    mThermalFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mEnvironmentContact(0, Element::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
    mEnvironmentContact3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),

    mSUPGConvDiff2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mSUPGConvDiff3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mSUPGConv3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mSUPGConvDiffPhaseChange2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mSUPGConvDiffPhaseChange3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mSUPGConv2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
   // mPoisson3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))))
    
    



{}
//
void KratosThermoMechanicalApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosThermoMechanicalApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE(NODE_PROPERTY_ID)
    KRATOS_REGISTER_VARIABLE(HTC)
    KRATOS_REGISTER_VARIABLE(REF_ID)
    KRATOS_REGISTER_VARIABLE(PARTICLE_RADIUS)
    KRATOS_REGISTER_VARIABLE(POSETIVE_DISTANCE)
    KRATOS_REGISTER_VARIABLE(NAGATIVE_DISTANCE)    
    KRATOS_REGISTER_VARIABLE(IS_ESCAPED)  
    KRATOS_REGISTER_VARIABLE(IS_SOLIDIFIED) 
    KRATOS_REGISTER_VARIABLE(SOLIDFRACTION)
    KRATOS_REGISTER_VARIABLE(SOLIDIF_TIME)
    KRATOS_REGISTER_VARIABLE(SOLIDIF_MODULUS)
    KRATOS_REGISTER_VARIABLE(MACRO_POROSITY)
    KRATOS_REGISTER_VARIABLE(FILLTIME)
    KRATOS_REGISTER_VARIABLE(MAX_VEL)
    KRATOS_REGISTER_VARIABLE(IS_GRAVITY_FILLING)
    KRATOS_REGISTER_VARIABLE(VOLUME_FRACTION ) 
    KRATOS_REGISTER_VARIABLE(KAPPA ) 
    KRATOS_REGISTER_VARIABLE(EPSILON )     

// 		KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);
//
// 		KRATOS_REGISTER_ELEMENT("Elem2D", mElem2D);
// 		KRATOS_REGISTER_ELEMENT("Elemt3D", mElem3D);
    KRATOS_REGISTER_CONDITION("HeatContact2D", mHeatContact2D);
    KRATOS_REGISTER_CONDITION("HeatContact3D", mHeatContact3D);
    KRATOS_REGISTER_CONDITION("ThermalFace2D", mThermalFace2D);
    KRATOS_REGISTER_CONDITION("ThermalFace3D", mThermalFace3D);
    KRATOS_REGISTER_CONDITION("EnvironmentContact", mEnvironmentContact);
    KRATOS_REGISTER_CONDITION("EnvironmentContact3D", mEnvironmentContact3D);    
    


    KRATOS_REGISTER_ELEMENT("SUPGConvDiff2D", mSUPGConvDiff2D);
    KRATOS_REGISTER_ELEMENT("SUPGConvDiff3D", mSUPGConvDiff3D);
    KRATOS_REGISTER_ELEMENT("SUPGConv3D", mSUPGConv3D);
    KRATOS_REGISTER_ELEMENT("SUPGConvDiffPhaseChange2D", mSUPGConvDiffPhaseChange2D);
    KRATOS_REGISTER_ELEMENT("SUPGConvDiffPhaseChange3D", mSUPGConvDiffPhaseChange3D); 
    KRATOS_REGISTER_ELEMENT("SUPGConv2D", mSUPGConv2D);
    //KRATOS_REGISTER_ELEMENT("Poisson3D", mPoisson3D);    
}

}  // namespace Kratos.


