//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "compressible_potential_flow_application.h"
#include "compressible_potential_flow_application_variables.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos {

KratosCompressiblePotentialFlowApplication::KratosCompressiblePotentialFlowApplication():
    KratosApplication("CompressiblePotentialFlowApplication"),
    mCompressiblePotentialFlowElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mCompressiblePotentialFlowElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mPotentialWallCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2)))),
    mPotentialWallCondition3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3))))
  {}

void KratosCompressiblePotentialFlowApplication::Register() 
{
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
 	std::cout << "Initializing KratosCompressiblePotentialFlowApplication... " << std::endl;

        // Register Variables (defined in compressible_potential_flow_application_variables.h)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_INFINITY);
        KRATOS_REGISTER_VARIABLE( WAKE_ELEMENTAL_DISTANCES );
        KRATOS_REGISTER_VARIABLE( POSITIVE_POTENTIAL );
        KRATOS_REGISTER_VARIABLE( NEGATIVE_POTENTIAL );
        KRATOS_REGISTER_VARIABLE( WAKE_DISTANCE );

        //Register elements
        KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement2D3N",mCompressiblePotentialFlowElement2D3N); //this is the name the element should have according to the naming convention
        KRATOS_REGISTER_ELEMENT("CompressiblePotentialFlowElement3D4N",mCompressiblePotentialFlowElement3D4N); //this is the name the element should have according to the naming convention

        //Register conditions
        KRATOS_REGISTER_CONDITION("PotentialWallCondition2D2N",mPotentialWallCondition2D2N); //this is the name the element should have according to the naming convention
        KRATOS_REGISTER_CONDITION("PotentialWallCondition3D3N",mPotentialWallCondition3D3N); //this is the name the element should have according to the naming convention

        

        
}

}  // namespace Kratos.
