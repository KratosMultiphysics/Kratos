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
#include "geometries/line_3d_2.h"
#include "geometries/point_3d.h"
#include "blood_flow_application.h"
#include "includes/variables.h"
#include "custom_elements/artery_element.h"

namespace Kratos
{
	//Example
    KRATOS_CREATE_VARIABLE(double, FLOW)
    KRATOS_CREATE_VARIABLE(double, TERMINAL_RESISTANCE)
    KRATOS_CREATE_VARIABLE(double, PRESSURE_VENOUS)
    KRATOS_CREATE_VARIABLE(double, BETA)
    KRATOS_CREATE_VARIABLE(double, C0)
    KRATOS_CREATE_VARIABLE(double, SYSTOLIC_PRESSURE)
    KRATOS_CREATE_VARIABLE(double, DYASTOLIC_PRESSURE)
    KRATOS_CREATE_VARIABLE(double, AVERAGE_PRESSURE)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WORK )
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosBloodFlowApplication::KratosBloodFlowApplication(): 
        mArteryElement(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mArtery11Condition(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mArtery1Dto3DCondition(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mArtery3Dto1DCondition(0, Element::GeometryType::Pointer(new Point3D<Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
        //mArtery3Dto1DCondition(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mArtery12Condition(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mArteryInletCondition(0, Element::GeometryType::Pointer(new Point3D<Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
        mArteryInletConditionPressure(0, Element::GeometryType::Pointer(new Point3D<Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
        mArteryOutletCondition(0, Element::GeometryType::Pointer(new Point3D<Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
        mArteryOutletFreeCondition(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>()))))
    {}
 	
 	void KratosBloodFlowApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosBloodFlowApplication.... " << std::endl;
 
        KRATOS_REGISTER_VARIABLE( FLOW )
        KRATOS_REGISTER_VARIABLE( TERMINAL_RESISTANCE )
        KRATOS_REGISTER_VARIABLE( PRESSURE_VENOUS )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WORK )
        KRATOS_REGISTER_VARIABLE( BETA )
        KRATOS_REGISTER_VARIABLE( C0 )
        KRATOS_REGISTER_VARIABLE( SYSTOLIC_PRESSURE )
        KRATOS_REGISTER_VARIABLE( DYASTOLIC_PRESSURE )
        KRATOS_REGISTER_VARIABLE( AVERAGE_PRESSURE )


// 		KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

        KRATOS_REGISTER_ELEMENT("ArteryElement", mArteryElement);
        KRATOS_REGISTER_CONDITION("Artery1Dto3DCondition", mArtery1Dto3DCondition);
        KRATOS_REGISTER_CONDITION("Artery3Dto1DCondition", mArtery3Dto1DCondition);
        KRATOS_REGISTER_CONDITION("Artery11Condition", mArtery11Condition);
        KRATOS_REGISTER_CONDITION("Artery12Condition", mArtery12Condition);
        KRATOS_REGISTER_CONDITION("ArteryInletCondition", mArteryInletCondition);
        KRATOS_REGISTER_CONDITION("ArteryInletConditionPressure", mArteryInletConditionPressure);
        KRATOS_REGISTER_CONDITION("ArteryOutletCondition", mArteryOutletCondition);
        KRATOS_REGISTER_CONDITION("ArteryOutletFreeCondition", mArteryOutletFreeCondition);



 	}

}  // namespace Kratos.


