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
#include "geometries/line_2d_2.h"
#include "geometries/point_2d.h"
#include "fluid_rve_lagrange_multipliers_application.h"
#include "includes/variables.h"


namespace Kratos
{
	//Example
 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_VELOCITY)
 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS)
 	
 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_X_COMPONENT )
 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_Y_COMPONENT )
 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_Z_COMPONENT )
 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_PRESSURE )

 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_X_COMPONENT_ANTIPERIODIC )
 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_Y_COMPONENT_ANTIPERIODIC )
 	KRATOS_CREATE_VARIABLE(int, NODE_PAIR_Z_COMPONENT_ANTIPERIODIC )
 	KRATOS_CREATE_VARIABLE(double, LAGRANGE_MULTIPLIER_TANGENT_VELOCITY )
 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_NORMAL_VELOCITY)
 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(REYNOLDS_STRESS_2D)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_TILDA_TILDA_STRESS_2D)


//	KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosFluidRveLagrangeMultipliersApplication::KratosFluidRveLagrangeMultipliersApplication():
		mMeanVelocityLagrangeMultiplierCondition2D ( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >( Condition::GeometryType::PointsArrayType (1) ) ) ),
		mVelocityGradientsLagrangeMultiplierCondition2D ( 0, Condition::GeometryType::Pointer( new Triangle2D3<Node<3> >( Condition::GeometryType::PointsArrayType (3) ) ) ),
		mInverseTangentVelocityPeriodicCondition2D2N ( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType (2) ) ) ),
		mInverseNormalVelocityPeriodicCondition2D2N ( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType (2) ) ) )

 	{}
 	
 	void KratosFluidRveLagrangeMultipliersApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosFluidRveLagrangeMultipliersApplication... " << std::endl;
 
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_VELOCITY);
		
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS);
		
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_X_COMPONENT)
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_Y_COMPONENT)
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_Z_COMPONENT)
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_PRESSURE)		
	
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_X_COMPONENT_ANTIPERIODIC)
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_Y_COMPONENT_ANTIPERIODIC)
		KRATOS_REGISTER_VARIABLE(NODE_PAIR_Z_COMPONENT_ANTIPERIODIC)
		KRATOS_REGISTER_VARIABLE(LAGRANGE_MULTIPLIER_TANGENT_VELOCITY)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_NORMAL_VELOCITY);
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(REYNOLDS_STRESS_2D);
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_TILDA_TILDA_STRESS_2D);

 		KRATOS_REGISTER_CONDITION( "MeanVelocityLagrangeMultiplierCondition2D", mMeanVelocityLagrangeMultiplierCondition2D ) //and our condition
 		KRATOS_REGISTER_CONDITION( "VelocityGradientsLagrangeMultiplierCondition2D", mVelocityGradientsLagrangeMultiplierCondition2D ) //and our condition
 		KRATOS_REGISTER_CONDITION( "InverseTangentVelocityPeriodicCondition2D2N", mInverseTangentVelocityPeriodicCondition2D2N ) //and our condition
 		KRATOS_REGISTER_CONDITION( "InverseNormalVelocityPeriodicCondition2D2N", mInverseNormalVelocityPeriodicCondition2D2N ) //and our condition

 	}

}  // namespace Kratos.


