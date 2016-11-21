//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand@{KRATOS_APP_AUTHOR}
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "iga_structural_mechanics_application.h"
#include "iga_structural_mechanics_application_variables.h"

#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "includes/condition.h"  

#include "geometries/geometry.h"

//#include "custom_geometries/meshless_geometry.h"


namespace Kratos {

KratosIGAStructuralMechanicsApplication::KratosIGAStructuralMechanicsApplication():
	mMeshlessElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
	mMeshlessMembraneElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(16)))),
	mMeshlessLaplaceElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
	mMeshlessShellElement(0, Element::GeometryType::Pointer(new Geometry<Node<3>>(Element::GeometryType::PointsArrayType(1)))),
	//mContinuityConditionLagrange(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	mMeshlessLagrangeCouplingCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	//mContinuityConditionPenalty(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	//mLoadCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	//mSupportCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	mMeshlessLoadCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	mMeshlessSupportRotationCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1)))),
	mMeshlessPenaltyCouplingRotationCondition(0, Condition::GeometryType::Pointer(new Geometry<Node<3> >(Element::GeometryType::PointsArrayType(1))))
{}

void KratosIGAStructuralMechanicsApplication::Register() {
 	// calling base class register to register Kratos components
 	KratosApplication::Register();
	std::cout << "     KRATOS .___  ________    _____     " << std::endl;
	std::cout << "            |   |/  _____/   /  _  \\    " << std::endl;
	std::cout << "            |   /   \\  ___  /  /_\\  \\   " << std::endl;
	std::cout << "            |   \\    \\_\\  \\/    |    \\  " << std::endl;
	std::cout << "            |___|\\______  /\\____|__  /  " << std::endl;
	std::cout << "                        \\/         \\/ STRUCTURAL MECHANICS" << std::endl;
 	std::cout << "     Initializing KratosIGAStructuralMechanicsApplication... " << std::endl;


	KRATOS_REGISTER_VARIABLE( INTEGRATION_WEIGHT)
	KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_VALUES)
	KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_LOCAL_DERIVATIVES)
	KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
	KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
	
	// edge integration
	KRATOS_REGISTER_VARIABLE( TANGENTS)
	
	// penalty factor
	KRATOS_REGISTER_VARIABLE( PENALTY_FACTOR)

	// coupling and support
	KRATOS_REGISTER_VARIABLE(DISPLACEMENT_ROTATION_FIX)
	// for load condition
	KRATOS_REGISTER_VARIABLE(LOAD_TYPE)
	KRATOS_REGISTER_VARIABLE( DISTRIBUTED_LOAD_FACTOR)

	// Register the meshless element
	KRATOS_REGISTER_ELEMENT("MeshlessElement", mMeshlessElement)
	KRATOS_REGISTER_ELEMENT("MeshlessMembraneElement", mMeshlessMembraneElement)
	KRATOS_REGISTER_ELEMENT("MeshlessLaplaceElement", mMeshlessLaplaceElement)
	KRATOS_REGISTER_ELEMENT("MeshlessShellElement", mMeshlessShellElement)
	
	// Register meshless condition
	KRATOS_REGISTER_CONDITION("MeshlessSupportRotationCondition", mMeshlessSupportRotationCondition)
	KRATOS_REGISTER_CONDITION("MeshlessLoadCondition", mMeshlessLoadCondition)
	KRATOS_REGISTER_CONDITION("MeshlessLagrangeCouplingCondition", mMeshlessLagrangeCouplingCondition)
	KRATOS_REGISTER_CONDITION("MeshlessPenaltyCouplingRotationCondition", mMeshlessPenaltyCouplingRotationCondition)

	//// Register outdated conditions
	//KRATOS_REGISTER_CONDITION("LoadCondition", mLoadCondition)
	//KRATOS_REGISTER_CONDITION("SupportCondition", mSupportCondition)
	//KRATOS_REGISTER_CONDITION("ContinuityConditionPenalty", mContinuityConditionPenalty)
	//KRATOS_REGISTER_CONDITION("ContinuityConditionLagrange", mContinuityConditionLagrange)
}
}  // namespace Kratos.
