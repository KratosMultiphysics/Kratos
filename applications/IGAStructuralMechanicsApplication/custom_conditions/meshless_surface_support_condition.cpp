//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/meshless_surface_support_condition.h"
#include "utilities/math_utils.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
MeshlessSurfaceSupportCondition::MeshlessSurfaceSupportCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : MeshlessBaseCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
MeshlessSurfaceSupportCondition::MeshlessSurfaceSupportCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : MeshlessBaseCondition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer MeshlessSurfaceSupportCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return MeshlessBaseCondition::Pointer(new MeshlessSurfaceSupportCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// Destructor
MeshlessSurfaceSupportCondition::~MeshlessSurfaceSupportCondition()
{
}
//************************************************************************************
//************************************************************************************
void MeshlessSurfaceSupportCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

    const unsigned int number_of_points = GetGeometry().size();

    //resizing the system in case it does not have the right size
    if(rLeftHandSideMatrix.size1() != number_of_points*3)
        rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points*3,number_of_points*3); //resetting LHS

    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points*3)
        rRightHandSideVector.resize(number_of_points*3,false);
    rRightHandSideVector = ZeroVector(number_of_points*3); //resetting RHS

    //Read in Data
  const Vector& ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
	//const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

	const array_1d<double, 3>& support = this->GetValue(DISPLACEMENT);

	double Penalty = this->GetValue(PENALTY_FACTOR);
  //KRATOS_WATCH(Penalty)
	const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
	//Vector Tangents = this->GetValue(TANGENTS);

	const int displacement_rotation_fix = this->GetValue(DISPLACEMENT_ROTATION_FIX);

	// Read out information of which elements are fixed
	// int cheaper to store than 4 bool
	int rot = displacement_rotation_fix / 1000;
	int dispX = (displacement_rotation_fix % 1000) / 100;
	int dispY = (displacement_rotation_fix % 100) / 10;
	int dispZ = (displacement_rotation_fix % 10) / 1;

	// DISPLACEMENTS
	Matrix Stiffness =  ZeroMatrix(3, ShapeFunctionsN.size()*3);
	for (unsigned int i = 0; i < ShapeFunctionsN.size(); i++)
	{
		if (dispX == 1)
			Stiffness(0, 3 * i) = ShapeFunctionsN[i];

		if (dispY == 1)
			Stiffness(1, 3 * i + 1) = ShapeFunctionsN[i];

		if (dispZ == 1)
			Stiffness(2, 3 * i + 2) = ShapeFunctionsN[i];
	}

  Vector TDisplacements(number_of_points*3);
	for (unsigned int i = 0; i < number_of_points; i++)
	{
		const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
		int index = 3*i;
		TDisplacements[index]     = (disp[0] - support[0]);
		TDisplacements[index + 1] = (disp[1] - support[1]);
		TDisplacements[index + 2] = (disp[2] - support[2]);
	}

  noalias(rLeftHandSideMatrix) += prod(trans(Stiffness), Stiffness);
	noalias(rRightHandSideVector) -= prod(prod(trans(Stiffness), Stiffness), TDisplacements);


  Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
  array_1d<double, 3> g1, g2, g3;
  GetBasisVectors(DN_De, g1, g2, g3);

  CrossProduct(g3, g1, g2);

  double dArea = norm_2(g3);

  KRATOS_WATCH(rLeftHandSideMatrix)
    KRATOS_WATCH(rRightHandSideVector)

    KRATOS_WATCH(support)
    KRATOS_WATCH(dArea)

  std::cout << "Meshless Surface Support Condition" << std::endl;
        KRATOS_WATCH(Penalty)
          KRATOS_WATCH(integration_weight)
	//MAPPING
	rLeftHandSideMatrix *= integration_weight * Penalty * dArea;
	rRightHandSideVector *= integration_weight * Penalty *dArea;

  KRATOS_WATCH(rLeftHandSideMatrix)
    KRATOS_WATCH(rRightHandSideVector)

  //KRATOS_WATCH(rLeftHandSideMatrix)
  KRATOS_CATCH("")
} // MeshlessSurfaceSupportCondition::MeshlessSurfaceSupportCondition


//************************************************************************************
//************************************************************************************
void MeshlessSurfaceSupportCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}


//************************************************************************************
//************************************************************************************
void MeshlessSurfaceSupportCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
	KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
	unsigned int dim = number_of_nodes * 3;

	if (rResult.size() != dim)
		rResult.resize(dim);

	for (unsigned int i = 0; i < number_of_nodes; i++)
	{
		int index = i * 3;
		rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
		rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
		rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	}

	KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void MeshlessSurfaceSupportCondition::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{

	ElementalDofList.resize(0);

	for (unsigned int i = 0; i < GetGeometry().size(); i++)
	{
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
		ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	}
}


} // Namespace Kratos


