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
//                   Michael Breitenberger
//                   Riccardo Rossi
//

// System includes

// External includes
#include "includes/define.h"
#include "custom_conditions/meshless_penalty_coupling_crack_condition.h"
#include "utilities/math_utils.h"

// Project includes
#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	MeshlessPenaltyCouplingCrackCondition::MeshlessPenaltyCouplingCrackCondition(IndexType NewId, GeometryType::Pointer pGeometry)
		: MeshlessBaseCouplingCondition(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}


	//************************************************************************************
	//************************************************************************************
	MeshlessPenaltyCouplingCrackCondition::MeshlessPenaltyCouplingCrackCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: MeshlessBaseCouplingCondition(NewId, pGeometry, pProperties)
	{
	}


	Condition::Pointer MeshlessPenaltyCouplingCrackCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return MeshlessBaseCouplingCondition::Pointer(new MeshlessPenaltyCouplingCrackCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}


	// Destructor
	MeshlessPenaltyCouplingCrackCondition::~MeshlessPenaltyCouplingCrackCondition()
	{
	}
	//************************************************************************************
	//************************************************************************************
	void MeshlessPenaltyCouplingCrackCondition::Initialize()
	{
		MeshlessBaseCouplingCondition::Initialize();

		Vector ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
		Matrix DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

		//Vector ShapeFunctionAndDerivatives = ZeroVector(3 * ShapeFunctionsN.size());
		//for (unsigned int i = 0; i < ShapeFunctionsN.size(); i++)
		//{
		//	ShapeFunctionAndDerivatives(i) = ShapeFunctionsN(i);
		//	ShapeFunctionAndDerivatives(i + ShapeFunctionsN.size()) = DN_De(i, 0);
		//	ShapeFunctionAndDerivatives(i + ShapeFunctionsN.size() * 2) = DN_De(i, 1);
		//}

		// Initialize Material
		mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
		ProcessInfo emptyProcessInfo = ProcessInfo();
		mConstitutiveLawVector->SetValue(INTEGRATION_WEIGHT, 1.0, emptyProcessInfo);
		mConstitutiveLawVector->SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, DN_De, emptyProcessInfo);
		mConstitutiveLawVector->InitializeMaterial(GetProperties(), GetGeometry(), ShapeFunctionsN);
	}
	//************************************************************************************
	//************************************************************************************
	void MeshlessPenaltyCouplingCrackCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY


		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

		Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
		Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
		Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		const unsigned int number_of_points = GetGeometry().size();

		//resizing the system in case it does not have the right size
		if (rLeftHandSideMatrix.size1() != number_of_points * 3)
			rLeftHandSideMatrix.resize(number_of_points * 3, number_of_points * 3, false);

		noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points * 3, number_of_points * 3); //resetting LHS

		//resizing as needed the RHS
		if (rRightHandSideVector.size() != number_of_points * 3)
			rRightHandSideVector.resize(number_of_points * 3, false);

		rRightHandSideVector = ZeroVector(number_of_points * 3); //resetting RHS

		const Vector& ShapeFunctionsN = this->GetValue(SHAPE_FUNCTION_VALUES);
		//KRATOS_WATCH(ShapeFunctionsN)
		const Vector& NSlave = this->GetValue(SHAPE_FUNCTION_VALUES_SLAVE);
		//KRATOS_WATCH(NSlave)
		Vector ShapeFunctions = ZeroVector(ShapeFunctionsN.size() + NSlave.size());
		for (unsigned int i = 0; i < ShapeFunctionsN.size(); i++)
		{
		  ShapeFunctions[i] = ShapeFunctionsN[i];
		}
		for (unsigned int i = 0; i < NSlave.size(); i++)
		{
		  ShapeFunctions[i + ShapeFunctionsN.size()] = -NSlave[i];
		}
		double Penalty = this->GetValue(PENALTY_FACTOR);
		const double Weighting = this->GetValue(INTEGRATION_WEIGHT);
		const Vector& localTrimTangents = this->GetValue(TANGENTS);
		const Matrix& ShapeFunctionDerivativesMaster = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

		array_1d<double, 2> localTrimTangentsMaster;
		localTrimTangentsMaster[0] = localTrimTangents[0];
		localTrimTangentsMaster[1] = localTrimTangents[1];

		Matrix DN_De_Master = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
		Matrix JMaster = ZeroMatrix(3, 2);
		JacobianElement(DN_De_Master, JMaster, true);

		array_1d<double, 3> g1;
		array_1d<double, 3> g2, g3;

		g1[0] = JMaster(0, 0);
		g2[0] = JMaster(0, 1);
		g1[1] = JMaster(1, 0);
		g2[1] = JMaster(1, 1);
		g1[2] = JMaster(2, 0);
		g2[2] = JMaster(2, 1);

		CrossProduct(g3, g1, g2);

		Vector gab = ZeroVector(3);
		Vector gab0 = ZeroVector(3);

		array_1d<double, 3> t2 = localTrimTangentsMaster(0)*g1 + localTrimTangentsMaster(1)*g2;  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		array_1d<double, 3> t1;
		CrossProduct(t1, t2, g3);
		//t2 = t2 / norm_2(t2);
		//t1 = t1 / norm_2(t1);

		array_1d<double, 3> T2 = localTrimTangentsMaster(0)*mg1_0_master + localTrimTangentsMaster(1)*mg2_0_master;  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		array_1d<double, 3> T1;
		CrossProduct(T1, T2, mg3_0_master);
		//T2 = T2 / norm_2(T2);
		//T1 = T1 / norm_2(T1);


		gab[0] = pow(t1[0], 2) + pow(t1[1], 2) + pow(t1[2], 2);
		gab[1] = pow(t2[0], 2) + pow(t2[1], 2) + pow(t2[2], 2);
		gab[2] = t1[0] * t2[0] + t1[1] * t2[1] + t1[2] * t2[2];

		gab0[0] = pow(T1[0], 2) + pow(T1[1], 2) + pow(T1[2], 2);
		gab0[1] = pow(T2[0], 2) + pow(T2[1], 2) + pow(T2[2], 2);
		gab0[2] = T1[0] * T2[0] + T1[1] * T2[1] + T1[2] * T2[2];

		Vector StrainVector = ZeroVector(3);

		StrainVector[0] = 0.5 * (gab[0] - gab0[0]);
		StrainVector[1] = 0.5 * (gab[1] - gab0[1]);
		StrainVector[2] = 0.5 * (gab[2] - gab0[2]);


		Values.SetStrainVector(StrainVector); //this is the input parameter
		Vector StressVector = ZeroVector(3);
		Values.SetStressVector(StressVector);    //this is an ouput parameter
		Matrix DMembrane = ZeroMatrix(3, 3);
		Values.SetConstitutiveMatrix(DMembrane); //this is an ouput parameter
		mConstitutiveLawVector->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);
		mConstitutiveLawVector->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);
		double damage_t = 0.0;
		mConstitutiveLawVector->GetValue(DAMAGE_T, damage_t);
		double rest_percentage = pow(1.0 - damage_t,0.1);
		double blablabla = 1.0;
		if (damage_t > 0.1)
		{
			blablabla = 0.0;
			this->SetValue(PENALTY_FACTOR, 0.0);
		}
		KRATOS_WATCH(Penalty)
			KRATOS_WATCH(blablabla)
		const int displacement_rotation_fix = this->GetValue(DISPLACEMENT_ROTATION_FIX);
    //KRATOS_WATCH(displacement_rotation_fix)

		// Read out information of which elements are fixed
		// int cheaper to store than 4 bool
		int rot = displacement_rotation_fix / 1000;
		int dispX = (displacement_rotation_fix % 1000) / 100;
		int dispY = (displacement_rotation_fix % 100) / 10;
		int dispZ = (displacement_rotation_fix % 10) / 1;

		if (rot == 1)
		{
			Vector Phi_r = ZeroVector(number_of_points * 3);
			Vector Phi_r_Lambda = ZeroVector(number_of_points * 3);
			Matrix Phi_rs = ZeroMatrix(number_of_points * 3, number_of_points * 3);
			array_1d<double, 2> Diff_Phi;
			Diff_Phi.clear();

			CaculateRotationalShapeFunctions(Phi_r, Phi_r_Lambda, Phi_rs, Diff_Phi);

			for (unsigned int i = 0; i < number_of_points * 3; i++)
			{
				for (unsigned int j = 0; j < number_of_points * 3; j++)
				{
					rLeftHandSideMatrix(i, j) = Phi_r(i)*Phi_r(j) + Diff_Phi(0)*Phi_rs(i, j);
				}
				rRightHandSideVector[i] = Diff_Phi(0)*Phi_r(i);
			}
		}

		KRATOS_WATCH(DMembrane)

		//FOR DISPLACEMENTS
		Matrix Hcomplete = ZeroMatrix(3, number_of_points * 3);
		for (unsigned int i = 0; i < number_of_points; i++)
		{
			if (dispX == 1)
				Hcomplete(0, 3 * i) = ShapeFunctions[i] * Penalty;// *DMembrane(0, 0);

			if (dispY == 1)
				Hcomplete(1, 3 * i + 1) = ShapeFunctions[i] * Penalty * blablabla;//  * DMembrane(1, 1);

			if (dispZ == 1)
				Hcomplete(2, 3 * i + 2) = ShapeFunctions[i] * Penalty;//  * DMembrane(2, 2);
		}
		//KRATOS_WATCH(ShapeFunctions)
		Vector TDisplacements(number_of_points * 3);
		for (unsigned int i = 0; i < number_of_points; i++)
		{
			const array_1d<double, 3> disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
			int index = 3 * i;
			TDisplacements[index]     = disp[0];
			TDisplacements[index + 1] = disp[1];
			TDisplacements[index + 2] = disp[2];
		}

		double JGeometrictoParameter;
		MappingGeometricToParameterMasterElement(ShapeFunctionDerivativesMaster, localTrimTangentsMaster, JGeometrictoParameter);
		
		noalias(rLeftHandSideMatrix) += prod(trans(Hcomplete), Hcomplete);
		noalias(rRightHandSideVector) += prod(prod(trans(Hcomplete), Hcomplete), TDisplacements);

		//Mapping:
		rLeftHandSideMatrix *= (Weighting * JGeometrictoParameter);
		rRightHandSideVector *= (Weighting * JGeometrictoParameter);

		KRATOS_CATCH("")
	} // MeshlessPenaltyCouplingCrackCondition::MeshlessPenaltyCouplingCrackCondition

	  //************************************************************************************
	  //************************************************************************************
	void MeshlessPenaltyCouplingCrackCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		MatrixType temp(0, 0);
		CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
	}


	//************************************************************************************
	//************************************************************************************
	void MeshlessPenaltyCouplingCrackCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void MeshlessPenaltyCouplingCrackCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
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


