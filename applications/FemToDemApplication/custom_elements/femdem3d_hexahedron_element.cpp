//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velázquez
//

#include "includes/define.h"
#include "femdem3d_hexahedron_element.hpp"
#include "includes/element.h"
#include "includes/node.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{
	//***********************DEFAULT CONSTRUCTOR******************************************
	//************************************************************************************

	FemDem3DHexahedronElement::FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: FemDem3DElement(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}
	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	FemDem3DHexahedronElement::FemDem3DHexahedronElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: FemDem3DElement(NewId, pGeometry, pProperties)
	{
		//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

		// Each component == Each edge
		mNumberOfEdges = 12;
		mF_sigmas = ZeroVector(mNumberOfEdges);   // Equivalent stress
		mThresholds = ZeroVector(mNumberOfEdges); // Stress mThreshold on edge
		mDamages = ZeroVector(mNumberOfEdges); // Converged mDamage on each edge
		mNonConvergedDamages = ZeroVector(mNumberOfEdges); // mDamages on edges of "i" iteration
		mNonConvergedFsigmas = ZeroVector(mNumberOfEdges); // Equivalent stress of "i" iteration
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	FemDem3DHexahedronElement::FemDem3DHexahedronElement(FemDem3DHexahedronElement const &rOther)
		: FemDem3DElement(rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
		//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
	}

	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	FemDem3DHexahedronElement &FemDem3DHexahedronElement::operator=(FemDem3DHexahedronElement const &rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE
		FemDem3DElement::operator=(rOther);
		return *this;
	}

	//*********************************OPERATIONS*****************************************
	//************************************************************************************

	Element::Pointer FemDem3DHexahedronElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
	{
		//NEEDED TO CREATE AN ELEMENT
		return Element::Pointer(new FemDem3DHexahedronElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
	}

	//************************************CLONE*******************************************
	//************************************************************************************

	Element::Pointer FemDem3DHexahedronElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
	{

		//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
		//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE
		FemDem3DHexahedronElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
		return Element::Pointer(new FemDem3DHexahedronElement(NewElement));
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	FemDem3DHexahedronElement::~FemDem3DHexahedronElement()
	{
	}

	void FemDem3DHexahedronElement::CalculateLocalSystem(
		MatrixType &rLeftHandSideMatrix,
		VectorType &rRightHandSideVector,
		ProcessInfo &rCurrentProcessInfo)
	{
		// provisional elastic
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) / 2;

		const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
		unsigned int system_size = number_of_nodes * dimension;
		if (rLeftHandSideMatrix.size1() != system_size)
			rLeftHandSideMatrix.resize(system_size, system_size, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);

		if (rRightHandSideVector.size() != system_size)
			rRightHandSideVector.resize(system_size, false);
		noalias(rRightHandSideVector) = ZeroVector(system_size);

		Matrix DeltaPosition(number_of_nodes, dimension);
		noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

		Matrix InvJ(dimension, dimension);;
		Matrix B = ZeroMatrix(voigt_size, dimension * number_of_nodes);
		Matrix DN_DX(number_of_nodes, dimension);

		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0].resize(dimension, dimension, false);
		noalias(J[0]) = ZeroMatrix(dimension, dimension);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

		Vector damages_edges = ZeroVector(mNumberOfEdges);
		Vector average_stress_edge = ZeroVector(voigt_size);
		Vector average_strain_edge = ZeroVector(voigt_size);
		double characteristic_length;

		// Loop over edges to compute the elemental damage
		for (unsigned int edge = 0; edge < mNumberOfEdges; edge++) {
			this->CalculateAverageStressOnEdge(average_stress_edge, edge);
			this->CalculateAverageStrainOnEdge(average_strain_edge, edge);
			this->CalculateCharacteristicLength(characteristic_length, edge);
			Vector integrated_edge_stress;

			this->IntegrateStressDamageMechanics(
				integrated_edge_stress, damages_edges(edge),
				average_strain_edge, average_stress_edge, edge,
				characteristic_length);

			this->SetNonConvergedDamages(damages_edges(edge), edge);
		}
		double damage_element = this->CalculateElementalDamage(damages_edges);
		if (damage_element >= 0.999)
			damage_element = 0.999;
		this->SetNonConvergedDamages(damage_element);

		//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
		const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

		for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
			const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
			Vector N = row(Ncontainer, PointNumber);

			double detJ = 0;
			Matrix InvJ(dimension, dimension);
			noalias(InvJ) = ZeroMatrix(dimension, dimension);
			MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

			//compute cartesian derivatives for this integration point  [dN/dx_n]
			noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

			double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;
			this->CalculateDeformationMatrix(B, DN_DX);

			Matrix ConstitutiveMatrix = ZeroMatrix(voigt_size, voigt_size);
			const double E = this->GetProperties()[YOUNG_MODULUS];
			const double nu = this->GetProperties()[POISSON_RATIO];
			this->CalculateConstitutiveMatrix(ConstitutiveMatrix, E, nu);

			Vector strain_vector, predictive_stress_vector;
			this->CalculateInfinitesimalStrain(strain_vector, DN_DX);
			this->CalculateStressVector(predictive_stress_vector, ConstitutiveMatrix, strain_vector);
			const Vector integrated_stress_vector = (1.0 - damage_element) * predictive_stress_vector;

			noalias(rLeftHandSideMatrix) += prod(
				trans(B), IntegrationWeight * (1.0 - damage_element) * Matrix(prod(ConstitutiveMatrix, B)));  // LHS

            Vector VolumeForce = ZeroVector(dimension);
			VolumeForce = this->CalculateVolumeForce(VolumeForce, N);

			// RHS
			for (unsigned int i = 0; i < number_of_nodes; i++) {
				int index = dimension * i;
				for (unsigned int j = 0; j < dimension; j++) {
					rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
				}
			}

			//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
			noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(B), integrated_stress_vector);
		}
		KRATOS_CATCH("")
	}

	void FemDem3DHexahedronElement::CalculateAverageStressOnEdge(
		Vector& rEdgeStressVector,
		const unsigned int edge
		)
	{
		const Geometry<Node<3>>& nodes_element = this->GetGeometry();
		const Matrix node_edges = this->GetEdgeNodeNumbering();
		const Vector stress_node1 = nodes_element[node_edges(edge, 0)].GetValue(STRESS_VECTOR);
		const Vector stress_node2 = nodes_element[node_edges(edge, 1)].GetValue(STRESS_VECTOR);

		// Here, the average edge stress is the average between the
		// extrapolated stress on the nodes
		noalias(rEdgeStressVector) = 0.5 * (stress_node1 + stress_node2);
	}

	void FemDem3DHexahedronElement::CalculateAverageStrainOnEdge(
		Vector& rEdgeStrainVector,
		const unsigned int edge
		)
	{
		const Geometry<Node<3>>& nodes_element = this->GetGeometry();
		const Matrix node_edges = this->GetEdgeNodeNumbering();
		const Vector strain_node1 = nodes_element[node_edges(edge, 0)].GetValue(STRAIN_VECTOR);
		const Vector strain_node2 = nodes_element[node_edges(edge, 1)].GetValue(STRAIN_VECTOR);

		// Here, the average edge strain is the average between the
		// extrapolated strain on the nodes
		noalias(rEdgeStrainVector) = 0.5 * (strain_node1 + strain_node2);
	}

	void FemDem3DHexahedronElement::CalculateCharacteristicLength(double& rcharacteristic_length, const int Edge)
	{
		const Geometry<Node<3>> &nodes_element = this->GetGeometry();
		const Matrix node_edges = this->GetEdgeNodeNumbering();

		const double X1 = nodes_element[node_edges(Edge, 0)].X();
		const double X2 = nodes_element[node_edges(Edge, 1)].X();
		const double Y1 = nodes_element[node_edges(Edge, 0)].Y();
		const double Y2 = nodes_element[node_edges(Edge, 1)].Y();
		const double Z1 = nodes_element[node_edges(Edge, 0)].Z();
		const double Z2 = nodes_element[node_edges(Edge, 1)].Z();

		rcharacteristic_length = std::sqrt(
			(X1 - X2) * (X1 - X2) + (Y1 - Y2) * (Y1 - Y2) + (Z1 - Z2) * (Z1 - Z2));
	}

	// Computes the damage of the element considering different fracture modes
	double FemDem3DHexahedronElement::CalculateElementalDamage(const Vector &EdgeDamages)
	{
		// 11 modes of fracture of the hexaedron
		Vector DamageModeFracture = ZeroVector(11);
		const double one_third = 1.0 / 3.0;

		DamageModeFracture[0]  = one_third * (EdgeDamages[0] + EdgeDamages[3] + EdgeDamages[8]);
		DamageModeFracture[1]  = one_third * (EdgeDamages[0] + EdgeDamages[1] + EdgeDamages[9]);
		DamageModeFracture[2]  = one_third * (EdgeDamages[1] + EdgeDamages[2] + EdgeDamages[10]);
		DamageModeFracture[3]  = one_third * (EdgeDamages[2] + EdgeDamages[3] + EdgeDamages[11]);
		DamageModeFracture[4]  = one_third * (EdgeDamages[4] + EdgeDamages[7] + EdgeDamages[8]);
		DamageModeFracture[5]  = one_third * (EdgeDamages[4] + EdgeDamages[5] + EdgeDamages[9]);
		DamageModeFracture[6]  = one_third * (EdgeDamages[6] + EdgeDamages[5] + EdgeDamages[10]);
		DamageModeFracture[7]  = one_third * (EdgeDamages[6] + EdgeDamages[7] + EdgeDamages[11]);
		DamageModeFracture[8]  = 0.25 * (EdgeDamages[8] + EdgeDamages[9] + EdgeDamages[10] + EdgeDamages[11]);
		DamageModeFracture[9]  = 0.25 * (EdgeDamages[0] + EdgeDamages[2] + EdgeDamages[6] + EdgeDamages[4]);
		DamageModeFracture[10] = 0.25 * (EdgeDamages[1] + EdgeDamages[3] + EdgeDamages[7] + EdgeDamages[5]);

		return this->GetMaxValue(DamageModeFracture);
	}


	void FemDem3DHexahedronElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
	{
		//Loop over edges
		for (unsigned int cont = 0; cont < this->GetNumberOfEdges(); cont++) {
			this->SetConvergedDamages(this->GetNonConvergedDamages(cont), cont);
			this->SetConvergedEquivalentStress(this->GetNonConvergedEquivalentStress(cont), cont);
			const double current_equivalent_stress = this->GetConvergedEquivalentStress(cont);
			
			if (current_equivalent_stress > this->GetThreshold(cont)) {
				this->SetThreshold(current_equivalent_stress, cont);
			}

		} // End Loop over edges

		const double damage_element = this->GetNonConvergedDamage();
		this->SetConvergedDamage(damage_element);

		if (damage_element >= 0.98) {
			this->Set(ACTIVE, false);
		}
		this->ResetNonConvergedVars();
		this->SetValue(DAMAGE_ELEMENT, damage_element);
	}


	// Double values
	void FemDem3DHexahedronElement::GetValueOnIntegrationPoints(
		const Variable<double> &rVariable,
		std::vector<double> &rValues,
		const ProcessInfo &rCurrentProcessInfo)
	{
		if (rVariable == DAMAGE_ELEMENT || rVariable == IS_DAMAGED || rVariable == STRESS_THRESHOLD) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}

	// Vector Values
	void FemDem3DHexahedronElement::GetValueOnIntegrationPoints(
		const Variable<Vector> &rVariable,
		std::vector<Vector> &rValues,
		const ProcessInfo &rCurrentProcessInfo)
	{
		if (rVariable == STRAIN_VECTOR) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		} else if (rVariable == STRESS_VECTOR) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		} else if (rVariable == STRESS_VECTOR_INTEGRATED) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}

	// Tensor variables
	void FemDem3DHexahedronElement::GetValueOnIntegrationPoints(
		const Variable<Matrix> &rVariable,
		std::vector<Matrix> &rValues,
		const ProcessInfo &rCurrentProcessInfo)
	{
		if (rVariable == STRAIN_TENSOR) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		} else if (rVariable == STRESS_TENSOR) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		} else if (rVariable == STRESS_TENSOR_INTEGRATED) {
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}

	// double variables
	void FemDem3DHexahedronElement::CalculateOnIntegrationPoints(
		const Variable<double> &rVariable,
		std::vector<double> &rOutput,
		const ProcessInfo &rCurrentProcessInfo)
	{
		ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
		const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

		if (rOutput.size() != integration_points.size())
			rOutput.resize(integration_points.size());

		if (rVariable == DAMAGE_ELEMENT) {
			for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
				rOutput[PointNumber] = this->GetDamage();
			}
		} else if (rVariable == IS_DAMAGED) {
			for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
				rOutput[PointNumber] = this->GetValue(IS_DAMAGED);
			}
		} else if (rVariable == STRESS_THRESHOLD) {
			for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
				rOutput[PointNumber] = this->GetValue(STRESS_THRESHOLD);
			}
		}
	}

	// 	VECTOR VARIABLES
	void FemDem3DHexahedronElement::CalculateOnIntegrationPoints(
		const Variable<Vector> &rVariable,
		std::vector<Vector> &rOutput,
		const ProcessInfo &rCurrentProcessInfo)
	{
		const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
		if (rOutput.size() != integration_points.size())
			rOutput.resize(integration_points.size());

		if (rVariable == STRESS_VECTOR) {
			//create and initialize element variables:
			ElementDataType Variables;
			this->InitializeElementData(Variables,rCurrentProcessInfo);

			//create constitutive law parameters:
			ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

			//set constitutive law flags:
			Flags &ConstitutiveLawOptions=Values.GetOptions();

			ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

			//reading integration points
			for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ ) {
				//compute element kinematic variables B, F, DN_DX ...
				this->CalculateKinematics(Variables,PointNumber);

				//set general variables to constitutivelaw parameters
				this->SetElementData(Variables,Values,PointNumber);

				//call the constitutive law to update material variables
				mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

				if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
						rOutput[PointNumber].resize( Variables.StressVector.size(), false );
				rOutput[PointNumber] = Variables.StressVector;
			}
		} else if (rVariable == STRAIN_VECTOR) {
			//create and initialize element variables:
			ElementDataType Variables;
			this->InitializeElementData(Variables,rCurrentProcessInfo);

			//create constitutive law parameters:
			ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

			//set constitutive law flags:
			Flags &ConstitutiveLawOptions=Values.GetOptions();

			ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

			//reading integration points
			for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ ) {
				//compute element kinematic variables B, F, DN_DX ...
				this->CalculateKinematics(Variables,PointNumber);

				//set general variables to constitutivelaw parameters
				this->SetElementData(Variables,Values,PointNumber);

				if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
						rOutput[PointNumber].resize( Variables.StrainVector.size(), false );
				rOutput[PointNumber] = Variables.StrainVector;
			}
		} else if (rVariable == STRESS_VECTOR_INTEGRATED) {
			//create and initialize element variables:
			ElementDataType Variables;
			this->InitializeElementData(Variables,rCurrentProcessInfo);

			//create constitutive law parameters:
			ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

			//set constitutive law flags:
			Flags &ConstitutiveLawOptions=Values.GetOptions();

			ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

			//reading integration points
			for (unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++) {
				//compute element kinematic variables B, F, DN_DX ...
				this->CalculateKinematics(Variables,PointNumber);

				//set general variables to constitutivelaw parameters
				this->SetElementData(Variables,Values,PointNumber);

				//call the constitutive law to update material variables
				mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

				if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
						rOutput[PointNumber].resize( Variables.StressVector.size(), false );

				const double damage = this->GetDamage();
				rOutput[PointNumber] = (1.0 - damage) * Variables.StressVector;
			}
		}
	}

	// 	TENSOR VARIABLES
	void FemDem3DHexahedronElement::CalculateOnIntegrationPoints(
		const Variable<Matrix> &rVariable,
		std::vector<Matrix> &rOutput,
		const ProcessInfo &rCurrentProcessInfo)
	{
		const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
		if (rOutput.size() != integration_points.size())
			rOutput.resize(integration_points.size());
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		if (rVariable == STRESS_TENSOR) {
			std::vector<Vector> StressVector;
			this->CalculateOnIntegrationPoints( STRESS_VECTOR, StressVector, rCurrentProcessInfo );

			//loop integration points
			for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ ) {
				if ( rOutput[PointNumber].size2() != dimension )
					rOutput[PointNumber].resize( dimension, dimension, false );
				rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
			}

		} else if (rVariable == STRAIN_TENSOR) {
			std::vector<Vector> StrainVector;
			this->CalculateOnIntegrationPoints( STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

			//loop integration points
			for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ ) {
				if ( rOutput[PointNumber].size2() != dimension )
					rOutput[PointNumber].resize( dimension, dimension, false );
				rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
			}

		} else if (rVariable == STRESS_TENSOR_INTEGRATED) {
			std::vector<Vector> IntegratedStressVector;
			this->CalculateOnIntegrationPoints( STRESS_VECTOR_INTEGRATED, IntegratedStressVector, rCurrentProcessInfo );

			//loop integration points
			for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ ) {
				if ( rOutput[PointNumber].size2() != dimension )
					rOutput[PointNumber].resize( dimension, dimension, false );
				rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(IntegratedStressVector[PointNumber]);
			}

		}
	}

} // namespace Kratos