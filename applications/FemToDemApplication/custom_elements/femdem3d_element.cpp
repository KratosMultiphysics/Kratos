//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//


#include "femdem3d_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"


namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

FemDem3DElement::FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: SmallDisplacementElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

FemDem3DElement::FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: SmallDisplacementElement(NewId, pGeometry, pProperties)
{
	//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

	// Each component == Each edge
	mNumberOfEdges = 6;
	mNonConvergedThresholds = ZeroVector(mNumberOfEdges);   // Equivalent stress
	mThresholds = ZeroVector(mNumberOfEdges); // Stress mThreshold on edge
	mDamages = ZeroVector(mNumberOfEdges); // Converged mDamage on each edge
	mNonConvergedDamages = ZeroVector(mNumberOfEdges); // mDamages on edges of "i" iteration
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

FemDem3DElement::FemDem3DElement(FemDem3DElement const &rOther)
	: SmallDisplacementElement(rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
	//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

FemDem3DElement &FemDem3DElement::operator=(FemDem3DElement const &rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

	SmallDisplacementElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer FemDem3DElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	//NEEDED TO CREATE AN ELEMENT
	return Element::Pointer(new FemDem3DElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer FemDem3DElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

	//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
	//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

	FemDem3DElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

	return Element::Pointer(new FemDem3DElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

FemDem3DElement::~FemDem3DElement()
{
}

void FemDem3DElement::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
	this->ComputeEdgeNeighbours(rCurrentProcessInfo);
	this->InitializeInternalVariablesAfterMapping();
}

void FemDem3DElement::InitializeInternalVariablesAfterMapping()
{
	// After the mapping, the thresholds of the edges ( are equal to 0.0) are imposed equal to the IP threshold
	const double element_threhsold = mThreshold;
	if (mThresholds[0] + mThresholds[1] + mThresholds[2] + mThresholds[3] +
	    mThresholds[4] + mThresholds[5] < std::numeric_limits<double>::epsilon()) {
		for (unsigned int edge = 0; edge < this->GetNumberOfEdges(); edge++) {
			mThresholds[edge] = element_threhsold;
		}
	}

	// IDEM with the edge damages
	const double damage_element = mDamage;
	if (mDamages[0] + mDamages[1] + mDamages[2] + mDamages[3] + mDamages[4] +
		mDamages[5] < std::numeric_limits<double>::epsilon()) {
		for (unsigned int edge = 0; edge < this->GetNumberOfEdges(); edge++) {
			mDamages[edge] = damage_element;
		}
	}
}

void FemDem3DElement::ComputeEdgeNeighbours(ProcessInfo &rCurrentProcessInfo)
{
	std::vector<std::vector<Element *>> EdgeNeighboursContainer;
	Geometry<Node<3>> &NodesCurrentElement = this->GetGeometry();

	Node<3> &pNode0 = NodesCurrentElement[0];
	Node<3> &pNode1 = NodesCurrentElement[1];
	Node<3> &pNode2 = NodesCurrentElement[2];
	Node<3> &pNode3 = NodesCurrentElement[3];

	// Neighbour elements of each node of the current element
	WeakPointerVector<Element> &NeighNode0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
	WeakPointerVector<Element> &NeighNode1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
	WeakPointerVector<Element> &NeighNode2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
	WeakPointerVector<Element> &NeighNode3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);

	// Nodal neighbours container
	std::vector<WeakPointerVector<Element>> NodalNeighbours;
	NodalNeighbours.push_back(NeighNode0);
	NodalNeighbours.push_back(NeighNode1);
	NodalNeighbours.push_back(NeighNode2);
	NodalNeighbours.push_back(NeighNode3);

	// Aux indexes
	Matrix nodes_indexes = ZeroMatrix(6, 2);
	this->SetNodeIndexes(nodes_indexes);

	// Loop over EDGES to assign the elements that share that edge -> Fill mEdgeNeighboursContainer
	for (unsigned int edge = 0; edge < 6; edge++) {
		const int NodeIndex1 = nodes_indexes(edge, 0);
		const int NodeIndex2 = nodes_indexes(edge, 1);

		// Neigh elements of local node 1 and 2  //
		WeakPointerVector<Element> &neigh_of_node_1 = NodalNeighbours[NodeIndex1];
		WeakPointerVector<Element> &neigh_of_node_2 = NodalNeighbours[NodeIndex2];

		const int NodeId1 = NodesCurrentElement[NodeIndex1].Id();
		const int NodeId2 = NodesCurrentElement[NodeIndex2].Id();

		std::vector<Element *> edge_shared_elements_node_1;
		// Loop over neigh elements of the node 1
		for (unsigned int neigh_elem = 0; neigh_elem < neigh_of_node_1.size(); neigh_elem++) {
			// Nodes of the neigh element
			Geometry<Node<3>> &NodesNeighElem = neigh_of_node_1[neigh_elem].GetGeometry();

			// Loop over the nodes of the neigh element
			for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
				const int NeighElementNodeId = NodesNeighElem[neigh_elem_node].Id();

				if (NeighElementNodeId == NodeId2 && this->Id() != neigh_of_node_1[neigh_elem].Id()) {
					edge_shared_elements_node_1.push_back(&neigh_of_node_1[neigh_elem]); // ( [] returns an Element object!!)
				}
			}
		}

		std::vector<Element *> edge_shared_elements_node_2;
		// Loop over neigh elements of the node 2
		for (unsigned int neigh_elem = 0; neigh_elem < neigh_of_node_2.size(); neigh_elem++) {
			// Nodes of the neigh element
			Geometry<Node<3>> &NodesNeighElem = neigh_of_node_2[neigh_elem].GetGeometry();

			// Loop over the nodes of the neigh element
			for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
				const int NeighElementNodeId = NodesNeighElem[neigh_elem_node].Id();

				if (NeighElementNodeId == NodeId1 && this->Id() != neigh_of_node_2[neigh_elem].Id()) {
					edge_shared_elements_node_2.push_back(&neigh_of_node_2[neigh_elem]);
				}
			}
		}
		// Let's create the vector of neighbour elements for this edge
		std::vector<Element *> edge_shared_elements = edge_shared_elements_node_1;
		// Add the neigh elements from the node 2
		for (unsigned int i = 0; i < edge_shared_elements_node_2.size(); i++) {
			int aux = 0;

			for (unsigned int j = 0; j < edge_shared_elements.size(); j++) {
				if (edge_shared_elements_node_2[i]->Id() == edge_shared_elements[j]->Id())
					aux++;
			}
			if (aux == 0)
				edge_shared_elements.push_back(edge_shared_elements_node_2[i]);
		}
		EdgeNeighboursContainer.push_back(edge_shared_elements);
	} // End loop edges

	// Storages the information inside the element
	this->SaveEdgeNeighboursContainer(EdgeNeighboursContainer);

} // End finding edge neighbour elements

void FemDem3DElement::UpdateDataBase()
{	
	for (unsigned int edge = 0; edge < mNumberOfEdges; edge++) {
		mDamages[edge] = mNonConvergedDamages[edge];
		mThresholds[edge] = mNonConvergedThresholds[edge];
	}

	double converged_damage, converged_threshold;
	converged_damage = this->CalculateElementalDamage(mDamages);
	if (converged_damage > mDamage) mDamage = converged_damage;
	converged_threshold = this->CalculateElementalDamage(mThresholds);
	if (converged_threshold > mThreshold) mThreshold = converged_threshold;
}

void FemDem3DElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
	this->UpdateDataBase();

	if (mDamage >= 0.98) {
		this->Set(ACTIVE, false);
	}
}

void FemDem3DElement::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
	//1.-Initialize sizes for the system components:
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	unsigned int voigt_size = dimension * (dimension + 1) / 2;

	Vector StrainVector(voigt_size);
	noalias(StrainVector) = ZeroVector(voigt_size);
	Vector StressVector(voigt_size);
	noalias(StressVector) = ZeroVector(voigt_size);
	Matrix ConstitutiveMatrix(voigt_size, voigt_size);
	noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
	Matrix B(voigt_size, dimension * number_of_nodes);
	noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
	Matrix DN_DX(number_of_nodes, dimension);
	noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

	//deffault values for the infinitessimal theory
	double detF = 1;
	Matrix F(dimension, dimension);
	noalias(F) = identity_matrix<double>(dimension);

	//3.-Calculate elemental system:
	//reading integration points
	const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

	//get the shape functions [N] (for the order of the default integration method)
	const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

	//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
	const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

	//calculate delta position (here coincides with the current displacement)
	Matrix DeltaPosition(number_of_nodes, dimension);
	noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

	//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
	GeometryType::JacobiansType J;
	J.resize(1, false);
	J[0].resize(dimension, dimension, false);
	noalias(J[0]) = ZeroMatrix(dimension, dimension);
	J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

	// Loop Over Integration Points
	for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
		Matrix InvJ(dimension, dimension);
		noalias(InvJ) = ZeroMatrix(dimension, dimension);
		double detJ = 0;
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

		KRATOS_ERROR_IF(detJ < 0) << " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 )" << std::endl;
		
		//compute cartesian derivatives for this integration point  [dN/dx_n]
		noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

		//set shape functions for this integration point
		Vector N = row(Ncontainer, PointNumber);

		//b.-compute infinitessimal strain
		this->CalculateInfinitesimalStrain(StrainVector, DN_DX);
		this->SetValue(STRAIN_VECTOR, StrainVector);

		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

		//set constitutive law variables: (it passes only references to this local variables)
		Values.SetStrainVector(StrainVector);
		Values.SetStressVector(StressVector);
		Values.SetConstitutiveMatrix(ConstitutiveMatrix);
		Values.SetShapeFunctionsDerivatives(DN_DX);
		Values.SetShapeFunctionsValues(N);
		//values to be set:
		Values.SetDeterminantF(detF);
		Values.SetDeformationGradientF(F);

		//set constitutive law flags:
		Flags& ConstitutiveLawOptions = Values.GetOptions();

		//compute stress and constitutive matrix
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		//CALL THE CONSTITUTIVE LAW (for this integration point)
		//(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
		mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
		this->SetValue(STRESS_VECTOR, Values.GetStressVector());
	}
}

void FemDem3DElement::CalculateLocalSystem(
	MatrixType &rLeftHandSideMatrix,
	VectorType &rRightHandSideVector,
	ProcessInfo &rCurrentProcessInfo)
{
	//*****************************
	KRATOS_TRY

	//1.-Initialize sizes for the system components:
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	const unsigned int voigt_size = dimension * (dimension + 1) / 2;
	const unsigned int system_size = number_of_nodes * dimension;

	if (rLeftHandSideMatrix.size1() != system_size)
		rLeftHandSideMatrix.resize(system_size, system_size, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);

	if (rRightHandSideVector.size() != system_size)
		rRightHandSideVector.resize(system_size, false);
	noalias(rRightHandSideVector) = ZeroVector(system_size);

	Vector StrainVector(voigt_size);
	noalias(StrainVector) = ZeroVector(voigt_size);
	Vector StressVector(voigt_size);
	noalias(StressVector) = ZeroVector(voigt_size);
	Matrix ConstitutiveMatrix(voigt_size, voigt_size);
	noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
	Matrix B(voigt_size, dimension * number_of_nodes);
	noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
	Matrix DN_DX(number_of_nodes, dimension);
	noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

	//deffault values for the infinitessimal theory
	double detF = 1;
	Matrix F(dimension, dimension);
	noalias(F) = identity_matrix<double>(dimension);

	//3.-Calculate elemental system:

	//reading integration points
	const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

	//get the shape functions [N] (for the order of the default integration method)
	const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

	//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
	const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

	//calculate delta position (here coincides with the current displacement)
	Matrix DeltaPosition(number_of_nodes, dimension);
	noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

	//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
	GeometryType::JacobiansType J;
	J.resize(1, false);
	J[0].resize(dimension, dimension, false);
	noalias(J[0]) = ZeroMatrix(dimension, dimension);
	J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

	for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
		Matrix InvJ(dimension, dimension);
		noalias(InvJ) = ZeroMatrix(dimension, dimension);
		double detJ = 0;
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

		double integration_weight = integration_points[PointNumber].Weight() * detJ;

		KRATOS_ERROR_IF(detJ < 0) << " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0" << std::endl;
		//compute cartesian derivatives for this integration point  [dN/dx_n]
		noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

		//set shape functions for this integration point
		Vector N = row(Ncontainer, PointNumber);

		//b.-compute infinitessimal strain
		this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

		//set constitutive law variables: (it passes only references to this local variables)
		Values.SetStrainVector(StrainVector);
		Values.SetStressVector(StressVector);
		Values.SetConstitutiveMatrix(ConstitutiveMatrix);
		Values.SetShapeFunctionsDerivatives(DN_DX);
		Values.SetShapeFunctionsValues(N);
		//values to be set:
		Values.SetDeterminantF(detF);
		Values.SetDeformationGradientF(F);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();

		//compute stress and constitutive matrix
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		//CALL THE CONSTITUTIVE LAW (for this integration point)
		//(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
		mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
		const Vector& r_characteristic_lengths = this->CalculateCharacteristicLengths();
		bool is_damaging = false;

		// Loop over edges of the element
		for (unsigned int edge = 0; edge < mNumberOfEdges; edge++) {
			std::vector<Element*> EdgeNeighbours = this->GetEdgeNeighbourElements(edge);
			Vector average_stress_edge, average_strain_edge;

			this->CalculateAverageStressOnEdge(average_stress_edge, EdgeNeighbours);
			this->CalculateAverageStrainOnEdge(average_strain_edge, EdgeNeighbours);

			double damage_edge = mDamages[edge];
			double threshold = mThresholds[edge];

			this->IntegrateStressDamageMechanics(threshold, 
												 damage_edge,
												 average_strain_edge, 
												 average_stress_edge, 
												 edge, 
												 r_characteristic_lengths[edge],
												 is_damaging);
			mNonConvergedDamages[edge] = damage_edge;
			mNonConvergedThresholds[edge] = threshold;
		} // End loop over edges

		const double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);
		const Vector& stress_vector = this->GetValue(STRESS_VECTOR);
		const Vector& integrated_stress_vector = (1.0 - damage_element) * stress_vector;

		Matrix constitutive_matrix = Values.GetConstitutiveMatrix();
		this->CalculateDeformationMatrix(B, DN_DX);

		Matrix tangent_tensor;
		if (is_damaging == true && std::abs(StrainVector[0] + StrainVector[1] + StrainVector[2]
			+ StrainVector[3]+ StrainVector[4] + StrainVector[5]) > tolerance) {

			this->CalculateTangentTensor(tangent_tensor, StrainVector, integrated_stress_vector, constitutive_matrix);
			noalias(rLeftHandSideMatrix) += prod(trans(B), integration_weight * Matrix(prod(tangent_tensor, B)));
		} else {
			noalias(rLeftHandSideMatrix) += prod(trans(B), integration_weight * (1.0 - damage_element) * Matrix(prod(constitutive_matrix, B)));
		}

		Vector VolumeForce = ZeroVector(dimension);
		VolumeForce = this->CalculateVolumeForce(VolumeForce, N);
		for (unsigned int i = 0; i < number_of_nodes; i++) {
			const int index = dimension * i;
			for (unsigned int j = 0; j < dimension; j++) {
				rRightHandSideVector[index + j] += integration_weight * N[i] * VolumeForce[j];
			}
		}
		//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
		noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), integrated_stress_vector);
	}
	KRATOS_CATCH("")
}

void FemDem3DElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
	//1.-Initialize sizes for the system components:
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	const unsigned int voigt_size = dimension * (dimension + 1) / 2;
	const unsigned int system_size = number_of_nodes * dimension;

	if (rLeftHandSideMatrix.size1() != system_size)
		rLeftHandSideMatrix.resize(system_size, system_size, false);
	noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size);

	Vector StrainVector(voigt_size);
	noalias(StrainVector) = ZeroVector(voigt_size);
	Vector StressVector(voigt_size);
	noalias(StressVector) = ZeroVector(voigt_size);
	Matrix ConstitutiveMatrix(voigt_size, voigt_size);
	noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
	Matrix B(voigt_size, dimension * number_of_nodes);
	noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
	Matrix DN_DX(number_of_nodes, dimension);
	noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

	//default values for the infinitessimal theory
	double detF = 1;
	Matrix F(dimension, dimension);
	noalias(F) = identity_matrix<double>(dimension);

	//3.-Calculate elemental system:

	//reading integration points
	const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

	//get the shape functions [N] (for the order of the default integration method)
	const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

	//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
	const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

	//calculate delta position (here coincides with the current displacement)
	Matrix DeltaPosition(number_of_nodes, dimension);
	noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

	//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
	GeometryType::JacobiansType J;
	J.resize(1, false);
	J[0].resize(dimension, dimension, false);
	noalias(J[0]) = ZeroMatrix(dimension, dimension);
	J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

	for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {
		Matrix InvJ(dimension, dimension);
		noalias(InvJ) = ZeroMatrix(dimension, dimension);
		double detJ = 0;
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

		double integration_weight = integration_points[PointNumber].Weight() * detJ;
		KRATOS_ERROR_IF(detJ < 0) << " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 " << std::endl;

		//compute cartesian derivatives for this integration point  [dN/dx_n]
		noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

		//set shape functions for this integration point
		Vector N = row(Ncontainer, PointNumber);

		//b.-compute infinitessimal strain
		this->CalculateInfinitesimalStrain(StrainVector, DN_DX);

		ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

		//set constitutive law variables: (it passes only references to this local variables)
		Values.SetStrainVector(StrainVector);
		Values.SetStressVector(StressVector);
		Values.SetConstitutiveMatrix(ConstitutiveMatrix);
		Values.SetShapeFunctionsDerivatives(DN_DX);
		Values.SetShapeFunctionsValues(N);
		//values to be set:
		Values.SetDeterminantF(detF);
		Values.SetDeformationGradientF(F);

		//set constitutive law flags:
		Flags &ConstitutiveLawOptions = Values.GetOptions();

		//compute stress and constitutive matrix
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
		ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

		//CALL THE CONSTITUTIVE LAW (for this integration point)
		//(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
		mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
		Matrix constitutive_matrix = Values.GetConstitutiveMatrix();
		this->CalculateDeformationMatrix(B, DN_DX);
		const double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);
		noalias(rLeftHandSideMatrix) += prod(trans(B), integration_weight * (1.0 - damage_element) * Matrix(prod(constitutive_matrix, B)));
	}
}

void FemDem3DElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	const unsigned int voigt_size = dimension * (dimension + 1) / 2;
	const unsigned int system_size = number_of_nodes * dimension;

	if (rRightHandSideVector.size() != system_size)
		rRightHandSideVector.resize(system_size, false);
	noalias(rRightHandSideVector) = ZeroVector(system_size);

	Matrix B(voigt_size, dimension * number_of_nodes);
	noalias(B) = ZeroMatrix(voigt_size, dimension * number_of_nodes);
	Matrix DN_DX(number_of_nodes, dimension);
	noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);

	const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
	const GeometryType::IntegrationPointsArrayType &integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
	const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

	Matrix DeltaPosition(number_of_nodes, dimension);
	noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

	GeometryType::JacobiansType J;
	J.resize(1, false);
	J[0].resize(dimension, dimension, false);
	noalias(J[0]) = ZeroMatrix(dimension, dimension);
	J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

	for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++) {

		Matrix InvJ(dimension, dimension);
		noalias(InvJ) = ZeroMatrix(dimension, dimension);
		double detJ = 0;
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);
		noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

		Vector N = row(Ncontainer, PointNumber);
		double integration_weight = integration_points[PointNumber].Weight() * detJ;
		Vector VolumeForce = ZeroVector(dimension);
		VolumeForce = this->CalculateVolumeForce(VolumeForce, N);
		for (unsigned int i = 0; i < number_of_nodes; i++) {
			const int index = dimension * i;
			for (unsigned int j = 0; j < dimension; j++) {
				rRightHandSideVector[index + j] += integration_weight * N[i] * VolumeForce[j];
			}
		}
		
		const double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);
		const Vector& stress_vector = this->GetValue(STRESS_VECTOR);
		const Vector& integrated_stress_vector = (1.0 - damage_element) * stress_vector;

		this->CalculateDeformationMatrix(B, DN_DX);
		noalias(rRightHandSideVector) -= integration_weight * prod(trans(B), integrated_stress_vector);
	}
}



void FemDem3DElement::CalculateDeformationMatrix(Matrix &rB, const Matrix &rDN_DX)
{
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	for (unsigned int i = 0; i < number_of_nodes; i++) {
		unsigned int index = 3 * i;

		rB(0, index + 0) = rDN_DX(i, 0);
		rB(1, index + 1) = rDN_DX(i, 1);
		rB(2, index + 2) = rDN_DX(i, 2);

		rB(3, index + 0) = rDN_DX(i, 1);
		rB(3, index + 1) = rDN_DX(i, 0);

		rB(4, index + 1) = rDN_DX(i, 2);
		rB(4, index + 2) = rDN_DX(i, 1);

		rB(5, index + 0) = rDN_DX(i, 2);
		rB(5, index + 2) = rDN_DX(i, 0);
	}
}

void FemDem3DElement::CalculateConstitutiveMatrix(
	Matrix &rConstitutiveMatrix,
	const double rYoungModulus,
	const double rPoissonCoefficient)
{
	rConstitutiveMatrix.clear();

	// 3D linear elastic constitutive matrix
	rConstitutiveMatrix(0, 0) = (rYoungModulus * (1.0 - rPoissonCoefficient) / ((1.0 + rPoissonCoefficient) * (1.0 - 2.0 * rPoissonCoefficient)));
	rConstitutiveMatrix(1, 1) = rConstitutiveMatrix(0, 0);
	rConstitutiveMatrix(2, 2) = rConstitutiveMatrix(0, 0);

	rConstitutiveMatrix(3, 3) = rConstitutiveMatrix(0, 0) * (1.0 - 2.0 * rPoissonCoefficient) / (2.0 * (1.0 - rPoissonCoefficient));
	rConstitutiveMatrix(4, 4) = rConstitutiveMatrix(3, 3);
	rConstitutiveMatrix(5, 5) = rConstitutiveMatrix(3, 3);

	rConstitutiveMatrix(0, 1) = rConstitutiveMatrix(0, 0) * rPoissonCoefficient / (1.0 - rPoissonCoefficient);
	rConstitutiveMatrix(1, 0) = rConstitutiveMatrix(0, 1);

	rConstitutiveMatrix(0, 2) = rConstitutiveMatrix(0, 1);
	rConstitutiveMatrix(2, 0) = rConstitutiveMatrix(0, 1);

	rConstitutiveMatrix(1, 2) = rConstitutiveMatrix(0, 1);
	rConstitutiveMatrix(2, 1) = rConstitutiveMatrix(0, 1);
}

void FemDem3DElement::CalculateDN_DX(Matrix &rDN_DX, int PointNumber)
{
	const unsigned int number_of_nodes = GetGeometry().size();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

	//get the shape functions parent coordinates derivative [dN/d�] (for the order of the default integration method)
	const GeometryType::ShapeFunctionsGradientsType &DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
	
	//calculate delta position (here coincides with the current displacement)
	Matrix DeltaPosition = ZeroMatrix(number_of_nodes, dimension);
	DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

	//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
	GeometryType::JacobiansType J;
	J.resize(1, false);
	J[0] = ZeroMatrix(1, 1);
	J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

	//calculating the inverse of the jacobian for this integration point[d�/dx_n]
	Matrix InvJ = ZeroMatrix(dimension, dimension);
	double detJ;
	MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

	KRATOS_ERROR_IF(detJ < 0) << " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = "<< detJ << std::endl;

	//compute cartesian derivatives for this integration point  [dN/dx_n]
	rDN_DX = prod(DN_De[PointNumber], InvJ);
}

void FemDem3DElement::CalculateInfinitesimalStrain(Vector &rStrainVector, const Matrix &rDN_DX)
{
	KRATOS_TRY

	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

	Matrix H = zero_matrix<double>(dimension); //[dU/dx_n]

	for (unsigned int i = 0; i < number_of_nodes; i++) {
		array_1d<double, 3> &Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

		H(0, 0) += Displacement[0] * rDN_DX(i, 0);
		H(0, 1) += Displacement[0] * rDN_DX(i, 1);
		H(0, 2) += Displacement[0] * rDN_DX(i, 2);
		H(1, 0) += Displacement[1] * rDN_DX(i, 0);
		H(1, 1) += Displacement[1] * rDN_DX(i, 1);
		H(1, 2) += Displacement[1] * rDN_DX(i, 2);
		H(2, 0) += Displacement[2] * rDN_DX(i, 0);
		H(2, 1) += Displacement[2] * rDN_DX(i, 1);
		H(2, 2) += Displacement[2] * rDN_DX(i, 2);
	}

	//Infinitesimal Strain Calculation
	if (rStrainVector.size() != 6)
		rStrainVector.resize(6, false);

	rStrainVector[0] = H(0, 0);
	rStrainVector[1] = H(1, 1);
	rStrainVector[2] = H(2, 2);
	rStrainVector[3] = (H(0, 1) + H(1, 0)); // xy
	rStrainVector[4] = (H(1, 2) + H(2, 1)); // yz
	rStrainVector[5] = (H(0, 2) + H(2, 0)); // xz

	KRATOS_CATCH("")
}

void FemDem3DElement::CalculateStressVector(
	Vector &rStressVector,
	const Matrix &rConstitutiveMAtrix,
	const Vector &rInfinitesimalStrainVector
	)
{
	//Infinitesimal Strain Calculation
	if (rStressVector.size() != 6)
		rStressVector.resize(6, false);

	noalias(rStressVector) = prod(rConstitutiveMAtrix, rInfinitesimalStrainVector);
}

void FemDem3DElement::CalculatePrincipalStresses(
	Vector &rPrincipalStressVector, 
	const Vector &StressVector
	)
{
	rPrincipalStressVector.resize(3);
	const double I1 = this->CalculateI1Invariant(StressVector);
	const double I2 = this->CalculateI2Invariant(StressVector);
	const double I3 = this->CalculateI3Invariant(StressVector);
	const double II1 = I1 * I1;

	const double Num = (2.0 * II1 - 9.0 * I2) * I1 + 27.0 * I3;
	const double Denom = (II1 - 3.0 * I2);

	if (Denom != 0.0) {
		double phi = Num / (2.0 * Denom * std::sqrt(Denom));

		if (std::abs(phi) > 1.0) {
			if (phi > 0.0)
				phi = 1.0;
			else
				phi = -1.0;
		}

		const double acosphi = std::acos(phi);
		phi = acosphi / 3.0;

		const double aux1 = 2.0 / 3.0 * std::sqrt(II1 - 3.0 * I2);
		const double aux2 = I1 / 3.0;

		rPrincipalStressVector[0] = aux2 + aux1 * std::cos(phi);
		rPrincipalStressVector[1] = aux2 + aux1 * std::cos(phi - 2.09439510239);
		rPrincipalStressVector[2] = aux2 + aux1 * std::cos(phi - 4.18879020478);
	} else {
		rPrincipalStressVector = ZeroVector(3);
	}
}

void FemDem3DElement::FinalizeNonLinearIteration(ProcessInfo &CurrentProcessInfo)
{
}

void FemDem3DElement::AverageVector(Vector &rAverageVector, const Vector &v, const Vector &w)
{
	int n = v.size();
	int m = w.size();
	if (n != m)
		KRATOS_ERROR << "The dimension of the vectors are different or null";
	rAverageVector.resize(n);

	for (unsigned int cont = 0; cont < n; cont++) {
		rAverageVector[cont] = (v[cont] + w[cont]) * 0.5;
	}
}

void FemDem3DElement::CalculateAverageStressOnEdge(
	Vector &rAverageVector,
	const std::vector<Element *>& VectorOfElems)
{
	Vector CurrentElementStress = this->GetValue(STRESS_VECTOR);
	rAverageVector = CurrentElementStress;
	int counter = 0;

	for (unsigned int elem = 0; elem < VectorOfElems.size(); elem++) {
		rAverageVector += VectorOfElems[elem]->GetValue(STRESS_VECTOR);
		counter++;
	}
	rAverageVector /= (counter + 1);
}

void FemDem3DElement::CalculateAverageStrainOnEdge(
	Vector &rAverageVector,
	const std::vector<Element*>& VectorOfElems)
{
	Vector CurrentElementStress = this->GetValue(STRAIN_VECTOR);
	rAverageVector = CurrentElementStress;
	int counter = 0;

	for (unsigned int elem = 0; elem < VectorOfElems.size(); elem++) {
		rAverageVector += VectorOfElems[elem]->GetValue(STRAIN_VECTOR);
		counter++;
	}
	rAverageVector /= (counter + 1);
}

// Double values
void FemDem3DElement::GetValueOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rValues,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT || rVariable == IS_DAMAGED || rVariable == STRESS_THRESHOLD || rVariable == EQUIVALENT_STRESS_VM) {
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
}

// Vector Values
void FemDem3DElement::GetValueOnIntegrationPoints(
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
void FemDem3DElement::GetValueOnIntegrationPoints(
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

// DOUBLE VARIABLES
void FemDem3DElement::CalculateOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rOutput,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT) {
		rOutput.resize(1);
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
			rOutput[PointNumber] = double(this->GetValue(DAMAGE_ELEMENT));
		}
	} else if (rVariable == IS_DAMAGED) {
		rOutput.resize(1);
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
			rOutput[PointNumber] = double(this->GetValue(IS_DAMAGED));
		}
	} else if (rVariable == STRESS_THRESHOLD) {
		rOutput.resize(1);
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
			rOutput[PointNumber] = double(this->GetValue(STRESS_THRESHOLD));
		}
	} else if (rVariable == EQUIVALENT_STRESS_VM) {
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
			const Vector stress_vector = this->GetValue(STRESS_VECTOR);
			const double I1 = this->CalculateI1Invariant(stress_vector);
			Vector deviator;
			this->CalculateDeviatorVector(deviator, stress_vector, I1);
			const double J2 = this->CalculateJ2Invariant(deviator);
			rOutput[PointNumber] = double(std::sqrt(3.0 * J2));
		}
	}
}

// 	VECTOR VARIABLES
void FemDem3DElement::CalculateOnIntegrationPoints(
	const Variable<Vector> &rVariable,
	std::vector<Vector> &rOutput,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == STRESS_VECTOR) {
		rOutput[0] = this->GetValue(STRESS_VECTOR);
	} else if (rVariable == STRAIN_VECTOR) {
		rOutput[0] = this->GetValue(STRAIN_VECTOR);
	} else if (rVariable == STRESS_VECTOR_INTEGRATED) {
		rOutput[0] = (1.0 - mDamage) * (this->GetValue(STRESS_VECTOR));
	}
}

// 	TENSOR VARIABLES
void FemDem3DElement::CalculateOnIntegrationPoints(
	const Variable<Matrix> &rVariable,
	std::vector<Matrix> &rOutput,
	const ProcessInfo &rCurrentProcessInfo)
{
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

	if (rOutput[0].size2() != dimension)
		rOutput[0].resize(dimension, dimension, false);

	if (rVariable == STRESS_TENSOR) {
		rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(STRESS_VECTOR));
	} else if (rVariable == STRAIN_TENSOR) {
		rOutput[0] = MathUtils<double>::StrainVectorToTensor(this->GetValue(STRAIN_VECTOR));
	} else if (rVariable == STRESS_TENSOR_INTEGRATED) {
		rOutput[0] = MathUtils<double>::StressVectorToTensor((1.0 - mDamage) * (this->GetValue(STRESS_VECTOR)));
	}
}

// Fills the array of characteristic lengths of the element
Vector FemDem3DElement::CalculateCharacteristicLengths()
{
	Vector lengths = ZeroVector(mNumberOfEdges);
	Geometry<Node<3>> &NodesElem = this->GetGeometry();
	Matrix Indexes;
	this->SetNodeIndexes(Indexes);

	for (unsigned int edge = 0; edge < mNumberOfEdges; edge++) {
		const double X1 = NodesElem[Indexes(edge, 0)].X0();
		const double X2 = NodesElem[Indexes(edge, 1)].X0();
		const double Y1 = NodesElem[Indexes(edge, 0)].Y0();
		const double Y2 = NodesElem[Indexes(edge, 1)].Y0();
		const double Z1 = NodesElem[Indexes(edge, 0)].Z0();
		const double Z2 = NodesElem[Indexes(edge, 1)].Z0();

		lengths[edge] = std::sqrt(std::pow((X1 - X2), 2.0) + std::pow((Y1 - Y2), 2.0) + std::pow((Z1 - Z2), 2.0));
	}
	return lengths;
}

void FemDem3DElement::Get2MaxValues(Vector &MaxValues, double a, double b, double c)
{
	MaxValues.resize(2);
	Vector V;
	V.resize(3);
	V[0] = a;
	V[1] = b;
	V[2] = c;
	int n = 3;

	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n - 1; j++) {
			if (V[j] > V[j + 1]) {
				double aux = V[j];
				V[j] = V[j + 1];
				V[j + 1] = aux;
			}
		}
	}
	MaxValues[0] = V[2];
	MaxValues[1] = V[1];
}

void FemDem3DElement::Get2MinValues(Vector &MaxValues, double a, double b, double c)
{
	MaxValues.resize(2);
	Vector V;
	V.resize(3);
	V[0] = a;
	V[1] = b;
	V[2] = c;
	int n = 3;

	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n - 1; j++) {
			if (V[j] > V[j + 1]) {
				double aux = V[j];
				V[j] = V[j + 1];
				V[j + 1] = aux;
			}
		}
	}
	MaxValues[0] = V[1];
	MaxValues[1] = V[0];
}

double FemDem3DElement::CalculateI1Invariant(const Vector& rStressVector)
{
	return rStressVector[0] + rStressVector[1] + rStressVector[2];
}

double FemDem3DElement::CalculateI2Invariant(const Vector& rStressVector)
{
	return (rStressVector[0] + rStressVector[2]) * rStressVector[1] + rStressVector[0] * rStressVector[2] +
		   -rStressVector[3] * rStressVector[3] - rStressVector[4] * rStressVector[4] - rStressVector[5] * rStressVector[5];
}

double FemDem3DElement::CalculateI3Invariant(const Vector& rStressVector)
{
	return (rStressVector[1] * rStressVector[2] - rStressVector[4] * rStressVector[4]) * rStressVector[0] -
		   rStressVector[1] * rStressVector[5] * rStressVector[5] - rStressVector[2] * rStressVector[3] * rStressVector[3] +
		   2.0 * rStressVector[3] * rStressVector[4] * rStressVector[5];
}

void FemDem3DElement::CalculateDeviatorVector(Vector& rDeviator, const Vector& rStressVector, const double I1)
{
	rDeviator.resize(6);
	rDeviator = rStressVector;
	const double Pmean = I1 / 3.0;

	rDeviator[0] -= Pmean;
	rDeviator[1] -= Pmean;
	rDeviator[2] -= Pmean;
}

double FemDem3DElement::CalculateJ2Invariant(const Vector& rDeviator)
{
	return 0.5 * (rDeviator[0] * rDeviator[0] + rDeviator[1] * rDeviator[1] + rDeviator[2] * rDeviator[2]) +
		   (rDeviator[3] * rDeviator[3] + rDeviator[4] * rDeviator[4] + rDeviator[5] * rDeviator[5]);
}

double FemDem3DElement::CalculateJ3Invariant(const Vector& rDeviator)
{
	return rDeviator[0] * (rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4]) +
		   rDeviator[3] * (-rDeviator[3] * rDeviator[2] + rDeviator[5] * rDeviator[4]) +
		   rDeviator[5] * (rDeviator[3] * rDeviator[4] - rDeviator[5] * rDeviator[1]);
}

double FemDem3DElement::CalculateLodeAngle(double J2, double J3)
{
	if (std::abs(J2) > tolerance) {
		double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
		if (sint3 < -0.95) {
			sint3 = -1.0;
		} else if (sint3 > 0.95) {
			sint3 = 1.0;
		}
		return std::asin(sint3) / 3.0;
	} else {
		return 0.0;
	}
}

void FemDem3DElement::CalculateMassMatrix(MatrixType &rMassMatrix, ProcessInfo &rCurrentProcessInfo)
{
	KRATOS_TRY

	bool ComputeLumpedMassMatrix = false;
	if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX))
		if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == true)
			ComputeLumpedMassMatrix = true;
	if (ComputeLumpedMassMatrix == false) {
		//create local system components
		LocalSystemComponents LocalSystem;

		//calculation flags
		LocalSystem.CalculationFlags.Set(SmallDisplacementElement::COMPUTE_LHS_MATRIX);

		VectorType RightHandSideVector = Vector();

		//Initialize sizes for the system components:
		this->InitializeSystemMatrices(rMassMatrix, RightHandSideVector, LocalSystem.CalculationFlags);

		//Set Variables to Local system components
		LocalSystem.SetLeftHandSideMatrix(rMassMatrix);
		LocalSystem.SetRightHandSideVector(RightHandSideVector);

		//Calculate elemental system
		CalculateDynamicSystem(LocalSystem, rCurrentProcessInfo);
	} else {
		//lumped
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int MatSize = dimension * number_of_nodes;

		if (rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize, MatSize, false);

		noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);

		double TotalMass = 0;
		TotalMass = this->CalculateTotalMass(TotalMass, rCurrentProcessInfo);

		Vector LumpFact(number_of_nodes);
		noalias(LumpFact) = ZeroVector(number_of_nodes);

		LumpFact = GetGeometry().LumpingFactors(LumpFact);

		for (unsigned int i = 0; i < number_of_nodes; i++) {
			const double temp = LumpFact[i] * TotalMass;
			for (unsigned int j = 0; j < dimension; j++) {
				unsigned int index = i * dimension + j;
				rMassMatrix(index, index) = temp;
			}
		}
	}

	KRATOS_CATCH("")
}
Vector &FemDem3DElement::CalculateVolumeForce(Vector &rVolumeForce, const Vector &rN)
{
	KRATOS_TRY

	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

	if (rVolumeForce.size() != dimension)
		rVolumeForce.resize(dimension, false);

	noalias(rVolumeForce) = ZeroVector(dimension);

	for (unsigned int j = 0; j < number_of_nodes; j++) {
		if (GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION)) { // it must be checked once at the begining only
			array_1d<double, 3> &VolumeAcceleration = GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
			for (unsigned int i = 0; i < dimension; i++)
				rVolumeForce[i] += rN[j] * VolumeAcceleration[i];
		}
	}

	rVolumeForce *= GetProperties()[DENSITY];
	return rVolumeForce;

	KRATOS_CATCH("")
}

double FemDem3DElement::GetMaxValue(Vector Strain)
{
	Vector V;
	int n = Strain.size();
	V.resize(n);

	for (int cont = 0; cont < n; cont++) {
		V[cont] = Strain[cont];
	}
	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < n - 1; j++) {
			if (V[j] > V[j + 1]) {
				double aux = V[j];
				V[j] = V[j + 1];
				V[j + 1] = aux;
			}
		}
	}
	return V[n - 1];
}

double FemDem3DElement::GetMaxAbsValue(
	const Vector& rArrayValues
	)
{
	const SizeType dimension = rArrayValues.size();

	IndexType counter = 0;
	double aux = 0.0;
	for (IndexType i = 0; i < dimension; ++i) {
		if (std::abs(rArrayValues[i]) > aux) {
			aux = std::abs(rArrayValues[i]);
			++counter;
		}
	}
	return aux;
}

double FemDem3DElement::GetMinAbsValue(
	const Vector& rArrayValues)
{
	const SizeType dimension = rArrayValues.size();

	IndexType counter = 0;
	double aux = std::numeric_limits<double>::max();
	for (IndexType i = 0; i < dimension; ++i) {
		if (std::abs(rArrayValues[i]) < aux) {
			aux = std::abs(rArrayValues[i]);
			++counter;
		}
	}
	return aux;
}

// ******* DAMAGE MECHANICS YIELD SURFACES AND EXPONENTIAL SOFTENING ********
void FemDem3DElement::IntegrateStressDamageMechanics(
	double& rThreshold,
	double &rDamage,
	const Vector &rStrainVector,
	const Vector &rStressVector,
	int Edge,
	double Length,
	bool& rIsDamaging
	)
{
	const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];
	if (yield_surface == "ModifiedMohrCoulomb") {
		this->ModifiedMohrCoulombCriterion(rThreshold, rDamage, 
			rStressVector, Edge, Length, rIsDamaging);
	} else if (yield_surface == "SimoJu") {
		this->SimoJuCriterion(rThreshold, rDamage, 
			rStrainVector, rStressVector, Edge, Length, rIsDamaging);
	} else if (yield_surface == "Rankine") {
		this->RankineCriterion(rThreshold, rDamage, 
			rStressVector, Edge, Length, rIsDamaging);
	} else if (yield_surface == "DruckerPrager") {
		this->DruckerPragerCriterion(rThreshold, rDamage, 
			rStressVector, Edge, Length, rIsDamaging);
	} else if (yield_surface == "RankineFragile") {
		this->RankineFragileLaw(rThreshold, rDamage, 
			rStressVector, Edge, Length, rIsDamaging);
	} else if (yield_surface == "Elastic") {
		this->ElasticLaw(rThreshold, rDamage, 
			rStressVector, Edge, Length, rIsDamaging);
	} else {
		KRATOS_ERROR << "Yield Surface not defined "<< std::endl;
	}
}

void FemDem3DElement::ModifiedMohrCoulombCriterion(
	double& rThreshold,
	double &rDamage, 
	const Vector &rStressVector, 
	const int Edge, 
	const double Length,
	bool& rIsDamaging
	)
{
	const auto& properties = this->GetProperties();
	const double sigma_c = properties[YIELD_STRESS_C];
	const double sigma_t = properties[YIELD_STRESS_T];
	double friction_angle = properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180.0; // In radians!
	const double E = properties[YOUNG_MODULUS];
	const double Gt = properties[FRAC_ENERGY_T];

	KRATOS_WARNING_IF("friction_angle", friction_angle < tolerance) << "Friction Angle not defined, assumed equal to 32deg" << std::endl;
	KRATOS_ERROR_IF(sigma_c < tolerance) << "Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa " << std::endl;
	KRATOS_ERROR_IF(sigma_t < tolerance) << "Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa" << std::endl;
	KRATOS_ERROR_IF(Gt < tolerance) << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa " << std::endl;
	// Check input variables
	if (friction_angle < tolerance) {
		friction_angle = 32.0 * Globals::Pi / 180;
		KRATOS_WARNING("friction_angle") << " Friction Angle not defined, assumed equal to 32 deg " << std::endl;
	}

	const double R = std::abs(sigma_c / sigma_t);
	const double Rmorh = std::pow(std::tan((Globals::Pi / 4.0) + friction_angle / 2.0), 2.0);
	const double alpha_r = R / Rmorh;
	const double c_max = std::abs(sigma_c);
	const double sinphi = std::sin(friction_angle);

	const double I1 = this->CalculateI1Invariant(rStressVector);
	Vector Deviator = ZeroVector(6);
	this->CalculateDeviatorVector(Deviator, rStressVector, I1);
	const double J2 = this->CalculateJ2Invariant(Deviator);
	const double J3 = this->CalculateJ3Invariant(Deviator);
	const double K1 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) * sinphi;
	const double K2 = 0.5 * (1.0 + alpha_r) - 0.5 * (1.0 - alpha_r) / sinphi;
	const double K3 = 0.5 * (1.0 + alpha_r) * sinphi - 0.5 * (1.0 - alpha_r);
	const double n = sigma_c / sigma_t;
	const double A = 1.00 / (n * n * Gt * E / (Length * std::pow(sigma_c, 2)) - 0.5);

	// Check Modified Mohr-Coulomb criterion
	double uniaxial_stress;
	if (std::abs(I1) < tolerance) {
		uniaxial_stress = 0.0;
	} else {
		const double theta = CalculateLodeAngle(J2, J3);
		uniaxial_stress = (2.00 * std::tan(Globals::Pi * 0.25 + friction_angle * 0.5) / std::cos(friction_angle)) *
						  ((I1 * K3 / 3.0) + std::sqrt(J2) * (K1 * std::cos(theta) - K2 * std::sin(theta) * sinphi /
			              std::sqrt(3.0)));
	}

	if (rThreshold < tolerance) {
		rThreshold = c_max; // 1st iteration sets threshold as c_max
	}  else if (c_max > rThreshold) { // remeshing stuff
		rThreshold = c_max;
	}

	const double F = uniaxial_stress - rThreshold;
	if (F <= 0.0) { // Elastic region --> Damage is constant
		rDamage = this->GetConvergedDamages(Edge);
	} else {
		this->CalculateExponentialDamage(rDamage, A, uniaxial_stress, c_max);
		rThreshold = uniaxial_stress;
		rIsDamaging = true;
	}
}

void FemDem3DElement::RankineCriterion(
	double& rThreshold,
	double &rDamage, 
	const Vector &rStressVector, 
	const int Edge, 
	const double Length,
	bool& rIsDamaging
	)
{
	Vector PrincipalStressVector = ZeroVector(3);
	this->CalculatePrincipalStresses(PrincipalStressVector, rStressVector);

	const auto& properties = this->GetProperties();
	const double sigma_c = properties[YIELD_STRESS_C];
	const double sigma_t = properties[YIELD_STRESS_T];
	const double E = properties[YOUNG_MODULUS];
	const double Gt = properties[FRAC_ENERGY_T];
	const double c_max = std::abs(sigma_t);
	KRATOS_ERROR_IF(sigma_c < tolerance) << "Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa " << std::endl;
	KRATOS_ERROR_IF(sigma_t < tolerance) << "Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa" << std::endl;
	KRATOS_ERROR_IF(Gt < tolerance) << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa " << std::endl;

	const double A = 1.00 / (Gt * E / (Length * std::pow(sigma_c, 2)) - 0.5);
	KRATOS_ERROR_IF(A < tolerance) << " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T" << std::endl;

	double uniaxial_stress = GetMaxValue(PrincipalStressVector);

	if (rThreshold < tolerance) {
		rThreshold = c_max; // 1st iteration sets threshold as c_max
	}  else if (c_max > rThreshold) { // remeshing stuff
		rThreshold = c_max;
	}

	const double F = uniaxial_stress - rThreshold;

	if (F <= 0) { // Elastic region --> Damage is constant
		rDamage = this->GetConvergedDamages(Edge);
	} else {
		this->CalculateExponentialDamage(rDamage, A, uniaxial_stress, c_max);
		rThreshold = uniaxial_stress;
		rIsDamaging = true;
	}
}

void FemDem3DElement::DruckerPragerCriterion(
	double& rThreshold,
	double &rDamage, 
	const Vector &rStressVector, 
	const int Edge, 
	const double Length,
	bool& rIsDamaging
	)
{
	const auto& properties = this->GetProperties();
	const double sigma_c = properties[YIELD_STRESS_C];
	const double sigma_t = properties[YIELD_STRESS_T];
	double friction_angle = properties[INTERNAL_FRICTION_ANGLE] * Globals::Pi / 180; // In radians!
	const double E = properties[YOUNG_MODULUS];
	const double Gt = properties[FRAC_ENERGY_T];

	KRATOS_WARNING_IF("friction_angle", friction_angle < tolerance) << "Friction Angle not defined, assumed equal to 32deg" << std::endl;
	KRATOS_ERROR_IF(sigma_c < tolerance) << "Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa " << std::endl;
	KRATOS_ERROR_IF(sigma_t < tolerance) << "Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa" << std::endl;
	KRATOS_ERROR_IF(Gt < tolerance) << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa " << std::endl;
	// Check input variables
	if (friction_angle < tolerance) {
		friction_angle = 32.0 * Globals::Pi / 180;
		std::cout << "Friction Angle not defined, assumed equal to 32deg " << std::endl;
	}
	const double c_max = std::abs(sigma_t * (3.0 + std::sin(friction_angle)) / (3.0 * std::sin(friction_angle) - 3.0));
	const double I1 = CalculateI1Invariant(rStressVector);
	Vector Deviator = ZeroVector(6);
	this->CalculateDeviatorVector(Deviator, rStressVector, I1);
	const double J2 = CalculateJ2Invariant(Deviator);
	const double A = 1.00 / (Gt * E / (Length * std::pow(sigma_c, 2)) - 0.5);
	KRATOS_ERROR_IF(A < tolerance) << " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T" << std::endl;

	// Check DruckerPrager criterion
	double uniaxial_stress;
	if (std::abs(I1) < tolerance) {
		uniaxial_stress = 0.0;
	} else {
		const double CFL = -std::sqrt(3.0) * (3.0 - std::sin(friction_angle)) / (3.0 * std::sin(friction_angle) - 3.0);
		const double TEN0 = 2.0 * I1 * std::sin(friction_angle) / (std::sqrt(3.0) * (3.0 - std::sin(friction_angle))) + std::sqrt(J2);
		uniaxial_stress = std::abs(CFL * TEN0);
	}

	if (rThreshold < tolerance) {
		rThreshold = c_max; // 1st iteration sets threshold as c_max
	}  else if (c_max > rThreshold) { // remeshing stuff
		rThreshold = c_max;
	}

	const double F = uniaxial_stress - rThreshold;

	if (F <= 0) { // Elastic region --> Damage is constant
		rDamage = this->GetConvergedDamages(Edge);
	} else {
		this->CalculateExponentialDamage(rDamage, A, uniaxial_stress, c_max);
		rThreshold = uniaxial_stress;
		rIsDamaging = true;
	}
}

void FemDem3DElement::SimoJuCriterion(
	double& rThreshold,
	double& rDamage,
	const Vector &rStrainVector,
	const Vector &rStressVector,
	const int Edge,
	const double Length,
	bool& rIsDamaging
	)
{
	Vector PrincipalStressVector = ZeroVector(3);
	this->CalculatePrincipalStresses(PrincipalStressVector, rStressVector);

	const auto& properties = this->GetProperties();
	const double sigma_t = properties[YIELD_STRESS_T];
	const double sigma_c = properties[YIELD_STRESS_C];
	const double E = properties[YOUNG_MODULUS];
	const double Gt = properties[FRAC_ENERGY_T];
	const double n = std::abs(sigma_c / sigma_t);
	const double c_max = std::abs(sigma_c) / std::sqrt(E);

	double SumA = 0.0, SumB = 0.0, SumC = 0.0;
	for (unsigned int cont = 0; cont < 3; cont++) {
		SumA += std::abs(PrincipalStressVector[cont]);
		SumB += 0.5 * (PrincipalStressVector[cont] + std::abs(PrincipalStressVector[cont]));
		SumC += 0.5 * (-PrincipalStressVector[cont] + std::abs(PrincipalStressVector[cont]));
	}
	const double ere0 = SumB / SumA;
	const double ere1 = SumC / SumA;

	// Check SimoJu criterion
	double uniaxial_stress;
	if (rStrainVector[0] + rStrainVector[1] + rStrainVector[2] < tolerance) {
		uniaxial_stress = 0.0;
	} else {
		double auxf = 0.0;
		for (unsigned int cont = 0; cont < 6; cont++) {
			auxf += rStrainVector[cont] * rStressVector[cont]; // E*S
		}
		uniaxial_stress = std::sqrt(auxf);
		uniaxial_stress *= (ere0 * n + ere1);
	}

	if (rThreshold < tolerance) {
		rThreshold = c_max; // 1st iteration sets threshold as c_max
	}  else if (c_max > rThreshold) { // remeshing stuff
		rThreshold = c_max;
	}
	const double F = uniaxial_stress - rThreshold;

	const double A = 1.00 / (Gt * n * n * E / (Length * std::pow(sigma_c, 2)) - 0.5);
	KRATOS_ERROR_IF(A < tolerance) << " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T" << std::endl;

	if (F <= 0.0) { // Elastic region --> Damage is constant
		rDamage = this->GetConvergedDamages(Edge);
	} else {
		this->CalculateExponentialDamage(rDamage, A, uniaxial_stress, c_max);
		rThreshold = uniaxial_stress;
		rIsDamaging = true;
	}
}

void FemDem3DElement::RankineFragileLaw(
	double& rThreshold,
	double &rDamage, 
	const Vector &rStressVector, 
	const int Edge, 
	const double Length,
	bool& rIsDamaging
	)
{
	Vector PrincipalStressVector = ZeroVector(3);
	this->CalculatePrincipalStresses(PrincipalStressVector, rStressVector);

	const auto& properties = this->GetProperties();
	const double sigma_c = properties[YIELD_STRESS_C];
	const double sigma_t = properties[YIELD_STRESS_T];
	const double E = properties[YOUNG_MODULUS];
	const double Gt = properties[FRAC_ENERGY_T];
	const double c_max = std::abs(sigma_t);

	const double  A = 1.00 / (Gt * E / (Length * std::pow(sigma_c, 2)) - 0.5);
	KRATOS_ERROR_IF(A < tolerance) << " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T" << std::endl;

	// F = f-c = 0 classical definition of yield surface
	const double uniaxial_stress = GetMaxValue(PrincipalStressVector);

	if (rThreshold < tolerance) {
		rThreshold = c_max; // 1st iteration sets threshold as c_max
	}  else if (c_max > rThreshold) { // remeshing stuff
		rThreshold = c_max;
	}
	const double F = uniaxial_stress - rThreshold;

	if (F <= 0.0) { // Elastic region --> Damage is constant
		rDamage = this->GetConvergedDamages(Edge);
	} else {
		rDamage = 0.98; // Fragile  law
		rThreshold = uniaxial_stress;
		rIsDamaging = true;
	}
}

void FemDem3DElement::ElasticLaw(
	double& rThreshold,
	double &rDamage, 
	const Vector &rStressVector, 
	const int Edge, 
	const double Length,
	bool& rIsDamaging
	)
{
	const auto& properties = this->GetProperties();
	const double sigma_t = properties[YIELD_STRESS_T];
	const double c_max = std::abs(sigma_t);
	rDamage = 0.0;
	if (rThreshold < tolerance) {
		rThreshold = c_max;
	} // 1st iteration sets threshold as c_max
}

void FemDem3DElement::CalculateExponentialDamage(
	double& rDamage,
	const double DamageParameter,
	const double UniaxialStress,
	const double InitialThrehsold
	)
{
	rDamage = 1.0 - (InitialThrehsold / UniaxialStress) * std::exp(DamageParameter *
			 (1.0 - UniaxialStress / InitialThrehsold)); // Exponential softening law
	if (rDamage > 0.99) rDamage = 0.99;
}

// Computes the damage of the element considering different fracture modes
double FemDem3DElement::CalculateElementalDamage(const Vector &EdgeDamages)
{
	// 7 modes of fracture of the tetrahedron
	Vector DamageModeFracture = ZeroVector(7);
	const double one_third = 1.0 / 3.0;
	DamageModeFracture[0] = one_third * (EdgeDamages[0] + EdgeDamages[1] + EdgeDamages[2]);
	DamageModeFracture[1] = one_third * (EdgeDamages[0] + EdgeDamages[3] + EdgeDamages[4]);
	DamageModeFracture[2] = one_third * (EdgeDamages[1] + EdgeDamages[3] + EdgeDamages[5]);
	DamageModeFracture[3] = 0.25 * (EdgeDamages[1] + EdgeDamages[2] + EdgeDamages[3] + EdgeDamages[4]);
	DamageModeFracture[4] = 0.25 * (EdgeDamages[0] + EdgeDamages[1] + EdgeDamages[4] + EdgeDamages[5]);
	DamageModeFracture[5] = one_third * (EdgeDamages[2] + EdgeDamages[4] + EdgeDamages[5]);
	DamageModeFracture[6] = 0.25 * (EdgeDamages[0] + EdgeDamages[2] + EdgeDamages[3] + EdgeDamages[5]);
	return this->GetMaxValue(DamageModeFracture);
}

void FemDem3DElement::SetValueOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rValues,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT) {
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
			mDamage = rValues[PointNumber];
		}
	} else if (rVariable == STRESS_THRESHOLD) {
		for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
			mThreshold = rValues[PointNumber];
		}
	} else {
		for (unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints().size(); ++point_number) {
			this->SetValue(rVariable, rValues[point_number]);
		}
	}
}

void FemDem3DElement::SetValueOnIntegrationPoints(
	const Variable<Vector> &rVariable,
	std::vector<Vector> &rValues,
	const ProcessInfo &rCurrentProcessInfo)
{
	for (unsigned int point_number = 0; point_number < GetGeometry().IntegrationPoints().size(); ++point_number) {
		this->SetValue(rVariable, rValues[point_number]);
	}
}


// Methods to compute the tangent tensor by numerical derivation
void FemDem3DElement::CalculateTangentTensor(
	Matrix& TangentTensor,
	const Vector& rStrainVectorGP,
	const Vector& rStressVectorGP,
	const Matrix& rElasticMatrix
	)
{
	const double number_components = rStrainVectorGP.size();
	TangentTensor.resize(number_components, number_components);
	Vector perturbed_stress, perturbed_strain;
	perturbed_strain.resize(number_components);
	perturbed_stress.resize(number_components);
	
	for (unsigned int component = 0; component < number_components; component++) {
		double perturbation;
		this->CalculatePerturbation(rStrainVectorGP, perturbation, component);
		this->PerturbateStrainVector(perturbed_strain, rStrainVectorGP, perturbation, component);
		this->IntegratePerturbedStrain(perturbed_stress, perturbed_strain, rElasticMatrix);
		const Vector& delta_stress = perturbed_stress - rStressVectorGP;
		this->AssignComponentsToTangentTensor(TangentTensor, delta_stress, perturbation, component);
	}
}

void FemDem3DElement::CalculatePerturbation(
	const Vector& rStrainVectorGP,
	double& rPerturbation,
	const int Component
	)
{
	double perturbation_1, perturbation_2;
	if (std::abs(rStrainVectorGP[Component]) > tolerance) {
		perturbation_1 = 1.0e-5 * rStrainVectorGP[Component];
	} else {
		double min_strain_component = this->GetMinAbsValue(rStrainVectorGP);
		perturbation_1 = 1.0e-5 * min_strain_component;
	}
	const double max_strain_component = this->GetMaxAbsValue(rStrainVectorGP);
	perturbation_2 = 1.0e-10 * max_strain_component;
	rPerturbation = std::max(perturbation_1, perturbation_2);
}

void FemDem3DElement::PerturbateStrainVector(
	Vector& rPerturbedStrainVector,
	const Vector& rStrainVectorGP,
	const double Perturbation,
	const int Component
	)
{
    noalias(rPerturbedStrainVector) = rStrainVectorGP;
    rPerturbedStrainVector[Component] += Perturbation;
}

void FemDem3DElement::IntegratePerturbedStrain(
	Vector& rPerturbedStressVector,
	const Vector& rPerturbedStrainVector,
	const Matrix& rElasticMatrix
	)
{
	const Vector& perturbed_predictive_stress = prod(rElasticMatrix, rPerturbedStrainVector);
	Vector damages_edges = ZeroVector(mNumberOfEdges);
	const Vector& r_characteristic_lengths = this->CalculateCharacteristicLengths();

	for (unsigned int edge = 0; edge < mNumberOfEdges; edge++) {
		std::vector<Element*> EdgeNeighbours = this->GetEdgeNeighbourElements(edge);

		// We compute the average stress/strain on edge to integrate
		Vector average_stress_vector = perturbed_predictive_stress;
		Vector average_strain_vector = rPerturbedStrainVector;
		int counter = 0;
		for (unsigned int elem = 0; elem < EdgeNeighbours.size(); elem++) {
			average_stress_vector += EdgeNeighbours[elem]->GetValue(STRESS_VECTOR);
			average_strain_vector += EdgeNeighbours[elem]->GetValue(STRAIN_VECTOR);
			counter++;
		}
		average_stress_vector /= (counter + 1);
		average_strain_vector /= (counter + 1);

		bool dummy = false; 

		Vector perturbed_integrated_stress;
		double damage_edge = mDamages[edge];
		double threshold = mThresholds[edge];

		this->IntegrateStressDamageMechanics(threshold,
											 damage_edge,
											 average_strain_vector, 
											 average_stress_vector, 
											 edge, 
											 r_characteristic_lengths[edge],
											 dummy);
		damages_edges[edge] = damage_edge;
	} // Loop edges
	const double damage_element = this->CalculateElementalDamage(damages_edges);
	rPerturbedStressVector = (1.0 - damage_element) * perturbed_predictive_stress;
}

void FemDem3DElement::AssignComponentsToTangentTensor(
	Matrix& rTangentTensor,
	const Vector& rDeltaStress,
	const double Perturbation,
	const int Component
	)
{
	const int voigt_size = rDeltaStress.size();
	for (IndexType row = 0; row < voigt_size; ++row) {
		rTangentTensor(row, Component) = rDeltaStress[row] / Perturbation;
	}
}




} // namespace Kratos