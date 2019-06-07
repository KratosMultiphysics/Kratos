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

// System includes

// External includes

// Project includes
#include "generic_small_strain_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: SmallDisplacementElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: SmallDisplacementElement(NewId, pGeometry, pProperties)
{
	// BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

	// Each component == Each edge
	if (mNonConvergedThresholds.size() != NumberOfEdges)
    	mNonConvergedThresholds.resize(NumberOfEdges);
	noalias(mNonConvergedThresholds) = ZeroVector(NumberOfEdges);   // Equivalent stress

	if (mThresholds.size() != NumberOfEdges)
    	mThresholds.resize(NumberOfEdges);
	noalias(mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

	if (mDamages.size() != NumberOfEdges)
    	mDamages.resize(NumberOfEdges);
	noalias(mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge

	if (mNonConvergedDamages.size() != NumberOfEdges)
    	mNonConvergedDamages.resize(NumberOfEdges);
	noalias(mNonConvergedDamages) = ZeroVector(NumberOfEdges); // mDamages on edges of "i" iteration
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim,TyieldSurf>::GenericSmallStrainFemDemElement(GenericSmallStrainFemDemElement const& rOther)
	: SmallDisplacementElement(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf> &GenericSmallStrainFemDemElement<TDim,TyieldSurf>::operator=(GenericSmallStrainFemDemElement const& rOther)
{
	SmallDisplacementElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
{
	return Element::Pointer(new GenericSmallStrainFemDemElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
	GenericSmallStrainFemDemElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
	return Element::Pointer(new GenericSmallStrainFemDemElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim,TyieldSurf>::~GenericSmallStrainFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::InitializeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
	this->ComputeEdgeNeighbours(rCurrentProcessInfo);
	this->InitializeInternalVariablesAfterMapping();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{
	KRATOS_TRY

	//create and initialize element variables:
	ElementDataType variables;
	this->InitializeElementData(variables,rCurrentProcessInfo);

	//create constitutive law parameters:
	ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

	//set constitutive law flags:
	Flags &ConstitutiveLawOptions = values.GetOptions();

	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	//reading integration points
	const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	for (SizeType point_number = 0; point_number < r_integration_points.size(); point_number++) {
		//compute element kinematic variables B, F, DN_DX ...
		this->CalculateKinematics(variables,point_number);

		//calculate material response
		this->CalculateMaterialResponse(variables, values, point_number);
	}
	this->SetValue(STRESS_VECTOR, values.GetStressVector());
	this->SetValue(STRAIN_VECTOR, values.GetStrainVector());

	KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::FinalizeNonLinearIteration(
    ProcessInfo& CurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateLocalSystem(
	MatrixType& rLeftHandSideMatrix,
	VectorType& rRightHandSideVector,
	ProcessInfo& rCurrentProcessInfo
    )
{
	KRATOS_TRY
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

	//resizing as needed the LHS
	unsigned int mat_size = number_of_nodes * (dimension); 
	if (rLeftHandSideMatrix.size1() != mat_size)
	rLeftHandSideMatrix.resize(mat_size, mat_size, false);

	noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
	if (rRightHandSideVector.size() != mat_size )
	rRightHandSideVector.resize(mat_size, false);

	noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS

	//create and initialize element variables:
	ElementDataType variables;
	this->InitializeElementData(variables, rCurrentProcessInfo);
	//create constitutive law parameters:
	ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

	//set constitutive law flags:
	Flags& ConstitutiveLawOptions = values.GetOptions();

	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	// reading integration points
	const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	//auxiliary terms
	Vector volume_force(dimension);
	const double characteristic_length = this->CalculateCharacteristicLength(this);

	for (SizeType point_number = 0; point_number < r_integration_points.size(); point_number++) {

		//compute element kinematic variables B, F, DN_DX ...
		this->CalculateKinematics(variables,point_number);

		//calculate material response
		this->CalculateMaterialResponse(variables, values, point_number);

		//calculating weights for integration on the "reference configuration"
		variables.IntegrationWeight = r_integration_points[point_number].Weight() * variables.detJ;
		variables.IntegrationWeight = this->CalculateIntegrationWeight( variables.IntegrationWeight );

        bool is_damaging = false;
        if (yield_surface != "Elastic") {
			// Loop over edges of the element...
			Vector average_stress_edge(VoigtSize);
			Vector average_strain_edge(VoigtSize);
			noalias(average_stress_edge) = this->GetValue(STRESS_VECTOR);
			noalias(average_strain_edge) = this->GetValue(STRAIN_VECTOR);

			for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {

				this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
				this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

				double damage_edge = mDamages[edge];
				double threshold = mThresholds[edge];
				
				this->IntegrateStressDamageMechanics(threshold, damage_edge, average_strain_edge, 
					average_stress_edge, edge, characteristic_length, values, is_damaging);

				mNonConvergedDamages[edge] = damage_edge;
				mNonConvergedThresholds[edge] = threshold;
			} // Loop over edges
		} else {
			noalias(mNonConvergedDamages) = ZeroVector(NumberOfEdges);
		}

		// Calculate the elemental Damage...
		const double damage_element = this->CalculateElementalDamage(mNonConvergedDamages);
		const Vector& r_predictive_stress_vector = values.GetStressVector();
		const Vector& r_integrated_stress_vector = (1.0 - damage_element) * r_predictive_stress_vector;
		const Vector& r_strain_vector = values.GetStrainVector();

		//contributions to stiffness matrix calculated on the reference config
		const Matrix& r_elastic_matrix =  values.GetConstitutiveMatrix();
		Matrix tangent_tensor;

		if (is_damaging == true && norm_2(r_strain_vector) > tolerance) {
			this->CalculateTangentTensor(tangent_tensor, r_strain_vector, r_integrated_stress_vector, r_elastic_matrix, values);
			noalias(rLeftHandSideMatrix) += variables.IntegrationWeight * prod(trans(variables.B), Matrix(prod(tangent_tensor, variables.B)));
		} else {
			noalias(rLeftHandSideMatrix) += variables.IntegrationWeight * (1.0 - damage_element) * prod(trans(variables.B), Matrix(prod(r_elastic_matrix, variables.B)));
		}
		
		//contribution to external forces
		volume_force = this->CalculateVolumeForce(volume_force, variables.N);

		// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
		this->CalculateAndAddExternalForces(rRightHandSideVector, variables, volume_force, variables.IntegrationWeight );

		// operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
		rRightHandSideVector -= variables.IntegrationWeight * prod(trans(variables.B), r_integrated_stress_vector);
	}
	KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, 
    ProcessInfo& rCurrentProcessInfo
    )
{
	KRATOS_TRY
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

	//resizing as needed the LHS
	unsigned int mat_size = number_of_nodes * (dimension); 
	rLeftHandSideMatrix = ZeroMatrix(mat_size, mat_size);

	//create and initialize element variables:
	ElementDataType variables;
	this->InitializeElementData(variables,rCurrentProcessInfo);

	//create constitutive law parameters:
	ConstitutiveLaw::Parameters values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

	//set constitutive law flags:
	Flags& r_constitutive_law_options = values.GetOptions();
	r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	//reading integration points
	const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	for (SizeType point_number = 0; point_number < r_integration_points.size(); point_number++ ) {
		
		//compute element kinematic variables B, F, DN_DX ...
		this->CalculateKinematics(variables, point_number);

		//calculate material response
		this->CalculateMaterialResponse(variables, values, point_number);

		//calculating weights for integration on the "reference configuration"
		variables.IntegrationWeight = r_integration_points[point_number].Weight() * variables.detJ;
		variables.IntegrationWeight = this->CalculateIntegrationWeight( variables.IntegrationWeight );

		// Calculate the elemental Damage...
		const double damage_element = this->CalculateElementalDamage(mDamages);
		const Matrix& r_elastic_matrix =  values.GetConstitutiveMatrix();

		noalias(rLeftHandSideMatrix) += variables.IntegrationWeight * (1.0 - damage_element) * prod(trans(variables.B), Matrix(prod(r_elastic_matrix, variables.B)));
	}
	KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
	KRATOS_TRY
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

	//resizing as needed the LHS
	unsigned int mat_size = number_of_nodes * (dimension); 
	rRightHandSideVector = ZeroVector(mat_size);
	//create and initialize element variables:
	ElementDataType variables;
	this->InitializeElementData(variables, rCurrentProcessInfo);

	//create constitutive law parameters:
	ConstitutiveLaw::Parameters values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

	//set constitutive law flags:
	Flags& ConstitutiveLawOptions = values.GetOptions();

	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

	//reading integration points
	const GeometryType::IntegrationPointsArrayType& r_integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	//auxiliary terms
	Vector volume_force(dimension);
	noalias(volume_force) = ZeroVector(dimension);
	const double characteristic_length = this->CalculateCharacteristicLength(this);

	for (SizeType point_number = 0; point_number < r_integration_points.size(); point_number++) {
		
		//compute element kinematic variables B, F, DN_DX ...
		this->CalculateKinematics(variables, point_number);

		//calculate material response
		this->CalculateMaterialResponse(variables,values,point_number);

		//calculating weights for integration on the "reference configuration"
		variables.IntegrationWeight = r_integration_points[point_number].Weight() * variables.detJ;
		variables.IntegrationWeight = this->CalculateIntegrationWeight(variables.IntegrationWeight);

		bool is_damaging = false;
		Vector damages(NumberOfEdges);
		if (yield_surface != "Elastic") {
			// Loop over edges of the element...
			Vector average_stress_edge(VoigtSize);
			Vector average_strain_edge(VoigtSize);
			noalias(average_stress_edge) = this->GetValue(STRESS_VECTOR);
			noalias(average_strain_edge) = this->GetValue(STRAIN_VECTOR);

			for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
				this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
				this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

				double damage_edge = mDamages[edge];
				double threshold = mThresholds[edge];
				
				this->IntegrateStressDamageMechanics(threshold, damage_edge, average_strain_edge, 
					average_stress_edge, edge, characteristic_length, values, is_damaging);

				damages[edge] = damage_edge;
			} // Loop over edges
		} else {
			noalias(damages) = ZeroVector(NumberOfEdges);
		}

		// Calculate the elemental Damage...
		const double damage_element = this->CalculateElementalDamage(damages);
		const Vector& r_predictive_stress_vector = values.GetStressVector();
		const Vector& r_integrated_stress_vector = (1.0 - damage_element) * r_predictive_stress_vector;
		const Vector& r_strain_vector = values.GetStrainVector();

		this->SetValue(STRESS_VECTOR, r_predictive_stress_vector);
		this->SetValue(STRESS_VECTOR_INTEGRATED, r_integrated_stress_vector);
		this->SetValue(STRAIN_VECTOR, r_strain_vector);
		
		//contribution to external forces
		volume_force = this->CalculateVolumeForce(volume_force, variables.N);

		// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
		this->CalculateAndAddExternalForces(rRightHandSideVector, variables, volume_force, variables.IntegrationWeight );

		// operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
		rRightHandSideVector -= variables.IntegrationWeight * prod(trans(variables.B), r_integrated_stress_vector);
	}
	KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::IntegrateStressDamageMechanics(
	double& rThreshold,
	double& rDamage,
	const Vector& rStrainVector,
	const Vector& rStressVector,
	const int Edge,
	const double CharacteristicLength,
    ConstitutiveLaw::Parameters& rValues,
	bool& rIsDamaging
	)
{
    double uniaxial_stress;
    this->CalculateEquivalentStress(rStressVector, rStrainVector, uniaxial_stress, rValues);

    double initial_threshold;
    this->GetInitialUniaxialThreshold(rValues, initial_threshold);

    double damage_parameter; // A parameter
    this->CalculateDamageParameter(rValues, damage_parameter, CharacteristicLength);

	if (rThreshold < tolerance) {
		rThreshold = initial_threshold; // 1st iteration sets threshold as c_max
	}  else if (initial_threshold > rThreshold) { // remeshing stuff
		rThreshold = initial_threshold;
	}

    const double F = uniaxial_stress - rThreshold;
	if (F <= 0.0) { // Elastic region --> Damage is constant
		rDamage = mDamages[Edge];
	} else {
		this->CalculateExponentialDamage(rDamage, damage_parameter, uniaxial_stress, initial_threshold);
		rThreshold = uniaxial_stress;
		rIsDamaging = true;
	}
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<3,0>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge3D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<2,0>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<2,1>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<2,2>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<2,3>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<2,4>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}
template<>
void GenericSmallStrainFemDemElement<2,5>::CalculateAverageVariableOnEdge(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector,
    const int edge
    )
{
    this->CalculateAverageVariableOnEdge2D(pCurrentElement, ThisVariable, rAverageVector, edge);
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateAverageVariableOnEdge2D(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector, // must be assigned to the current element outside this method
    const int edge
    )
{
    auto& r_elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_elem_neigb.size() == 0) << " Neighbour Elements not calculated" << std::endl;
    auto& neighbour = r_elem_neigb[edge];

    rAverageVector += pCurrentElement->GetValue(ThisVariable);
    rAverageVector *= 0.5;
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateAverageVariableOnEdge3D(
    const Element* pCurrentElement,
    const Variable<Vector> ThisVariable,
    Vector& rAverageVector, // must be assigned to the current element outside this method
    const int edge
    )
{
    std::vector<Element*> p_edge_neighbours = this->GetEdgeNeighbourElements(edge);
	int counter = 0;

	for (unsigned int elem = 0; elem < p_edge_neighbours.size(); elem++) {
		rAverageVector += p_edge_neighbours[elem]->GetValue(ThisVariable);
		counter++;
	}
	rAverageVector /= (counter + 1);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double GenericSmallStrainFemDemElement<3,0>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<3,1>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<3,2>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<3,3>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<3,4>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<3,5>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
    return this->CalculateElementalDamage3D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<2,0>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<2,1>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<2,2>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<2,3>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<2,4>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}
template<>
double GenericSmallStrainFemDemElement<2,5>::CalculateElementalDamage(const Vector& rEdgeDamages)
{
   return this->CalculateElementalDamage2D(rEdgeDamages);
}

template<unsigned int TDim, unsigned int TyieldSurf>
double GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateElementalDamage3D(const Vector& rEdgeDamages)
{
	// 7 modes of fracture of the tetrahedron
	Vector damage_mode_fracture = ZeroVector(7);
	const double one_third = 1.0 / 3.0;
	damage_mode_fracture[0] = one_third * (rEdgeDamages[0] + rEdgeDamages[1] + rEdgeDamages[2]);
	damage_mode_fracture[1] = one_third * (rEdgeDamages[0] + rEdgeDamages[3] + rEdgeDamages[4]);
	damage_mode_fracture[2] = one_third * (rEdgeDamages[1] + rEdgeDamages[3] + rEdgeDamages[5]);
	damage_mode_fracture[3] = 0.25 * (rEdgeDamages[1] + rEdgeDamages[2] + rEdgeDamages[3] + rEdgeDamages[4]);
	damage_mode_fracture[4] = 0.25 * (rEdgeDamages[0] + rEdgeDamages[1] + rEdgeDamages[4] + rEdgeDamages[5]);
	damage_mode_fracture[5] = one_third * (rEdgeDamages[2] + rEdgeDamages[4] + rEdgeDamages[5]);
	damage_mode_fracture[6] = 0.25 * (rEdgeDamages[0] + rEdgeDamages[2] + rEdgeDamages[3] + rEdgeDamages[5]);
	return ConstitutiveLawUtilities<VoigtSize>::GetMaxValue(damage_mode_fracture);
}

template<unsigned int TDim, unsigned int TyieldSurf>
double GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateElementalDamage2D(const Vector& rEdgeDamages)
{
	Vector two_max_values;
	ConstitutiveLawUtilities<VoigtSize>::Get2MaxValues(two_max_values, rEdgeDamages[0], rEdgeDamages[1], rEdgeDamages[2]);
	return 0.5*(two_max_values[0] + two_max_values[1]);
}


/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::UpdateDataBase()
{	
	for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
		mDamages[edge] = mNonConvergedDamages[edge];
		mThresholds[edge] = mNonConvergedThresholds[edge];
	}

	const double converged_damage = this->CalculateElementalDamage(mDamages);
	if (converged_damage > mDamage) mDamage = converged_damage;
	const double converged_threshold = this->CalculateElementalDamage(mThresholds);
	if (converged_threshold > mThreshold) mThreshold = converged_threshold;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
	this->UpdateDataBase();

	if (mDamage >= 0.98) {
		this->Set(ACTIVE, false);
		// We set a "flag" to generate the DEM 
		rCurrentProcessInfo[GENERATE_DEM] = true;
	}
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::InitializeInternalVariablesAfterMapping()
{
    // After the mapping, the thresholds of the edges (are equal to 0.0) are imposed equal to the IP threshold
    const double element_threhsold = mThreshold;
    if (norm_2(mThresholds) < std::numeric_limits<double>::epsilon()) {
        for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
            mThresholds[edge] = element_threhsold;
        }
    }

    // IDEM with the edge damages
    const double damage_element = mDamage;
    if (norm_2(mDamages) < std::numeric_limits<double>::epsilon()) {
        for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
            mDamages[edge] = damage_element;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<3,0>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
{
    this->AuxComputeEdgeNeighbours(rCurrentProcessInfo);
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<2,0>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,1>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,2>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,3>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,4>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}
template<> 
void GenericSmallStrainFemDemElement<2,5>::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo) {}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::AuxComputeEdgeNeighbours(
    ProcessInfo& rCurrentProcessInfo
    )
{
	std::vector<std::vector<Element*>> edge_neighbours_container;
	auto& r_nodes_current_element = this->GetGeometry();

	auto& pNode0 = r_nodes_current_element[0];
	auto& pNode1 = r_nodes_current_element[1];
	auto& pNode2 = r_nodes_current_element[2];
	auto& pNode3 = r_nodes_current_element[3];

	// Neighbour elements of each node of the current element
	GlobalPointersVector<Element>& r_neigh_node_0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
	GlobalPointersVector<Element>& r_neigh_node_1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
	GlobalPointersVector<Element>& r_neigh_node_2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
	GlobalPointersVector<Element>& r_neigh_node_3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);

	// Nodal neighbours container
	std::vector<GlobalPointersVector<Element>> nodal_neighbours;
	nodal_neighbours.push_back(r_neigh_node_0);
	nodal_neighbours.push_back(r_neigh_node_1);
	nodal_neighbours.push_back(r_neigh_node_2);
	nodal_neighbours.push_back(r_neigh_node_3);

	// Aux indexes
	Matrix nodes_indexes = ZeroMatrix(6, 2);
	this->SetNodeIndexes(nodes_indexes);

	// Loop over EDGES to assign the elements that share that edge -> Fill mEdgeNeighboursContainer
	for (unsigned int edge = 0; edge < 6; edge++) {
		const int node_index_1 = nodes_indexes(edge, 0);
		const int node_index_2 = nodes_indexes(edge, 1);

		// Neigh elements of local node 1 and 2  //
		auto& r_neigh_of_node_1 = nodal_neighbours[node_index_1];
		auto& r_neigh_of_node_2 = nodal_neighbours[node_index_2];

		const int node_id_1 = r_nodes_current_element[node_index_1].Id();
		const int node_id_2 = r_nodes_current_element[node_index_2].Id();

		std::vector<Element*> edge_shared_elements_node_1;
		// Loop over neigh elements of the node 1
		for (unsigned int neigh_elem = 0; neigh_elem < r_neigh_of_node_1.size(); neigh_elem++) {
			// Nodes of the neigh element
			auto& r_nodes_neigh_elem = r_neigh_of_node_1[neigh_elem].GetGeometry();

			// Loop over the nodes of the neigh element
			for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
				const int neigh_element_node_id = r_nodes_neigh_elem[neigh_elem_node].Id();

				if (neigh_element_node_id == node_id_2 && this->Id() != r_neigh_of_node_1[neigh_elem].Id()) {
					edge_shared_elements_node_1.push_back(&r_neigh_of_node_1[neigh_elem]); // ( [] returns an Element object!!)
				}
			}
		}

		std::vector<Element *> edge_shared_elements_node_2;
		// Loop over neigh elements of the node 2
		for (unsigned int neigh_elem = 0; neigh_elem < r_neigh_of_node_2.size(); neigh_elem++) {
			// Nodes of the neigh element
			auto &r_nodes_neigh_elem = r_neigh_of_node_2[neigh_elem].GetGeometry();

			// Loop over the nodes of the neigh element
			for (unsigned int neigh_elem_node = 0; neigh_elem_node < 4; neigh_elem_node++) {
				const int neigh_element_node_id = r_nodes_neigh_elem[neigh_elem_node].Id();

				if (neigh_element_node_id == node_id_1 && this->Id() != r_neigh_of_node_2[neigh_elem].Id()) {
					edge_shared_elements_node_2.push_back(&r_neigh_of_node_2[neigh_elem]);
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
		edge_neighbours_container.push_back(edge_shared_elements);
	} // End loop edges

	// Storages the information inside the element
	this->SaveEdgeNeighboursContainer(edge_neighbours_container);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<2,0>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressModifiedMohrCoulomb(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<3,0>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressModifiedMohrCoulomb(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<2,1>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressRankine(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressRankine(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<2,2>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressSimoJu(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressSimoJu(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<2,3>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressDruckerPrager(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressDruckerPrager(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<2,4>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressHuberVonMises(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressHuberVonMises(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<2,5>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressTresca(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::CalculateEquivalentStress(
    const array_1d<double,VoigtSize>& rPredictiveStressVector,
    const Vector& rStrainVector,
    double& rEquivalentStress,
    ConstitutiveLaw::Parameters& rValues
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateEquivalentStressTresca(rPredictiveStressVector, rStrainVector, rEquivalentStress, rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<2,0>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdModifiedMohrCoulomb(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<3,0>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdModifiedMohrCoulomb(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<2,1>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdRankine(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdRankine(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<2,2>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdSimoJu(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdSimoJu(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<2,3>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdDruckerPrager(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdDruckerPrager(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<2,4>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdHuberVonMises(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdHuberVonMises(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<2,5>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdTresca(rValues, rThreshold);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::GetInitialUniaxialThreshold(
    ConstitutiveLaw::Parameters& rValues,
    double& rThreshold
    )
{
    ConstitutiveLawUtilities<VoigtSize>::GetInitialUniaxialThresholdTresca(rValues, rThreshold);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<2,0>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterModifiedMohrCoulomb(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<3,0>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterModifiedMohrCoulomb(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<2,1>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterRankine(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterRankine(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<2,2>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterSimoJu(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterSimoJu(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<2,3>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterDruckerPrager(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterDruckerPrager(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<2,4>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterHuberVonMises(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterHuberVonMises(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<2,5>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterTresca(rValues, rAParameter, CharacteristicLength);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::CalculateDamageParameter(
    ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength
    )
{
    ConstitutiveLawUtilities<VoigtSize>::CalculateDamageParameterTresca(rValues, rAParameter, CharacteristicLength);
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
double GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateCharacteristicLength(
    GenericSmallStrainFemDemElement<TDim,TyieldSurf> *pCurrentElement
    )
{
    return pCurrentElement->GetGeometry().Length();
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Vector& GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateVolumeForce(
    Vector& rVolumeForce, 
    const Vector& rN
    )
{
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
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateExponentialDamage(
	double& rDamage,
	const double DamageParameter,
	const double UniaxialStress,
	const double InitialThrehsold
	)
{
	rDamage = 1.0 - (InitialThrehsold / UniaxialStress) * std::exp(DamageParameter *
			 (1.0 - UniaxialStress / InitialThrehsold)); // Exponential softening law
	if (rDamage > 0.999) rDamage = 0.999;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateTangentTensor(
	Matrix& rTangentTensor,
	const Vector& rStrainVectorGP,
	const Vector& rStressVectorGP,
	const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
	)
{
	const double number_components = rStrainVectorGP.size();
	rTangentTensor.resize(number_components, number_components);
	Vector perturbed_stress, perturbed_strain;
	perturbed_strain.resize(number_components);
	perturbed_stress.resize(number_components);
	
	for (unsigned int component = 0; component < number_components; component++) {
		double perturbation;
		this->CalculatePerturbation(rStrainVectorGP, perturbation, component);
		this->PerturbateStrainVector(perturbed_strain, rStrainVectorGP, perturbation, component);
		this->IntegratePerturbedStrain(perturbed_stress, perturbed_strain, rElasticMatrix, rValues);
		const Vector& r_delta_stress = perturbed_stress - rStressVectorGP;
		this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress, perturbation, component);
	}
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculatePerturbation(
	const Vector& rStrainVectorGP,
	double& rPerturbation,
	const int Component
	)
{
	double perturbation_1, perturbation_2;
	if (std::abs(rStrainVectorGP[Component]) > tolerance) {
		perturbation_1 = 1.0e-5 * rStrainVectorGP[Component];
	} else {
		double min_strain_component = ConstitutiveLawUtilities<VoigtSize>::GetMinAbsValue(rStrainVectorGP);
		perturbation_1 = 1.0e-5 * min_strain_component;
	}
	const double max_strain_component = ConstitutiveLawUtilities<VoigtSize>::GetMaxAbsValue(rStrainVectorGP);
	perturbation_2 = 1.0e-10 * max_strain_component;
	rPerturbation = std::max(perturbation_1, perturbation_2);
	if (rPerturbation < 1e-8) rPerturbation = 1e-8;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::PerturbateStrainVector(
	Vector& rPerturbedStrainVector,
	const Vector& rStrainVectorGP,
	const double Perturbation,
	const int Component
	)
{
    noalias(rPerturbedStrainVector) = rStrainVectorGP;
    rPerturbedStrainVector[Component] += Perturbation;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::IntegratePerturbedStrain(
	Vector& rPerturbedStressVector,
	const Vector& rPerturbedStrainVector,
	const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
	)
{
	const Vector& r_perturbed_predictive_stress = prod(rElasticMatrix, rPerturbedStrainVector);
	Vector damages_edges = ZeroVector(NumberOfEdges);
	const double characteristic_length = this->CalculateCharacteristicLength(this);
    bool dummy = false;

    Vector average_stress_edge = r_perturbed_predictive_stress;
    Vector average_strain_edge = rPerturbedStrainVector;

	for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
        this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
        this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

        double damage_edge = mDamages[edge];
        double threshold = mThresholds[edge];
        
        this->IntegrateStressDamageMechanics(threshold, damage_edge, average_strain_edge, 
            average_stress_edge, edge, characteristic_length, rValues, dummy);

        damages_edges[edge] = damage_edge;
	} // Loop edges
	const double damage_element = this->CalculateElementalDamage(damages_edges);
	rPerturbedStressVector = (1.0 - damage_element) * r_perturbed_predictive_stress;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::AssignComponentsToTangentTensor(
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

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::GetValueOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rValues,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT || rVariable == IS_DAMAGED || rVariable == STRESS_THRESHOLD || rVariable == EQUIVALENT_STRESS_VM) {
		CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
	}
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::GetValueOnIntegrationPoints(
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

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::GetValueOnIntegrationPoints(
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

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
	const Variable<double> &rVariable,
	std::vector<double> &rOutput,
	const ProcessInfo &rCurrentProcessInfo)
{
	if (rVariable == DAMAGE_ELEMENT) {
		rOutput.resize(1);
		for (unsigned int point_number = 0; point_number < 1; point_number++) {
			rOutput[point_number] = mDamage;
		}
	} else if (rVariable == IS_DAMAGED) {
		rOutput.resize(1);
		for (unsigned int point_number = 0; point_number < 1; point_number++) {
			rOutput[point_number] = double(this->GetValue(IS_DAMAGED));
		}
	} else if (rVariable == STRESS_THRESHOLD) {
		rOutput.resize(1);
		for (unsigned int point_number = 0; point_number < 1; point_number++) {
			rOutput[point_number] = mThreshold;
		}
	}
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
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

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
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

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainFemDemElement<2,0>;
template class GenericSmallStrainFemDemElement<2,1>;
template class GenericSmallStrainFemDemElement<2,2>;
template class GenericSmallStrainFemDemElement<2,3>;
template class GenericSmallStrainFemDemElement<2,4>;
template class GenericSmallStrainFemDemElement<2,5>;
template class GenericSmallStrainFemDemElement<3,0>;
template class GenericSmallStrainFemDemElement<3,1>;
template class GenericSmallStrainFemDemElement<3,2>;
template class GenericSmallStrainFemDemElement<3,3>;
template class GenericSmallStrainFemDemElement<3,4>;
template class GenericSmallStrainFemDemElement<3,5>;
} // namespace Kratos