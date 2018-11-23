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

#include "custom_elements/femdem3d_large_displacement_element.hpp"
#include "utilities/geometry_utilities.h"
#include "fem_to_dem_application_variables.h"


namespace Kratos
{
//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::FemDem3DLargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry)
	: FemDem3DElement(NewId, pGeometry)
{
	//DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::FemDem3DLargeDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
	: FemDem3DElement(NewId, pGeometry, pProperties)
{
	//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
	mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::FemDem3DLargeDisplacementElement(FemDem3DLargeDisplacementElement const &rOther)
	: FemDem3DElement(rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
	//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

FemDem3DLargeDisplacementElement &FemDem3DLargeDisplacementElement::operator=(FemDem3DLargeDisplacementElement const &rOther)
{
	//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

	FemDem3DElement::operator=(rOther);
	return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer FemDem3DLargeDisplacementElement::Create(IndexType NewId, NodesArrayType const &rThisNodes, PropertiesType::Pointer pProperties) const
{
	//NEEDED TO CREATE AN ELEMENT
	return Element::Pointer(new FemDem3DLargeDisplacementElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer FemDem3DLargeDisplacementElement::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{

	//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
	//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

	FemDem3DLargeDisplacementElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

	return Element::Pointer(new FemDem3DLargeDisplacementElement(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

FemDem3DLargeDisplacementElement::~FemDem3DLargeDisplacementElement()
{
}

void FemDem3DLargeDisplacementElement::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;

    const SizeType mat_size = number_of_nodes * dimension;
    B.resize(strain_size, dimension * number_of_nodes);

    Matrix constitutive_matrix = ZeroMatrix(strain_size, strain_size);
    const double E = this->GetProperties()[YOUNG_MODULUS];
    const double nu = this->GetProperties()[POISSON_RATIO];
    this->CalculateConstitutiveMatrix(constitutive_matrix, E, nu);

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    Matrix DeltaPosition(number_of_nodes, dimension);
    noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
    DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

    // Loop over Gauss Points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {

        J = this->GetGeometry().Jacobian(J, point_number, mThisIntegrationMethod);
        const double detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, mThisIntegrationMethod);

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);
        
        Vector stress_vector, strain_vector;
        stress_vector.resize(strain_size);
        strain_vector.resize(strain_size);
        this->CalculateGreenLagrangeStrainVector(strain_vector, F);

        this->SetValue(STRAIN_VECTOR, strain_vector);
        
        // S = C:E -> Assume small deformations
        this->CalculateStressVectorPredictor(stress_vector, constitutive_matrix, strain_vector);
        this->SetValue(STRESS_VECTOR, stress_vector);
    }
}

void FemDem3DLargeDisplacementElement::CalculateLocalSystem(
	MatrixType &rLeftHandSideMatrix,
	VectorType &rRightHandSideVector,
	ProcessInfo &rCurrentProcessInfo)
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size = GetStrainSize();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    double detJ0;

    const SizeType mat_size = number_of_nodes * dimension;
    B.resize(strain_size, dimension * number_of_nodes);

    Matrix constitutive_matrix = ZeroMatrix(strain_size, strain_size);
    const double E = this->GetProperties()[YOUNG_MODULUS];
    const double nu = this->GetProperties()[POISSON_RATIO];
    this->CalculateConstitutiveMatrix(constitutive_matrix, E, nu);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);

    rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    Matrix DeltaPosition(number_of_nodes, dimension);
    noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
    DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

    // Loop over Gauss Points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

		J = this->GetGeometry().Jacobian(J, point_number, mThisIntegrationMethod);
        detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, mThisIntegrationMethod);

        const double integration_weigth = integration_points[point_number].Weight() * detJ0;
        const Matrix &Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
		Vector N = row(Ncontainer, point_number);

        Vector VolumeForce = ZeroVector(dimension);
		VolumeForce = this->CalculateVolumeForce(VolumeForce, N);
		// Taking into account Volume Force into de RHS
		for (unsigned int i = 0; i < number_of_nodes; i++) {
			int index = dimension * i;
			for (unsigned int j = 0; j < dimension; j++) {
				rRightHandSideVector[index + j] += integration_weigth * N[i] * VolumeForce[j];
			}
		}

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);

        Vector IntegratedStressVector = ZeroVector(strain_size);
		Vector DamagesOnEdges = ZeroVector(6);

		// Loop over edges of the element
		for (unsigned int edge = 0; edge < 6; edge++) {
			std::vector<Element *> EdgeNeighbours = this->GetEdgeNeighbourElements(edge);
			Vector AverageStressVector, AverageStrainVector, IntegratedStressVectorOnEdge;

			this->CalculateAverageStressOnEdge(AverageStressVector, EdgeNeighbours);
			this->CalculateAverageStrainOnEdge(AverageStrainVector, EdgeNeighbours);

			double DamageEdge;
			const double Lchar = this->Get_l_char(edge);
			this->IntegrateStressDamageMechanics(IntegratedStressVectorOnEdge, DamageEdge,
												 AverageStrainVector, AverageStressVector, edge, Lchar);

			this->SetNonConvergedDamages(DamageEdge, edge);
			DamagesOnEdges[edge] = DamageEdge;
		} // End loop over edges

        double damage_element = this->CalculateElementalDamage(DamagesOnEdges);
		if (damage_element >= 0.999)
			damage_element = 0.999;
		this->SetNonConvergedDamages(damage_element);

		const Vector &stress_vector = this->GetValue(STRESS_VECTOR);
        Vector integrated_stress_vector = ZeroVector(strain_size);
		integrated_stress_vector = (1.0 - damage_element) * stress_vector;
		this->SetIntegratedStressVector(integrated_stress_vector);

        this->CalculateAndAddMaterialK(rLeftHandSideMatrix, B, constitutive_matrix, integration_weigth);
        this->CalculateGeometricK(rLeftHandSideMatrix, DN_DX, integrated_stress_vector, integration_weigth);
        this->CalculateAndAddInternalForcesVector(rRightHandSideVector, B, integrated_stress_vector, integration_weigth);
    }
} // CalculateLocalSystem

double FemDem3DLargeDisplacementElement::CalculateDerivativesOnReferenceConfiguration(
    Matrix &rJ0,
    Matrix &rInvJ0,
    Matrix &rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod)
{
    GeometryType &r_geom = GetGeometry();
    GeometryUtils::JacobianOnInitialConfiguration(r_geom, r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);

    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);

    const Matrix &rDN_De = GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);

    return detJ0;
}

void FemDem3DLargeDisplacementElement::FinalizeNonLinearIteration(ProcessInfo &CurrentProcessInfo)
{
}

void FemDem3DLargeDisplacementElement::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    KRATOS_TRY

	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (int i = 0; i < number_of_nodes; ++i) {
        const auto index = dimension * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
    }
    KRATOS_CATCH("")
}

void FemDem3DLargeDisplacementElement::CalculateAndAddMaterialK(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    )
{
    KRATOS_TRY

    double damage = this->GetNonConvergedDamage();

    // Secant Constitutive Tensor
    noalias(rLeftHandSideMatrix) += (1.0 - damage) * IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));

    KRATOS_CATCH("")
}

void FemDem3DLargeDisplacementElement::CalculateGeometricK(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& DN_DX,
    const Vector& StressVector,
    const double IntegrationWeight
    )
{
    KRATOS_TRY

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(StressVector);
    Matrix reduced_Kg = prod(DN_DX, IntegrationWeight * Matrix(prod(stress_tensor, trans(DN_DX))));
    MathUtils<double>::ExpandAndAddReducedMatrix(rLeftHandSideMatrix, reduced_Kg, dimension);

    KRATOS_CATCH("")
}

void FemDem3DLargeDisplacementElement::CalculateGreenLagrangeStrainVector(
    Vector &rStrainVector,
    const Matrix &rF
    )
{
    KRATOS_TRY

    Matrix strain_tensor;
    strain_tensor.resize(3, 3);
    Matrix identity = identity_matrix<double>(3);
    noalias(strain_tensor) = 0.5 * (prod(trans(rF), rF) - identity);
    rStrainVector = MathUtils<double>::StrainTensorToVector(strain_tensor, rStrainVector.size());

    KRATOS_CATCH("")
}

void FemDem3DLargeDisplacementElement::CalculateStressVectorPredictor(
    Vector& rStressVector, 
    const Matrix& rConstitutiveMAtrix, 
    const Vector& rStrainVector)
{
    noalias(rStressVector) = prod(rConstitutiveMAtrix, rStrainVector);
}

void FemDem3DLargeDisplacementElement::CalculateAndAddInternalForcesVector(
    Vector& rRightHandSideVector, 
    const Matrix& rB, 
    const Vector& rStressVector, 
    const double IntegrationWeight
    )
{
    noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(rB), rStressVector);
}


} // namespace Kratos