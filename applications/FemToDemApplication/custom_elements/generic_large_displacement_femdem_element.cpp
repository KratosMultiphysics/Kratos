//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//

// System includes

// External includes

// Project includes
#include "custom_elements/generic_large_displacement_femdem_element.hpp"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GenericSmallStrainFemDemElement<TDim, TyieldSurf>(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::GenericLargeDisplacementFemDemElement(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    )
    : GenericSmallStrainFemDemElement<TDim, TyieldSurf>(NewId, pGeometry, pProperties)
{
    // BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
    this->mThisIntegrationMethod = this->GetGeometry().GetDefaultIntegrationMethod();

    // Each component == Each edge
    if (this->mNonConvergedThresholds.size() != NumberOfEdges)
        this->mNonConvergedThresholds.resize(NumberOfEdges);
    noalias(this->mNonConvergedThresholds) = ZeroVector(NumberOfEdges);   // Equivalent stress

    if (this->mThresholds.size() != NumberOfEdges)
        this->mThresholds.resize(NumberOfEdges);
    noalias(this->mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (this->mDamages.size() != NumberOfEdges)
        this->mDamages.resize(NumberOfEdges);
    noalias(this->mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge

    if (this->mNonConvergedDamages.size() != NumberOfEdges)
        this->mNonConvergedDamages.resize(NumberOfEdges);
    noalias(this->mNonConvergedDamages) = ZeroVector(NumberOfEdges); // mDamages on edges of "i" iteration
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::GenericLargeDisplacementFemDemElement(
    GenericLargeDisplacementFemDemElement const& rOther
    )
    : GenericSmallStrainFemDemElement<TDim,TyieldSurf>(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>& GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>::operator=(
    GenericLargeDisplacementFemDemElement<TDim, TyieldSurf> const& rOther
    )
{
    GenericSmallStrainFemDemElement<TDim, TyieldSurf>::operator=(rOther);
    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Create(
    IndexType NewId, 
    NodesArrayType const& rThisNodes, 
    PropertiesType::Pointer pProperties
    ) const 
{
    return Element::Pointer(new GenericLargeDisplacementFemDemElement(NewId, this->GetGeometry().Create(rThisNodes), pProperties));
}

//************************************CLONE*******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
    GenericLargeDisplacementFemDemElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    return Element::Pointer(new GenericLargeDisplacementFemDemElement<TDim, TyieldSurf>(NewElement));
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::~GenericLargeDisplacementFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    B.resize(VoigtSize, dimension * number_of_nodes);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    //set constitutive law flags:
    Flags& r_constitutive_law_options = values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //create and initialize element variables:
    ElementDataType variables;
    this->InitializeElementData(variables, rCurrentProcessInfo);
    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Loop over Gauss Points
    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        //calculate the elastic matrix
        this->CalculateMaterialResponse(variables, values, point_number);
        Matrix& r_constitutive_matrix = values.GetConstitutiveMatrix();
        J = this->GetGeometry().Jacobian(J, point_number, this->mThisIntegrationMethod);
        const double detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, this->mThisIntegrationMethod);

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);
        
        Vector stress_vector, strain_vector;
        stress_vector.resize(VoigtSize);
        strain_vector.resize(VoigtSize);
        this->CalculateGreenLagrangeStrainVector(strain_vector, F);

        this->SetValue(STRAIN_VECTOR, strain_vector);
        
        // S = C:E -> Assume small deformations
        this->CalculateStressVectorPredictor(stress_vector, r_constitutive_matrix, strain_vector);
        this->SetValue(STRESS_VECTOR, stress_vector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    double detJ0;

    const SizeType mat_size = number_of_nodes * dimension;
    B.resize(VoigtSize, dimension * number_of_nodes);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    //set constitutive law flags:
    Flags& r_constitutive_law_options = values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //create and initialize element variables:
    ElementDataType variables;
    this->InitializeElementData(variables, rCurrentProcessInfo);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);
    rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const double characteristic_length = this->CalculateCharacteristicLength(this);

    // Loop over Gauss Points
    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
        //calculate the elastic matrix
        this->CalculateMaterialResponse(variables, values, point_number);
        Matrix& r_constitutive_matrix = values.GetConstitutiveMatrix();

        J = this->GetGeometry().Jacobian(J, point_number, this->mThisIntegrationMethod);
        detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, this->mThisIntegrationMethod);

        double integration_weigth = r_integration_points[point_number].Weight() * detJ0;
        integration_weigth = this->CalculateIntegrationWeight(integration_weigth);
        const Matrix& Ncontainer = this->GetGeometry().ShapeFunctionsValues(this->mThisIntegrationMethod);
        Vector N = row(Ncontainer, point_number);

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);
        Vector strain_vector;
        this->CalculateGreenLagrangeStrainVector(strain_vector, F);
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

                double damage_edge = this->mDamages[edge];
                double threshold = this->mThresholds[edge];
                
                this->IntegrateStressDamageMechanics(threshold, damage_edge, average_strain_edge, 
                    average_stress_edge, edge, characteristic_length, values, is_damaging);

                this->mNonConvergedDamages[edge] = damage_edge;
                this->mNonConvergedThresholds[edge] = threshold;
            } // Loop over edges
        }  else {
            noalias(this->mNonConvergedDamages) = ZeroVector(NumberOfEdges);
        }

        const double damage_element = this->CalculateElementalDamage(this->mNonConvergedDamages);

        Vector stress_vector = ZeroVector(VoigtSize);
        this->CalculateStressVectorPredictor(stress_vector, r_constitutive_matrix, strain_vector);
        const Vector& r_integrated_stress_vector = (1.0 - damage_element) * stress_vector;
        Matrix tangent_tensor;

        if (is_damaging == true) { // Tangent Tensor
            this->CalculateTangentTensor(tangent_tensor, strain_vector, r_integrated_stress_vector, F, r_constitutive_matrix, values);
            rLeftHandSideMatrix += integration_weigth * prod(trans(B), Matrix(prod(tangent_tensor, B)));
        } else { // Secant
            this->CalculateAndAddMaterialK(rLeftHandSideMatrix, B, r_constitutive_matrix, integration_weigth, damage_element);
        }

        this->CalculateGeometricK(rLeftHandSideMatrix, DN_DX, r_integrated_stress_vector, integration_weigth);
        Vector volume_force = ZeroVector(dimension);
        volume_force = this->CalculateVolumeForce(volume_force, N);
        this->CalculateAndAddExternalForces(rRightHandSideVector, variables, volume_force, integration_weigth);
        this->CalculateAndAddInternalForcesVector(rRightHandSideVector, B, r_integrated_stress_vector, integration_weigth);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, 
    ProcessInfo& rCurrentProcessInfo
    )
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    double detJ0;

    const SizeType mat_size = number_of_nodes * dimension;
    B.resize(VoigtSize, dimension * number_of_nodes);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    //set constitutive law flags:
    Flags& r_constitutive_law_options = values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //create and initialize element variables:
    ElementDataType variables;
    this->InitializeElementData(variables, rCurrentProcessInfo);

    if (rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    // Loop over Gauss Points
    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number ) {
        //calculate the elastic matrix
        this->CalculateMaterialResponse(variables, values, point_number);
        Matrix& r_constitutive_matrix = values.GetConstitutiveMatrix();

        J = this->GetGeometry().Jacobian(J, point_number, this->mThisIntegrationMethod);
        detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, this->mThisIntegrationMethod);

        double integration_weigth = r_integration_points[point_number].Weight() * detJ0;
        integration_weigth = this->CalculateIntegrationWeight(integration_weigth);
        const Matrix& r_Ncontainer = this->GetGeometry().ShapeFunctionsValues(this->mThisIntegrationMethod);
        Vector N = row(r_Ncontainer, point_number);

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);

        const double damage_element = this->CalculateElementalDamage(this->mDamages);

        const Vector& r_stress_vector = this->GetValue(STRESS_VECTOR);
        const Vector& r_integrated_stress_vector = (1.0 - damage_element) * r_stress_vector;

        this->CalculateAndAddMaterialK(rLeftHandSideMatrix, B, r_constitutive_matrix, integration_weigth, damage_element);
        this->CalculateGeometricK(rLeftHandSideMatrix, DN_DX, r_integrated_stress_vector, integration_weigth);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    // Kinematic variables
    Matrix B, F, DN_DX, InvJ0, J, J0;
    double detJ0;

    const SizeType mat_size = number_of_nodes * dimension;
    B.resize(VoigtSize, dimension * number_of_nodes);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    //set constitutive law flags:
    Flags& r_constitutive_law_options = values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //create and initialize element variables:
    ElementDataType variables;
    this->InitializeElementData(variables, rCurrentProcessInfo);

    // Resizing as needed the RHS
    if (rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size, false);
    rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const double characteristic_length = this->CalculateCharacteristicLength(this);

    // Loop over Gauss Points
    for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
        //calculate the elastic matrix
        this->CalculateMaterialResponse(variables, values, point_number);
        Matrix& r_constitutive_matrix = values.GetConstitutiveMatrix();

        J = this->GetGeometry().Jacobian(J, point_number, this->mThisIntegrationMethod);
        detJ0 = this->CalculateDerivativesOnReferenceConfiguration(J0, InvJ0, DN_DX, point_number, this->mThisIntegrationMethod);

        double integration_weigth = r_integration_points[point_number].Weight() * detJ0;
        integration_weigth = this->CalculateIntegrationWeight(integration_weigth);
        const Matrix &Ncontainer = this->GetGeometry().ShapeFunctionsValues(this->mThisIntegrationMethod);
        Vector N = row(Ncontainer, point_number);

        GeometryUtils::DeformationGradient(J, InvJ0, F);
        this->CalculateB(B, F, DN_DX);
        Vector strain_vector;
        this->CalculateGreenLagrangeStrainVector(strain_vector, F);
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

                double damage_edge = this->mDamages[edge];
                double threshold = this->mThresholds[edge];
                
                this->IntegrateStressDamageMechanics(threshold, damage_edge, average_strain_edge, 
                    average_stress_edge, edge, characteristic_length, values, is_damaging);

                damages[edge] = damage_edge;
            } // Loop over edges
        }  else {
            noalias(damages) = ZeroVector(NumberOfEdges);
        }
        const double damage_element = this->CalculateElementalDamage(damages);

        Vector stress_vector = ZeroVector(VoigtSize);
        this->CalculateStressVectorPredictor(stress_vector, r_constitutive_matrix, strain_vector);
        const Vector& r_integrated_stress_vector = (1.0 - damage_element) * stress_vector;

        this->SetValue(STRESS_VECTOR, stress_vector);
        this->SetValue(STRESS_VECTOR_INTEGRATED, r_integrated_stress_vector);
        this->SetValue(STRAIN_VECTOR, strain_vector);

        Vector volume_force = ZeroVector(dimension);
        volume_force = this->CalculateVolumeForce(volume_force, N);
        this->CalculateAndAddExternalForces(rRightHandSideVector, variables, volume_force, integration_weigth);
        this->CalculateAndAddInternalForcesVector(rRightHandSideVector, B, r_integrated_stress_vector, integration_weigth);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericLargeDisplacementFemDemElement<2,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,0>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,1>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,2>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,3>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,4>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,5>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<2,6>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB2D(rB, rF, rDN_DX);
}
template<>
void GenericLargeDisplacementFemDemElement<3,6>::CalculateB(Matrix& rB, const Matrix& rF, const Matrix& rDN_DX)
{
    this->CalculateB3D(rB, rF, rDN_DX);
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateB2D(
    Matrix& rB, 
    const Matrix& rF, 
    const Matrix& rDN_DX
    )
{
    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < number_of_nodes; i++) {
        unsigned int index = 2 * i;
        rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
        rB(2, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
        rB(2, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);

    }
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateB3D(
    Matrix& rB, 
    const Matrix& rF, 
    const Matrix& rDN_DX
    )
{
    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < number_of_nodes; ++i) {
        const auto index = 3 * i;
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
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateGreenLagrangeStrainVector(
    Vector& rStrainVector,
    const Matrix& rF
    )
{
    Matrix strain_tensor;
    if (strain_tensor.size1() != TDim)
        strain_tensor.resize(TDim, TDim);
    
    if (rStrainVector.size() != VoigtSize)
        rStrainVector.resize(VoigtSize);

    Matrix identity = identity_matrix<double>(TDim);
    noalias(strain_tensor) = 0.5 * (prod(trans(rF), rF) - identity);
    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(strain_tensor, rStrainVector.size());
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateStressVectorPredictor(
    Vector& rStressVector, 
    const Matrix& rConstitutiveMAtrix, 
    const Vector& rStrainVector)
{
    if (rStressVector.size() != VoigtSize)
        rStressVector.resize(VoigtSize);
    noalias(rStressVector) = prod(rConstitutiveMAtrix, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateAndAddInternalForcesVector(
    Vector& rRightHandSideVector, 
    const Matrix& rB, 
    const Vector& rStressVector, 
    const double IntegrationWeight
    )
{
    noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(rB), rStressVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateGeometricK(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& rDN_DX,
    const Vector& rStressVector,
    const double IntegrationWeight
    )
{
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(rStressVector);
    Matrix reduced_Kg = prod(rDN_DX, IntegrationWeight * Matrix(prod(stress_tensor, trans(rDN_DX))));
    MathUtils<double>::ExpandAndAddReducedMatrix(rLeftHandSideMatrix, reduced_Kg, dimension);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateAndAddMaterialK(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight,
    const double Damage
    )
{
    // Secant Constitutive Tensor
    noalias(rLeftHandSideMatrix) += (1.0 - Damage) * IntegrationWeight * prod(trans(B), Matrix(prod(D, B)));
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
double GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateDerivativesOnReferenceConfiguration(
    Matrix& rJ0,
    Matrix& rInvJ0,
    Matrix& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    )
{
    GeometryType& r_geom = this->GetGeometry();
    GeometryUtils::JacobianOnInitialConfiguration(r_geom, r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    const Matrix& rDN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    return detJ0;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void  GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateTangentTensor(
    Matrix& rTangentTensor,
    const Vector& rStrainVectorGP,
    const Vector& rStressVectorGP,
    const Matrix& rDeformationGradientGP,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double number_components = rStrainVectorGP.size();
    rTangentTensor.resize(number_components, number_components);
    Vector perturbed_stress, perturbed_strain;
    perturbed_strain.resize(number_components);
    perturbed_stress.resize(number_components);
    Matrix perturbed_deformation_gradient;
    perturbed_deformation_gradient.resize(number_components, number_components);
    const double size_1 = rDeformationGradientGP.size1();
    const double size_2 = rDeformationGradientGP.size2();
    
    for (unsigned int i_component = 0; i_component < size_1; i_component++) {
        for (unsigned int j_component = i_component; j_component < size_2; j_component++) {
            double perturbation;
            const int component_voigt_index = this->CalculateVoigtIndex(number_components, i_component, j_component);
            this->CalculatePerturbation(rStrainVectorGP, perturbation, component_voigt_index);
            this->PerturbateDeformationGradient(perturbed_deformation_gradient, rDeformationGradientGP, perturbation, i_component, j_component);
            this->CalculateGreenLagrangeStrainVector(perturbed_strain, perturbed_deformation_gradient);
            this->IntegratePerturbedStrain(perturbed_stress, perturbed_strain, rElasticMatrix, rValues);
            const Vector& r_delta_stress = perturbed_stress - rStressVectorGP;
            this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress, perturbation, component_voigt_index);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::PerturbateDeformationGradient(
    Matrix& rPerturbedDeformationGradient,
    const Matrix& rDeformationGradientGP,
    const double Perturbation,
    const int ComponentI,
    const int ComponentJ
    )
{
    rPerturbedDeformationGradient = rDeformationGradientGP;
    if (ComponentI == ComponentJ) {
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += Perturbation;
    } else {
        rPerturbedDeformationGradient(ComponentI, ComponentJ) += 0.5 * Perturbation;
        rPerturbedDeformationGradient(ComponentJ, ComponentI) += 0.5 * Perturbation;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
int GenericLargeDisplacementFemDemElement<TDim,TyieldSurf>::CalculateVoigtIndex(
    const SizeType VoigtSize,
    const int ComponentI,
    const int ComponentJ
    )
{
    if (VoigtSize == 6) {
        switch(ComponentI) {
            case 0:
                switch(ComponentJ) {
                    case 0:
                        return 0;
                    case 1:
                        return 3;
                    case 2:
                        return 5;
                    default:
                        return 0;
                }
            case 1:
                switch(ComponentJ) {
                    case 0:
                        return 3;
                    case 1:
                        return 1;
                    case 2:
                        return 4;
                    default:
                        return 0;
                }
            case 2:
                switch(ComponentJ) {
                    case 0:
                        return 5;
                    case 1:
                        return 4;
                    case 2:
                        return 2;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    } else {
        switch(ComponentI) {
            case 0:
                switch(ComponentJ) {
                    case 0:
                        return 0;
                    case 1:
                        return 2;
                    default:
                        return 0;
                }
            case 1:
                switch(ComponentJ) {
                    case 0:
                        return 2;
                    case 1:
                        return 1;
                    default:
                        return 0;
                }
            default:
                return 0;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericLargeDisplacementFemDemElement<2,0>;
template class GenericLargeDisplacementFemDemElement<2,1>;
template class GenericLargeDisplacementFemDemElement<2,2>;
template class GenericLargeDisplacementFemDemElement<2,3>;
template class GenericLargeDisplacementFemDemElement<2,4>;
template class GenericLargeDisplacementFemDemElement<2,5>;
template class GenericLargeDisplacementFemDemElement<2,6>;
template class GenericLargeDisplacementFemDemElement<3,0>;
template class GenericLargeDisplacementFemDemElement<3,1>;
template class GenericLargeDisplacementFemDemElement<3,2>;
template class GenericLargeDisplacementFemDemElement<3,3>;
template class GenericLargeDisplacementFemDemElement<3,4>;
template class GenericLargeDisplacementFemDemElement<3,5>;
template class GenericLargeDisplacementFemDemElement<3,6>;

} // namespace Kratos