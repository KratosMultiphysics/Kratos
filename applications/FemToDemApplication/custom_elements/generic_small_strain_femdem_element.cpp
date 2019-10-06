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
#include "custom_elements/generic_small_strain_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry)
{
    // DO NOT ADD DOFS HERE!!!
    if (this->mThresholds.size() != NumberOfEdges)
        this->mThresholds.resize(NumberOfEdges);
    noalias(this->mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (this->mDamages.size() != NumberOfEdges)
        this->mDamages.resize(NumberOfEdges);
    noalias(this->mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim, TyieldSurf>::GenericSmallStrainFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry, pProperties)
{
    // DO NOT ADD DOFS HERE!!!
    if (this->mThresholds.size() != NumberOfEdges)
        this->mThresholds.resize(NumberOfEdges);
    noalias(this->mThresholds) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (this->mDamages.size() != NumberOfEdges)
        this->mDamages.resize(NumberOfEdges);
    noalias(this->mDamages) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericSmallStrainFemDemElement>( NewId, this->GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericSmallStrainFemDemElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
GenericSmallStrainFemDemElement<TDim,TyieldSurf>::~GenericSmallStrainFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Pointer p_new_elem = Kratos::make_intrusive<GenericSmallStrainFemDemElement>(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::InitializeNonLinearIteration(
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = this->mConstitutiveLawVector[0]->GetStrainSize();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, this->UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // Reading integration points
    const Properties& r_properties = this->GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(this->mThisIntegrationMethod);

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->mThisIntegrationMethod);

    for ( IndexType point_number = 0; point_number < this->mConstitutiveLawVector.size(); ++point_number ) {
        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->mThisIntegrationMethod);

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this->GetStressMeasure());
    }
    this->SetValue(STRESS_VECTOR, this_constitutive_variables.StressVector);
    this->SetValue(STRAIN_VECTOR, this_constitutive_variables.StrainVector);


    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    // todo
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
    // TODO
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
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    const auto& r_geometry = this->GetGeometry();

    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    rThisKinematicVariables.detJ0 = this->CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    this->CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.DN_DX);

    // Compute equivalent F
    this->GetValuesVector(rThisKinematicVariables.Displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    // ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    // rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

// template<>
// void GenericSmallStrainFemDemElement<2,0>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate2DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<3,0>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }
template<>
void GenericSmallStrainFemDemElement<2,1>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
// template<>
// void GenericSmallStrainFemDemElement<3,1>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<2,2>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate2DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<3,2>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<2,3>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate2DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<3,3>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<2,4>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate2DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<3,4>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<2,5>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate2DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<3,5>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<2,6>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate2DB(rB, rDN_DX);
// }
// template<>
// void GenericSmallStrainFemDemElement<3,6>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
// {
//     this->Calculate3DB(rB, rDN_DX);
// }

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Calculate2DB(
    Matrix& rB, 
    const Matrix& rDN_DX
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for (SizeType i = 0; i < number_of_nodes; ++i ) {
        rB(0, i*2    ) = rDN_DX(i, 0);
        rB(1, i*2 + 1) = rDN_DX(i, 1);
        rB(2, i*2    ) = rDN_DX(i, 1);
        rB(2, i*2 + 1) = rDN_DX(i, 0);
    }

    KRATOS_CATCH( "" )
}
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Calculate3DB(
    Matrix& rB, 
    const Matrix& rDN_DX
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for ( SizeType i = 0; i < number_of_nodes; ++i ) {
        rB(0, i*3    ) = rDN_DX(i, 0);
        rB(1, i*3 + 1) = rDN_DX(i, 1);
        rB(2, i*3 + 2) = rDN_DX(i, 2);
        rB(3, i*3    ) = rDN_DX(i, 1);
        rB(3, i*3 + 1) = rDN_DX(i, 0);
        rB(4, i*3 + 1) = rDN_DX(i, 2);
        rB(4, i*3 + 2) = rDN_DX(i, 1);
        rB(5, i*3    ) = rDN_DX(i, 2);
        rB(5, i*3 + 2) = rDN_DX(i, 0);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

// template<unsigned int TDim, unsigned int TyieldSurf>
// void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::ComputeEquivalentF(
//     Matrix& rF,
//     const Vector& rStrainTensor
//     )
// {
//     const SizeType dim = GetGeometry().WorkingSpaceDimension();

//     if(dim == 2) {
//         rF(0,0) = 1.0+rStrainTensor(0);
//         rF(0,1) = 0.5*rStrainTensor(2);
//         rF(1,0) = 0.5*rStrainTensor(2);
//         rF(1,1) = 1.0+rStrainTensor(1);
//     } else {
//         rF(0,0) = 1.0+rStrainTensor(0);
//         rF(0,1) = 0.5*rStrainTensor(3);
//         rF(0,2) = 0.5*rStrainTensor(5);
//         rF(1,0) = 0.5*rStrainTensor(3);
//         rF(1,1) = 1.0+rStrainTensor(1);
//         rF(1,2) = 0.5*rStrainTensor(4);
//         rF(2,0) = 0.5*rStrainTensor(5);
//         rF(2,1) = 0.5*rStrainTensor(4);
//         rF(2,2) = 1.0+rStrainTensor(2);
//     }
// }

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    this->mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Displacements vector
    Vector displacements(mat_size);
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& displacement = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k) {
            rValues[index + k] = displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/
// template class GenericSmallStrainFemDemElement<2,0>;
template class GenericSmallStrainFemDemElement<2,1>;
// template class GenericSmallStrainFemDemElement<2,2>;
// template class GenericSmallStrainFemDemElement<2,3>;
// template class GenericSmallStrainFemDemElement<2,4>;
// template class GenericSmallStrainFemDemElement<2,5>;
// template class GenericSmallStrainFemDemElement<2,6>;
// template class GenericSmallStrainFemDemElement<3,0>;
// template class GenericSmallStrainFemDemElement<3,1>;
// template class GenericSmallStrainFemDemElement<3,2>;
// template class GenericSmallStrainFemDemElement<3,3>;
// template class GenericSmallStrainFemDemElement<3,4>;
// template class GenericSmallStrainFemDemElement<3,5>;
// template class GenericSmallStrainFemDemElement<3,6>;
} // namespace Kratos