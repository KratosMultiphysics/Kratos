//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
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
    const auto& r_geometry         = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension       = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size     = this->mConstitutiveLawVector[0]->GetStrainSize();

    BaseSolidElement::KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    BaseSolidElement::ConstitutiveVariables this_constitutive_variables(strain_size);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters cl_values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options = cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_values.SetStrainVector(this_constitutive_variables.StrainVector);
    cl_values.SetStressVector(this_constitutive_variables.StressVector);
    cl_values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // Reading integration points
    const Properties& r_properties = this->GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(this->mThisIntegrationMethod);

    // Reading integration points
    const auto& integration_points = r_geometry.IntegrationPoints(this->mThisIntegrationMethod);

    for (IndexType point_number = 0; point_number < this->mConstitutiveLawVector.size(); ++point_number) {
        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->mThisIntegrationMethod);

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this->GetStressMeasure());
        this->SetValue(STRESS_VECTOR, this_constitutive_variables.StressVector);
        this->SetValue(STRAIN_VECTOR, this_constitutive_variables.StrainVector);
    }


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
    KRATOS_TRY;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];

    BaseSolidElement::KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    BaseSolidElement::ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag) { // Calculation of the matrix is required
        if (rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false);

        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    // Reading integration points and local gradients
    const auto& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());
    ConstitutiveLaw::Parameters cl_values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options=cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if (CalculateStiffnessMatrixFlag) {
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    cl_values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;
    const double characteristic_length = this->CalculateCharacteristicLength(this);

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this->GetStressMeasure());

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = this->GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if (dimension == 2 && this->GetProperties().Has(THICKNESS))
            int_to_reference_weight *= this->GetProperties()[THICKNESS];

        bool is_damaging = false;
        Vector damages_edges = ZeroVector(NumberOfEdges);
        if (yield_surface != "Elastic") {
            // Loop over edges of the element...
            Vector average_stress_edge(VoigtSize);
            Vector average_strain_edge(VoigtSize);

            for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
                noalias(average_stress_edge) = this_constitutive_variables.StressVector;
                noalias(average_strain_edge) = this_constitutive_variables.StrainVector;
                this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
                this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);
 
                damages_edges[edge] = this->mDamages[edge];
                double threshold    = this->mThresholds[edge];
                
                this->IntegrateStressDamageMechanics(threshold, damages_edges[edge], average_strain_edge, 
                                                     average_stress_edge, edge, characteristic_length, cl_values, 
                                                     is_damaging);
            } // Loop over edges
        }
        // Calculate the elemental Damage...
        const double damage_element = this->CalculateElementalDamage(damages_edges);

        const Vector& r_strain_vector = this_constitutive_variables.StrainVector;
        const Vector& r_stress_vector = this_constitutive_variables.StressVector;
        const Vector& r_integrated_stress_vector = (1.0 - damage_element)*r_stress_vector;

        if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            if (is_damaging == true && norm_2(r_strain_vector) > tolerance) {
                Matrix tangent_tensor;
                if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 0) {
                    tangent_tensor = this_constitutive_variables.D;
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 1) {
                    tangent_tensor = (1.0 - damage_element) * this_constitutive_variables.D;
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 2) {
                    this->CalculateTangentTensor(tangent_tensor, r_strain_vector, r_integrated_stress_vector, this_constitutive_variables.D, cl_values);
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 3) {
                    this->CalculateTangentTensorSecondOrder(tangent_tensor, r_strain_vector, r_integrated_stress_vector, this_constitutive_variables.D, cl_values);
                }
                this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, tangent_tensor, int_to_reference_weight);
            } else {
                this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, (1.0-damage_element)*this_constitutive_variables.D, int_to_reference_weight);
            }
        }
        if (CalculateResidualVectorFlag) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, r_integrated_stress_vector, int_to_reference_weight);
        }
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::FinalizeSolutionStep(
    ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const SizeType number_of_nodes   = this->GetGeometry().size();
    const SizeType dimension         = this->GetGeometry().WorkingSpaceDimension();
    const auto strain_size           = this->GetStrainSize();
    const std::string& yield_surface = this->GetProperties()[YIELD_SURFACE];

    BaseSolidElement::KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    BaseSolidElement::ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options = cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // If strain has to be computed inside of the constitutive law with PK2
    cl_values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;
    const double characteristic_length = this->CalculateCharacteristicLength(this);

    const GeometryType& r_geometry = this->GetGeometry();
    const Properties& r_properties = this->GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(this->mThisIntegrationMethod);

    // Computing in all integrations points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this->GetStressMeasure());

        // Call the constitutive law to update material variables
        this->mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(cl_values, this->GetStressMeasure());

        // TODO: Deprecated, remove this
        this->mConstitutiveLawVector[point_number]->FinalizeSolutionStep(r_properties, r_geometry, row( N_values, point_number ), rCurrentProcessInfo);

        bool is_damaging = false;
        if (yield_surface != "Elastic") {
            // Loop over edges of the element...
            Vector average_stress_edge(VoigtSize);
            Vector average_strain_edge(VoigtSize);

            for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
                noalias(average_stress_edge) = this_constitutive_variables.StressVector;
                noalias(average_strain_edge) = this_constitutive_variables.StrainVector;
                this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
                this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);
                
                this->IntegrateStressDamageMechanics(this->mThresholds[edge], this->mDamages[edge], 
                                                     average_strain_edge, average_stress_edge, edge, 
                                                     characteristic_length, cl_values, is_damaging);
            } // Loop over edges
        }
        // Calculate the elemental Damage...
        this->mDamage = this->CalculateElementalDamage(this->mDamages);

        if (this->mDamage >= 0.98) {
            this->Set(ACTIVE, false);
            this->mDamage = 0.98;
            // We set a "flag" to generate the DEM 
            rCurrentProcessInfo[GENERATE_DEM] = true;
        }
    }
    KRATOS_CATCH( "" )
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
void  GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateTangentTensorSecondOrder(
    Matrix& rTangentTensor,
    const Vector& rStrainVectorGP,
    const Vector& rStressVectorGP,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const double number_components = rStrainVectorGP.size();
    rTangentTensor.resize(number_components, number_components);
    Vector perturbed_stress_plus, perturbed_strain_plus;
    Vector perturbed_stress_minus, perturbed_strain_minus;
    perturbed_strain_minus.resize(number_components);
    perturbed_strain_plus.resize(number_components);
    perturbed_stress_minus.resize(number_components);
    perturbed_stress_plus.resize(number_components);

    for (unsigned int component = 0; component < number_components; component++) {
        double perturbation;
        this->CalculatePerturbation(rStrainVectorGP, perturbation, component);

        this->PerturbateStrainVector(perturbed_strain_plus, rStrainVectorGP, perturbation, component);
        this->IntegratePerturbedStrain(perturbed_stress_plus, perturbed_strain_plus, rElasticMatrix, rValues);

        this->PerturbateStrainVector(perturbed_strain_minus, rStrainVectorGP, -perturbation, component);
        this->IntegratePerturbedStrain(perturbed_stress_minus, perturbed_strain_minus, rElasticMatrix, rValues);

        const Vector& r_delta_stress = perturbed_stress_plus - perturbed_stress_minus;
        this->AssignComponentsToTangentTensor(rTangentTensor, r_delta_stress, 2.0*perturbation, component);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateKinematicVariables(
    BaseSolidElement::KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    const auto& r_geometry = this->GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());
    rThisKinematicVariables.detJ0 = this->CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    this->CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.DN_DX);

    // Compute equivalent F
    this->GetValuesVector(rThisKinematicVariables.Displacements);
    Vector strain_vector = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void GenericSmallStrainFemDemElement<2,0>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,0>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<2,1>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,1>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<2,2>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,2>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<2,3>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,3>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<2,4>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,4>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<2,5>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,5>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<2,6>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate2DB(rB, rDN_DX);
}
template<>
void GenericSmallStrainFemDemElement<3,6>::CalculateB(Matrix& rB, const Matrix& rDN_DX)
{
    this->Calculate3DB(rB, rDN_DX);
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::Calculate2DB(
    Matrix& rB, 
    const Matrix& rDN_DX
    )
{
    KRATOS_TRY

    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();

    for (SizeType i = 0; i < number_of_nodes; ++i) {
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

    for (SizeType i = 0; i < number_of_nodes; ++i) {
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

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    )
{
    const SizeType dim = this->GetGeometry().WorkingSpaceDimension();

    if (dim == 2) {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(2);
        rF(1,0) = 0.5*rStrainTensor(2);
        rF(1,1) = 1.0+rStrainTensor(1);
    } else {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(3);
        rF(0,2) = 0.5*rStrainTensor(5);
        rF(1,0) = 0.5*rStrainTensor(3);
        rF(1,1) = 1.0+rStrainTensor(1);
        rF(1,2) = 0.5*rStrainTensor(4);
        rF(2,0) = 0.5*rStrainTensor(5);
        rF(2,1) = 0.5*rStrainTensor(4);
        rF(2,2) = 1.0+rStrainTensor(2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateConstitutiveVariables(
    BaseSolidElement::KinematicVariables& rThisKinematicVariables,
    BaseSolidElement::ConstitutiveVariables& rThisConstitutiveVariables,
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
    BaseSolidElement::KinematicVariables& rThisKinematicVariables,
    BaseSolidElement::ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    const auto& r_geometry         = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension       = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size        = number_of_nodes * dimension;

    // Displacements vector
    Vector displacements(mat_size);
    this->GetValuesVector(displacements);

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

template<unsigned int TDim, unsigned int TyieldSurf>
bool GenericSmallStrainFemDemElement<TDim,TyieldSurf>::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericSmallStrainFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );


    if ( rVariable == INSITU_STRESS ) {
        const SizeType strain_size = this->mConstitutiveLawVector[0]->GetStrainSize();
        Vector strain_vector( strain_size );

        for ( IndexType point_number = 0; point_number < this->mConstitutiveLawVector.size(); ++point_number ) {
            if ( rOutput[point_number].size() != strain_vector.size() )
                rOutput[point_number].resize( strain_vector.size(), false );

            rOutput[point_number] = this->mConstitutiveLawVector[point_number]->GetValue( INSITU_STRESS, rOutput[point_number] );
        }
    } else if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
        // Create and initialize element variables:
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = this->mConstitutiveLawVector[0]->GetStrainSize();

        BaseSolidElement::KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        BaseSolidElement::ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& cl_options=cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            //call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR) {
                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            } else {
                // Compute material reponse
                this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points,ConstitutiveLaw::StressMeasure_PK2);
            }

            if (rOutput[point_number].size() != strain_size )
                rOutput[point_number].resize( strain_size, false);

            rOutput[point_number] = this_constitutive_variables.StressVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR) {
        // Create and initialize element variables:
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = this->mConstitutiveLawVector[0]->GetStrainSize();

        BaseSolidElement::KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        BaseSolidElement::ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &cl_options = cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        const ConstitutiveLaw::StressMeasure this_stress_measure = rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ? ConstitutiveLaw::StressMeasure_PK2 : ConstitutiveLaw::StressMeasure_Kirchhoff;

        //reading integration points
        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            // Compute material reponse
            this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, this_stress_measure);

            if ( rOutput[point_number].size() != strain_size)
                rOutput[point_number].resize( strain_size, false );

            rOutput[point_number] = this_constitutive_variables.StrainVector;
        }
    } else if (rVariable == STRESS_VECTOR_INTEGRATED) {
        const GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod() );
        const SizeType number_of_integration_points = integration_points.size();
        if ( rOutput.size() != number_of_integration_points )
            rOutput.resize( number_of_integration_points );
        // Create and initialize element variables:
        const SizeType number_of_nodes = this->GetGeometry().size();
        const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
        const SizeType strain_size = this->mConstitutiveLawVector[0]->GetStrainSize();

        BaseSolidElement::KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
        BaseSolidElement::ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags& cl_options = cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            //call the constitutive law to update material variables
            // Compute material reponse
            CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);

            if ( rOutput[point_number].size() != strain_size )
                rOutput[point_number].resize(strain_size, false);
            rOutput[point_number] = (1.0 - this->mDamage)*this_constitutive_variables.StressVector;
        }
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
template class GenericSmallStrainFemDemElement<2,6>;
template class GenericSmallStrainFemDemElement<3,0>;
template class GenericSmallStrainFemDemElement<3,1>;
template class GenericSmallStrainFemDemElement<3,2>;
template class GenericSmallStrainFemDemElement<3,3>;
template class GenericSmallStrainFemDemElement<3,4>;
template class GenericSmallStrainFemDemElement<3,5>;
template class GenericSmallStrainFemDemElement<3,6>;
} // namespace Kratos