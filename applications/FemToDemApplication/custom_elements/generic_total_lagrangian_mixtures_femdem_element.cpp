//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//
// System includes

// External includes

// Project includes
#include "custom_elements/generic_total_lagrangian_mixtures_femdem_element.hpp"
#include "fem_to_dem_application_variables.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry)
{
    GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::BaseType(NewId, pGeometry);

    if (mAcumulatedPlasticStrains.size() != NumberOfEdges)
        mAcumulatedPlasticStrains.resize(NumberOfEdges);
    noalias(mAcumulatedPlasticStrains) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (mPlasticityThresholds.size() != NumberOfEdges)
        mPlasticityThresholds.resize(NumberOfEdges);
    noalias(mPlasticityThresholds) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge

    if (mPlasticStrains.size() != NumberOfEdges)
        mPlasticStrains.resize(NumberOfEdges);
    for (IndexType i = 0; i < NumberOfEdges; ++i) {
        mPlasticStrains[i].resize(VoigtSize);
        noalias(mPlasticStrains[i]) = ZeroVector(VoigtSize);
    }

}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry, pProperties)
{
    GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::BaseType(NewId, pGeometry, pProperties);

    if (mAcumulatedPlasticStrains.size() != NumberOfEdges)
        mAcumulatedPlasticStrains.resize(NumberOfEdges);
    noalias(mAcumulatedPlasticStrains) = ZeroVector(NumberOfEdges); // Stress mThreshold on edge

    if (this->mPlasticityThresholds.size() != NumberOfEdges)
        this->mPlasticityThresholds.resize(NumberOfEdges);
    noalias(this->mPlasticityThresholds) = ZeroVector(NumberOfEdges); // Converged mDamage on each edge

    if (this->mPlasticStrains.size() != NumberOfEdges)
        this->mPlasticStrains.resize(NumberOfEdges);
    for (IndexType i = 0; i < NumberOfEdges; ++i) {
        this->mPlasticStrains[i].resize(VoigtSize);
        noalias(mPlasticStrains[i]) = ZeroVector(VoigtSize);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>(NewId, this->GetGeometry().Create( ThisNodes ), pProperties);
}

//************************************************************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::~GenericTotalLagrangianMixturesFemDemElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Pointer p_new_elem = Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
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
Vector GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::IntegrateSmoothedConstitutiveLaw(
    const std::string& rYieldSurface,
    ConstitutiveLaw::Parameters& rValues,
    const BaseSolidElement::ConstitutiveVariables& rThisConstVars,
    const BaseSolidElement::KinematicVariables& rKinVariables,
    Vector& rStrainVector,
    double& rDamageElement,
    bool& rIsDamaging,
    const double CharacteristicLength,
    const bool SaveIntVars
    )
{
    Vector damages_edges = ZeroVector(NumberOfEdges);
    Vector plastic_strain_edges_sum(VoigtSize);
    noalias(plastic_strain_edges_sum) = ZeroVector(VoigtSize);
    double plastic_uniaxial_stress = 0.0;

    if (rYieldSurface != "Elastic") {
        // Loop over edges of the element...
        Vector average_stress_edge(VoigtSize);
        Vector average_strain_edge(VoigtSize);

        for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
            noalias(average_stress_edge) = rThisConstVars.StressVector;
            noalias(average_strain_edge) = rThisConstVars.StrainVector;
            this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
            this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

            if (!SaveIntVars) {
                damages_edges[edge] = this->mDamages[edge];
                double threshold    = this->mThresholds[edge];

                // Matrix constitutive damage model
                this->IntegrateStressDamageMechanics(threshold, damages_edges[edge], average_strain_edge,
                                                     average_stress_edge, edge, CharacteristicLength, rValues,
                                                     rIsDamaging);
                rDamageElement = this->CalculateElementalDamage(damages_edges);

                // Fibre plasticity constitutive model
                Vector plastic_strain            = this->mPlasticStrains[edge];
                double acumulated_plastic_strain = this->mAcumulatedPlasticStrains[edge];
                double plasticity_threshold      = this->mPlasticityThresholds[edge];


                Vector r_fiber_stress = this->IntegrateStressPlasticity(rValues, average_strain_edge,
                                                                        plastic_strain, acumulated_plastic_strain,
                                                                        plasticity_threshold, plastic_uniaxial_stress,
                                                                        rIsDamaging);
                plastic_strain_edges_sum += plastic_strain; // we smooth the plastic strain
            } else {
                this->IntegrateStressDamageMechanics(this->mThresholds[edge], this->mDamages[edge], average_strain_edge,
                                                     average_stress_edge, edge, CharacteristicLength, rValues,
                                                     rIsDamaging);
                this->mDamage = this->CalculateElementalDamage(this->mDamages);
                rDamageElement = this->mDamage;

                Vector r_fiber_stress = this->IntegrateStressPlasticity(rValues, average_strain_edge,
                                                                        this->mPlasticStrains[edge], this->mAcumulatedPlasticStrains[edge],
                                                                        this->mPlasticityThresholds[edge], this->mUniaxialStress,
                                                                        rIsDamaging);
                plastic_strain_edges_sum += this->mPlasticStrains[edge]; // we smooth the plastic strain
            }
        } // Loop over edges
    }

    plastic_strain_edges_sum /= NumberOfEdges;

    this->CalculateGreenLagrangeStrainVector(rStrainVector, rKinVariables.F);
    const Vector& r_stress_vector = rThisConstVars.StressVector;

    Matrix fiber_constitutive_matrix(VoigtSize, VoigtSize);
    noalias(fiber_constitutive_matrix) = ZeroMatrix(VoigtSize, VoigtSize);
    auto &r_mat_props    = rValues.GetMaterialProperties();
    const double young   = r_mat_props[YOUNG_MODULUS_FIBER];
    const double poisson = r_mat_props[POISSON_RATIO_FIBER];
    ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(fiber_constitutive_matrix, young, poisson);

    const Vector &r_stress_vector_fiber = prod(fiber_constitutive_matrix, rStrainVector - plastic_strain_edges_sum);

    const double fiber_vol_part = r_mat_props[FIBER_VOLUMETRIC_PART];

    return (1.0 - fiber_vol_part)*(1.0 - rDamageElement) * r_stress_vector + fiber_vol_part * r_stress_vector_fiber;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Vector GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::IntegrateStressPlasticity(
    ConstitutiveLaw::Parameters& rValues,
    const Vector& rStrainVector,
    Vector& rPlasticStrainVector,
    double& rAcumulatedPlasticStrain,
    double& rThreshold,
    double& rUniaxialStress,
    bool& rIsPlastifying
    )
{
    auto &r_mat_props = rValues.GetMaterialProperties();

    // E^e = E - E^p
    const Vector &r_elastic_strain = rStrainVector - rPlasticStrainVector;
    Vector volumetric_elastic_strain(VoigtSize), deviatoric_elastic_strain(VoigtSize);

    // Now we perform the splitting E = Ed + Evol
    ConstitutiveLawUtilities<VoigtSize>::CalculateVolumetricStrainVector(
        r_elastic_strain, volumetric_elastic_strain);
    ConstitutiveLawUtilities<VoigtSize>::CalculateDeviatoricStrainVector(
        r_elastic_strain, volumetric_elastic_strain, deviatoric_elastic_strain);

    const double young         = r_mat_props[YOUNG_MODULUS_FIBER];
    const double poisson       = r_mat_props[POISSON_RATIO_FIBER];
    const double bulk_modulus  = ConstitutiveLawUtilities<VoigtSize>::CalculateBulkModulus(young, poisson);
    const double shear_modulus = ConstitutiveLawUtilities<VoigtSize>::CalculateShearModulus(young, poisson);

    const Vector &pressure_stress = bulk_modulus * volumetric_elastic_strain;
    Vector r_deviator_stress      = 2.0 * shear_modulus * deviatoric_elastic_strain;

    rUniaxialStress = std::sqrt(1.5 * MathUtils<double>::Dot(r_deviator_stress, r_deviator_stress));
    this->ComputePlasticThreshold(rAcumulatedPlasticStrain, rThreshold, rValues);
    const double stress_excess = rUniaxialStress - rThreshold;

    // Check if plasticity occurs
    if (stress_excess > tolerance) { // We must compute plasticity
        rIsPlastifying = true;
        const double hardening_modulus  = r_mat_props.Has(HARDENING_MODULUS) ? r_mat_props[HARDENING_MODULUS] : 0.0;
        const double plastic_multiplier = stress_excess / (3.0 * shear_modulus + hardening_modulus);
        rAcumulatedPlasticStrain += plastic_multiplier;
        this->ComputePlasticThreshold(rAcumulatedPlasticStrain, rThreshold, rValues);
        rPlasticStrainVector += plastic_multiplier * std::sqrt(1.5) * r_deviator_stress / MathUtils<double>::Norm(r_deviator_stress);
        r_deviator_stress = (1.0 - 3.0 * plastic_multiplier * shear_modulus / rUniaxialStress) * r_deviator_stress;
    }
    return pressure_stress + r_deviator_stress; // stress integrated
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::ComputePlasticMultiplier(
    const double UniaxialStress,
    const double Threshold,
    double& rPlasticMultiplier,
    ConstitutiveLaw::Parameters& rValues
    )
{
    auto &r_mat_props              = rValues.GetMaterialProperties();
    const double young_fiber       = r_mat_props[YOUNG_MODULUS_FIBER];
    const double poisson_fiber     = r_mat_props[POISSON_RATIO_FIBER];
    const double shear_modulus     = ConstitutiveLawUtilities<VoigtSize>::CalculateShearModulus(young_fiber, poisson_fiber);
    const double hardening_modulus = r_mat_props.Has(HARDENING_MODULUS) ? r_mat_props[HARDENING_MODULUS] : 0.0;

    // It assumes perfect plasticity of linear hardening, hence no iteration is required
    rPlasticMultiplier = (UniaxialStress - Threshold) / (3.0 * shear_modulus + hardening_modulus);
}

/***********************************************************************************/
/***********************************************************************************/
template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::ComputePlasticThreshold(
    const double AcumulatedPlasticStrain,
    double& rThreshold,
    ConstitutiveLaw::Parameters& rValues
    )
{
    auto &r_mat_props = rValues.GetMaterialProperties();
    const double yield_stress      = r_mat_props[YIELD_STRESS_T_FIBER];
    const double hardening_modulus = r_mat_props.Has(HARDENING_MODULUS) ? r_mat_props[HARDENING_MODULUS] : 0.0;
    rThreshold = yield_stress + hardening_modulus * AcumulatedPlasticStrain;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::IntegratePerturbedStrain(
    Vector& rPerturbedStressVector,
    const Vector& rPerturbedStrainVector,
    const Matrix& rElasticMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    Vector damages_edges = ZeroVector(NumberOfEdges);
    Vector plastic_strain_edges_sum(VoigtSize);
    noalias(plastic_strain_edges_sum) = ZeroVector(VoigtSize);

    Vector average_stress_edge(VoigtSize);
    Vector average_strain_edge(VoigtSize);

    const Vector& r_perturbed_predictive_stress = prod(rElasticMatrix, rPerturbedStrainVector);
    const double characteristic_length = this->CalculateCharacteristicLength(this);

    for (unsigned int edge = 0; edge < NumberOfEdges; edge++) {
        noalias(average_stress_edge) = r_perturbed_predictive_stress;
        noalias(average_strain_edge) = rPerturbedStrainVector;
        this->CalculateAverageVariableOnEdge(this, STRESS_VECTOR, average_stress_edge, edge);
        this->CalculateAverageVariableOnEdge(this, STRAIN_VECTOR, average_strain_edge, edge);

        damages_edges[edge] = this->mDamages[edge];
        double threshold    = this->mThresholds[edge];
        bool dummy = false;

        // Matrix constitutive damage model
        this->IntegrateStressDamageMechanics(threshold, damages_edges[edge], average_strain_edge,
                                                average_stress_edge, edge, characteristic_length, rValues,
                                                dummy);

        // Fibre plasticity constitutive model
        Vector plastic_strain            = this->mPlasticStrains[edge];
        double acumulated_plastic_strain = this->mAcumulatedPlasticStrains[edge];
        double plasticity_threshold      = this->mPlasticityThresholds[edge];
        double plastic_uniaxial_stress = 0.0;

        Vector r_fiber_stress = this->IntegrateStressPlasticity(rValues, average_strain_edge,
                                                                plastic_strain, acumulated_plastic_strain,
                                                                plasticity_threshold, plastic_uniaxial_stress,
                                                                dummy);
        plastic_strain_edges_sum += plastic_strain; // we smooth the plastic strain

    } // Loop over edges

    plastic_strain_edges_sum /= NumberOfEdges;
    const double damage_element = this->CalculateElementalDamage(damages_edges);

    Matrix fiber_constitutive_matrix(VoigtSize, VoigtSize);
    noalias(fiber_constitutive_matrix) = ZeroMatrix(VoigtSize, VoigtSize);
    auto &r_mat_props    = rValues.GetMaterialProperties();
    const double young   = r_mat_props[YOUNG_MODULUS_FIBER];
    const double poisson = r_mat_props[POISSON_RATIO_FIBER];
    ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(fiber_constitutive_matrix, young, poisson);

    const Vector &r_stress_vector_fiber = prod(fiber_constitutive_matrix, rPerturbedStrainVector - plastic_strain_edges_sum);

    const double fiber_vol_part = r_mat_props[FIBER_VOLUMETRIC_PART];

    rPerturbedStressVector = (1.0 - fiber_vol_part)*(1.0 - damage_element) * r_perturbed_predictive_stress +
        fiber_vol_part * r_stress_vector_fiber;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<double> &rVariable,
    std::vector<double> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    if (rVariable == ACUMULATED_PLASTIC_STRAIN) {
        rOutput.resize(1);
        for (unsigned int point_number = 0; point_number < 1; point_number++) {
            rOutput[point_number] = this->CalculateElementalDamage(this->mAcumulatedPlasticStrains);
        }
    } else if (rVariable == PLASTIC_UNIAXIAL_STRESS) {
        rOutput.resize(1);
        for (unsigned int point_number = 0; point_number < 1; point_number++) {
            rOutput[point_number] = this->mUniaxialStress;
        }
    }
}

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

    if (rVariable == FIBER_STRESS_VECTOR) {
        const GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
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
        Flags& cl_options=cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, this->UseElementProvidedStrain());
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);

        // Reading integration points
        for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
            // Compute element kinematics B, F, DN_DX ...
            this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

            //call the constitutive law to update material variables
            // Compute material reponse
            this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);

            if ( rOutput[point_number].size() != strain_size ) {}
                rOutput[point_number].resize(strain_size, false);

            Matrix fiber_constitutive_matrix(VoigtSize, VoigtSize);
            noalias(fiber_constitutive_matrix) = ZeroMatrix(VoigtSize, VoigtSize);
            auto &r_mat_props    = cl_values.GetMaterialProperties();
            const double young   = r_mat_props[YOUNG_MODULUS_FIBER];
            const double poisson = r_mat_props[POISSON_RATIO_FIBER];
            ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(fiber_constitutive_matrix, young, poisson);

            rOutput[point_number] = prod(fiber_constitutive_matrix,
                                         this_constitutive_variables.StrainVector - this->CalculateAveragePlasticStrain());
        }
    } else if (rVariable == PLASTIC_STRAIN_VECTOR) {
        rOutput[0] = this->CalculateAveragePlasticStrain();
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints( this->GetIntegrationMethod() );
    const SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
    BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    std::vector<Vector> stress_vector;

    if (rVariable == FIBER_STRESS_TENSOR) {
        this->CalculateOnIntegrationPoints(FIBER_STRESS_VECTOR, stress_vector, rCurrentProcessInfo);

        // Loop integration points
        for (IndexType point_number = 0; point_number < this->mConstitutiveLawVector.size(); ++point_number) {
            if (rOutput[point_number].size2() != dimension)
                rOutput[point_number].resize(dimension, dimension, false);

            rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
        }
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
void GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
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

    if (CalculateStiffnessMatrixFlag == true) { // Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) { // Calculation of the matrix is required
        if (rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize(mat_size, false);
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    // Reading integration points
    const auto& integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& cl_options = cl_values.GetOptions();
    cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, this->UseElementProvidedStrain());
    cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
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
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
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

        Vector r_strain_vector;
        double damage_element = 0.0;
        bool is_damaging = false;
        const Vector &r_integrated_stress_vector = this->IntegrateSmoothedConstitutiveLaw(
                                                        yield_surface, cl_values, this_constitutive_variables,
                                                        this_kinematic_variables, r_strain_vector, damage_element,
                                                        is_damaging, characteristic_length, false);

        if (CalculateStiffnessMatrixFlag == true) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            /* Material stiffness matrix */
            if (is_damaging == true && norm_2(r_strain_vector) > tolerance) {
                 Matrix tangent_tensor(VoigtSize, VoigtSize);
                if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 0) {
                    noalias(tangent_tensor) = this_constitutive_variables.D;
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 1) {
                    noalias(tangent_tensor) = (1.0 - damage_element) * this_constitutive_variables.D;
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 2) {
                    this->CalculateTangentTensor(tangent_tensor, r_strain_vector, r_integrated_stress_vector, this_kinematic_variables.F, this_constitutive_variables.D, cl_values);
                } else if (rCurrentProcessInfo[TANGENT_CONSTITUTIVE_TENSOR] == 3) {
                    this->CalculateTangentTensorSecondOrder(tangent_tensor, r_strain_vector, r_integrated_stress_vector, this_kinematic_variables.F, this_constitutive_variables.D, cl_values);
                }
                this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, tangent_tensor, int_to_reference_weight);
            } else {
                Matrix tangent_tensor(VoigtSize, VoigtSize);
                const double kf = this->GetProperties()[FIBER_VOLUMETRIC_PART];
                const double km = 1.0 - kf;

                Matrix fiber_constitutive_matrix(VoigtSize, VoigtSize);
                noalias(fiber_constitutive_matrix) = ZeroMatrix(VoigtSize, VoigtSize);
                auto &r_mat_props    = this->GetProperties();
                const double young   = r_mat_props[YOUNG_MODULUS_FIBER];
                const double poisson = r_mat_props[POISSON_RATIO_FIBER];
                ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(fiber_constitutive_matrix, young, poisson);

                noalias(tangent_tensor) = km * (1.0 - damage_element) * this_constitutive_variables.D + kf * fiber_constitutive_matrix;
                this->CalculateAndAddKm(rLeftHandSideMatrix, this_kinematic_variables.B, tangent_tensor, int_to_reference_weight);
            }
            /* Geometric stiffness matrix */
            this->CalculateAndAddKg(rLeftHandSideMatrix, this_kinematic_variables.DN_DX, this_constitutive_variables.StressVector, int_to_reference_weight);
        }

        if (CalculateResidualVectorFlag == true) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, r_integrated_stress_vector, int_to_reference_weight);
        }
    }
    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericTotalLagrangianMixturesFemDemElement<2,0>;
template class GenericTotalLagrangianMixturesFemDemElement<2,1>;
template class GenericTotalLagrangianMixturesFemDemElement<2,2>;
template class GenericTotalLagrangianMixturesFemDemElement<2,3>;
template class GenericTotalLagrangianMixturesFemDemElement<2,4>;
template class GenericTotalLagrangianMixturesFemDemElement<2,5>;
template class GenericTotalLagrangianMixturesFemDemElement<2,6>;
template class GenericTotalLagrangianMixturesFemDemElement<3,0>;
template class GenericTotalLagrangianMixturesFemDemElement<3,1>;
template class GenericTotalLagrangianMixturesFemDemElement<3,2>;
template class GenericTotalLagrangianMixturesFemDemElement<3,3>;
template class GenericTotalLagrangianMixturesFemDemElement<3,4>;
template class GenericTotalLagrangianMixturesFemDemElement<3,5>;
template class GenericTotalLagrangianMixturesFemDemElement<3,6>;
} // namespace Kratos