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
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************DEFAULT CONSTRUCTOR******************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
GenericTotalLagrangianMixturesFemDemElement<TDim, TyieldSurf>::GenericTotalLagrangianMixturesFemDemElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GenericTotalLagrangianFemDemElement<TDim,TyieldSurf>(NewId, pGeometry)
{
    BaseType::BaseType(NewId, pGeometry);

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
    BaseType::BaseType(NewId, pGeometry, pProperties);

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

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>( NewId, this->GetGeometry().Create( ThisNodes ), pProperties );
}

//************************************************************************************
//************************************************************************************

template<unsigned int TDim, unsigned int TyieldSurf>
Element::Pointer GenericTotalLagrangianMixturesFemDemElement<TDim,TyieldSurf>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<GenericTotalLagrangianMixturesFemDemElement>( NewId, pGeom, pProperties );
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
    const ConstitutiveVariables& rThisConstVars,
    const KinematicVariables& rKinVariables,
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
                damages_edges[edge] = mDamages[edge];
                double threshold = mThresholds[edge];

                // Matrix constitutive damage model
                this->IntegrateStressDamageMechanics(threshold, damages_edges[edge], average_strain_edge,
                                                     average_stress_edge, edge, CharacteristicLength, rValues,
                                                     rIsDamaging);
                rDamageElement = this->CalculateElementalDamage(damages_edges); 

                // Fibre plasticity constitutive model
                Vector plastic_strain            = mPlasticStrains[edge];
                double acumulated_plastic_strain = mAcumulatedPlasticStrains[edge];
                double plasticity_threshold      = mPlasticityThresholds[edge];

                Vector r_fiber_stress = this->IntegrateStressPlasticity(rValues, average_strain_edge,
                                                                        plastic_strain, acumulated_plastic_strain,
                                                                        plasticity_threshold, rIsDamaging);
                plastic_strain_edges_sum += plastic_strain; // we smooth the plastic strain
            } else {
                this->IntegrateStressDamageMechanics(mThresholds[edge], mDamages[edge], average_strain_edge,
                                                     average_stress_edge, edge, CharacteristicLength, rValues,
                                                     rIsDamaging);
                mDamage = this->CalculateElementalDamage(mDamages);
                rDamageElement = mDamage;

                Vector r_fiber_stress = this->IntegrateStressPlasticity(rValues, average_strain_edge,
                                                                        mPlasticStrains[edge], mAcumulatedPlasticStrains[edge],
                                                                        mPlasticityThresholds[edge], rIsDamaging);
                plastic_strain_edges_sum += mPlasticStrains[edge]; // we smooth the plastic strain
            }
        } // Loop over edges
    }

    plastic_strain_edges_sum /= NumberOfEdges;

    this->CalculateGreenLagrangeStrainVector(rStrainVector, rKinVariables.F);
    const Vector& r_stress_vector = rThisConstVars.StressVector;

    Matrix fiber_constitutive_matrix(VoigtSize, VoigtSize);
    auto &r_mat_props    = rValues.GetMaterialProperties();
    const double young   = r_mat_props[YOUNG_MODULUS_FIBER];
    const double poisson = r_mat_props[POISSON_RATIO_FIBER];
    ConstitutiveLawUtilities<VoigtSize>::CalculateElasticMatrix(fiber_constitutive_matrix, young, poisson);

    const Vector &r_stress_vector_fiber = prod(fiber_constitutive_matrix, rStrainVector - plastic_strain_edges_sum);

    const double fiber_vol_part = r_mat_props[FIBER_VOLUMETRIC_PART];
    
    return (1.0 - fiber_vol_part)*(1.0 - rDamageElement) * r_stress_vector + fiber_vol_part * r_stress_vector_fiber;;
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

    const double young             = r_mat_props[YOUNG_MODULUS_FIBER];
    const double poisson           = r_mat_props[POISSON_RATIO_FIBER];
    const double bulk_modulus      = ConstitutiveLawUtilities<VoigtSize>::CalculateBulkModulus(young, poisson);
    const double shear_modulus     = ConstitutiveLawUtilities<VoigtSize>::CalculateShearModulus(young, poisson);

    const Vector &pressure_stress = bulk_modulus * volumetric_elastic_strain;
    Vector r_deviator_stress      = 2.0 * shear_modulus * deviatoric_elastic_strain;

    const double uniaxial_stress = std::sqrt(1.5 * MathUtils<double>::Dot(r_deviator_stress, r_deviator_stress));
    this->ComputePlasticThreshold(rAcumulatedPlasticStrain, rThreshold, rValues);
    const double stress_excess = uniaxial_stress - rThreshold;

    // Check if plasticity occurs
    if (stress_excess > tolerance) { // We must compute plasticity
        rIsPlastifying = true;
        const double hardening_modulus  = r_mat_props.Has(HARDENING_MODULUS) ? r_mat_props[HARDENING_MODULUS] : 0.0;
        const double plastic_multiplier = stress_excess / (3.0 * shear_modulus + hardening_modulus);
        rAcumulatedPlasticStrain       += plastic_multiplier;
        this->ComputePlasticThreshold(rAcumulatedPlasticStrain, rThreshold, rValues);
        rPlasticStrainVector           += plastic_multiplier * std::sqrt(1.5) * r_deviator_stress / MathUtils<double>::Norm(r_deviator_stress);
        r_deviator_stress = (1.0 - 3.0 * plastic_multiplier * shear_modulus / uniaxial_stress) * r_deviator_stress;
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
    const double yield_stress      = r_mat_props[YIELD_STRESS];
    const double hardening_modulus = r_mat_props.Has(HARDENING_MODULUS) ? r_mat_props[HARDENING_MODULUS] : 0.0;
    rThreshold = yield_stress + hardening_modulus * AcumulatedPlasticStrain;
}



/***********************************************************************************/
/***********************************************************************************/




/***********************************************************************************/
/***********************************************************************************/




/***********************************************************************************/
/***********************************************************************************/









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