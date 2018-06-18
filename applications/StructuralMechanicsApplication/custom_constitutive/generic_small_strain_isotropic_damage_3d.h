// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//

#if !defined (KRATOS_GENERIC_SMALL_STRAIN_ISOTROPIC_DAMAGE_3D_H_INCLUDED)
#define  KRATOS_GENERIC_SMALL_STRAIN_ISOTROPIC_DAMAGE_3D_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{
/**
 * @class GenericConstitutiveLawIntegrator
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */
template <class ConstLawIntegratorType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainIsotropicDamage3D
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainIsotropicDamage3D);

    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainIsotropicDamage3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>::Pointer p_clone
            (new GenericSmallStrainIsotropicDamage3D<ConstLawIntegratorType>(*this));
        return p_clone;
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainIsotropicDamage3D (const GenericSmallStrainIsotropicDamage3D& rOther)
    : ConstitutiveLaw(rOther)
    {
    }
    /**
    * Destructor.
    */
    ~GenericSmallStrainIsotropicDamage3D() override
    {
    }

    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override
    {
        const std::string& yield = NewParameters["yield_surface"].GetString();
        const std::string& potential = NewParameters["plastic_potential"].GetString();
        ConstitutiveLaw::Pointer constitutive_law;

        if (yield == "VonMises") {
            if (potential == "VonMises") {
                // constitutive_law = Kratos::make_shared
                //     <GenericConstitutiveLawIntegratorDamage
                //         <VonMisesYieldSurface
                //             <VonMisesPlasticPotential>>>()
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <VonMisesYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <VonMisesYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <VonMisesYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <VonMisesYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong" << std::endl;
            }
        } else if (yield == "ModifiedMohrCoulomb") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <ModifiedMohrCoulombYieldSurface
                            <VonMisesPlasticPotential>>>
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <ModifiedMohrCoulombYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <ModifiedMohrCoulombYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <ModifiedMohrCoulombYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong" << std::endl;
            }
        } else if (yield == "Tresca") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <TrescaYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <TrescaYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <TrescaYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone() 
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <TrescaYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong" << std::endl;
            }
        } else if (yield == "DruckerPrager") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <DruckerPragerYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <DruckerPragerYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <DruckerPragerYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <DruckerPragerYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong" << std::endl;
            }
        } else if (yield == "Rankine") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <RankineYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <RankineYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <RankineYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <RankineYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong" << std::endl;
            }
        } else if (yield == "SimoJu") {
                        if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <SimoJuYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <SimoJuYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <SimoJuYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicDamage3D 
                    <GenericConstitutiveLawIntegratorDamage
                        <SimoJuYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong" << std::endl;
            }
        } else {
            KRATOS_ERROR << "Yield Surface not defined or wrong" << std::endl;
        }
        return constitutive_law;
    }
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    int GetVoigtSize() {return 6;}
    int GetWorkingSpaceDimension() {return 3;}

    double GetThreshold() {return mThreshold;}
    double GetDamage() {return mDamage;}
    double GetNonConvThreshold() {return mNonConvThreshold;}
    double GetNonConvDamage() {return mNonConvDamage;}

    void SetThreshold(const double& toThreshold) {mThreshold = toThreshold;}
    void SetDamage(const double& toDamage) {mDamage = toDamage;}
    void SetNonConvThreshold(const double& toThreshold) {mNonConvThreshold = toThreshold;}
    void SetNonConvDamage(const double& toDamage) {mNonConvDamage = toDamage;}


    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }

    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override
    {
        // Integrate Stress Damage
        const Properties& rMaterialProperties = rValues.GetMaterialProperties();
        const int VoigtSize = this->GetVoigtSize();
        Vector& IntegratedStressVector = rValues.GetStressVector();
        Matrix& TangentTensor = rValues.GetConstitutiveMatrix(); // todo modify after integration

        // Elastic Matrix
        Matrix C;
        this->CalculateElasticMatrix(C, rMaterialProperties);

        double Threshold, Damage;
        // In the 1st step must be set
        if (std::abs(this->GetThreshold()) < tolerance)
        {
            ConstLawIntegratorType::YieldSurfaceType::GetInitialUniaxialThreshold(rMaterialProperties, Threshold);
            this->SetThreshold(Threshold);
        }

        // Converged values
        Threshold = this->GetThreshold();
        Damage    = this->GetDamage();

        // S0 = C:(E-Ep)
        Vector PredictiveStressVector = prod(C, rValues.GetStrainVector());

        // Initialize Plastic Parameters
        double UniaxialStress;
        ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(PredictiveStressVector, 
            rValues.GetStrainVector(), UniaxialStress, rMaterialProperties);

        const double F = UniaxialStress - Threshold; 

        if (F <= 0.0) {   // Elastic case
            noalias(IntegratedStressVector) = PredictiveStressVector;
            this->SetNonConvDamage(Damage);
            this->SetNonConvThreshold(Threshold);
            
            noalias(TangentTensor) = (1.0 - Damage)*C;
        } else { // Damage case
            const double CharacteristicLength = rValues.GetElementGeometry().Length();

            // This routine updates the PredictiveStress to verify the yield surf
            ConstLawIntegratorType::IntegrateStressVector(PredictiveStressVector, UniaxialStress,
                Damage, Threshold, rMaterialProperties, CharacteristicLength);

            // Updated Values
            noalias(IntegratedStressVector) = PredictiveStressVector; 
            this->SetNonConvDamage(Damage);
            this->SetNonConvThreshold(Threshold);

            //noalias(TangentTensor) = (1 - Damage)*C; // Secant Tensor
            this->CalculateTangentTensor(rValues); 
            TangentTensor = rValues.GetConstitutiveMatrix();
        }
        
    } // End CalculateMaterialResponseCauchy

    void CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues) 
    {
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this);
    }

    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        this->SetDamage(this->GetNonConvDamage());
        this->SetThreshold(this->GetNonConvThreshold());
    }

    void CalculateElasticMatrix(Matrix &rElasticityTensor,
        const Properties &rMaterialProperties
    )
    {
        const double E = rMaterialProperties[YOUNG_MODULUS];
        const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
        const double lambda =
            E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
        const double mu = E / (2.0 + 2.0 * poisson_ratio);

        if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
            rElasticityTensor.resize(6, 6, false);
        rElasticityTensor.clear();

        rElasticityTensor(0, 0) = lambda + 2.0 * mu;
        rElasticityTensor(0, 1) = lambda;
        rElasticityTensor(0, 2) = lambda;
        rElasticityTensor(1, 0) = lambda;
        rElasticityTensor(1, 1) = lambda + 2.0 * mu;
        rElasticityTensor(1, 2) = lambda;
        rElasticityTensor(2, 0) = lambda;
        rElasticityTensor(2, 1) = lambda;
        rElasticityTensor(2, 2) = lambda + 2.0 * mu;
        rElasticityTensor(3, 3) = mu;
        rElasticityTensor(4, 4) = mu;
        rElasticityTensor(5, 5) = mu;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    // Converged values
    double mDamage = 0.0;
    double mThreshold = 0.0;

    // Non Converged values
    double mNonConvDamage = 0.0;
    double mNonConvThreshold = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
        rSerializer.save("Damage", mDamage);
        rSerializer.save("Threshold", mThreshold);
        rSerializer.save("NonConvDamage", mNonConvDamage);
        rSerializer.save("mNonConvThreshold", mNonConvThreshold);

    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("Damage", mDamage);
        rSerializer.load("Threshold", mThreshold);
        rSerializer.load("NonConvDamage", mDamage);
        rSerializer.load("mNonConvThreshold", mNonConvThreshold);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace kratos
#endif
