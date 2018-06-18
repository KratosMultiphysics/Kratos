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

#if !defined (KRATOS_SMALL_STRAIN_ISOTROPIC_PLASTICITY_FACTORY_3D_H_INCLUDED)
#define  KRATOS_SMALL_STRAIN_ISOTROPIC_PLASTICITY_FACTORY_3D_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"

#include "includes/constitutive_law.h"
#include "structural_mechanics_application_variables.h"

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"

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
 * @class SmallStrainIsotropicPlasticityFactory3D
 * @ingroup StructuralMechanicsApplication
 * @brief: dummy class to register, only implements create()
 * @details
 * @author Alejandro Cornejo & Lucia Barbu
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallStrainIsotropicPlasticityFactory3D
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(SmallStrainIsotropicPlasticityFactory3D);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    SmallStrainIsotropicPlasticityFactory3D()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        SmallStrainIsotropicPlasticityFactory3D::Pointer p_clone
            (new SmallStrainIsotropicPlasticityFactory3D(*this));
        return p_clone;
    }

    /**
    * Copy constructor.
    */
    SmallStrainIsotropicPlasticityFactory3D (const SmallStrainIsotropicPlasticityFactory3D& rOther)
    : ConstitutiveLaw(rOther)
    {
    }
    /**
    * Destructor.
    */
    ~SmallStrainIsotropicPlasticityFactory3D() override
    {
    }

    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override
    {
        const std::string& yield = NewParameters["yield_surface"].GetString();
        const std::string& potential = NewParameters["plastic_potential"].GetString();
        ConstitutiveLaw::Pointer constitutive_law;

        KRATOS_ERROR_IF(yield == "SimoJu")  << "SimoJu yield surface not available in plasticity "  << yield << std::endl;
        KRATOS_ERROR_IF(yield == "Rankine") << "Rankine yield surface not available in plasticity " << yield << std::endl;

       if (yield == "VonMises") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <VonMisesYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <VonMisesYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <VonMisesYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <VonMisesYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong: " << potential <<std::endl;
            }
        } else if (yield == "ModifiedMohrCoulomb") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <ModifiedMohrCoulombYieldSurface
                            <VonMisesPlasticPotential>>>
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <ModifiedMohrCoulombYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <ModifiedMohrCoulombYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <ModifiedMohrCoulombYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong: " << potential <<std::endl;
            }
        } else if (yield == "Tresca") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <TrescaYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <TrescaYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <TrescaYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone() 
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <TrescaYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong: " << potential <<std::endl;
            }
        } else if (yield == "DruckerPrager") {
            if (potential == "VonMises") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <DruckerPragerYieldSurface
                            <VonMisesPlasticPotential>>>().Clone()
                ;
            } else if (potential == "ModifiedMohrCoulomb") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <DruckerPragerYieldSurface
                            <ModifiedMohrCoulombPlasticPotential>>>().Clone()
                ;
            } else if (potential == "DruckerPrager") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <DruckerPragerYieldSurface
                            <DruckerPragerPlasticPotential>>>().Clone()
                ;
            } else if (potential == "Tresca") {
                constitutive_law = GenericSmallStrainIsotropicPlasticity3D 
                    <GenericConstitutiveLawIntegratorPlasticity
                        <DruckerPragerYieldSurface
                            <TrescaPlasticPotential>>>().Clone()
                ;
            } else {
                KRATOS_ERROR << "Plastic Potential not defined or wrong: " << potential <<std::endl;
            }
        } else {
            KRATOS_ERROR << "Yield Surface not defined or wrong: " << yield << std::endl;
        }
        return constitutive_law;
    }
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
    {
        this->CalculateMaterialResponseCauchy(rValues);
    }

    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
    {
        KRATOS_ERROR << "The base dummy class has been called, defina a correct constitutive law" << std::endl;
    } 


    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
    ) override
    {
        KRATOS_ERROR << "The base dummy class has been called, defina a correct constitutive law" << std::endl;
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


    ///@}

}; // Class GenericYieldSurface

} // namespace kratos
#endif
