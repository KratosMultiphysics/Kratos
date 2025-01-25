// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// External includes

// Project includes
#include "timoshenko_beam_elastic_constitutive_law.h"

namespace Kratos
{
/**
 * @class TimoshenkoBeamPlaneStrainElasticConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class is used by the Timoshenko beam elements in such a way that the CL computes the axial, bending and shear forces according to elastic relations but assuming plane strain conditions: 
 * ^ y'
 * |
 * O---------O -> x'
 * 
 * The following elastic relations come from the genralized Hookean elastic CL, assuming Sigma_y = 0.

 * Nx = EA * E_l / (1-nu^2)
 * Mz = EI * Kappa / (1-nu^2)
 * Vxy = GAs * Gamma_xy
 * Nz = nu * Nx
 * Mx = nu * Mz

 * being:
 * E: Young's modulus
 * A: Cross section area
 * As: Reduced shear cross section area
 * I: Inertia
 * G: Shear modulus
 * nu: Poisson ratio

 * As well as its derivatives dNx_dE_l, dMz_dKappa, dVxy_dGamma_xy, dNz_dE_l and dMx_dKappa
 * @details This means that the input strain is size 5: (E_l, kappa, gamma_xy, 0, 0); output stress vector is size 5: (Nx, Mz, Vxy, Nz, Mx); The constitutive law also retrieves the derivatives via the CalculateValue methods.
 * The beam is supposed to be in the x-y plane. Rotations along the z-axis
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TimoshenkoBeamPlaneStrainElasticConstitutiveLaw : public TimoshenkoBeamElasticConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    using BaseType = TimoshenkoBeamElasticConstitutiveLaw;
    using SizeType = std::size_t;
    /**
     * Counted pointer of TimoshenkoBeamPlaneStrainElasticConstitutiveLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(TimoshenkoBeamPlaneStrainElasticConstitutiveLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    TimoshenkoBeamPlaneStrainElasticConstitutiveLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    TimoshenkoBeamPlaneStrainElasticConstitutiveLaw(const TimoshenkoBeamPlaneStrainElasticConstitutiveLaw &rOther);

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() const override
    {
        return 5;
    }

    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

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


    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }

}; // Class TimoshenkoBeamPlaneStrainElasticConstitutiveLaw
}  // namespace Kratos.
