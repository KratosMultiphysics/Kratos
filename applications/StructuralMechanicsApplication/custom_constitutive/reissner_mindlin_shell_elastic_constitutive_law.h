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
#include "beam_constitutive_law.h"

namespace Kratos
{
/**
 * @class ReissnerMindlinShellElasticConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class is used by the Reissner-Mindlin shell elements in such a way that the CL computes the membrane, bending and shear generalized stresses according to elastic relations: 
 * Sm = Dm · Em (3 components)
 * Sb = Db · Kappa (3 components)
 * Ss = Ds · Gamma (2 components)
 * As well as the generalized constitutive matrices Dm, Db and Ds.
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ReissnerMindlinShellElasticConstitutiveLaw : public BeamConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    using ProcessInfoType = ProcessInfo;
    using BaseType = BeamConstitutiveLaw;
    using SizeType = std::size_t;

    /**
     * Counted pointer of ReissnerMindlinShellElasticConstitutiveLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(ReissnerMindlinShellElasticConstitutiveLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    ReissnerMindlinShellElasticConstitutiveLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    ReissnerMindlinShellElasticConstitutiveLaw(const ReissnerMindlinShellElasticConstitutiveLaw &rOther);

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
        // 3 membrane + 3 bending + 2 shear
        return 8;
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
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BeamConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BeamConstitutiveLaw)
    }

}; // Class ReissnerMindlinShellElasticConstitutiveLaw
}  // namespace Kratos.
