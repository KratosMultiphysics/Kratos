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
 * @class TimoshenkoBeamElasticConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class is used by the Timoshenko beam elements in such a way that the CL computes the axial, bending and shear forces according to elastic relations: 
 * N = EA E_l
 * M = EI Kappa
 * V = GAs Gamma_xy
 * As well as its derivatives dN_dE_l, dM_dKappa and dV_dGamma_xy
 * @details This means that the input strain is size 3: (E_l, kappa, gamma_xy); output stress vector is size 3: (N, M, V); The constitutive law also retrieves the derivatives via the CalculateValue methods.
 * @author Alejandro Cornejo
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TimoshenkoBeamElasticConstitutiveLaw : public BeamConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    using ProcessInfoType = ProcessInfo;
    using BaseType = BeamConstitutiveLaw;
    using SizeType = std::size_t;
    /**
     * Counted pointer of TimoshenkoBeamElasticConstitutiveLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(TimoshenkoBeamElasticConstitutiveLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    TimoshenkoBeamElasticConstitutiveLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    TimoshenkoBeamElasticConstitutiveLaw(const TimoshenkoBeamElasticConstitutiveLaw &rOther);

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
        return 3;
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

}; // Class TimoshenkoBeamElasticConstitutiveLaw
}  // namespace Kratos.
