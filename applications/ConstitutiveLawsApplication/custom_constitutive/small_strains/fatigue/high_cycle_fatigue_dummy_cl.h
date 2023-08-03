// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/small_strains/fatigue/hcf_data_container.h"

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
 * @class HighCycleFatigueDummyCl
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for plane strain cases
 * @details This class derives from the linear elastic case on 3D
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) HighCycleFatigueDummyCl
    : public ElasticIsotropic3D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of HighCycleFatigueDummyCl
    KRATOS_CLASS_POINTER_DEFINITION( HighCycleFatigueDummyCl );

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HighCycleFatigueDummyCl();

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HighCycleFatigueDummyCl (const HighCycleFatigueDummyCl& rOther);


    /**
     * @brief Destructor.
     */
    ~HighCycleFatigueDummyCl() override;

    // /**
    //  * @brief It calculates the stress vector
    //  * @param rStrainVector The strain vector in Voigt notation
    //  * @param rStressVector The stress vector in Voigt notation
    //  * @param rValues Parameters of the constitutive law
    //  */
    // void CalculatePK2Stress(
    //     const ConstitutiveLaw::StrainVectorType& rStrainVector,
    //     ConstitutiveLaw::StressVectorType& rStressVector,
    //     ConstitutiveLaw::Parameters& rValues
    //     ) override;

     /**
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

     /**
     * @brief Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief If the CL requires to finalize the material response, called by the element in FinalizeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }



    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    HCFDataContainer mFatigueData = HCFDataContainer();
    Vector mPreviousStresses = ZeroVector(2); // [S_t-2, S_t-1]
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    bool mMaxDetected = false;
    bool mMinDetected = false;
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

}; // Class HighCycleFatigueDummyCl
}  // namespace Kratos.
