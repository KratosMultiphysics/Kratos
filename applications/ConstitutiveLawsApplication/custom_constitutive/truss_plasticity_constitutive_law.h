// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#if !defined (KRATOS_TRUSS_PLASTICITY_LAW_H_INCLUDED)
#define  KRATOS_TRUSS_PLASTICITY_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "includes/checks.h"

namespace Kratos
{

/**
 * @namespace TrussPlasticityConstitutiveLaw
 *
 * @brief This constitutive law represents a linear hardening plasticity 1D law
 *
 * @author Klaus B Sautter
 */


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TrussPlasticityConstitutiveLaw : public ConstitutiveLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of TrussPlasticityConstitutiveLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( TrussPlasticityConstitutiveLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    TrussPlasticityConstitutiveLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    TrussPlasticityConstitutiveLaw (const TrussPlasticityConstitutiveLaw& rOther);


    /**
     * Destructor.
     */
    ~TrussPlasticityConstitutiveLaw() override;

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

    void SetValue(const Variable<double>& rThisVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;

    double& GetValue(const Variable<double>& rThisVariable,
                    double& rValue) override;

    array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                    const Variable<double>& rThisVariable,
                    double& rValue) override;

    /**
     * @brief This function calculates the yield function for the plastic model
     * @param rMaterialProperties Material Properties of the problem
     * @param rCurrentStress Current stress value
     */
    double TrialYieldFunction(const Properties& rMaterialProperties,
     const double& rCurrentStress);

    /**
     * @brief This function checks if the current stress is in the plastic regime
     * @param rMaterialProperties Material Properties of the problem
     * @param rCurrentStress Current stress value
     */
    bool CheckIfIsPlasticRegime(Parameters& rValues,const double& rCurrentStress);

    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    void CalculateMaterialResponsePK2Custom(Parameters& rValues, double& rCurrentAccumulatedPlasticStrain, double& rCurrentPlasticAlpha);

    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 1;
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
    ) override;

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
    //the following members are
    bool mCurrentInElasticFlag = false;/// This flags tells if we are in a elastic or ineslastic regime
    double mPlasticAlpha = 0.0;
    double mAccumulatedPlasticStrain = 0.0;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.save("PlasticAlpha", mPlasticAlpha);
        rSerializer.save("AccumulatedPlasticStrain", mAccumulatedPlasticStrain);
        rSerializer.save("CurrentInElasticFlag", mCurrentInElasticFlag);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.load("PlasticAlpha", mPlasticAlpha);
        rSerializer.load("AccumulatedPlasticStrain", mAccumulatedPlasticStrain);
        rSerializer.load("CurrentInElasticFlag", mCurrentInElasticFlag);
    }


}; // Class TrussPlasticityConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_TRUSS_PLASTICITY_LAW_H_INCLUDED  defined
