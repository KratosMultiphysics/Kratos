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

    void SetValue(const Variable<bool>& rVariable,
                          const bool& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;

    void SetValue(const Variable<double>& rThisVariable,
                          const double& rValue,
                          const ProcessInfo& rCurrentProcessInfo) override;

    double& GetValue(const Variable<double>& rThisVariable,
                    double& rValue) override;

    bool& GetValue(const Variable<bool>& rThisVariable,bool& rValue) override;

    array_1d<double, 3 > & GetValue(const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue) override;

    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                    const Variable<double>& rThisVariable,
                    double& rValue) override;

    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                const Variable<Vector>& rThisVariable,
                Vector& rValue) override;

    array_1d<double, 3 > & CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue) override;

    void CalculateMaterialResponse(const Vector& rStrainVector,
                                    const Matrix& rDeformationGradient,
                                    Vector& rStressVector,
                                    Matrix& rAlgorithmicTangent,
                                    const ProcessInfo& rCurrentProcessInfo,
                                    const Properties& rMaterialProperties,
                                    const GeometryType& rElementGeometry,
                                    const Vector& rShapeFunctionsValues,
                                    bool CalculateStresses = true,
                                    int CalculateTangent = true,
                                    bool SaveInternalVariables = true) override;

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
    bool CheckIfIsPlasticRegime(const Properties& rMaterialProperties,
        const double& rCurrentStress);

    void FinalizeNonLinearIteration(const Properties& rMaterialProperties,
					    const GeometryType& rElementGeometry,
					    const Vector& rShapeFunctionsValues,
					    const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const Properties& rMaterialProperties,
                        const GeometryType& rElementGeometry,
                        const Vector& rShapeFunctionsValues,
                        const ProcessInfo& rCurrentProcessInfo) override;


    /**
     * @brief This function checks if the predicted stress state is in the elastic regime
     * but was in the plastic regime in the previous non_linear iteration step
     */
    bool CheckPlasticIterationHistory() const
    {
        bool check_flag = false;
        if(this->mInElasticFlagVector[1] && !this->mInElasticFlagVector[0])
        {
            check_flag = true;
        }
        return check_flag;
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
    double mStressState = 0.0; // The current stress state

    //the following members are vectors where
    //[0] is the current non_linear iteration step
    //[1] is the previous non_linear iteration step
    BoundedVector<bool, 2> mInElasticFlagVector = ZeroVector(2); /// This flags tells if we are in a elastic or ineslastic regime
    BoundedVector<double, 2> mPlasticAlphaVector = ZeroVector(2); /// The current plastic increment
    BoundedVector<double, 2> mAccumulatedPlasticStrainVector = ZeroVector(2); /// The current accumulated plastic strain


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
        rSerializer.save("StressState", this->mStressState);
        rSerializer.save("PlasticAlpha", this->mPlasticAlphaVector);
        rSerializer.save("AccumulatedPlasticStrain", this->mAccumulatedPlasticStrainVector);
        rSerializer.save("InelasticFlag", this->mInElasticFlagVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw);
        rSerializer.load("StressState", this->mStressState);
        rSerializer.load("PlasticAlpha", this->mPlasticAlphaVector);
        rSerializer.load("AccumulatedPlasticStrain", this->mAccumulatedPlasticStrainVector);
        rSerializer.load("InelasticFlag", this->mInElasticFlagVector);
    }


}; // Class TrussPlasticityConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_DUMMY_TRUSS_LAW_H_INCLUDED  defined
