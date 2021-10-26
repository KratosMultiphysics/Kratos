// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Sergio Jimenez/Alejandro Cornejo/Lucia Barbu
//  Collaborator:
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_HIGH_CYCLE_FATIGUE_LAW_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_HIGH_CYCLE_FATIGUE_LAW_H_INCLUDED
// System includes

// External includes

// Project includes
#include "custom_constitutive/generic_small_strain_isotropic_damage.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

// The size type definition
typedef std::size_t SizeType;

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
 * @class GenericSmallStrainHighCycleFatigueLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the base class which defines the constitutive law used for high cycle fatigue (HCF) in small deformation
 * @details This class uses the GenericSmallStrainIsotropicDamage class once the load is applied Nf cycles. The code has been written following the approach proposed by S. Oller et al. in A continuum mechanics model for mechanical fatigue analysis (2005)
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author Sergio Jiménez/Alejandro Cornejo/Lucia Barbu
 */
template <class TConstLawIntegratorType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) GenericSmallStrainHighCycleFatigueLaw
    : public GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>
{
public:
    ///@name Type Definitions
    ///@{

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TConstLawIntegratorType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TConstLawIntegratorType::VoigtSize;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainHighCycleFatigueLaw);

    /// The node definition
    typedef Node<3> NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    /// Definition of the base class
    typedef GenericSmallStrainIsotropicDamage<TConstLawIntegratorType> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainHighCycleFatigueLaw()
    {
    }

    GenericSmallStrainHighCycleFatigueLaw(  const double FatigueReductionFactor,
                                            const double PreviousStress0,
                                            const double PreviousStress1,
                                            const double MaxStress,
                                            const double MinStress,
                                            const unsigned int NumberOfCyclesGlobal,
                                            const double FatigueReductionParameter)
    {
        mFatigueReductionFactor = FatigueReductionFactor;
        Vector PreviousStresses = ZeroVector(2);
        PreviousStresses[0] = PreviousStress0;
        PreviousStresses[1] = PreviousStress1;
        mPreviousStresses = PreviousStresses;
        mMaxStress = MaxStress;
        mMinStress = MinStress;
        mNumberOfCyclesGlobal = NumberOfCyclesGlobal;
        mFatigueReductionParameter = FatigueReductionParameter;
    }
    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainHighCycleFatigueLaw<TConstLawIntegratorType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainHighCycleFatigueLaw(const GenericSmallStrainHighCycleFatigueLaw &rOther)
        : GenericSmallStrainIsotropicDamage<TConstLawIntegratorType>(rOther),
            mFatigueReductionFactor(rOther.mFatigueReductionFactor),
            mPreviousStresses(rOther.mPreviousStresses),
            mMaxStress(rOther.mMaxStress),
            mMinStress(rOther.mMinStress),
            mPreviousMaxStress(rOther.mPreviousMaxStress),
            mPreviousMinStress(rOther.mPreviousMinStress),
            mNumberOfCyclesGlobal(rOther.mNumberOfCyclesGlobal),
            mNumberOfCyclesLocal(rOther.mNumberOfCyclesLocal),
            mFatigueReductionParameter(rOther.mFatigueReductionParameter),
            mStressVector(rOther.mStressVector),
            mMaxDetected(rOther.mMaxDetected),
            mMinDetected(rOther.mMinDetected),
            mWohlerStress(rOther.mWohlerStress)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainHighCycleFatigueLaw() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

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
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * returns whether this constitutive Law has specified variable
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    using ConstitutiveLaw::Has;
    bool Has(const Variable<int>& rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (bool)
     * @param rThisVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<bool>& rThisVariable,
        const bool& Value,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    using ConstitutiveLaw::SetValue;
    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(
        const Variable<bool>& rThisVariable,
        bool& rValue) override;

    /**
     * @brief Returns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    using ConstitutiveLaw::GetValue;
    int& GetValue(
        const Variable<int>& rThisVariable,
        int& rValue) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return true;
    }

        /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

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
    ///@name Access
    ///@{
    Vector GetStressVector() {return mStressVector;}
    void SetStressVector(const Vector toStressVector) {mStressVector = toStressVector;}
    ///@}
    ///@name Member Variables
    ///@{
    double mFatigueReductionFactor = 1.0;
    Vector mPreviousStresses = ZeroVector(2); // [S_t-2, S_t-1]
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    unsigned int mNumberOfCyclesGlobal = 1; // Total number of cycles in the whole analysis
    unsigned int mNumberOfCyclesLocal = 1; // Equivalent number of cycles for the current cyclic load
    double mFatigueReductionParameter = 0.0; // B0
    Vector mStressVector = ZeroVector(VoigtSize);
    bool mMaxDetected = false; // Maximum's indicator in the current cycle
    bool mMinDetected = false; // Minimum's indicator in the current cycle
    double mWohlerStress = 1.0; // Normalised Wohler stress required for building the life prediction curves (SN curves)
    double mThresholdStress = 0.0; // Endurance limit of the fatigue model.
    double mReversionFactorRelativeError = 0.0; // Relative error of the R = Smin / Smax between cycles inducing recalculation of Nlocal and advanciing process.
    double mMaxStressRelativeError = 0.0; // Relative error of Smax between cycles inducing recalculation of Nlocal and advanciing process.
    bool mNewCycleIndicator = false; // New cycle identifier required for the advancing process.
    double mCyclesToFailure = 0.0; // Nf. Required for the advanciing process.
    double mPreviousCycleTime = 0.0; // Instanced variable used in the advanciing process for the conversion between time and number of cycles.
    double mPeriod = 0.0; // Instanced variable used in the advanciing process for the conversion between time and number of cycles.

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}

    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("FatigueReductionFactor", mFatigueReductionFactor);
        rSerializer.save("PreviousStresses", mPreviousStresses);
        rSerializer.save("MaxStress", mMaxStress);
        rSerializer.save("MinStress", mMinStress);
        rSerializer.save("PreviousMaxStress", mPreviousMaxStress);
        rSerializer.save("PreviousMinStress", mPreviousMinStress);
        rSerializer.save("NumberOfCyclesGlobal", mNumberOfCyclesGlobal);
        rSerializer.save("NumberOfCyclesLocal", mNumberOfCyclesLocal);
        rSerializer.save("FatigueReductionParameter", mFatigueReductionParameter);
        rSerializer.save("StressVector", mStressVector);
        rSerializer.save("MaxDetected", mMaxDetected);
        rSerializer.save("MinDetected", mMinDetected);
        rSerializer.save("WohlerStress", mWohlerStress);
        rSerializer.save("ThresholdStress", mThresholdStress);
        rSerializer.save("ReversionFactorRelativeError", mReversionFactorRelativeError);
        rSerializer.save("MaxStressRelativeError", mMaxStressRelativeError);
        rSerializer.save("NewCycleIndicator", mNewCycleIndicator);
        rSerializer.save("CyclesToFailure", mCyclesToFailure);
        rSerializer.save("PreviousCycleTime", mPreviousCycleTime);
        rSerializer.save("Period", mPeriod);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("FatigueReductionFactor", mFatigueReductionFactor);
        rSerializer.load("PreviousStresses", mPreviousStresses);
        rSerializer.load("MaxStress", mMaxStress);
        rSerializer.load("MinStress", mMinStress);
        rSerializer.load("PreviousMaxStress", mPreviousMaxStress);
        rSerializer.load("PreviousMinStress", mPreviousMinStress);
        rSerializer.load("NumberOfCyclesGlobal", mNumberOfCyclesGlobal);
        rSerializer.load("NumberOfCyclesLocal", mNumberOfCyclesLocal);
        rSerializer.load("FatigueReductionParameter", mFatigueReductionParameter);
        rSerializer.load("StressVector", mStressVector);
        rSerializer.load("MaxDetected", mMaxDetected);
        rSerializer.load("MinDetected", mMinDetected);
        rSerializer.load("WohlerStress", mWohlerStress);
        rSerializer.load("ThresholdStress", mThresholdStress);
        rSerializer.load("ReversionFactorRelativeError", mReversionFactorRelativeError);
        rSerializer.load("MaxStressRelativeError", mMaxStressRelativeError);
        rSerializer.load("NewCycleIndicator", mNewCycleIndicator);
        rSerializer.load("CyclesToFailure", mCyclesToFailure);
        rSerializer.load("PreviousCycleTime", mPreviousCycleTime);
        rSerializer.load("Period", mPeriod);
    }
    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif