// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainHighCycleFatigueLaw
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
    //typedef typename GenericSmallStrainIsotropicDamage<TConstLawIntegratorType> BaseType;

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
                                            const double PreviousMaxStress,
                                            const double PreviousMinStress,
                                            const unsigned int NumberOfCycles,
                                            const double FatigueReductionParameter,
                                            const bool HasCountedCycle)
    {
        mFatigueReductionFactor = FatigueReductionFactor;
        Vector PreviousStresses = ZeroVector(2);
        PreviousStresses[0] = PreviousStress0;
        PreviousStresses[1] = PreviousStress1;
        mPreviousStresses = PreviousStresses;
        mMaxStress = MaxStress;
        mMinStress = MinStress;
        mPreviousMinStress = PreviousMinStress;
        mNumberOfCycles = NumberOfCycles;
        mFatigueReductionParameter = FatigueReductionParameter;
        mHasCountedCycle = HasCountedCycle;
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
            mNumberOfCycles(rOther.mNumberOfCycles), 
            mFatigueReductionParameter(rOther.mFatigueReductionParameter), 
            mStressVector(rOther.mStressVector), 
            mHasCountedCycle(rOther.mHasCountedCycle)
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
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

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
    bool Has(const Variable<int>& rThisVariable) override;

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
     * @brief Sets the value of a specified variable (integer)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     */
    void SetValue(
        const Variable<int>& rThisVariable, 
        const int& rValue);

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
     * @brief Returns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(
        const Variable<int>& rThisVariable, 
        int& rValue) override;

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
    double GetFatigueReductionFactor() {return mFatigueReductionFactor;}
    void SetFatigueReductionFactor(const double toFred) {mFatigueReductionFactor = toFred;}

    Vector GetPreviousStresses() {return mPreviousStresses;}
    void SetPreviousStresses(const Vector toPreviousStresses) {mPreviousStresses = toPreviousStresses;}

    double GetMaxStress() {return mMaxStress;}
    void SetMaxStress(const double toMaxStress) {mMaxStress = toMaxStress;}

    double GetMinStress() {return mMinStress;}
    void SetMinStress(const double toMinStress) {mMinStress = toMinStress;}

    double GetPreviousMaxStress() {return mPreviousMaxStress;}
    void SetPreviousMaxStress(const double toPreviousMaxStress) {mPreviousMaxStress = toPreviousMaxStress;}

    double GetPreviousMinStress() {return mPreviousMinStress;}
    void SetPreviousMinStress(const double toPreviousMinStress) {mPreviousMinStress = toPreviousMinStress;}

    unsigned int GetNumberOfCycles() {return mNumberOfCycles;}
    void AddCycle() {mNumberOfCycles++;}
    void SetNumberOfCycles(const unsigned int toCycles) {mNumberOfCycles = toCycles;}

    double GetFatigueReductionParameter() {return mFatigueReductionParameter;}
    void SetFatigueReductionParameter(const double toFatigueReductionParameter) {mFatigueReductionParameter = toFatigueReductionParameter;}

    Vector GetStressVector() {return mStressVector;}
    void SetStressVector(const Vector toStressVector) {mStressVector = toStressVector;}

    void ResetCycleCounter(){mHasCountedCycle = false;}
    void SetCycleCounter(const bool tocycle){mHasCountedCycle = tocycle;}
    bool GetCycleCounter() {return mHasCountedCycle;}
    
    ///@}
    ///@name Member Variables
    ///@{
    double mFatigueReductionFactor = 1.0;
    Vector mPreviousStresses = ZeroVector(2); // [S_t-2, S_t-1]
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    unsigned int mNumberOfCycles = 1;
    double mFatigueReductionParameter = 0.0; // B0
    Vector mStressVector = ZeroVector(VoigtSize);
    bool mHasCountedCycle = false;

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
        rSerializer.save("NumberOfCycles", mNumberOfCycles);
        rSerializer.save("FatigueReductionParameter", mFatigueReductionParameter);
        rSerializer.save("StressVector", mStressVector);
        rSerializer.save("HasCountedCycle", mHasCountedCycle);
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
        rSerializer.load("NumberOfCycles", mNumberOfCycles);
        rSerializer.load("FatigueReductionParameter", mFatigueReductionParameter);
        rSerializer.load("StressVector", mStressVector);
        rSerializer.load("HasCountedCycle", mHasCountedCycle);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
