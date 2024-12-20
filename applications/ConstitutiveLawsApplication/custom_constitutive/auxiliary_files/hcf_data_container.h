//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo
//
//

# pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "constitutive_laws_application_variables.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class HCFDataContainer
 * @ingroup ConstitutiveLawsApplication
 * @brief Defining all the methods and variables required in fatigue simulations
 * @details Give access to methods and restores all member variables in fatigue simulations
 * @author Alireza Taherzadeh-Fard, Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) HCFDataContainer
{
public:

    struct FatigueVariables {
        double MaxStress = 0.0;
        double MinStress = 0.0;
        bool MaxIndicator = false;
        bool MinIndicator = false;
        Vector PreviousStresses = ZeroVector(2);
        double FatigueReductionFactor = 1.0;
        double ReversionFactorRelativeError = 0.0;
        double MaxStressRelativeError = 0.0;
        unsigned int GlobalNumberOfCycles = 1;
        unsigned int LocalNumberOfCycles = 1;
        double B0 = 0.0;
        double PreviousMaxStress = 0.0;
        double PreviousMinStress = 0.0;
        double WohlerStress = 1.0;
        double Sth = 0.0;
        double ReferenceDamage = 0.0;
        double CyclesToFailure = 0.0;
        bool NewCycle = false;
        double Alphat = 0.0;
        double PreviousReversionFactor = 0.0;
        double ReversionFactor = 0.0;
        bool AdvanceStrategyApplied = false;
        bool DamageActivation = false;
        bool current_load_type = false;
        bool new_model_part = false;
        bool is_initiated = false;
        double nc_initiation = 0.0;
        double fred_initiation = 1.0;
        Vector hcf_coefficients = ZeroVector(10);
    };

    static constexpr double tolerance = 1.0e-3;

    typedef std::size_t SizeType;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType VoigtSize = 6;

    HCFDataContainer()
    {};

    // Defining fatigue methods

    /**
     * @brief This method sets maximum and minimum stresses
     * for fatigue calculations
     */
    void CalculateSminAndSmax(const double CurrentStress,
                            HCFDataContainer::FatigueVariables &rFatigueVariables);

    /**
     * @brief This method identifies if the overall load state
     * is tension or compression
     */
    double CalculateTensionOrCompressionIdentifier(const Vector& rStressVector);

    /**
     * @brief This method calculates the reversion factor
     * based on the maximum and minimum stresses
     */
    double CalculateReversionFactor(const double MaxStress, const double MinStress);

    /**
     * @brief This method sets the variables required
     * for calculating fatigue reduction factor and Wohler stress
     */
    void CalculateFatigueParameters(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables, const Variable<double>& rVariable = YIELD_STRESS);

    /**
     * @brief This method calculates fatigue reduction factor
     * and Wohler stress
     */
    void CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables, const Variable<double>& rVariable = YIELD_STRESS);


    /**
     * @brief This method initializes all the values
     * in the FatigueVariables
     */
    void InitializeFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables, ConstitutiveLaw::Parameters& rValues);


    /**
     * @brief This method updates all member variables
     */
    void UpdateFatigueVariables(HCFDataContainer::FatigueVariables &rFatigueVariables);


    /**
     * @brief This method computes fatigue-related quantities
     */
    void FinalizeSolutionStep(HCFDataContainer::FatigueVariables &rFatigueVariables,
                            const Properties& rMaterialProperties,
                            const ProcessInfo& rCurrentProcessInfo,
                            ConstitutiveLaw::StressVectorType stress_vector,
                            double uniaxial_stress,
                            double damage,
                            double threshold,
                            const Variable<double>& rVariable = YIELD_STRESS,
                            const Variable<bool>& rMethodVariable = DEFAULT_METHOD);

    /**
     * @brief This method returns fatigue reduction factor
     */
    double GetFatigueReductionFactor() {
        return mFatigueReductionFactor;
    }

    /**
     * @brief This method returns fatigue reduction factor
     */
    double GetWohlerStress() {
        return mWohlerStress;
    }

    /**
     * @brief This method returns fatigue reduction factor
     */
    int GetLocalNumberOfCycles() {
        return mNumberOfCyclesLocal;
    }

    /**
     * @brief This method sets local number of cycles
     */
    void SetLocalNumberOfCycles(int rValue) {
        mNumberOfCyclesLocal = rValue;
    }

    /**
     * @brief This method returns the previous cycle time
     */
    double GetPreviousCycleTime() {
        return mPreviousCycleTime;
    }

    /**
     * @brief This method sets the previous cycle time
     */
    void SetPreviousCycleTime(double rValue) {
        mPreviousCycleTime = rValue;
    }

    /**
     * @brief This method sets the previous cycle damage
     */
    void SetPreviousCycleDamage(double rValue) {
        mPreviousCycleDamage = rValue;
    }

    /**
     * @brief This method returns the previous cycle damage
     */
    double GetPreviousCycleDamage() {
        return mPreviousCycleDamage;
    }

    /**
     * @brief This method returns the cycle period
     */
    double GetCyclePeriod() {
        return mPeriod;
    }

    /**
     * @brief This method sets the cycle period
     */
    void SetCyclePeriod(double rValue) {
        mPeriod = rValue;
    }

    /**
     * @brief This method returns maximum stress relative error
     */
    double GetMaxStressRelativeError() {
        return mMaxStressRelativeError;
    }

    /**
     * @brief This method returns reversion factor relative error
     */
    double GetReversionFactorRelativeError() {
        return mReversionFactorRelativeError;
    }

    /**
     * @brief This method returns the threshold stress
     */
    double GetThresholdStress() {
        return mThresholdStress;
    }

    /**
     * @brief This method returns the threshold stress
     */
    double GetReferenceDamage() {
        return mReferenceDamage;
    }

    /**
     * @brief This method returns the maximum stress
     */
    double GetMaximumStress() {
        return mMaxStress;
    }

    /**
     * @brief This method returns the maximum stress
     */
    double GetMinimumStress() {
        return mMinStress;
    }

    /**
     * @brief This method returns the maximum stress
     */
    double GetReversionFactor() {
        return mReversionFactor;
    }

    /**
     * @brief This method returns the maximum stress
     */
    double GetPreviousMaximumStress() {
        return mPreviousMaxStress;
    }

    /**
     * @brief This method returns cycles to failure
     */
    double GetCyclesToFailure() {
        return mCyclesToFailure;
    }

    /**
     * @brief This method returns AIT control counter
     */
    int GetAITControlCounter() {
        return mAITControlParameter;
    }

    /**
     * @brief This method returns the global number of cycles
     */
    int GetGlobalNumberOfCycles() {
        return mNumberOfCyclesGlobal;
    }

    /**
     * @brief This method returns the increment in number of cycles
     */
    int GetIncrementInNumberOfCycles() {
        return mIncrementInNumberOfCycles;
    }

    /**
     * @brief This method sets global number of cycles
     */
    void SetIncrementInNumberOfCycles(int rValue) {
        mIncrementInNumberOfCycles = rValue;
    }

    /**
     * @brief This method sets global number of cycles
     */
    void SetGlobalNumberOfCycles(int rValue) {
        mNumberOfCyclesGlobal = rValue;
    }

    /**
     * @brief This method returns the new cycle indicator
     */
    bool GetNewCycleIndicator() {
        return mNewCycleIndicator;
    }


    /**
     * @brief This method checks the fatigue inout properties
     */
    int Check(const Properties& rMaterialProperties);



private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{
    double mMaxStress = 0.0;
    double mMinStress = 0.0;
    bool mMaxDetected = false;
    bool mMinDetected = false;
    Vector mPreviousStresses = ZeroVector(2);
    double mFatigueReductionFactor = 1.0;
    double mReversionFactorRelativeError = 0.0;
    double mMaxStressRelativeError = 0.0;
    unsigned int mNumberOfCyclesGlobal = 1;
    unsigned int mNumberOfCyclesLocal = 1;
    unsigned int mIncrementInNumberOfCycles = 0;
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    double mWohlerStress = 1.0;
    double mThresholdStress = 0.0;
    double mReferenceDamage = 0.0;
    double mCyclesToFailure = std::numeric_limits<double>::infinity();
    bool mNewCycleIndicator = false;
    double mPreviousCycleTime = 0.0; // Instanced variable used in the advanciing process for the conversion between time and number of cycles.
    double mPeriod = 0.0; // Instanced variable used in the advanciing process for the conversion between time and number of cycles.
    double mReversionFactor = 0.0;
    unsigned int mAITControlParameter = 0;
    bool mIsInitiated = false;
    double mNcInitiation = 0.0;
    double mFredInitiation = 1.0;
    double mPreviousCycleDamage = 0.0;
    ///@}

friend class Serializer;

void save(Serializer& rSerializer) const
{
    rSerializer.save("MaxStress",mMaxStress);
    rSerializer.save("MinStress",mMinStress);
    rSerializer.save("MaxDetected",mMaxDetected);
    rSerializer.save("MinDetected",mMinDetected);
    rSerializer.save("PreviousStresses",mPreviousStresses);
    rSerializer.save("FatigueReductionFactor",mFatigueReductionFactor);
    rSerializer.save("ReversionFactorRelativeError",mReversionFactorRelativeError);
    rSerializer.save("MaxStressRelativeError",mMaxStressRelativeError);
    rSerializer.save("NumberOfCyclesGlobal",mNumberOfCyclesGlobal);
    rSerializer.save("NumberOfCyclesLocal",mNumberOfCyclesLocal);
    rSerializer.save("PreviousMaxStress",mPreviousMaxStress);
    rSerializer.save("PreviousMinStress",mPreviousMinStress);
    rSerializer.save("WohlerStress",mWohlerStress);
    rSerializer.save("ThresholdStress",mThresholdStress);
    rSerializer.save("CyclesToFailure",mCyclesToFailure);
    rSerializer.save("NewCycleIndicator",mNewCycleIndicator);
}

void load(Serializer& rSerializer)
{
    rSerializer.load("MaxStress",mMaxStress);
    rSerializer.load("MinStress",mMinStress);
    rSerializer.load("MaxDetected",mMaxDetected);
    rSerializer.load("MinDetected",mMinDetected);
    rSerializer.load("PreviousStresses",mPreviousStresses);
    rSerializer.load("FatigueReductionFactor",mFatigueReductionFactor);
    rSerializer.load("ReversionFactorRelativeError",mReversionFactorRelativeError);
    rSerializer.load("MaxStressRelativeError",mMaxStressRelativeError);
    rSerializer.load("NumberOfCyclesGlobal",mNumberOfCyclesGlobal);
    rSerializer.load("NumberOfCyclesLocal",mNumberOfCyclesLocal);
    rSerializer.load("PreviousMaxStress",mPreviousMaxStress);
    rSerializer.load("PreviousMinStress",mPreviousMinStress);
    rSerializer.load("WohlerStress",mWohlerStress);
    rSerializer.load("ThresholdStress",mThresholdStress);
    rSerializer.load("CyclesToFailure",mCyclesToFailure);
    rSerializer.load("NewCycleIndicator",mNewCycleIndicator);
}

}; // class
///@}

///@} addtogroup block

} // namespace Kratos
