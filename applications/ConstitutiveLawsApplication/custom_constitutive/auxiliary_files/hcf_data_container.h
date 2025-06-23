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
        double CyclesToFailure = 0.0;
        bool NewCycle = false;
        double Alphat = 0.0;
        double PreviousReversionFactor = 0.0;
        double ReversionFactor = 0.0;
        bool AdvanceStrategyApplied;
        bool DamageActivation;
        double CFactor = 1.0;
        double UltimateStress = 1.0;
    };

    static constexpr double tolerance = 1.0e-3;

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
     * @brief This method calculates the ultimate stress
     * depending on the used softening law
     */
    double UltimateStressDamage(const Properties& rMaterialParameters);

    /**
     * @brief This method sets the variables required
     * for calculating fatigue reduction factor and Wohler stress
     */
    void CalculateFatigueParameters(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables);

    /**
     * @brief This method calculates fatigue reduction factor
     * and Wohler stress
     */
    void CalculateFatigueReductionFactorAndWohlerStress(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables);


    /**
     * @brief This method initializes all the values
     * in the FatigueVariables
     */
    void InitializeFatigueVariables(const Properties& rMaterialParameters, HCFDataContainer::FatigueVariables &rFatigueVariables);


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
                            double uniaxial_stress);

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
    double mPreviousMaxStress = 0.0;
    double mPreviousMinStress = 0.0;
    double mWohlerStress = 1.0;
    double mThresholdStress = 0.0;
    double mCyclesToFailure = std::numeric_limits<double>::infinity();
    bool mNewCycleIndicator = false;
    double mCFactor = 1.0; // Fatigue reduction factor smoothness term
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
    rSerializer.save("CFactor",mCFactor);
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
    rSerializer.load("CFactor",mCFactor);
}

}; // class
///@}

///@} addtogroup block

} // namespace Kratos
