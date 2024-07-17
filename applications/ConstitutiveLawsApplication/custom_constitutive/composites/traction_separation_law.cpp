// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:                 BSD License
//                               license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo Velazquez
//                   Sergio Jimenez Reyes
//                   Lucia Gratiela Barbu
//
// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_constitutive/composites/traction_separation_law.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
// #include <Eigen/Dense>

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D()
    : BaseType()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D(const std::vector<double>& rCombinationFactors) : BaseType(rCombinationFactors)
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::TractionSeparationLaw3D(const TractionSeparationLaw3D<TDim>& rOther)
    : BaseType(rOther),
      mDelaminationDamageModeOne(rOther.mDelaminationDamageModeOne),
      mDelaminationDamageModeTwo(rOther.mDelaminationDamageModeTwo),
      mThresholdModeOne(rOther.mThresholdModeOne),
      mThresholdModeTwo(rOther.mThresholdModeTwo),
      mFatigueDataContainersModeOne(rOther.mFatigueDataContainersModeOne),
      mFatigueDataContainersModeTwo(rOther.mFatigueDataContainersModeTwo)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::Pointer TractionSeparationLaw3D<TDim>::Clone() const
{
    return Kratos::make_shared<TractionSeparationLaw3D>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

template<unsigned int TDim>
ConstitutiveLaw::Pointer TractionSeparationLaw3D<TDim>::Create(Kratos::Parameters NewParameters) const
{
    // We do some checks
    KRATOS_ERROR_IF_NOT(NewParameters.Has("combination_factors")) << "TractionSeparationLaw3D: Please define combination_factors" << std::endl;

    const SizeType number_of_factors = NewParameters["combination_factors"].size();

    // We create the vectors
    std::vector<double> combination_factors(number_of_factors);

    for (IndexType i_layer = 0; i_layer < number_of_factors; ++i_layer) {
        combination_factors[i_layer] = NewParameters["combination_factors"][i_layer].GetDouble();
    }

    KRATOS_ERROR_IF(number_of_factors == 0) << "Please define the combination factors" << std::endl;

    // We create the law
    return Kratos::make_shared<TractionSeparationLaw3D>(combination_factors);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

template<unsigned int TDim>
TractionSeparationLaw3D<TDim>::~TractionSeparationLaw3D()
{
};

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::Has(const Variable<bool>& rThisVariable)
{
    bool has = false;

    if (rThisVariable == CYCLE_INDICATOR) {
        has = true;
    } else {
        BaseType::Has(rThisVariable);
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::Has(const Variable<int>& rThisVariable)
{
    bool has = false;

    if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {
        has = true;
    } else if (rThisVariable == NUMBER_OF_CYCLES) {
        has = true;
    } else if (rThisVariable == INCREMENT_IN_NUMBER_OF_CYCLES) {
        has = true;
    } else if (rThisVariable == AIT_CONTROL_COUNTER) {
        has = true;
    } else {
        BaseType::Has(rThisVariable);
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::Has(const Variable<double>& rThisVariable)
{
    bool has = false;

    if (rThisVariable == DAMAGE) {
        has = true;
    } else if (rThisVariable == REFERENCE_DAMAGE) {
        has = true;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {
        has = true;
    } else if (rThisVariable == PREVIOUS_CYCLE) {
        has = true;
    } else if (rThisVariable == CYCLE_PERIOD) {
        has = true;
    } else if (rThisVariable == MAX_STRESS_RELATIVE_ERROR) {
        has = true;
    } else if (rThisVariable == REVERSION_FACTOR_RELATIVE_ERROR) {
        has = true;
    } else if (rThisVariable == THRESHOLD_STRESS) {
        has = true;
    } else if (rThisVariable == MAX_STRESS) {
        has = true;
    } else if (rThisVariable == CYCLES_TO_FAILURE) {
        has = true;
    } else if (rThisVariable == WOHLER_STRESS) {
        has = true;
    } else {
        BaseType::Has(rThisVariable);
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::Has(const Variable<Vector>& rThisVariable)
{
    bool has = false;

    if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_ONE) {
        has = true;
    } else if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_TWO) {
        has = true;
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR_VECTOR_MODE_ONE) {
        has = true;
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR_VECTOR_MODE_TWO) {
        has = true;
    } else if (rThisVariable == WOHLER_STRESS_VECTOR_MODE_ONE) {
        has = true;
    } else if (rThisVariable == WOHLER_STRESS_VECTOR_MODE_TWO) {
        has = true;
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES_MODE_ONE) {
        has = true;
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES_MODE_TWO) {
        has = true;
    }
      else {
        BaseType::Has(rThisVariable);
    }

    return has;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if (rThisVariable == CYCLE_INDICATOR) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetNewCycleIndicator();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetNewCycleIndicator();
        }

        // bool rValueII = mFatigueDataContainersModeTwo[0].GetNewCycleIndicator();
        // bool rValueI = mFatigueDataContainersModeOne[0].GetNewCycleIndicator();    //Change the loading mode here
        // rValue = rValueII || rValueI;
        return rValue;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetLocalNumberOfCycles();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetLocalNumberOfCycles();
        }
        // rValue = mFatigueDataContainersModeTwo[0].GetLocalNumberOfCycles();
        // rValue = mFatigueDataContainersModeOne[0].GetLocalNumberOfCycles();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == NUMBER_OF_CYCLES) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetGlobalNumberOfCycles();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetGlobalNumberOfCycles();
        }
        // double rValueII = mFatigueDataContainersModeTwo[0].GetGlobalNumberOfCycles();
        // double rValueI = mFatigueDataContainersModeOne[0].GetGlobalNumberOfCycles();    //Change the loading mode here
        // if (rValueII > rValueI) {
        //     rValue = rValueI;
        // } else {
        //     rValue = rValueII;
        // }
        return rValue;
    } else if (rThisVariable == INCREMENT_IN_NUMBER_OF_CYCLES) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetIncrementInNumberOfCycles();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetIncrementInNumberOfCycles();
        }
        return rValue;
    } else if (rThisVariable == AIT_CONTROL_COUNTER) {

        rValue = mFatigueDataContainersModeTwo[0].GetAITControlCounter();
        // rValue = mFatigueDataContainersModeOne[0].GetAITControlCounter();    //Change the loading mode here
        return rValue;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // if (rThisVariable == DAMAGE) {

    //     if (mDelaminationDamageModeTwo[1] > 0.0) {
    //         rValue = mDelaminationDamageModeTwo[1];
    //     } else if (mDelaminationDamageModeOne[1] > 0.0) {
    //         rValue = mDelaminationDamageModeOne[1];
    //     } else {
    //         rValue = 0.0;
    //     }
    //     // rValue = mDelaminationDamageModeTwo[1];
    //     // rValue = mDelaminationDamageModeOne[1];    //Change the loading mode here
    //     return rValue;
    // } else 
    if (rThisVariable == REFERENCE_DAMAGE) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetReferenceDamage();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetReferenceDamage();
        }
        // if ((1.0 - mFatigueDataContainersModeTwo[0].GetReferenceDamage()) * mFatigueDataContainersModeTwo[0].GetMaximumStress() > mFatigueDataContainersModeTwo[0].GetThresholdStress()
        //     && (1.0 - mFatigueDataContainersModeOne[0].GetReferenceDamage()) * mFatigueDataContainersModeOne[0].GetMaximumStress() > mFatigueDataContainersModeOne[0].GetThresholdStress()) 
        // {
        //     rValue = mFatigueDataContainersModeTwo[0].GetReferenceDamage(); //It is optional. Could be mode one or mode two.
        // } else {
        //     rValue = 0.0;
        // }

        // rValue = mFatigueDataContainersModeTwo[0].GetReferenceDamage();
        // rValue = mFatigueDataContainersModeOne[0].GetReferenceDamage();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {

        rValue = mFatigueDataContainersModeTwo[0].GetPreviousCycleDamage();
        // rValue = mFatigueDataContainersModeOne[0].GetPreviousCycleDamage();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == PREVIOUS_CYCLE) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetPreviousCycleTime();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetPreviousCycleTime();
        }

        // rValue = mFatigueDataContainersModeTwo[0].GetPreviousCycleTime(); //It is optional. Could be mode one or mode two.
        // rValue = mFatigueDataContainersModeOne[0].GetPreviousCycleTime();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == CYCLE_PERIOD) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetCyclePeriod();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetCyclePeriod();
        }

        // rValue = mFatigueDataContainersModeTwo[0].GetCyclePeriod(); //It is optional. Could be mode one or mode two.
        // rValue = mFatigueDataContainersModeOne[0].GetCyclePeriod();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == MAX_STRESS_RELATIVE_ERROR) {

        if (mFatigueDataContainersModeTwo[0].GetMaxStressRelativeError() > mFatigueDataContainersModeOne[0].GetMaxStressRelativeError()) {
            rValue = mFatigueDataContainersModeTwo[0].GetMaxStressRelativeError();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetMaxStressRelativeError();
        }
        // rValue = mFatigueDataContainersModeTwo[0].GetMaxStressRelativeError();
        // rValue = mFatigueDataContainersModeOne[0].GetMaxStressRelativeError();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == REVERSION_FACTOR_RELATIVE_ERROR) {

        if (mFatigueDataContainersModeTwo[0].GetReversionFactorRelativeError() > mFatigueDataContainersModeOne[0].GetReversionFactorRelativeError()) {
            rValue = mFatigueDataContainersModeTwo[0].GetReversionFactorRelativeError();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetReversionFactorRelativeError();
        }
        // rValue = mFatigueDataContainersModeTwo[0].GetReversionFactorRelativeError();
        // rValue = mFatigueDataContainersModeOne[0].GetReversionFactorRelativeError();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == THRESHOLD_STRESS) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetThresholdStress();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetThresholdStress();
        }

        // if ((1.0 - mFatigueDataContainersModeTwo[0].GetReferenceDamage()) * mFatigueDataContainersModeTwo[0].GetMaximumStress() > mFatigueDataContainersModeTwo[0].GetThresholdStress()
        //     && (1.0 - mFatigueDataContainersModeOne[0].GetReferenceDamage()) * mFatigueDataContainersModeOne[0].GetMaximumStress() > mFatigueDataContainersModeOne[0].GetThresholdStress()) 
        // {
        //     rValue = mFatigueDataContainersModeTwo[0].GetThresholdStress(); //It is optional. Could be mode one or mode two.
        // } else {
        //     rValue = 0.0;
        // }

        // rValue = mFatigueDataContainersModeTwo[0].GetThresholdStress();
        // rValue = mFatigueDataContainersModeOne[0].GetThresholdStress();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == MAX_STRESS) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetMaximumStress();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetMaximumStress();
        }

        // if ((1.0 - mFatigueDataContainersModeTwo[0].GetReferenceDamage()) * mFatigueDataContainersModeTwo[0].GetMaximumStress() > mFatigueDataContainersModeTwo[0].GetThresholdStress()
        //     && (1.0 - mFatigueDataContainersModeOne[0].GetReferenceDamage()) * mFatigueDataContainersModeOne[0].GetMaximumStress() > mFatigueDataContainersModeOne[0].GetThresholdStress()) 
        // {
        //     rValue = mFatigueDataContainersModeTwo[0].GetMaximumStress(); //It is optional. Could be mode one or mode two.
        // } else {
        //     rValue = 0.0;
        // }


        // rValue = mFatigueDataContainersModeTwo[0].GetMaximumStress();
        // rValue = mFatigueDataContainersModeOne[0].GetMaximumStress();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == CYCLES_TO_FAILURE) {

        if (mFatigueLoadingStateParameter) {
            rValue = mFatigueDataContainersModeTwo[0].GetCyclesToFailure();
        } else {
            rValue = mFatigueDataContainersModeOne[0].GetCyclesToFailure();
        }

        // double N_residual_mode_one = mFatigueDataContainersModeOne[0].GetCyclesToFailure() - mFatigueDataContainersModeOne[0].GetLocalNumberOfCycles();
        // double N_residual_mode_two = mFatigueDataContainersModeTwo[0].GetCyclesToFailure() - mFatigueDataContainersModeTwo[0].GetLocalNumberOfCycles();

        // if (N_residual_mode_one > N_residual_mode_two) {
        //     rValue = mFatigueDataContainersModeTwo[0].GetCyclesToFailure();
        // } else {
        //     rValue = mFatigueDataContainersModeOne[0].GetCyclesToFailure();
        // }

        // rValue = mFatigueDataContainersModeTwo[0].GetCyclesToFailure();
        // rValue = mFatigueDataContainersModeOne[0].GetCyclesToFailure();    //Change the loading mode here
        return rValue;
    } else if (rThisVariable == WOHLER_STRESS) {

        rValue = mFatigueDataContainersModeTwo[0].GetReversionFactor();
        // rValue = mFatigueDataContainersModeOne[0].GetReversionFactor();    //Change the loading mode here
        return rValue;
    } else {
        return BaseType::GetValue(rThisVariable, rValue);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& TractionSeparationLaw3D<TDim>::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    const auto& r_combination_factors = this->GetCombinationFactors();

    rValue.clear();

    if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_ONE) {

        rValue.resize(6, false);
        Vector SortedmDelaminationDamageModeOne;
        SortedmDelaminationDamageModeOne.resize(r_combination_factors.size()+1, false);
        noalias(SortedmDelaminationDamageModeOne) = mDelaminationDamageModeOne;
        std::sort(SortedmDelaminationDamageModeOne.begin(), SortedmDelaminationDamageModeOne.end(),std::greater<>());
        for (int i = 0; i < 6; ++i) {
            rValue[i] = SortedmDelaminationDamageModeOne[i];
        }
        return rValue;
    } else if (rThisVariable == DELAMINATION_DAMAGE_VECTOR_MODE_TWO) {

        rValue.resize(6, false);
        Vector SortedmDelaminationDamageModeTwo;
        SortedmDelaminationDamageModeTwo.resize(r_combination_factors.size()+1, false);
        noalias(SortedmDelaminationDamageModeTwo) = mDelaminationDamageModeTwo;
        std::sort(SortedmDelaminationDamageModeTwo.begin(), SortedmDelaminationDamageModeTwo.end(),std::greater<>());
        for (int i = 0; i < 6; ++i) {
            rValue[i] = SortedmDelaminationDamageModeTwo[i];
        }
        // KRATOS_WATCH(rValue);
        return rValue;
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR_VECTOR_MODE_ONE) {

        rValue.resize(6, false); // Artificially changed to size 6 for printing purposes
        for (IndexType i = 0; i < r_combination_factors.size()-1; ++i) {
            rValue[i] = mFatigueDataContainersModeOne[i].GetFatigueReductionFactor();
        }
        return rValue;
    } else if (rThisVariable == FATIGUE_REDUCTION_FACTOR_VECTOR_MODE_TWO) {

        rValue.resize(6, false); // Artificially changed to size 6 for printing purposes
        for (IndexType i = 0; i < r_combination_factors.size()-1; ++i) {
            rValue[i] = mFatigueDataContainersModeTwo[i].GetFatigueReductionFactor();
        }
        return rValue;
    } else if (rThisVariable == WOHLER_STRESS_VECTOR_MODE_ONE) {

        rValue.resize(r_combination_factors.size()-1, false);
        for (IndexType i = 0; i < r_combination_factors.size()-1; ++i) {
            rValue[i] = mFatigueDataContainersModeOne[i].GetWohlerStress();
        }
        return rValue;
    } else if (rThisVariable == WOHLER_STRESS_VECTOR_MODE_TWO) {

        rValue.resize(r_combination_factors.size()-1, false);
        for (IndexType i = 0; i < r_combination_factors.size()-1; ++i) {
            rValue[i] = mFatigueDataContainersModeTwo[i].GetWohlerStress();
        }
        return rValue;
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES_MODE_ONE) {

        rValue.resize(r_combination_factors.size()-1, false);
        for (IndexType i = 0; i < r_combination_factors.size()-1; ++i) {
            rValue[i] = mFatigueDataContainersModeOne[i].GetLocalNumberOfCycles();
        }
        return rValue;
    } else if (rThisVariable == LOCAL_NUMBER_OF_CYCLES_MODE_TWO) {

        rValue.resize(r_combination_factors.size()-1, false);
        for (IndexType i = 0; i < r_combination_factors.size()-1; ++i) {
            rValue[i] = mFatigueDataContainersModeTwo[i].GetLocalNumberOfCycles();
        }
        return rValue;
    } else {
        return BaseType::GetValue(rThisVariable,rValue);
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == LOCAL_NUMBER_OF_CYCLES) {

        // if ((1.0 - mFatigueDataContainersModeTwo[0].GetReferenceDamage()) * mFatigueDataContainersModeTwo[0].GetMaximumStress() > mFatigueDataContainersModeTwo[0].GetThresholdStress()) {
        //     mFatigueDataContainersModeTwo[0].SetLocalNumberOfCycles(rValue);
        // }
        // if ((1.0 - mFatigueDataContainersModeOne[0].GetReferenceDamage()) * mFatigueDataContainersModeOne[0].GetMaximumStress() > mFatigueDataContainersModeOne[0].GetThresholdStress()) {
        //     mFatigueDataContainersModeOne[0].SetLocalNumberOfCycles(rValue);    //Change the loading mode here
        // }
    } else if (rThisVariable == NUMBER_OF_CYCLES) {

        // if ((1.0 - mFatigueDataContainersModeTwo[0].GetReferenceDamage()) * mFatigueDataContainersModeTwo[0].GetMaximumStress() > mFatigueDataContainersModeTwo[0].GetThresholdStress()) {
        //     mFatigueDataContainersModeTwo[0].SetGlobalNumberOfCycles(rValue);
        // }

        // if ((1.0 - mFatigueDataContainersModeOne[0].GetReferenceDamage()) * mFatigueDataContainersModeOne[0].GetMaximumStress() > mFatigueDataContainersModeOne[0].GetThresholdStress()) {
        //     mFatigueDataContainersModeOne[0].SetGlobalNumberOfCycles(rValue);    //Change the loading mode here
        // }
    } else if (rThisVariable == INCREMENT_IN_NUMBER_OF_CYCLES) {
        mFatigueDataContainersModeOne[0].SetIncrementInNumberOfCycles(rValue);
        mFatigueDataContainersModeTwo[0].SetIncrementInNumberOfCycles(rValue);
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rThisVariable == PREVIOUS_CYCLE) {

        mFatigueDataContainersModeTwo[0].SetPreviousCycleTime(rValue); //It is optional. Could be mode one or mode two.
        mFatigueDataContainersModeOne[0].SetPreviousCycleTime(rValue);    //Change the loading mode here
    } else if (rThisVariable == CYCLE_PERIOD) {

        if (mFatigueLoadingStateParameter) {
            mFatigueDataContainersModeTwo[0].SetCyclePeriod(rValue);
        } else {
            mFatigueDataContainersModeOne[0].SetCyclePeriod(rValue);
        }
        // mFatigueDataContainersModeTwo[0].SetCyclePeriod(rValue); //It is optional. Could be mode one or mode two.
        // mFatigueDataContainersModeOne[0].SetCyclePeriod(rValue);    //Change the loading mode here
    } else if (rThisVariable == PREVIOUS_CYCLE_DAMAGE) {

        mFatigueDataContainersModeTwo[0].SetPreviousCycleDamage(rValue);
        // mFatigueDataContainersModeOne[0].SetPreviousCycleDamage(rValue);    //Change the loading mode here
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double& TractionSeparationLaw3D<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if (rThisVariable == MODE_ONE_UNIAXIAL_STRESS || rThisVariable == MODE_TWO_UNIAXIAL_STRESS) {
        const int interface_identifier = 0;

        // Get Values to compute the constitutive law:
        Flags& r_flags = rValues.GetOptions();

        // Previous flags saved
        const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

        const Properties& r_material_properties = rValues.GetMaterialProperties();

        // The deformation gradient
        if (rValues.IsSetDeterminantF()) {
            const double determinant_f = rValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }

        // All the strains must be the same, therefore we can just simply compute the strain in the first layer
        if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateGreenLagrangeStrain(rValues);
            r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
        }

        // The global strain vector, constant
        const Vector strain_vector = rValues.GetStrainVector();

        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
        const auto& r_combination_factors = this->GetCombinationFactors();

        std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
        for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        // for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        //     rValue[i].resize(3, false);
        // }

        std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i].resize(6, false);
            for (IndexType j=0; j < 6; ++j) {
                interfacial_stress[i][j] = 0.0;
            }
        }

        for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {

            interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress
        }

        // double equivalent_stress_mode_one_fred = std::abs((layer_stress[interface_identifier][2] + layer_stress[interface_identifier + 1][2]) * 0.5);

        // Vector fatigue_interfacial_stress_vector_mode_one = ZeroVector(VoigtSize);
        // fatigue_interfacial_stress_vector_mode_one[2] = (layer_stress[interface_identifier][2] + layer_stress[interface_identifier + 1][2]) * 0.5;
        // const double sign_factor = mFatigueDataContainersModeOne[0].CalculateTensionOrCompressionIdentifier(fatigue_interfacial_stress_vector_mode_one);
        // equivalent_stress_mode_one_fred *= sign_factor;

        double sign_factor = 0.0;
        if (rThisVariable == MODE_ONE_UNIAXIAL_STRESS) {
            // rValue = interfacial_stress[interface_identifier][0];
            double equivalent_stress_mode_one_fred = std::abs((layer_stress[interface_identifier][2] + layer_stress[interface_identifier + 1][2]) * 0.5);

            Vector fatigue_interfacial_stress_vector_mode_one = ZeroVector(VoigtSize);
            fatigue_interfacial_stress_vector_mode_one[2] = (layer_stress[interface_identifier][2] + layer_stress[interface_identifier + 1][2]) * 0.5;
            sign_factor = mFatigueDataContainersModeOne[0].CalculateTensionOrCompressionIdentifier(fatigue_interfacial_stress_vector_mode_one);
            equivalent_stress_mode_one_fred *= sign_factor;

            rValue = equivalent_stress_mode_one_fred;
            return rValue;
        } else if (rThisVariable == MODE_TWO_UNIAXIAL_STRESS) {

            double aux_sum = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets(interfacial_stress[interface_identifier][1]) + AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets(interfacial_stress[interface_identifier][2]);
            double resultant_shear_stress = std::sqrt(std::pow(interfacial_stress[interface_identifier][1],2.0)+std::pow(interfacial_stress[interface_identifier][2],2.0));
            const double pre_indicator = aux_sum / resultant_shear_stress;

            if (pre_indicator < 0.7) {
                sign_factor = 1.0;
            } else {
                sign_factor = -1.0;
            }
            rValue = sign_factor * std::sqrt(std::pow(interfacial_stress[interface_identifier][1],2.0)+std::pow(interfacial_stress[interface_identifier][2],2.0));
            return rValue;
        }

    } else {
        return BaseType::CalculateValue(rValues, rThisVariable, rValue);
    }
    return rValue;

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
Vector& TractionSeparationLaw3D<TDim>::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == INTERFACIAL_STRESS_VECTOR) {
        const int interface_identifier = 0;

        // Get Values to compute the constitutive law:
        Flags& r_flags = rValues.GetOptions();

        // Previous flags saved
        const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

        const Properties& r_material_properties = rValues.GetMaterialProperties();

        // The deformation gradient
        if (rValues.IsSetDeterminantF()) {
            const double determinant_f = rValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }

        // All the strains must be the same, therefore we can just simply compute the strain in the first layer
        if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateGreenLagrangeStrain(rValues);
            r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
        }

        // The global strain vector, constant
        const Vector strain_vector = rValues.GetStrainVector();

        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
        const auto& r_combination_factors = this->GetCombinationFactors();

        rValue.resize(6, false);

        std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
        for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        // for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        //     rValue[i].resize(3, false);
        // }

        std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i].resize(6, false);
            for (IndexType j=0; j < 6; ++j) {
                interfacial_stress[i][j] = 0.0;
            }
        }

        for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {

            interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress
        }

        noalias(rValue) = interfacial_stress[interface_identifier];

        return rValue;
    } else if (rThisVariable == STRESS_VECTOR_COMP_1 ||
                rThisVariable == STRESS_VECTOR_COMP_2 ||
                rThisVariable == STRESS_VECTOR_COMP_3 ||
                rThisVariable == STRESS_VECTOR_COMP_4) {
        
        //

        // Get Values to compute the constitutive law:
        Flags& r_flags = rValues.GetOptions();

        // Previous flags saved
        const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
        const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

        const Properties& r_material_properties = rValues.GetMaterialProperties();

        // The deformation gradient
        if (rValues.IsSetDeterminantF()) {
            const double determinant_f = rValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }

        // All the strains must be the same, therefore we can just simply compute the strain in the first layer
        if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateGreenLagrangeStrain(rValues);
            r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
        }

        // The global strain vector, constant
        const Vector strain_vector = rValues.GetStrainVector();

        if (r_flags.Is(BaseType::COMPUTE_STRESS)) {
            // Set new flags
            r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
            r_flags.Set(BaseType::COMPUTE_STRESS, true);

            // Auxiliar stress vector
            const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();
            Vector auxiliar_stress_vector  = ZeroVector(VoigtSize);
            Vector delamination_damage_affected_stress_vector  = ZeroVector(VoigtSize);

            // The rotation matrix
            BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

            const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
            const auto& r_combination_factors = this->GetCombinationFactors();

            rValue.resize(6, false);

            std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
            for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
                layer_stress[i].resize(6, false);
            }

            std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
            for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
                interfacial_stress[i].resize(3, false);
            }

            std::vector<bool> negative_interfacial_stress_indicator(r_p_constitutive_law_vector.size()+1);
            for (IndexType i=0; i < r_p_constitutive_law_vector.size()+1; ++i) {
                negative_interfacial_stress_indicator[i] = false;
            }

            for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
                this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

                Properties& r_prop             = *(it_prop_begin + i_layer);
                ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

                // We rotate to local axes the strain
                noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

                rValues.SetMaterialProperties(r_prop);
                p_law->CalculateMaterialResponsePK2(rValues);

                // we return the stress and constitutive tensor to the global coordinates
                rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
                noalias(layer_stress[i_layer]) = rValues.GetStressVector();

                // we reset the properties and Strain
                rValues.SetMaterialProperties(r_material_properties);
                noalias(rValues.GetStrainVector()) = strain_vector;
            }

            const double tolerance = std::numeric_limits<double>::epsilon();
            Vector DelaminationDamageModeOne(r_p_constitutive_law_vector.size()+1);
            Vector DelaminationDamageModeTwo(r_p_constitutive_law_vector.size()+1);
            Vector ThresholdModeOne(r_p_constitutive_law_vector.size()-1);
            Vector ThresholdModeTwo(r_p_constitutive_law_vector.size()-1);

            noalias(DelaminationDamageModeOne) = mDelaminationDamageModeOne;
            noalias(DelaminationDamageModeTwo) = mDelaminationDamageModeTwo;
            noalias(ThresholdModeOne) = mThresholdModeOne;
            noalias(ThresholdModeTwo) = mThresholdModeTwo;

            for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {

                interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
                interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
                interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress

                double equivalent_stress_mode_one = interfacial_stress[i][0];
                double equivalent_stress_mode_two = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));

                if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < tolerance) {
                    negative_interfacial_stress_indicator[i+1] = true;
                }

                // Damage calculation

                const double T0n = r_material_properties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
                const double T0s = r_material_properties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
                const double GIc = r_material_properties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
                const double GIIc = r_material_properties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
                const double Ei = r_material_properties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? r_material_properties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface
                const double Gi = r_material_properties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? r_material_properties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface
                
                equivalent_stress_mode_one /= mFatigueDataContainersModeOne[i].GetFatigueReductionFactor();
                const double F_mode_one = equivalent_stress_mode_one - ThresholdModeOne[i];
                if (F_mode_one > tolerance) {

                    DelaminationDamageModeOne[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIc, Ei, T0n, equivalent_stress_mode_one);
                }

                equivalent_stress_mode_two /= mFatigueDataContainersModeTwo[i].GetFatigueReductionFactor();
                const double F_mode_two = equivalent_stress_mode_two - ThresholdModeTwo[i];
                if (F_mode_two > tolerance) {

                    DelaminationDamageModeTwo[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIIc, Gi, T0s, equivalent_stress_mode_two);
                }

                // End damage calculation
            }

            for(IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
                double layer_damage_variable_mode_one = 0.0;
                double layer_damage_variable_mode_two = 0.0;

                if (DelaminationDamageModeOne[i+1] > DelaminationDamageModeOne[i]) {
                    layer_damage_variable_mode_one = DelaminationDamageModeOne[i+1];
                } else {
                    layer_damage_variable_mode_one = DelaminationDamageModeOne[i];
                }

                if (DelaminationDamageModeTwo[i+1] > DelaminationDamageModeTwo[i]) {
                    layer_damage_variable_mode_two = DelaminationDamageModeTwo[i+1];
                } else {
                    layer_damage_variable_mode_two = DelaminationDamageModeTwo[i];
                }

                layer_stress[i][2] *= ((1.0-layer_damage_variable_mode_one));
                layer_stress[i][4] *= ((1.0-layer_damage_variable_mode_one) * (1.0-layer_damage_variable_mode_two));
                layer_stress[i][5] *= ((1.0-layer_damage_variable_mode_one) * (1.0-layer_damage_variable_mode_two));
            }

            // for(IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            //     const double factor = r_combination_factors[i];
            //     delamination_damage_affected_stress_vector += factor * layer_stress[i];
            // }

            auxiliar_stress_vector = delamination_damage_affected_stress_vector;

            if (rThisVariable == STRESS_VECTOR_COMP_1) {
                rValue = layer_stress[0];
            } else if (rThisVariable == STRESS_VECTOR_COMP_2) {
                rValue = layer_stress[1];
            } else if (rThisVariable == STRESS_VECTOR_COMP_3) {
                rValue = layer_stress[2];
            } else if (rThisVariable == STRESS_VECTOR_COMP_4) {
                rValue = layer_stress[3];
            }
            
            // noalias(rValues.GetStressVector()) = auxiliar_stress_vector;

            // if (flag_const_tensor) {
            //     this->CalculateTangentTensor(rValues, BaseType::StressMeasure_PK2);
            // }

            // Previous flags restored
            // r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
            // r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
            // r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        }

        //
    } else {
        return BaseType::CalculateValue(rValues, rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
bool TractionSeparationLaw3D<TDim>::ValidateInput(const Properties& rMaterialProperties)
{
    // We check it layer by layer
    bool valid_input = true;
    const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
    for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
        ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];
        Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
        if (p_law->ValidateInput(r_prop)) {
            valid_input = false;
            break;
        }
    }

    return valid_input;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::InitializeMaterial(
    const Properties& rMaterialProperties,
    const ConstitutiveLaw::GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();

    BaseType::InitializeMaterial(rMaterialProperties,rElementGeometry,rShapeFunctionsValues);

    mDelaminationDamageModeOne.resize(r_p_constitutive_law_vector.size()+1, false);
    noalias(mDelaminationDamageModeOne) = ZeroVector(r_p_constitutive_law_vector.size()+1);

    mDelaminationDamageModeTwo.resize(r_p_constitutive_law_vector.size()+1, false);
    noalias(mDelaminationDamageModeTwo) = ZeroVector(r_p_constitutive_law_vector.size()+1);

    mThresholdModeOne.resize(r_p_constitutive_law_vector.size()-1, false);
    for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        mThresholdModeOne[i] = rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH];
    }

    mThresholdModeTwo.resize(r_p_constitutive_law_vector.size()-1, false);
    for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        mThresholdModeTwo[i] = rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH];
    }

    mFatigueDataContainersModeOne.resize(r_p_constitutive_law_vector.size()-1);
    mFatigueDataContainersModeTwo.resize(r_p_constitutive_law_vector.size()-1);

    for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
        mFatigueDataContainersModeOne[i] = HCFDataContainer();
        mFatigueDataContainersModeTwo[i] = HCFDataContainer();
    }

    // Calculating A, B, and D matrices

    double Thickness = 0.015;
    double NumberOfLayers = r_p_constitutive_law_vector.size();

    std::vector<Vector> A(3);
    std::vector<Vector> A1(3);
    std::vector<Vector> A2(3);
    std::vector<Vector> B(3);
    std::vector<Vector> B1(3);
    std::vector<Vector> B2(3);
    std::vector<Vector> D(3);
    std::vector<Vector> D1(3);
    std::vector<Vector> D2(3);
    for (IndexType i=0; i < 3; ++i) {
        A[i].resize(3, false);
        A1[i].resize(3, false);
        A2[i].resize(3, false);
        B[i].resize(3, false);
        B1[i].resize(3, false);
        B2[i].resize(3, false);
        D[i].resize(3, false);
        D1[i].resize(3, false);
        D2[i].resize(3, false);
    }
    std::vector<Vector> a(NumberOfLayers);
    for (IndexType i=0; i < NumberOfLayers; ++i) {
        a[i].resize(3);
    }

    std::vector<Vector> StiffnessMatrix(6);
    std::vector<Vector> ComplianceMatrix(6);
    for (IndexType i=0; i < 6; ++i) {
        StiffnessMatrix[i].resize(6, false);
        ComplianceMatrix[i].resize(6, false);
    }

    double a11 = 0.0;
    double a22 = 0.0;
    double a33 = 0.0;


    // Calculating the [ABD] matrix of the original intact composite

    double LaminaThickness = Thickness / NumberOfLayers;
    Vector L(NumberOfLayers + 1);
    Vector T(NumberOfLayers);
    for(IndexType i=0; i < NumberOfLayers + 1; ++i) {
        L[i] = (Thickness / 2.0) - (LaminaThickness * i);
    }

    for(IndexType i=0; i < NumberOfLayers; ++i) {
        T[i] = rMaterialProperties[LAYER_EULER_ANGLES][3*i];
    }

    CalculateABDMatrices(NumberOfLayers, L, T, 0.0, A, B, D);
    KRATOS_WATCH(A);
    KRATOS_WATCH(B);
    KRATOS_WATCH(D);

    CalculateCompositeStiffnessAndComplianceMatrices(A, B, D, StiffnessMatrix, ComplianceMatrix, a11, a22, a33);
    KRATOS_WATCH(StiffnessMatrix);
    KRATOS_WATCH(ComplianceMatrix);

    a[0][0] = a11;
    a[0][1] = a22;
    a[0][2] = a33;

    // Calculating the [ABD] matrix of the original intact composite



    // Calculating the [ABD] matrix in delaminated case scenarios

    for(double i=0; i < NumberOfLayers - 1; ++i) {
        double n1 = i+1;
        double n2 = NumberOfLayers - (i+1);
        Vector L1(n1+1);
        Vector L2(n2+1);
        double S1 = 0.5 * (Thickness - ((1 + i) * LaminaThickness));
        double S2 = -0.5 * (1 + i) * LaminaThickness;
        Vector T1(n1);
        Vector T2(n2);
        for(IndexType j=0; j < n1; ++j) {
            T1[j] = T[j]; 
        }
        for(IndexType k=n1; k < (n1 + n2); ++k) {
            T2[k-n1] = T[k];
        }
        for(IndexType m=0; m < n1+1; ++m) {
            L1[m] = L[m] - S1;
        }
        for(IndexType h=n1; h < NumberOfLayers + 1; ++h) {
            L2[h-n1] = L[h] - S2;
        }
        
        CalculateABDMatrices(n1, L1, T1, S1, A1, B1, D1);
        CalculateABDMatrices(n2, L2, T2, S2, A2, B2, D2);

        for(IndexType k=0; k < 3; ++k) {
            for(IndexType h=0; h < 3; ++h) {
                A[k][h] = A1[k][h] + A2[k][h];
                B[k][h] = B1[k][h] + B2[k][h];
                D[k][h] = D1[k][h] + D2[k][h];
            }
        }

        CalculateCompositeStiffnessAndComplianceMatrices(A, B, D, StiffnessMatrix, ComplianceMatrix, a11, a22, a33);

        a[i+1][0] = a11;
        a[i+1][1] = a22;
        a[i+1][2] = a33;
    }

    std::vector<Vector> Dmax(NumberOfLayers + 1);
    for (IndexType i=0; i < NumberOfLayers + 1; ++i) {
        Dmax[i].resize(3);
    }

    for(IndexType i=0; i < NumberOfLayers - 1; ++i) {
        for(IndexType j=0; j < 3; ++j) {
            Dmax[i+1][j] = a[0][j] / a[i+1][j];
        } 
    }

    KRATOS_WATCH(a);
    KRATOS_WATCH(Dmax);

    // Calculating the [ABD] matrix in delaminated case scenarios

    // Vector L(4); 
    // L[0] = 0.0075; 
    // L[1] = 0.0025; 
    // L[2] = -0.0025; 
    // L[3] = -0.0075; 

    // Vector T(3);
    // T[0] = -45;
    // T[1] = 30;
    // T[2] = 0;

    // std::vector<Vector> A(3);
    // for (IndexType i=0; i < 3; ++i) {
    //     A[i].resize(3, false);
    // }

    // std::vector<Vector> B(3);
    // for (IndexType i=0; i < 3; ++i) {
    //     B[i].resize(3, false);
    // }

    // std::vector<Vector> D(3);
    // for (IndexType i=0; i < 3; ++i) {
    //     D[i].resize(3, false);
    // }

    // std::vector<Vector> K(6);
    // for (IndexType i=0; i < 6; ++i) {
    //     K[i].resize(6, false);
    // }

    // std::vector<Vector> ComplianceMatrix(6);
    // for (IndexType i=0; i < 6; ++i) {
    //     ComplianceMatrix[i].resize(6, false);
    // }

    // double S=0;


    // CalculateABDMatrices(3, L, T, S, A, B, D, K, ComplianceMatrix);
    // KRATOS_WATCH(A);
    // KRATOS_WATCH(B);
    // KRATOS_WATCH(D);
    // KRATOS_WATCH(StiffnessMatrix);
    // KRATOS_WATCH(ComplianceMatrix);
    // KRATOS_WATCH(a11);
    // KRATOS_WATCH(a22);
    // KRATOS_WATCH(a33);
   

    // Calculating A, B, and D matrices
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  TractionSeparationLaw3D<TDim>::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  TractionSeparationLaw3D<TDim>::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is(BaseType::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();
        Vector auxiliar_stress_vector  = ZeroVector(VoigtSize);
        Vector delamination_damage_affected_stress_vector  = ZeroVector(VoigtSize);

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();
        const auto& r_combination_factors = this->GetCombinationFactors();

        std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
        for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i].resize(3, false);
        }

        std::vector<bool> negative_interfacial_stress_indicator(r_p_constitutive_law_vector.size()+1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()+1; ++i) {
            negative_interfacial_stress_indicator[i] = false;
        }

        for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        const double tolerance = std::numeric_limits<double>::epsilon();
        Vector DelaminationDamageModeOne(r_p_constitutive_law_vector.size()+1);
        Vector DelaminationDamageModeTwo(r_p_constitutive_law_vector.size()+1);
        Vector ThresholdModeOne(r_p_constitutive_law_vector.size()-1);
        Vector ThresholdModeTwo(r_p_constitutive_law_vector.size()-1);

        noalias(DelaminationDamageModeOne) = mDelaminationDamageModeOne;
        noalias(DelaminationDamageModeTwo) = mDelaminationDamageModeTwo;
        noalias(ThresholdModeOne) = mThresholdModeOne;
        noalias(ThresholdModeTwo) = mThresholdModeTwo;

        for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {

            interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress

            double equivalent_stress_mode_one = interfacial_stress[i][0];
            double equivalent_stress_mode_two = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));

            if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < tolerance) {
                negative_interfacial_stress_indicator[i+1] = true;
            }

            // Damage calculation

            const double T0n = r_material_properties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
            const double T0s = r_material_properties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double GIc = r_material_properties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
            const double GIIc = r_material_properties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
            const double Ei = r_material_properties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? r_material_properties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface
            const double Gi = r_material_properties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? r_material_properties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface

            equivalent_stress_mode_one /= mFatigueDataContainersModeOne[i].GetFatigueReductionFactor();
            const double F_mode_one = equivalent_stress_mode_one - ThresholdModeOne[i];
            if (F_mode_one > tolerance) {

                DelaminationDamageModeOne[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIc, Ei, T0n, equivalent_stress_mode_one);
            }

            equivalent_stress_mode_two /= mFatigueDataContainersModeTwo[i].GetFatigueReductionFactor();
            const double F_mode_two = equivalent_stress_mode_two - ThresholdModeTwo[i];
            if (F_mode_two > tolerance) {

                DelaminationDamageModeTwo[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIIc, Gi, T0s, equivalent_stress_mode_two);
            }

            // End damage calculation
        }

        for(IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            double layer_damage_variable_mode_one = 0.0;
            double layer_damage_variable_mode_two = 0.0;

            if (DelaminationDamageModeOne[i+1] > DelaminationDamageModeOne[i]) {
                layer_damage_variable_mode_one = DelaminationDamageModeOne[i+1];
            } else {
                layer_damage_variable_mode_one = DelaminationDamageModeOne[i];
            }

            if (DelaminationDamageModeTwo[i+1] > DelaminationDamageModeTwo[i]) {
                layer_damage_variable_mode_two = DelaminationDamageModeTwo[i+1];
            } else {
                layer_damage_variable_mode_two = DelaminationDamageModeTwo[i];
            }

            layer_stress[i][2] *= ((1.0-layer_damage_variable_mode_one));
            layer_stress[i][4] *= ((1.0-layer_damage_variable_mode_one) * (1.0-layer_damage_variable_mode_two));
            layer_stress[i][5] *= ((1.0-layer_damage_variable_mode_one) * (1.0-layer_damage_variable_mode_two));
        }

        for(IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            const double factor = r_combination_factors[i];
            delamination_damage_affected_stress_vector += factor * layer_stress[i];
        }

        auxiliar_stress_vector = delamination_damage_affected_stress_vector;

        noalias(rValues.GetStressVector()) = auxiliar_stress_vector;

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, BaseType::StressMeasure_PK2);
        }

        // Previous flags restored
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void  TractionSeparationLaw3D<TDim>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain       = r_flags.Is(BaseType::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(BaseType::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress       = r_flags.Is(BaseType::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }

    // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    if (r_flags.IsNot(BaseType::USE_ELEMENT_PROVIDED_STRAIN)) {
        this->CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    // The global strain vector, constant
    const Vector strain_vector = rValues.GetStrainVector();

    if (r_flags.Is(BaseType::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(BaseType::COMPUTE_STRESS, true);

        // Auxiliar stress vector
        const auto it_prop_begin       = r_material_properties.GetSubProperties().begin();

        // The rotation matrix
        BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;

        const auto& r_p_constitutive_law_vector = this->GetConstitutiveLaws();

        std::vector<Vector> layer_stress(r_p_constitutive_law_vector.size());
        for (IndexType i=0; i < r_p_constitutive_law_vector.size(); ++i) {
            layer_stress[i].resize(6, false);
        }

        std::vector<Vector> interfacial_stress(r_p_constitutive_law_vector.size()-1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i].resize(3, false);
        }

        std::vector<bool> negative_interfacial_stress_indicator(r_p_constitutive_law_vector.size()+1);
        for (IndexType i=0; i < r_p_constitutive_law_vector.size()+1; ++i) {
            negative_interfacial_stress_indicator[i] = false;
        }

        for (IndexType i_layer = 0; i_layer < r_p_constitutive_law_vector.size(); ++i_layer) {
            this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);

            Properties& r_prop             = *(it_prop_begin + i_layer);
            ConstitutiveLaw::Pointer p_law = r_p_constitutive_law_vector[i_layer];

            // We rotate to local axes the strain
            noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);

            rValues.SetMaterialProperties(r_prop);
            p_law->CalculateMaterialResponsePK2(rValues);

            // we return the stress and constitutive tensor to the global coordinates
            rValues.GetStressVector()        = prod(trans(voigt_rotation_matrix), rValues.GetStressVector());
            noalias(layer_stress[i_layer]) = rValues.GetStressVector();

            p_law->FinalizeMaterialResponsePK2(rValues);

            // we reset the properties and Strain
            rValues.SetMaterialProperties(r_material_properties);
            noalias(rValues.GetStrainVector()) = strain_vector;
        }

        const double tolerance = std::numeric_limits<double>::epsilon();
        Vector DelaminationDamageModeOne(r_p_constitutive_law_vector.size()+1);
        Vector DelaminationDamageModeTwo(r_p_constitutive_law_vector.size()+1);
        Vector ThresholdModeOne(r_p_constitutive_law_vector.size()-1);
        Vector ThresholdModeTwo(r_p_constitutive_law_vector.size()-1);
        // vector <HCFDataContainer::FatigueVariables> HCFVariablesModeOne(r_p_constitutive_law_vector.size()-1);
        // vector <HCFDataContainer::FatigueVariables> HCFVariablesModeTwo(r_p_constitutive_law_vector.size()-1);
        // HCFDataContainer::FatigueVariables HCFVariablesModeOne = HCFDataContainer::FatigueVariables();
        // HCFDataContainer::FatigueVariables HCFVariablesModeTwo = HCFDataContainer::FatigueVariables();

        noalias(DelaminationDamageModeOne) = mDelaminationDamageModeOne;
        noalias(DelaminationDamageModeTwo) = mDelaminationDamageModeTwo;
        noalias(ThresholdModeOne) = mThresholdModeOne;
        noalias(ThresholdModeTwo) = mThresholdModeTwo;

        for(IndexType i=0; i < r_p_constitutive_law_vector.size()-1; ++i) {
            interfacial_stress[i][0] = AdvancedConstitutiveLawUtilities<VoigtSize>::MacaullyBrackets((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5); // interfacial normal stress
            interfacial_stress[i][1] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5; // interfacial shear stress
            interfacial_stress[i][2] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5; // interfacial shear stress

            double equivalent_stress_mode_one = interfacial_stress[i][0];
            double equivalent_stress_mode_one_fred = std::abs((layer_stress[i][2] + layer_stress[i+1][2]) * 0.5);
            double equivalent_stress_mode_two = std::sqrt(std::pow(interfacial_stress[i][1],2.0)+std::pow(interfacial_stress[i][2],2.0));

            // Fatigue calculations
            HCFDataContainer::FatigueVariables HCFVariablesModeOne = HCFDataContainer::FatigueVariables();
            HCFDataContainer::FatigueVariables HCFVariablesModeTwo = HCFDataContainer::FatigueVariables();

            mFatigueDataContainersModeOne[i].InitializeFatigueVariables(HCFVariablesModeOne, rValues);
            mFatigueDataContainersModeTwo[i].InitializeFatigueVariables(HCFVariablesModeTwo, rValues);
            HCFVariablesModeOne.hcf_coefficients = r_material_properties[HIGH_CYCLE_FATIGUE_COEFFICIENTS_MODE_ONE];
            HCFVariablesModeTwo.hcf_coefficients = r_material_properties[HIGH_CYCLE_FATIGUE_COEFFICIENTS_MODE_TWO];

            Vector fatigue_interfacial_stress_vector_mode_one = ZeroVector(VoigtSize);
            Vector fatigue_interfacial_stress_vector_mode_two = ZeroVector(VoigtSize);

            fatigue_interfacial_stress_vector_mode_one[2] = (layer_stress[i][2] + layer_stress[i+1][2]) * 0.5;
            fatigue_interfacial_stress_vector_mode_two[4] = (layer_stress[i][4] + layer_stress[i+1][4]) * 0.5;
            fatigue_interfacial_stress_vector_mode_two[5] = (layer_stress[i][5] + layer_stress[i+1][5]) * 0.5;

            mFatigueDataContainersModeOne[i].FinalizeSolutionStep(HCFVariablesModeOne,
                                            rValues.GetMaterialProperties(),
                                            rValues.GetProcessInfo(),
                                            fatigue_interfacial_stress_vector_mode_one,
                                            equivalent_stress_mode_one_fred,
                                            DelaminationDamageModeOne[i+1],
                                            ThresholdModeOne[i],
                                            INTERFACIAL_NORMAL_STRENGTH);

            mFatigueDataContainersModeTwo[i].FinalizeSolutionStep(HCFVariablesModeTwo,
                                            rValues.GetMaterialProperties(),
                                            rValues.GetProcessInfo(),
                                            fatigue_interfacial_stress_vector_mode_two,
                                            equivalent_stress_mode_two,
                                            DelaminationDamageModeTwo[i+1],
                                            ThresholdModeTwo[i],
                                            INTERFACIAL_SHEAR_STRENGTH,
                                            USER_DEFINED_METHOD);
            // End fatigue calculations

            // KRATOS_WATCH("------------------------------------------------------");
            // KRATOS_WATCH(mFatigueDataContainersModeTwo[0].GetMaxStressRelativeError());
            // KRATOS_WATCH(mFatigueDataContainersModeTwo[0].GetReversionFactorRelativeError());
            // KRATOS_WATCH(mFatigueDataContainersModeTwo[0].GetMaximumStress());
            // KRATOS_WATCH(mFatigueDataContainersModeTwo[0].GetPreviousMaximumStress());
            // KRATOS_WATCH("------------------------------------------------------");

            if ((layer_stress[i][2] + layer_stress[i+1][2] * 0.5) < tolerance) {
                negative_interfacial_stress_indicator[i+1] = true;
            }

            // Damage calculation
            const double T0n = r_material_properties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
            const double T0s = r_material_properties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? r_material_properties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : r_material_properties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
            const double GIc = r_material_properties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
            const double GIIc = r_material_properties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? r_material_properties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : r_material_properties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
            const double Ei = r_material_properties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? r_material_properties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface
            const double Gi = r_material_properties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? r_material_properties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : r_material_properties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface

            equivalent_stress_mode_one /= mFatigueDataContainersModeOne[i].GetFatigueReductionFactor();
            const double F_mode_one = equivalent_stress_mode_one - ThresholdModeOne[i];
            if (F_mode_one > tolerance) {

                DelaminationDamageModeOne[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIc, Ei, T0n, equivalent_stress_mode_one);

                mDelaminationDamageModeOne[i+1] = DelaminationDamageModeOne[i+1];
                mThresholdModeOne[i] = equivalent_stress_mode_one;
            }

            equivalent_stress_mode_two /= mFatigueDataContainersModeTwo[i].GetFatigueReductionFactor();
            const double F_mode_two = equivalent_stress_mode_two - ThresholdModeTwo[i];
            if (F_mode_two > tolerance) {

                DelaminationDamageModeTwo[i+1] = CalculateDelaminationDamageExponentialSoftening(rValues, GIIc, Gi, T0s, equivalent_stress_mode_two);

                mDelaminationDamageModeTwo[i+1] = DelaminationDamageModeTwo[i+1];
                mThresholdModeTwo[i] = equivalent_stress_mode_two;
            }

            mFatigueDataContainersModeOne[i].UpdateFatigueVariables(HCFVariablesModeOne);
            mFatigueDataContainersModeTwo[i].UpdateFatigueVariables(HCFVariablesModeTwo);

            // KRATOS_WATCH(HCFVariablesModeOne.LocalNumberOfCycles);
            // KRATOS_WATCH(HCFVariablesModeOne.Sth);
            // KRATOS_WATCH(HCFVariablesModeOne.MaxStress);
            // KRATOS_WATCH(HCFVariablesModeOne.MinStress);
            // KRATOS_WATCH(HCFVariablesModeOne.WohlerStress);

            // End damage calculation
        }

        // Previous flags restored
        r_flags.Set(BaseType::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(BaseType::COMPUTE_STRESS, flag_stress);
        r_flags.Set(BaseType::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }

    //Specifying the state of loading in the fatigue regime to be fed into the get/set methods. true->mode two / false->mode one
    double stress_mode_one = (1.0 - mFatigueDataContainersModeOne[0].GetReferenceDamage()) * mFatigueDataContainersModeOne[0].GetMaximumStress();
    double stress_mode_two = (1.0 - mFatigueDataContainersModeTwo[0].GetReferenceDamage()) * mFatigueDataContainersModeTwo[0].GetMaximumStress();
    double Sth_mode_one = mFatigueDataContainersModeOne[0].GetThresholdStress();
    double Sth_mode_two = mFatigueDataContainersModeTwo[0].GetThresholdStress();
    double Su_mode_one = r_material_properties[INTERFACIAL_NORMAL_STRENGTH];
    double Su_mode_two = r_material_properties[INTERFACIAL_SHEAR_STRENGTH];

    if ((stress_mode_one > Su_mode_one && stress_mode_two > Su_mode_two) || (stress_mode_one < Sth_mode_one && stress_mode_two < Sth_mode_two)) {
        mFatigueLoadingStateParameter = true; //It is optional and could be mode one (false) or mode two (true)
    } else if (stress_mode_one > Sth_mode_one && stress_mode_one < Su_mode_one && stress_mode_two > Sth_mode_two && stress_mode_two < Su_mode_two) {
    
        double N_residual_mode_one = mFatigueDataContainersModeOne[0].GetCyclesToFailure() - mFatigueDataContainersModeOne[0].GetLocalNumberOfCycles();
        double N_residual_mode_two = mFatigueDataContainersModeTwo[0].GetCyclesToFailure() - mFatigueDataContainersModeTwo[0].GetLocalNumberOfCycles();

        if (N_residual_mode_one > N_residual_mode_two) {
            mFatigueLoadingStateParameter = true;
        } else {
            mFatigueLoadingStateParameter = false;
        }
    } else if (stress_mode_one < Sth_mode_one && stress_mode_two > Sth_mode_two) {
        mFatigueLoadingStateParameter = true;
    } else if (stress_mode_two < Sth_mode_two && stress_mode_one > Sth_mode_one) {
        mFatigueLoadingStateParameter = false;
    } else if (stress_mode_one > Sth_mode_one && stress_mode_one < Su_mode_one && stress_mode_two > Su_mode_two) {
        mFatigueLoadingStateParameter = true;
    } else if (stress_mode_two > Sth_mode_two && stress_mode_two < Su_mode_two && stress_mode_one > Su_mode_one) {
        mFatigueLoadingStateParameter = false;
    }
    //Specifying the state of loading in the fatigue regime to be fed into the get/set methods
    KRATOS_CATCH("");

}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
double TractionSeparationLaw3D<TDim>::CalculateDelaminationDamageExponentialSoftening(
    ConstitutiveLaw::Parameters& rValues,
    const double GI,
    const double E,
    const double T0,
    const double equivalent_stress)
{
    const double characteristic_length = 0.6343 * (AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry()));
    const double AParameter = 1.00 / (GI * E / (characteristic_length * std::pow(T0, 2)) - 0.5); // Exponential

    KRATOS_ERROR_IF(AParameter < 0.0) << "AParameter is negative." << std::endl;

    double DelaminationDamage = 1.0 - (T0 / equivalent_stress) * std::exp(AParameter * (1.0 - equivalent_stress / T0)); // Exponential

    DelaminationDamage = (DelaminationDamage >= 0.99999) ? 0.99999 : DelaminationDamage;
    DelaminationDamage = (DelaminationDamage < 0.0) ? 0.0 : DelaminationDamage;
    return DelaminationDamage;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
std::vector<Vector>& TractionSeparationLaw3D<TDim>::CalculateABDMatrices(
    double n, 
    Vector L, 
    Vector T,
    double S,
    std::vector<Vector>& A,
    std::vector<Vector>& B,
    std::vector<Vector>& D)
{
    for(IndexType i=0; i < T.size(); ++i) {
        T[i] *= (Globals::Pi / 180.0);
    }
    

    double E1=181e9;
    double E2=10.3e9;
    double G12=7.17e9;
    double v12=0.28;
    double v21=(E2/E1)*v12;
    double Q11=E1/(1-(v12)*(v21));
    double Q22=E2/(1-(v12)*(21));
    double Q12=v12*E2/(1-(v12)*(v21));
    double Q66=G12;
    double U1=(3*Q11+3*Q22+2*Q12+4*Q66)/8;
    double U2=(Q11-Q22)/2;
    double U3=(Q11+Q22-2*Q12-4*Q66)/8;
    double U4=(Q11+Q22+6*Q12-4*Q66)/8;

    std::vector<std::vector<std::vector<double>>> QB(3, std::vector<std::vector<double>>(3, std::vector<double>(n)));
    
    for(IndexType k=0; k < n; ++k) {
        QB[0][0][k]=U1+U2*cos(2*T[k])+U3*cos(4*T[k]);
        QB[0][1][k]=U4-U3*cos(4*T[k]);
        QB[1][0][k]=QB[0][1][k];
        QB[1][1][k]=U1-U2*cos(2*T[k])+U3*cos(4*T[k]);
        QB[2][0][k]=0.5*U2*sin(2*T[k])+U3*sin(4*T[k]);
        QB[0][2][k]=0.5*U2*sin(2*T[k])+U3*sin(4*T[k]);
        QB[2][1][k]=0.5*U2*sin(2*T[k])-U3*sin(4*T[k]);
        QB[1][2][k]=0.5*U2*sin(2*T[k])-U3*sin(4*T[k]);
        QB[2][2][k]=0.5*(U1-U4)-U3*cos(4*T[k]);
    }

        // QB[0][0][k]=U1+U2*cos(2*T[k])+U3*cos(4*T[k]);
        // QB[0][1][k]=U4-U3*cos(4*T[k]);
        // QB[1][0][k]=QB[0][1][k];
        // QB[1][1][k]=U1-U2*cos(2*T[k])+U3*cos(4*T[k]);
        // QB[2][0][k]=0.5*U2*sin(2*T[k])+U3*sin(4*T[k]);
        // QB[0][2][k]=2*QB[2][0][k];
        // QB[2][1][k]=0.5*U2*sin(2*T[k])-U3*sin(4*T[k]);
        // QB[1][2][k]=2*QB[2][1][k];
        // QB[2][2][k]=(U1-U4)-2*U3*cos(4*T[k]);

    std::vector<std::vector<std::vector<double>>> AS(3, std::vector<std::vector<double>>(3, std::vector<double>(n)));
    std::vector<std::vector<std::vector<double>>> BS(3, std::vector<std::vector<double>>(3, std::vector<double>(n)));
    std::vector<std::vector<std::vector<double>>> DS(3, std::vector<std::vector<double>>(3, std::vector<double>(n)));

    for(IndexType k=0; k < n; ++k) {
        for(IndexType i=0; i < 3; ++i) {
            for(IndexType j=0; j < 3; ++j) {
                AS[i][j][k]=QB[i][j][k]*(L[k]-L[k+1]);
                BS[i][j][k]=QB[i][j][k]*(((L[k])*(L[k]))-(L[k+1])*(L[k+1]));
                DS[i][j][k]=QB[i][j][k]*((L[k])*(L[k])*(L[k])-(L[k+1])*(L[k+1])*(L[k+1]));

                // AS[i][j][k]=QB[i][j][k]*(L[k+1]-L[k]);
                // BS[i][j][k]=QB[i][j][k]*(((L[k+1])*(L[k+1]))-(L[k])*(L[k]));
                // DS[i][j][k]=QB[i][j][k]*((L[k+1])*(L[k+1])*(L[k+1])-(L[k])*(L[k])*(L[k]));
            }
        }
    }

    std::vector<Vector> ASF(3);
    std::vector<Vector> BSF(3);
    std::vector<Vector> DSF(3);
    for (IndexType i=0; i < 3; ++i) {
        ASF[i].resize(3, false);
        BSF[i].resize(3, false);
        DSF[i].resize(3, false);
        A[i].resize(3, false);
        B[i].resize(3, false);
        D[i].resize(3, false);
    }

    for(IndexType i=0; i < 3; ++i) {
        for(IndexType j=0; j < 3; ++j) {
            ASF[i][j] = 0.0;
            BSF[i][j] = 0.0;
            DSF[i][j] = 0.0;
            A[i][j] = 0.0;
            B[i][j] = 0.0;
            D[i][j] = 0.0;
        }
    }

    for(IndexType k=0; k < n; ++k) {
        for(IndexType i=0; i < 3; ++i) {
            for(IndexType j=0; j < 3; ++j) {
                ASF[i][j] += AS[i][j][k];
                BSF[i][j] += (0.5) * (BS[i][j][k]);
                DSF[i][j] += (0.33333333) * (DS[i][j][k]);
            }
        }
    }

   // Constructing transformed [ABD] matrices to the global mid-plane
   for(IndexType i=0; i < 3; ++i) {
        for(IndexType j=0; j < 3; ++j) {
            A[i][j] = ASF[i][j];
            B[i][j] = BSF[i][j] + (S * ASF[i][j]);
            D[i][j] = DSF[i][j] + (2 * S * BSF[i][j]) + (S * S * ASF[i][j]);
        }
   }

//    // Constructing stifness matrix

//    for (IndexType i=0; i < 6; ++i) {
//         K[i].resize(6, false);
//    }

//    for (int i = 0; i < A.size(); ++i) {
//         for (int j = 0; j < A[i].size(); ++j) {
//             K[i][j] = A[i][j];
//             K[i][j + A[i].size()] = B[i][j];
//             K[i + A.size()][j] = B[i][j];
//             K[i + A.size()][j + A[i].size()] = D[i][j];
//         }
//     }
   
//    // Calculating the inverse matrix

//     int matrix_size = 6;
//     std::vector<Vector> augmentedMatrix(matrix_size);
//     for (IndexType i=0; i < matrix_size; ++i) {
//         augmentedMatrix[i].resize(2*matrix_size, false);
//     }

//     for(IndexType i=0; i < matrix_size; ++i) {
//         for(IndexType j=0; j < matrix_size; ++j) {
//             augmentedMatrix[i][j] = K[i][j];
//         }
//     }

//      for(IndexType i=0; i < matrix_size; ++i) {
//         for(IndexType j=6; j < 12; ++j) {
//             augmentedMatrix[i][j] = 0.0;
//         }
//     }

//     for (IndexType i = 0; i < matrix_size; ++i) {
//         augmentedMatrix[i][i + matrix_size] = 1.0;
//     }

//     for (IndexType i = 0; i < matrix_size; ++i) {
//         double diagElem = augmentedMatrix[i][i];
//         if (diagElem == 0.0) throw std::runtime_error("Matrix is singular and cannot be inverted.");

//         for (IndexType j = 0; j < 2 * matrix_size; ++j) {
//             augmentedMatrix[i][j] /= diagElem;
//         }

//         for (IndexType k = 0; k < matrix_size; ++k) {
//             if (k != i) {
//                 double factor = augmentedMatrix[k][i];
//                 for (IndexType j = 0; j < 2 * matrix_size; ++j) {
//                     augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
//                 }
//             }
//         }
//     }

//     for (IndexType i=0; i < matrix_size; ++i) {
//         ComplianceMatrix[i].resize(matrix_size, false);
//     }

//     for (IndexType i = 0; i < matrix_size; ++i) {
//         for (IndexType j = 0; j < matrix_size; ++j) {
//             ComplianceMatrix[i][j] = augmentedMatrix[i][j + matrix_size];
//         }
//     }

   //

    return A;
    return B;
    return D;
    // return K;
    // return ComplianceMatrix;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateCompositeStiffnessAndComplianceMatrices(
    std::vector<Vector>& A,
    std::vector<Vector>& B,
    std::vector<Vector>& D,
    std::vector<Vector>& StiffnessMatrix,
    std::vector<Vector>& ComplianceMatrix,
    double& a11,
    double& a22,
    double& a33)
{
   // Constructing stifness matrix

   for (IndexType i=0; i < 6; ++i) {
        StiffnessMatrix[i].resize(6, false);
   }

   for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A[i].size(); ++j) {
            StiffnessMatrix[i][j] = A[i][j];
            StiffnessMatrix[i][j + A[i].size()] = B[i][j];
            StiffnessMatrix[i + A.size()][j] = B[i][j];
            StiffnessMatrix[i + A.size()][j + A[i].size()] = D[i][j];
        }
    }
   
   // Calculating the inverse matrix

    int matrix_size = 6;
    std::vector<Vector> augmentedMatrix(matrix_size);
    for (IndexType i=0; i < matrix_size; ++i) {
        augmentedMatrix[i].resize(2*matrix_size, false);
    }

    for(IndexType i=0; i < matrix_size; ++i) {
        for(IndexType j=0; j < matrix_size; ++j) {
            augmentedMatrix[i][j] = StiffnessMatrix[i][j];
        }
    }

     for(IndexType i=0; i < matrix_size; ++i) {
        for(IndexType j=6; j < 12; ++j) {
            augmentedMatrix[i][j] = 0.0;
        }
    }

    for (IndexType i = 0; i < matrix_size; ++i) {
        augmentedMatrix[i][i + matrix_size] = 1.0;
    }

    for (IndexType i = 0; i < matrix_size; ++i) {
        double diagElem = augmentedMatrix[i][i];
        if (diagElem == 0.0) throw std::runtime_error("Matrix is singular and cannot be inverted.");

        for (IndexType j = 0; j < 2 * matrix_size; ++j) {
            augmentedMatrix[i][j] /= diagElem;
        }

        for (IndexType k = 0; k < matrix_size; ++k) {
            if (k != i) {
                double factor = augmentedMatrix[k][i];
                for (IndexType j = 0; j < 2 * matrix_size; ++j) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    for (IndexType i=0; i < matrix_size; ++i) {
        ComplianceMatrix[i].resize(matrix_size, false);
    }

    for (IndexType i = 0; i < matrix_size; ++i) {
        for (IndexType j = 0; j < matrix_size; ++j) {
            ComplianceMatrix[i][j] = augmentedMatrix[i][j + matrix_size];
        }
    }

    a11 = ComplianceMatrix[0][0];
    a22 = ComplianceMatrix[1][1];
    a33 = ComplianceMatrix[2][2];
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
int TractionSeparationLaw3D<TDim>::Check(
    const Properties& rMaterialProperties,
    const ConstitutiveLaw::GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    BaseType::Check(rMaterialProperties,rElementGeometry,rCurrentProcessInfo);

    // Check if input parameters are completely defined
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH) || rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR)) << "INTERFACIAL_NORMAL_STRENGTH is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH) || rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR)) << "INTERFACIAL_SHEAR_STRENGTH is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY) || rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR)) << "MODE_ONE_FRACTURE_ENERGY is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY) || rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR)) << "MODE_TWO_FRACTURE_ENERGY is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS) || rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS_VECTOR)) << "TENSILE_INTERFACE_MODULUS is not a defined value" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS) || rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS_VECTOR)) << "SHEAR_INTERFACE_MODULUS is not a defined value" << std::endl;

    // Check the size of the vectors
    if (rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "INTERFACIAL_NORMAL_STRENGTH_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "INTERFACIAL_SHEAR_STRENGTH_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "MODE_ONE_FRACTURE_ENERGY_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "MODE_TWO_FRACTURE_ENERGY_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[TENSILE_INTERFACE_MODULUS_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "TENSILE_INTERFACE_MODULUS_VECTOR badly defined" << std::endl;
    }

    if (rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS_VECTOR)) {
        KRATOS_ERROR_IF(rMaterialProperties[SHEAR_INTERFACE_MODULUS_VECTOR].size() != (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1) << "SHEAR_INTERFACE_MODULUS_VECTOR badly defined" << std::endl;
    }

    // Check negative fracture energy
    if (rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR)) {
        const double SizeModeOne = rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR].size();
        for (IndexType i=0; i < SizeModeOne; ++i) {
            KRATOS_ERROR_IF(rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] < 0.0) << "MODE_ONE_FRACTURE_ENERGY is negative at interface " << i << std::endl;
        }
    } else {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_ONE_FRACTURE_ENERGY] < 0.0) << "MODE_ONE_FRACTURE_ENERGY is negative." << std::endl;
    }

    if (rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR)) {
        const double SizeModeTwo = rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR].size();
        for (IndexType i=0; i < SizeModeTwo; ++i) {
            KRATOS_ERROR_IF(rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] < 0.0) << "MODE_TWO_FRACTURE_ENERGY is negative at interface " << i << std::endl;
        }
    } else {
        KRATOS_ERROR_IF(rMaterialProperties[MODE_TWO_FRACTURE_ENERGY] < 0.0) << "MODE_TWO_FRACTURE_ENERGY is negative." << std::endl;
    }

    // Check fracture energy
    const double characteristic_length = 0.6343 * (AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rElementGeometry));

    for(IndexType i=0; i < (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1; ++i) {
        const double Tn = rMaterialProperties.Has(INTERFACIAL_NORMAL_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_NORMAL_STRENGTH]; // Interfacial Normal Strength
        const double GI = rMaterialProperties.Has(MODE_ONE_FRACTURE_ENERGY_VECTOR) ? rMaterialProperties[MODE_ONE_FRACTURE_ENERGY_VECTOR][i] : rMaterialProperties[MODE_ONE_FRACTURE_ENERGY]; // Mode I Energy Release Rate
        const double E = rMaterialProperties.Has(TENSILE_INTERFACE_MODULUS_VECTOR) ? rMaterialProperties[TENSILE_INTERFACE_MODULUS_VECTOR][i] : rMaterialProperties[TENSILE_INTERFACE_MODULUS]; // Tensile modulus of the interface

        const double AParameter_mode_one = 1.00 / (GI * E / (characteristic_length * std::pow(Tn, 2)) - 0.5); // Exponential
        KRATOS_ERROR_IF(AParameter_mode_one < 0.0) << "MODE_ONE_FRACTURE_ENERGY is too low at interface " << i << std::endl;
    }

    for(IndexType i=0; i < (rMaterialProperties[LAYER_EULER_ANGLES].size() / 3) - 1; ++i) {
        const double Ts = rMaterialProperties.Has(INTERFACIAL_SHEAR_STRENGTH_VECTOR) ? rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH_VECTOR][i] : rMaterialProperties[INTERFACIAL_SHEAR_STRENGTH]; // Interfacial Shear Strength
        const double GII = rMaterialProperties.Has(MODE_TWO_FRACTURE_ENERGY_VECTOR) ? rMaterialProperties[MODE_TWO_FRACTURE_ENERGY_VECTOR][i] : rMaterialProperties[MODE_TWO_FRACTURE_ENERGY]; // Mode II Energy Release Rate
        const double G = rMaterialProperties.Has(SHEAR_INTERFACE_MODULUS_VECTOR) ? rMaterialProperties[SHEAR_INTERFACE_MODULUS_VECTOR][i] : rMaterialProperties[SHEAR_INTERFACE_MODULUS]; // Shear modulus of the interface

        const double AParameter_mode_two = 1.00 / (GII * G / (characteristic_length * std::pow(Ts, 2)) - 0.5); // Exponential
        KRATOS_ERROR_IF(AParameter_mode_two < 0.0) << "MODE_TWO_FRACTURE_ENERGY is too low at interface " << i << std::endl;
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

template<unsigned int TDim>
void TractionSeparationLaw3D<TDim>::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues, const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 4);
    } else {
        BaseType::CalculateTangentTensor(rValues,rStressMeasure);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class TractionSeparationLaw3D<3>;

} // Namespace Kratos
