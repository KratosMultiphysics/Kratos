//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquín Irazábal González
//

// Application includes
#include "custom_constitutive/dam_joint_bilinear_cohesive_2D_law.hpp"

namespace Kratos
{

    double& DamJointBilinearCohesive2DLaw::GetValue(const Variable<double>& rThisVariable,
                                                    double& rValue )
    {
        KRATOS_TRY

        if(rThisVariable == STATE_VARIABLE)
        {
            rValue = mStateVariable;
        }
        else if (rThisVariable == UPLIFT_PRESSURE)
        {
            rValue = mUpliftPressure;
        }

        return rValue;

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------

    void DamJointBilinearCohesive2DLaw::SetValue(const Variable<double>& rThisVariable,
                                                 const double& rValue,
                                                 const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        if (rThisVariable == STATE_VARIABLE)
        {
            mStateVariable = rValue;
        }
        else if (rThisVariable == UPLIFT_PRESSURE)
        {
            mUpliftPressure = rValue;
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------

    void DamJointBilinearCohesive2DLaw::ComputeStressVector(Vector& rStressVector,
                                                            ConstitutiveLawVariables& rVariables,
                                                            Parameters& rValues)
    {
        BilinearCohesive2DLaw::ComputeStressVector(rStressVector, rVariables, rValues);

        // Add Uplift Pressure
        rStressVector[1] -= mUpliftPressure;
    }

} // Namespace Kratos
