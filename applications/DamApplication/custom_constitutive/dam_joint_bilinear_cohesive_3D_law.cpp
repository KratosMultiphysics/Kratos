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
#include "custom_constitutive/dam_joint_bilinear_cohesive_3D_law.hpp"

namespace Kratos
{

    double& DamJointBilinearCohesive3DLaw::GetValue(const Variable<double>& rThisVariable,
                                                    double& rValue)
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

    void DamJointBilinearCohesive3DLaw::SetValue(const Variable<double>& rThisVariable,
                                                 const double& rValue,
                                                 const ProcessInfo& rCurrentProcessInfo)
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

    void DamJointBilinearCohesive3DLaw::ComputeStressVector(Vector& rStressVector,
                                                            ConstitutiveLawVariables& rVariables,
                                                            Parameters& rValues)
    {
        BilinearCohesive3DLaw::ComputeStressVector(rStressVector, rVariables, rValues);

        // Add Uplift Pressure
        rStressVector[2] -= mUpliftPressure;
    }

} // Namespace Kratos
