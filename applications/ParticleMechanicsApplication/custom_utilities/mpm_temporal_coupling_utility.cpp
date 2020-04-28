//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Phillip Salatoskratos
//

// System includes

// External includes

// Project includes

#include "custom_utilities/mpm_temporal_coupling_utility.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

namespace Kratos
{

    void MPMTemporalCouplingUtility::CalculateCorrectiveLagrangianMultipliers(
        const SystemMatrixType& K_1, 
        const SystemMatrixType& K_2 )
    {
        Vector SubDomain1InterpolatedVelocities;

        // Store inverted mass matrix and coupling matrix for sub-domain 1 for the whole timestep
        if (mJ == 0) {
            // invert M1 to get mInvM1
            mInvM1 = ZeroMatrix(0);

            // Setup coupling matrix for domain 1;
            mCoupling1 = ZeroMatrix(0);
        }

        // Interpolate subdomain 1 velocities
        if (mJ == 0)
        {
            SubDomain1InterpolatedVelocities = mSubDomain1InitialVelocity;
        }
        else if (mJ == mTimeStepRatio)
        {
            SubDomain1InterpolatedVelocities = mSubDomain1FinalVelocity;
        }
        else
        {
            SubDomain1InterpolatedVelocities = 
                (1.0 - mJ / mTimeStepRatio) * mSubDomain1InitialVelocity
                + mJ / mTimeStepRatio * mSubDomain1FinalVelocity;
        }

        // Invert sub domain 2 mass matrix
        Matrix InvM2 = ZeroMatrix(0);

        // Establish sub domain 2 coupling matrix
        Matrix Coupling2 = ZeroMatrix(0);

        // Get sub domain 2 free velocities
        Vector subDomain2freeVelocities;

        // Assemble condensation operator H
        Matrix H = -1.0 * mSmallTimestep * (1.0 - mGamma[0]) * mInvM1 
            - mSmallTimestep * (1.0 - mGamma[1]) * InvM2;

        // Assemble condensed unbalanced interface velocities
        Vector interface_free_vel_1 = SubDomain1InterpolatedVelocities + mSubDomain1AccumulatedLinkVelocity;
        Vector b = prod(mCoupling1, interface_free_vel_1) + prod(Coupling2, subDomain2freeVelocities);

        // Calculate corrective Lagrangian multipliers
        Vector lamda = ZeroVector(0); // inv(H)*b;

        // Calculate link velocities
        Matrix temp1 = prod(mInvM1, trans(mCoupling1));
        Vector link_velocity_1 = (mGamma[0] * mSmallTimestep * prod(temp1, lamda));
        mSubDomain1AccumulatedLinkVelocity += link_velocity_1;

        Matrix temp2 = prod(InvM2, trans(Coupling2));
        Vector link_accel_2 = prod(temp2, lamda);
        Vector link_velocity_2 = mGamma[1] * mSmallTimestep * link_accel_2;

        const IndexType working_space_dimension = 2; // TODO update

        if (mJ == mTimeStepRatio)
        {
            const double time_step_1 = mSmallTimestep * mTimeStepRatio;
            const Vector link_accel_1 = mSubDomain1AccumulatedLinkVelocity / mGamma[0] / time_step_1;
            if (mGamma[0] == 0.0) // explicit correction
            {
                for (IndexType i = 0; i < link_velocity_1.size(); ++i)
                {
                    double nodal_mass = 0.0; // TODO update
                    for (IndexType dim = 0; dim < working_space_dimension; ++dim)
                    {
                        // domain_1_mom[i] += nodal_mass*mSubDomain1AccumulatedLinkVelocity[i+0];
                        // domain_1_force[i] += nodal_mass*link_accel_1[i+0];
                    }
                }
            }
            else // implicit correction
            {
                for (IndexType i = 0; i < link_accel_1.size(); ++i)
                {
                    for (IndexType dim = 0; dim < working_space_dimension; ++dim)
                    {
                        // domain_1_disp[i] += 0.25*time_step_1*time_step_1*link_accel_1[i+0];
                        // domain_1_accel[i] += link_accel_1[i+0];
                    }
                }
            }
        }

        if (mGamma[1] == 0.0) // explicit correction
        {
            for (IndexType i = 0; i < link_velocity_2.size(); ++i)
            {
                double nodal_mass = 0.0; // TODO update
                for (IndexType dim = 0; dim < working_space_dimension; ++dim)
                {
                    // domain_2_mom[i] += nodal_mass*link_velocity_2[i+0];
                    // domain_2_force[i] += nodal_mass*link_accel_2[i+0];
                }
            }
        }
        else // implicit correction
        {
            for (IndexType i = 0; i < link_accel_2.size(); ++i)
            {
                for (IndexType dim = 0; dim < working_space_dimension; ++dim)
                {
                    // domain_1_disp[i] += 0.25*mSmallTimestep*mSmallTimestep*link_accel_2[i+0];
                    // domain_1_accel[i] += link_accel_2[i+0];
                }
            }
        }
    }

 // end namespace MPMTemporalCouplingUtility
} // end namespace Kratos



