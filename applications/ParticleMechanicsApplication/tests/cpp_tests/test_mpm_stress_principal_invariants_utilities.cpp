//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/math_utils.h"
#include "custom_utilities/mpm_stress_principal_invariants_utility.h"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    // Generate hydrostatic principal stress of value 3
    Vector CreateHydrostaticStressVector()
    {
        Vector vector = ZeroVector(3);
        std::fill(vector.begin(), vector.end(), 3.0e6);
        return vector;
    }

    // Generate hydrostatic principal stress of value -8
    Vector CreateNegativeHydrostaticStressVector()
    {
        Vector vector = ZeroVector(3);
        std::fill(vector.begin(), vector.end(), -8.0e6);
        return vector;
    }

    // Generate testing principal stress with random but sorted values
    Vector CreateRandomSortedStressVector()
    {
        Vector vector = ZeroVector(3);
        vector[0] = 10.0e6;
        vector[1] = 2.0e6;
        vector[2] = -5.0e6;
        return vector;
    }

    // Generate testing principal stress with random and unsorted values
    Vector CreateRandomUnsortedStressVector()
    {
        Vector vector = ZeroVector(3);
        vector[0] = -3.0e6;
        vector[1] = 12.0e6;
        vector[2] = -1.0e6;
        return vector;
    }

    // Generate testing stress vector with random and unsorted values
    Vector CreateRandomUnsortedLargeStressVector()
    {
        Vector stress_vector = ZeroVector(6);
        stress_vector[0] = -19.0e6;
        stress_vector[1] = 4.6e6;
        stress_vector[2] = -8.3e6;
        stress_vector[3] = -4.7e6;
        stress_vector[4] = 11.8e6;
        stress_vector[5] = 6.45e6;
        return stress_vector;
    }

    // Generate dummy principal strain with random values
    Vector CreateRandomStrainVector()
    {
        Vector strain_dummy = ZeroVector(3);
        strain_dummy[0] = 1.0e-7;
        strain_dummy[1] =-2.0e-7;
        strain_dummy[2] = 3.0e-7;
        return strain_dummy;
    }

    // Generate dummy direction matrix with random values
    Matrix CreateRandomDirectionMatrix()
    {
        Matrix direction_dummy = ZeroMatrix(3);
        direction_dummy(0,0) = 1.0;
        direction_dummy(0,1) = 2.0;
        direction_dummy(0,2) = 3.0;

        direction_dummy(1,0) = 4.0;
        direction_dummy(1,1) = 5.0;
        direction_dummy(1,2) = 6.0;

        direction_dummy(2,0) = 7.0;
        direction_dummy(2,1) = 8.0;
        direction_dummy(2,2) = 9.0;

        return direction_dummy;
    }

    /**
    * Check whether the stress sorting is working as expected
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityStressSorting, KratosParticleMechanicsFastSuite)
    {
        // Stress dummy -- unsorted
        auto random_unsorted = CreateRandomUnsortedStressVector();

        // Strain dummy
        auto dummy_strain = CreateRandomStrainVector();

        // Direction dummy
        auto dummy_direction = CreateRandomDirectionMatrix();

        MPMStressPrincipalInvariantsUtility::SortPrincipalStress(random_unsorted, dummy_strain, dummy_direction);

        // KRATOS_CHECK_EQUAL(geomRegLen1->FacesNumber(), 6);
        // KRATOS_CHECK_EQUAL(geomRegLen2->FacesNumber(), 6);

    }

    /**
    * Check whether the tensor invariant calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityTensorInvariantsCalculation, KratosParticleMechanicsFastSuite)
    {
        auto positive_hydrostatic  = CreateHydrostaticStressVector();
        auto random_sorted         = CreateRandomSortedStressVector();
        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();

        // KRATOS_CHECK_LESS_EQUAL((0.0412844 - strain_vector[0])/strain_vector[0], tolerance);
        // KRATOS_CHECK_LESS_EQUAL((-0.00458716 - strain_vector[1])/strain_vector[1], tolerance);
        // KRATOS_CHECK_LESS_EQUAL((0.0 - strain_vector[2])/(strain_vector[2] + 1.0e-12), tolerance);
        // KRATOS_CHECK_LESS_EQUAL((-0.0458716 - strain_vector[3])/strain_vector[3], tolerance);
    }

    /**
    * Check whether the tensor invariant derivatives calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityTensorInvariantDerivativesCalculation, KratosParticleMechanicsFastSuite)
    {
        auto negative_hydrostatic  = CreateNegativeHydrostaticStressVector();
        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
    }

    /**
    * Check whether the tensor second invariant derivatives calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityTensorInvariantSecondDerivativesCalculation, KratosParticleMechanicsFastSuite)
    {
        auto random_sorted         = CreateRandomSortedStressVector();
        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
    }

    /**
    * Check whether the stress invariant calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityStressInvariantsCalculation, KratosParticleMechanicsFastSuite)
    {
        auto negative_hydrostatic  = CreateNegativeHydrostaticStressVector();
        auto random_unsorted       = CreateRandomUnsortedStressVector();
        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
    }

    /**
    * Check whether the stress invariant derivatives calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityStressInvariantDerivativesCalculation, KratosParticleMechanicsFastSuite)
    {
        auto random_unsorted       = CreateRandomUnsortedStressVector();
        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
    }

    /**
    * Check whether the stress second invariant derivatives calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityStressInvariantSecondDerivativesCalculation, KratosParticleMechanicsFastSuite)
    {
        auto positive_hydrostatic  = CreateHydrostaticStressVector();
        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
    }




} // namespace Testing
} // namespace Kratos
