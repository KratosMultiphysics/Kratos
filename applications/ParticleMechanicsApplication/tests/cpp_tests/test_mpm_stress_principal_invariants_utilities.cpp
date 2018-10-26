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
        stress_vector[0] = -19.00e6;
        stress_vector[1] =   4.60e6;
        stress_vector[2] =  -8.30e6;
        stress_vector[3] =  -4.70e6;
        stress_vector[4] =  11.80e6;
        stress_vector[5] =   6.45e6;
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

        // The function here just performed sorting, thus results need to be exact
        KRATOS_CHECK_EQUAL(random_unsorted[0], 12.0e6);
        KRATOS_CHECK_EQUAL(random_unsorted[1], -1.0e6);
        KRATOS_CHECK_EQUAL(random_unsorted[2], -3.0e6);

        KRATOS_CHECK_EQUAL(dummy_strain[0], -2.0e-7);
        KRATOS_CHECK_EQUAL(dummy_strain[1],  3.0e-7);
        KRATOS_CHECK_EQUAL(dummy_strain[2],  1.0e-7);

        KRATOS_CHECK_EQUAL(dummy_direction(0,0), 4.0);
        KRATOS_CHECK_EQUAL(dummy_direction(0,1), 7.0);
        KRATOS_CHECK_EQUAL(dummy_direction(0,2), 1.0);
        KRATOS_CHECK_EQUAL(dummy_direction(1,0), 5.0);
        KRATOS_CHECK_EQUAL(dummy_direction(1,1), 8.0);
        KRATOS_CHECK_EQUAL(dummy_direction(1,2), 2.0);
        KRATOS_CHECK_EQUAL(dummy_direction(2,0), 6.0);
        KRATOS_CHECK_EQUAL(dummy_direction(2,1), 9.0);
        KRATOS_CHECK_EQUAL(dummy_direction(2,2), 3.0);

    }

    /**
    * Check whether the tensor invariant calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityTensorInvariantsCalculation, KratosParticleMechanicsFastSuite)
    {
        double I1, J2, J3;
        auto positive_hydrostatic = CreateHydrostaticStressVector();
        MPMStressPrincipalInvariantsUtility::CalculateTensorInvariants(positive_hydrostatic, I1, J2, J3);

        KRATOS_CHECK_LESS_EQUAL((9.0e6 - I1)/I1, tolerance);
        KRATOS_CHECK_LESS_EQUAL((0.0 - J2)/J2, tolerance);
        KRATOS_CHECK_LESS_EQUAL((0.0 - J3)/J3, tolerance);

        auto random_sorted = CreateRandomSortedStressVector();
        MPMStressPrincipalInvariantsUtility::CalculateTensorInvariants(random_sorted, I1, J2, J3);

        KRATOS_CHECK_LESS_EQUAL((  7.0e6      - I1)/I1, tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 56.333333e6 - J2)/J2, tolerance);
        KRATOS_CHECK_LESS_EQUAL((-18.740722e6 - J3)/J3, tolerance);

        auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
        MPMStressPrincipalInvariantsUtility::CalculateTensorInvariants(random_unsorted_large, I1, J2, J3);

        KRATOS_CHECK_LESS_EQUAL((   -2.27e+07 - I1)/I1, tolerance);
        KRATOS_CHECK_LESS_EQUAL((-1.70812e+14 - J2)/J2, tolerance);
        KRATOS_CHECK_LESS_EQUAL((  4.8859e+20 - J3)/J3, tolerance);

    }

    /**
    * Check whether the tensor invariant derivatives calculation is performed correctly
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMStressPrincipalInvariantsUtilityTensorInvariantDerivativesCalculation, KratosParticleMechanicsFastSuite)
    {
        Vector derivative_I1, derivative_J2, derivative_J3;
        auto negative_hydrostatic  = CreateNegativeHydrostaticStressVector();
        MPMStressPrincipalInvariantsUtility::CalculateTensorInvariantsDerivatives(negative_hydrostatic, derivative_I1, derivative_J2, derivative_J3);

        KRATOS_CHECK_LESS_EQUAL(( 1.0 - derivative_I1[0])/derivative_I1[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 1.0 - derivative_I1[1])/derivative_I1[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 1.0 - derivative_I1[2])/derivative_I1[2], tolerance);

        KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J2[0])/derivative_J2[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J2[1])/derivative_J2[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J2[2])/derivative_J2[2], tolerance);

        KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J3[0])/derivative_J3[0], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J3[1])/derivative_J3[1], tolerance);
        KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J3[2])/derivative_J3[2], tolerance);

        // auto random_unsorted_large = CreateRandomUnsortedLargeStressVector();
        // MPMStressPrincipalInvariantsUtility::CalculateTensorInvariantsDerivatives(random_unsorted_large, derivative_I1, derivative_J2, derivative_J3);

        // KRATOS_CHECK_LESS_EQUAL(( 1.0 - derivative_I1[0])/derivative_I1[0], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 1.0 - derivative_I1[1])/derivative_I1[1], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 1.0 - derivative_I1[2])/derivative_I1[2], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_I1[3])/derivative_I1[3], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_I1[4])/derivative_I1[4], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_I1[5])/derivative_I1[5], tolerance);

        // KRATOS_CHECK_LESS_EQUAL(( -11433333.333333 - derivative_J2[0])/derivative_J2[0], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( -12166666.666667 - derivative_J2[1])/derivative_J2[1], tolerance);
        // KRATOS_CHECK_LESS_EQUAL((   -733333.333333 - derivative_J2[2])/derivative_J2[2], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( random_unsorted_large[3] - derivative_J2[3])/derivative_J2[3], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( random_unsorted_large[4] - derivative_J2[4])/derivative_J2[4], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( random_unsorted_large[5] - derivative_J2[5])/derivative_J2[5], tolerance);

        // std::cout << "TEST J3" << derivative_J3 << std::endl;

        // KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J3[0])/derivative_J3[0], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J3[1])/derivative_J3[1], tolerance);
        // KRATOS_CHECK_LESS_EQUAL(( 0.0 - derivative_J3[2])/derivative_J3[2], tolerance);
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
        double p, q, lode_angle;
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
