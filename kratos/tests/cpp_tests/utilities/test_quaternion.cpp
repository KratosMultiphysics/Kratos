//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquín Irazabal Gonzalez
//

// System includes
#include <limits>

// External includes


// Project includes
#include "testing/testing.h"
#include "includes/global_variables.h"

// Utility includes
#include "utilities/quaternion.h"

namespace Kratos
{
    namespace Testing
    {
        /// Tests

		/** Checks if the squared norm of a Quaternion is being calculated correctly.
         * x*x + y*y + z*z + w*w
         */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionSquaredNorm, KratosCoreFastSuite)
        {
            Quaternion<double> quaternion = Quaternion<double>(0.5, 0.4, 1.8, 2.4);
            const double squared_norm = quaternion.squaredNorm();

            KRATOS_CHECK_EQUAL(squared_norm, 9.41);
		}

		/** Checks if the norm of a Quaternion is being calculated correctly.
         * sqrt(x*x + y*y + z*z + w*w)
         */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionNorm, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Quaternion<double> quaternion = Quaternion<double>(0.7, 2.3, 0.9, 3.7);
            const double norm = quaternion.norm();

            KRATOS_CHECK_NEAR(norm, 4.5033321, tolerance);
        }

		/** Checks if the normalization of a Quaternion (make it a unit quaternion)
         * is being calculated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionNormalize, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Quaternion<double> quaternion = Quaternion<double>(5.2, 1.3, 4.5, 0.1);
            quaternion.normalize();

            KRATOS_CHECK_NEAR(quaternion.W(), 0.7429330, tolerance);
            KRATOS_CHECK_NEAR(quaternion.X(), 0.1857332, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Y(), 0.6429228, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Z(), 0.0142872, tolerance);
        }

		/** Checks if the conjugate of a Quaternion, which represents the opposite rotation,
         * is being calculated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionConjugate, KratosCoreFastSuite)
        {
            Quaternion<double> quaternion = Quaternion<double>(0.1, 0.2, 4.1, 3.5);
            Quaternion<double> quaternion_conjugate = quaternion.conjugate();

            KRATOS_CHECK_EQUAL(quaternion_conjugate.W(), 0.1);
            KRATOS_CHECK_EQUAL(quaternion_conjugate.X(),-0.2);
            KRATOS_CHECK_EQUAL(quaternion_conjugate.Y(),-4.1);
            KRATOS_CHECK_EQUAL(quaternion_conjugate.Z(),-3.5);
        }

		/** Checks if the transformation from Quaternion to Rotation Matrix is being
         * calculated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionToRotationMatrix, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Quaternion<double> quaternion = Quaternion<double>(0.4, 2.5, 0.3, 0.8);
            quaternion.normalize();

            BoundedMatrix<double, 3, 3> rotation_matrix;
            quaternion.ToRotationMatrix(rotation_matrix);

            BoundedMatrix<double, 3, 3> actual_rotation_matrix;
            actual_rotation_matrix(0,0) = 0.7955182; actual_rotation_matrix(0,1) = 0.1204482; actual_rotation_matrix(0,2) = 0.5938376;
            actual_rotation_matrix(1,0) = 0.2997199; actual_rotation_matrix(1,1) =-0.9299720; actual_rotation_matrix(1,2) =-0.2128852;
            actual_rotation_matrix(2,0) = 0.5266107; actual_rotation_matrix(2,1) = 0.3473389; actual_rotation_matrix(2,2) =-0.7759104;

            for (unsigned i = 0; i < 3; ++i)
                for (unsigned j = 0; j < 3; ++j)
                    KRATOS_CHECK_NEAR(rotation_matrix(i,j), actual_rotation_matrix(i,j), tolerance);
        }

		/** Checks if the transformation from Quaternion to Euler Angles is being
         * calculated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionToEulerAngles, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Quaternion<double> quaternion = Quaternion<double>(-1.3, 3.5,-0.1, 2.2);
            quaternion.normalize();

            array_1d<double, 3> euler_angles;
            quaternion.ToEulerAngles(euler_angles);

            KRATOS_CHECK_NEAR(euler_angles[0], 2.07594, tolerance);
            KRATOS_CHECK_NEAR(euler_angles[1],-1.88068, tolerance);
            KRATOS_CHECK_NEAR(euler_angles[2], 2.13307, tolerance);
        }

		/** Check if the transformation from Quaternion to Rotation Vector is being
         * calculated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionToRotationVector, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Quaternion<double> quaternion = Quaternion<double>(-3.7,-1.1, 4.2, 3.1);
            quaternion.normalize();

            Vector rotation_vector(3);
            quaternion.ToRotationVector(rotation_vector);

            Vector rotation_vector_comp(3);
            quaternion.ToRotationVector(rotation_vector_comp[0],
                                        rotation_vector_comp[1],
                                        rotation_vector_comp[2]);

            Vector actual_rotation_vector(3);
            actual_rotation_vector[0] = 0.397708; actual_rotation_vector[1] =-1.51852; actual_rotation_vector[2] =-1.12081;

            for (unsigned i = 0; i < 3; ++i) {
                KRATOS_CHECK_NEAR(rotation_vector[i], actual_rotation_vector[i], tolerance);
                KRATOS_CHECK_NEAR(rotation_vector_comp[i], actual_rotation_vector[i], tolerance);
            }
        }

		/** Check if the rotation of Vectors using Quaternions is being
         * calculated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionRotateVector3, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            Quaternion<double> quaternion = Quaternion<double>(-2.4, 0.3, 6.1,-1.1);
            quaternion.normalize();

            array_1d<double, 3> vector;
            vector[0] = 2.7; vector[1] =-5.1; vector[2] =-0.6;

            array_1d<double, 3> rotated_vector;

            quaternion.RotateVector3(vector, rotated_vector);
            quaternion.RotateVector3(vector);

            Vector actual_rotated_vector(3);
            actual_rotated_vector[0] =-1.39401; actual_rotated_vector[1] =-4.09286; actual_rotated_vector[2] = 3.86849;

            for (unsigned i = 0; i < 3; ++i) {
                KRATOS_CHECK_NEAR(rotated_vector[i], actual_rotated_vector[i], tolerance);
                KRATOS_CHECK_NEAR(vector[i], actual_rotated_vector[i], tolerance);

            }
        }

		/** Check if the Identity Quaternion (that represents a Zero rotation) is being
         * generated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionIdentity, KratosCoreFastSuite)
        {
            Quaternion<double> quaternion = Quaternion<double>::Identity();

            KRATOS_CHECK_EQUAL(quaternion.W(), 1.0);
            KRATOS_CHECK_EQUAL(quaternion.X(), 0.0);
            KRATOS_CHECK_EQUAL(quaternion.Y(), 0.0);
            KRATOS_CHECK_EQUAL(quaternion.Z(), 0.0);
        }


		/** Check if a Quaternion that represents a rotation of an angle 'radians' around the Axis (x, y, z)
         * is being generated correctly.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionFromAxisAngle, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            array_1d<double, 3> axis;
            axis[0] =-1.3; axis[1] =-2.1; axis[2] = 0.4;

            double radians = 4.8;

            Quaternion<double> quaternion = Quaternion<double>::FromAxisAngle(axis[0], axis[1], axis[2], radians);
            quaternion.normalize();

            KRATOS_CHECK_NEAR(quaternion.W(),-0.7373937, tolerance);
            KRATOS_CHECK_NEAR(quaternion.X(),-0.3509602, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Y(),-0.5669357, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Z(), 0.1079878, tolerance);
        }

		/** Check if a Quaternion is being generated correctly from a Rotation
         * Vector.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionFromRotationVector, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            array_1d<double, 3> rotation_vector;
            rotation_vector[0] = 0.3; rotation_vector[1] =-7.5; rotation_vector[2] =-3.9;

            Quaternion<double> quaternion = Quaternion<double>::FromRotationVector(rotation_vector);
            quaternion.normalize();

            Quaternion<double> quaternion_comp = Quaternion<double>::FromRotationVector(rotation_vector[0], rotation_vector[1], rotation_vector[2]);
            quaternion_comp.normalize();

            KRATOS_CHECK_NEAR(quaternion.W(),-0.4644623, tolerance);
            KRATOS_CHECK_NEAR(quaternion.X(),-0.0314087, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Y(), 0.7852186, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Z(), 0.4083137, tolerance);

            KRATOS_CHECK_NEAR(quaternion_comp.W(),-0.4644623, tolerance);
            KRATOS_CHECK_NEAR(quaternion_comp.X(),-0.0314087, tolerance);
            KRATOS_CHECK_NEAR(quaternion_comp.Y(), 0.7852186, tolerance);
            KRATOS_CHECK_NEAR(quaternion_comp.Z(), 0.4083137, tolerance);
        }

		/** Check if a Quaternion is being generated correctly from a Rotation
         * Matrix.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionFromRotationMatrix, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            BoundedMatrix<double, 3, 3> rotation_matrix;
            rotation_matrix(0,0) = 0.5; rotation_matrix(0,1) = 0.3; rotation_matrix(0,2) = 0.9;
            rotation_matrix(1,0) = 0.0; rotation_matrix(1,1) =-0.4; rotation_matrix(1,2) =-0.7;
            rotation_matrix(2,0) = 0.2; rotation_matrix(2,1) = 0.5; rotation_matrix(2,2) =-0.8;

            Quaternion<double> quaternion = Quaternion<double>::FromRotationMatrix(rotation_matrix);
            quaternion.normalize();

            KRATOS_CHECK_NEAR(quaternion.W(), 0.3789054, tolerance);
            KRATOS_CHECK_NEAR(quaternion.X(), 0.8525371, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Y(), 0.0947264, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Z(), 0.3473299, tolerance);
        }

		/** Check if a Quaternion is being generated correctly from a Euler
         * Angles.
		 */

        KRATOS_TEST_CASE_IN_SUITE(QuaternionFromEulerAngles, KratosCoreFastSuite)
        {
            constexpr double tolerance = 1e-5;

            array_1d<double, 3> euler_angles;
            euler_angles[0] = 4.9; euler_angles[1] =-7.2; euler_angles[2] = 0.3;

            Quaternion<double> quaternion = Quaternion<double>::FromEulerAngles(euler_angles);
            quaternion.normalize();

            KRATOS_CHECK_NEAR(quaternion.W(), 0.768422, tolerance);
            KRATOS_CHECK_NEAR(quaternion.X(), 0.294841, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Y(),-0.32999, tolerance);
            KRATOS_CHECK_NEAR(quaternion.Z(),-0.46228, tolerance);
        }
    } // namespace Testing
}  // namespace Kratos.
