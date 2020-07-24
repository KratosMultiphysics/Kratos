//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//

// Project includes
#include "custom_utilities/mpm_explicit_utilities.h"

namespace Kratos
{
    void MPMExplicitUtilities::CalculateAndAddExplicitInternalForce(
        Element& rElement,
        const Matrix& rDN_DX,
        const Vector& rMPStress,
        const double rMPVolume,
        const SizeType StrainSize,
        Vector& rRightHandSideVector)
    {
        KRATOS_TRY


        // Add in explicit internal force calculation (f_i = V * Sum_j [s_ij N_,j])
        // Refer to link for notation https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-Constitutive-Law-class
        GeometryType& rGeom = rElement.GetGeometry();
        const SizeType dimension = rGeom.WorkingSpaceDimension();
        const SizeType number_of_nodes = rGeom.PointsNumber();
        array_1d<double, 3> nodal_force_internal_normal = ZeroVector(3);
        for (IndexType i = 0; i < number_of_nodes; i++) 
        {
            if (dimension == 2 && StrainSize == 3) 
            {
                // StressVec = s00 s11  s01
                // Index        0   1   2 
                //f_x = V*(s_xx*dNdX + s_xy*dNdY)
                nodal_force_internal_normal[0] = rMPVolume *
                    (rMPStress[0] * rDN_DX(i, 0) +
                        rMPStress[2] * rDN_DX(i, 1));
                //f_y = V*(s_yy*dNdY + s_yx*dNdX)
                nodal_force_internal_normal[1] = rMPVolume *
                    (rMPStress[1] * rDN_DX(i, 1) +
                        rMPStress[2] * rDN_DX(i, 0));
            } 
            else if (dimension == 2 && StrainSize == 4) 
            {
                KRATOS_ERROR 
                    << "Call CalcuateAndAddAxisymmetricExplicitInternalForce instead of CalcuateAndAddExplicitInternalForce" 
                    << std::endl;
            } 
            else if (dimension == 3 && StrainSize == 6) 
            {
                // StressVec = s00 s11 s22   s01   s12   s02
                // Index        0   1   2     3     4      5
                //f_x = V*(s_xx*dNdX + s_xy*dNdY + s_xz*dNdZ)
                nodal_force_internal_normal[0] = rMPVolume *
                    (rMPStress[0] * rDN_DX(i, 0) +
                        rMPStress[3] * rDN_DX(i, 1) +
                        rMPStress[5] * rDN_DX(i, 2));
                //f_y = V*(s_yy*dNdY + s_yx*dNdX + s_yz*dNdZ)
                nodal_force_internal_normal[1] = rMPVolume *
                    (rMPStress[1] * rDN_DX(i, 1) +
                        rMPStress[3] * rDN_DX(i, 0) +
                        rMPStress[4] * rDN_DX(i, 2));
                //f_z = V*(s_zz*dNdZ + s_zx*dNdX + s_zy*dNdY)
                nodal_force_internal_normal[2] = rMPVolume *
                    (rMPStress[2] * rDN_DX(i, 2) +
                        rMPStress[5] * rDN_DX(i, 0) +
                        rMPStress[4] * rDN_DX(i, 1));
                rRightHandSideVector[dimension * i + 2] -= nodal_force_internal_normal[2]; //minus sign, internal forces
            } 
            else
            {
                KRATOS_ERROR << "Dimension = " << dimension << " and strain size = " << StrainSize
                    << " are invalid for MPM explicit internal force calculation." << std::endl;
            }
            rRightHandSideVector[dimension * i] -= nodal_force_internal_normal[0]; //minus sign, internal forces
            rRightHandSideVector[dimension * i + 1] -= nodal_force_internal_normal[1]; //minus sign, internal forces
        }

        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void MPMExplicitUtilities::CalculateAndAddAxisymmetricExplicitInternalForce(
        Element& rElement, 
        const Matrix& rDN_DX, 
        const Vector& rMPStress, 
        const double rMPVolume, 
        const SizeType StrainSize, 
        const double AxisymmetricRadius,
        Vector& rRightHandSideVector)
    {
        KRATOS_TRY

        // Add in explicit internal force calculation (Fint = Volume*divergence(sigma))
        // Refer to link for notation https://github.com/KratosMultiphysics/Kratos/wiki/How-to-use-the-Constitutive-Law-class
        GeometryType& rGeom = rElement.GetGeometry();
        const SizeType dimension = rGeom.WorkingSpaceDimension();
        const SizeType number_of_nodes = rGeom.PointsNumber();
        array_1d<double, 3> nodal_force_internal_normal = ZeroVector(3);
        const Matrix& r_N = rElement.GetGeometry().ShapeFunctionsValues();

        KRATOS_ERROR_IF_NOT(dimension == 2 && StrainSize == 4)
            << "Call CalcuateAndAddExplicitInternalForce instead of CalcuateAndAddAxisymmetricExplicitInternalForce"
            << std::endl;

        for (IndexType i = 0; i < number_of_nodes; i++) {
            // StressVec = srr szz  sthetatheta srz
            // Index        0   1   2           3

            nodal_force_internal_normal[0] = rMPVolume *
                (rMPStress[0] * rDN_DX(i, 0) +
                    rMPStress[2] * r_N(0,i) / AxisymmetricRadius +
                    rMPStress[3] * rDN_DX(i, 1));

            nodal_force_internal_normal[1] = rMPVolume *
                (rMPStress[1] * rDN_DX(i, 1) +
                    rMPStress[3] * rDN_DX(i, 0));
            
            rRightHandSideVector[dimension * i] -= nodal_force_internal_normal[0]; //minus sign, internal forces
            rRightHandSideVector[dimension * i + 1] -= nodal_force_internal_normal[1]; //minus sign, internal forces
        }
        KRATOS_CATCH("")
    }

    void MPMExplicitUtilities::UpdateGaussPointExplicit(
        const ProcessInfo& rCurrentProcessInfo,
        Element& rElement)
    {
        KRATOS_TRY

        const double& rDeltaTime = rCurrentProcessInfo[DELTA_TIME];
        const bool& isCentralDifference = rCurrentProcessInfo.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE);
        GeometryType& rGeom = rElement.GetGeometry();
        const SizeType number_of_nodes = rGeom.PointsNumber();
        const SizeType dimension = rGeom.WorkingSpaceDimension();
        const Matrix& r_N = rElement.GetGeometry().ShapeFunctionsValues();


        // False for stability (typical across research papers). 
        // One may set to true to reduces energy lost from kinematic aliasing
        bool isUpdateMPPositionFromUpdatedMPVelocity = false; 


        const ProcessInfo& rProcessInfo = ProcessInfo();
        std::vector<array_1d<double, 3 > > MP_PreviousVelocity;
        std::vector<array_1d<double, 3 > > MP_PreviousAcceleration; 
        array_1d<double, 3> MP_Velocity = ZeroVector(3);
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY, MP_PreviousVelocity,rProcessInfo);        
        rElement.CalculateOnIntegrationPoints(MP_ACCELERATION, MP_PreviousAcceleration, rProcessInfo);

        const double gamma = (isCentralDifference)
            ? 0.5
            : 1.0; // 0.5 for central difference, 1.0 for forward euler
        if (isCentralDifference) isUpdateMPPositionFromUpdatedMPVelocity = false;
        

        // Advance the material point predictor velocity
        for (IndexType i = 0; i < dimension; i++) {
            MP_Velocity[i] = MP_PreviousVelocity[0][i] + (1.0 - gamma) * rDeltaTime * MP_PreviousAcceleration[0][i];
        }


        // Calculate the MP displacement and acceleration for this timestep
        array_1d<double, 3> delta_xg = ZeroVector(3);
        array_1d<double, 3> MP_Acceleration = ZeroVector(3);
        for (IndexType i = 0; i < number_of_nodes; i++)
        {
            const double nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS);
            if (nodal_mass > std::numeric_limits<double>::epsilon())
            {
                const array_1d<double, 3>& r_nodal_momenta = rGeom[i].FastGetSolutionStepValue(NODAL_MOMENTUM);
                const array_1d<double, 3>& r_current_residual = rGeom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

                // Applicable to central difference only, the value in VELOCITY at the moment is actually the
                // predicted (middle) grid velocity
                const array_1d<double, 3>& r_middle_velocity = rGeom[i].FastGetSolutionStepValue(VELOCITY);

                for (IndexType j = 0; j < dimension; j++)
                {
                    // Update MP acceleration regardless of explicit method
                    MP_Acceleration[j] += r_N(0,i) * r_current_residual[j] / nodal_mass;

                    if (isCentralDifference)
                    {
                        delta_xg[j] += rDeltaTime * r_N(0, i) * r_middle_velocity[j];
                    }
                    else if (!isUpdateMPPositionFromUpdatedMPVelocity)
                    {
                        delta_xg[j] += rDeltaTime * r_N(0, i) * r_nodal_momenta[j] / nodal_mass;
                    }
                }
            }
        }

        // Update the MP Acceleration
        rElement.SetValuesOnIntegrationPoints(MP_ACCELERATION, { MP_Acceleration }, rProcessInfo);

        // Update the MP Velocity corrector
        for (IndexType j = 0; j < dimension; j++)
        {
            MP_Velocity[j] += gamma * rDeltaTime * MP_Acceleration[j];
        }
        rElement.SetValuesOnIntegrationPoints(MP_VELOCITY, { MP_Velocity }, rProcessInfo);

        // Update the MP Position
        if (isUpdateMPPositionFromUpdatedMPVelocity)
        {
            for (IndexType j = 0; j < dimension; j++)
            {
                delta_xg[j] = rDeltaTime * MP_Velocity[j];
            }
        }
        std::vector<array_1d<double, 3 > > xg;
        rElement.CalculateOnIntegrationPoints(MP_COORD, xg, rProcessInfo);
        const array_1d<double, 3>& new_xg = xg[0] + delta_xg;
        rElement.SetValuesOnIntegrationPoints(MP_COORD, { new_xg }, rProcessInfo);

        // Update the MP total displacement
        std::vector<array_1d<double, 3 > > MP_Displacement;
        rElement.CalculateOnIntegrationPoints(MP_DISPLACEMENT, MP_Displacement,rProcessInfo);
        MP_Displacement[0] += delta_xg;
        rElement.SetValuesOnIntegrationPoints(MP_DISPLACEMENT,MP_Displacement,rProcessInfo);

        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void MPMExplicitUtilities::CalculateMUSLGridVelocity(const ProcessInfo& rCurrentProcessInfo,
        Element& rElement)
    {
        KRATOS_TRY

        GeometryType& rGeom = rElement.GetGeometry();
        const SizeType dimension = rGeom.WorkingSpaceDimension();
        const SizeType number_of_nodes = rGeom.PointsNumber();
        const Matrix& r_N = rElement.GetGeometry().ShapeFunctionsValues();

        std::vector<array_1d<double, 3 > > MP_Velocity;
        std::vector<double> MP_Mass;
        rElement.CalculateOnIntegrationPoints(MP_VELOCITY, MP_Velocity, rCurrentProcessInfo);
        rElement.CalculateOnIntegrationPoints(MP_MASS, MP_Mass, rCurrentProcessInfo);

        for (IndexType i = 0; i < number_of_nodes; i++)
        {
            const double& r_nodal_mass = rGeom[i].FastGetSolutionStepValue(NODAL_MASS);

            if (r_nodal_mass > std::numeric_limits<double>::epsilon())
            {
                array_1d<double, 3>& r_current_velocity = rGeom[i].FastGetSolutionStepValue(VELOCITY);
                for (IndexType j = 0; j < dimension; j++)
                {
                    r_current_velocity[j] += r_N(0,i) * MP_Mass[0] * MP_Velocity[0][j] / r_nodal_mass;
                }
            }
        }

        KRATOS_CATCH("")
    }

    /***********************************************************************************/
    /***********************************************************************************/

    void MPMExplicitUtilities::CalculateExplicitKinematics(
        const ProcessInfo& rCurrentProcessInfo,
        Element& rElement,
        const Matrix& rDN_DX,
        Vector& rMPStrain,
        Matrix& rDeformationGradientIncrement,
        const SizeType StrainSize)
    {
        KRATOS_TRY

        GeometryType rGeom = rElement.GetGeometry();
        const double deltaTime = rCurrentProcessInfo[DELTA_TIME];
        const SizeType dimension = rGeom.WorkingSpaceDimension();
        const SizeType number_of_nodes = rGeom.PointsNumber();

        //Calculate velocity gradients
        Matrix velocityGradient = Matrix(dimension, dimension, 0.0);

        for (IndexType nodeIndex = 0; nodeIndex < number_of_nodes; nodeIndex++)
        {
            const array_1d<double, 3 >& nodal_velocity = rGeom[nodeIndex].FastGetSolutionStepValue(VELOCITY);

            for (IndexType i = 0; i < dimension; i++)
            {
                for (IndexType j = 0; j < dimension; j++)
                {
                    velocityGradient(i, j) += nodal_velocity[i] * rDN_DX(nodeIndex, j);
                }
            }
        }

        //Calculate rate of deformation and spin tensors
        const Matrix rateOfDeformation = 0.5 * (velocityGradient + trans(velocityGradient));
        const Matrix spinTensor = velocityGradient - rateOfDeformation;

        //Calculate objective Jaumann strain rate
        const Matrix jaumannRate = rateOfDeformation -
            (prod(spinTensor, rateOfDeformation)) * deltaTime +
            prod((rateOfDeformation * deltaTime), spinTensor);
        const Matrix strainIncrement = deltaTime * jaumannRate;

        // Apply strain increment to strain vector
        rMPStrain(0) += strainIncrement(0, 0); //e_xx
        rMPStrain(1) += strainIncrement(1, 1); //e_yy
        if (dimension == 2 && StrainSize == 3)
        {
            rMPStrain(2) += 2.0 * strainIncrement(0, 1); //e_xy
        }        
        else if((dimension == 3 && StrainSize == 6))
        {
            rMPStrain(2) += strainIncrement(2, 2); //e_zz

            rMPStrain(3) += 2.0 * strainIncrement(0, 1); //e_xy
            rMPStrain(4) += 2.0 * strainIncrement(1, 2); //e_yz
            rMPStrain(5) += 2.0 * strainIncrement(0, 2); //e_xz
        }  
        else
        {
            KRATOS_ERROR << "Dimension = " << dimension << " and strain size = " << StrainSize
                << " are invalid for MPM explicit kinematic calculation." << std::endl;
        }

        // Model compressibility
        rDeformationGradientIncrement = IdentityMatrix(dimension);
        if (rCurrentProcessInfo.GetValue(IS_COMPRESSIBLE)) rDeformationGradientIncrement += strainIncrement;

        KRATOS_CATCH("")
    }

    void MPMExplicitUtilities::CalculateExplicitAsymmetricKinematics(
        const ProcessInfo& rCurrentProcessInfo,
        Element& rElement,
        const Matrix& rDN_DX,
        Vector& rMPStrain,
        Matrix& rDeformationGradientIncrement,
        const SizeType StrainSize,
        const double AxisymmetricRadius)
    {
        KRATOS_TRY

        const GeometryType rGeom = rElement.GetGeometry();
        const double deltaTime = rCurrentProcessInfo[DELTA_TIME];
        const SizeType dimension = rGeom.WorkingSpaceDimension();
        const SizeType number_of_nodes = rGeom.PointsNumber();
        const Matrix& r_N = rElement.GetGeometry().ShapeFunctionsValues();

        //Calculate velocity gradients
        Matrix velocityGradient = Matrix(3, 3, 0.0); // for axisymmetric case

        for (IndexType nodeIndex = 0; nodeIndex < number_of_nodes; nodeIndex++)
        {
            const array_1d<double, 3 >& nodal_velocity = rGeom[nodeIndex].FastGetSolutionStepValue(VELOCITY);

            for (IndexType i = 0; i < dimension; i++)
            {
                for (IndexType j = 0; j < dimension; j++)
                {
                    velocityGradient(i, j) += nodal_velocity[i] * rDN_DX(nodeIndex, j);
                }
            }
        }
        if (dimension == 2 && StrainSize == 4) // axisymmetric case
        {
            for (IndexType nodeIndex = 0; nodeIndex < number_of_nodes; nodeIndex++)
            {
                const array_1d<double, 3 >& nodal_velocity = rGeom[nodeIndex].FastGetSolutionStepValue(VELOCITY);
                velocityGradient(2, 2) += nodal_velocity[0] * r_N(0,nodeIndex) / AxisymmetricRadius;
            }
        }

        //Calculate rate of deformation and spin tensors
        const Matrix rateOfDeformation = 0.5 * (velocityGradient + trans(velocityGradient));
        const Matrix spinTensor = velocityGradient - rateOfDeformation;

        //Calculate objective Jaumann strain rate
        const Matrix jaumannRate = rateOfDeformation -
            (prod(spinTensor, rateOfDeformation)) * deltaTime +
            prod((rateOfDeformation * deltaTime), spinTensor);
        const Matrix strainIncrement = deltaTime * jaumannRate;

        // Apply strain increment to strain vector
        rMPStrain(0) += strainIncrement(0, 0); //e_xx
        rMPStrain(1) += strainIncrement(1, 1); //e_yy
        if (dimension == 2 && StrainSize == 4)
        {
            rMPStrain(2) += strainIncrement(2, 2); //e_theta theta
            rMPStrain(3) += 2.0 * strainIncrement(0, 1); //e_xy
        }
        else
        {
            KRATOS_ERROR << "Dimension = " << dimension << " and strain size = " << StrainSize
                << " are invalid for MPM explicit asymmetric kinematic calculation." << std::endl;
        }

        // Model compressibility
        rDeformationGradientIncrement = IdentityMatrix(3);
        if (rCurrentProcessInfo.GetValue(IS_COMPRESSIBLE)) rDeformationGradientIncrement += strainIncrement;

        KRATOS_CATCH("")
    }
} // namespace Kratos