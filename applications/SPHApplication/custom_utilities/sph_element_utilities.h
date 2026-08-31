//  ____  ____  _   _                   _ _           _   _             
// / ___||  _ \| | | | __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __  
// \___ \| |_) | |_| |/ _` | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \ 
//  ___) |  __/|  _  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
// |____/|_|   |_| |_|\__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                         |_|   |_|                                    

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

#pragma once 

#include <cmath>
#include <algorithm>
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "sph_application_variables.h"

/**
 * @class SPHElementUtilities
 * @brief A collection of utility functions for SPH elements. Created to minimize code duplication.
 * @details The methods are static, so it can be called without constructing the class
 */

namespace Kratos
{
class SPHElementUtilities
{
public:

    using SizeType = std::size_t;
    using MatrixType = Matrix;
    using VectorType = Vector;

    /**
     * @brief Computes the body forces acting on the particle, the body is assumed to be constant in time  and equal to the gravity vector. 
     */
    static void GetLocalBodyForces(Element& rElement, VectorType& body_force, const int Step = 0);
    
    /**
     * @brief The flag is set to true if the lumped mass matrix must be computed.
     */
    static bool ComputeLumpedMassMatrix(const Properties& rProperties, const ProcessInfo& rProcessInfo);

    /**
    * @brief Computes the linear elastic acoustic tensor using the material elastic properties.
    * @param rNormal The propagation direction.
    */
    static void ComputeLinearElasticAcousticTensor(MatrixType& rAcousticTensor, const VectorType& rNormal, const Properties& rProperties);

    /**
    * @brief Computes the linear reconstructed displacement jump between two neighboring particles in the current configuration.
    */
    static void ComputeParticleJump(VectorType& rJumpVector, Element& rThisParticle, Element& rThisNeighbour, VectorType& rInitialDistance, const ProcessInfo& rProcessInfo);

    /**
    * @brief Computes the linear reconstructed velocity jump between two neighboring particles.
    */
    static void ComputeVelocityJump(VectorType& rJumpVector, Element& rThisParticle, Element& rThisNeighbour, VectorType& rInitialDistance, const int Step = 0);

    /**
    * @brief Computes pressure and shear wave speeds from material properties.
    * */
    static void ComputeWaveSpeed(double& PressureWaveSpeed, double& ShearWaveSpeed, const Properties& rProperties); 

    /**
     * @brief Computes pressure and shear wave speeds from material properties and the current deformation state.
     */
    static void ComputeDeformationDependentWaveSpeed(double& rPressureWaveSpeed, double& rShearWaveSpeed, const MatrixType& rDeformationGradient, const Properties& rProperties);
    
    /**
     * @brief Calculates the 2D and 3D deformation gradient.
     * @param rB The deformation gradient matrix.
     * @param rF The deformation gradient.
     * @param rDW_DX The kernel gradients.
     */
    static void Calculate2DB(MatrixType& rB, const MatrixType& rF, const MatrixType& rDW_DX, const SizeType NumberOfNeighbours);

    static void Calculate3DB(MatrixType& rB, const MatrixType& rF, const MatrixType& rDW_DX, const SizeType NumberOfNeighbours);

    /**
    * @brief Converts a non symmetric tensor to a vector.
    * @details This method extends the one in math_utils.h
    */
    template<std::size_t TDim>
    static VectorType NonSymmetricTensorToVector(const MatrixType& rTensor); 


};

}
