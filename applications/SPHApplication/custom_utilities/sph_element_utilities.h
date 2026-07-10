
#pragma once 

#include <cmath>
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "sph_application_variables.h"

/**
 * @class SPHElementUtilities
 * @brief 
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

    static void GetLocalBodyForces(Element& rElement, VectorType& body_force);
    
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
    * @brief Computes pressure and shear wave speeds from material properties.
    * */
    static void ComputeWaveSpeed(double& PressureWaveSpeed, double& ShearWaveSpeed, const Properties& rProperties); 


};

}