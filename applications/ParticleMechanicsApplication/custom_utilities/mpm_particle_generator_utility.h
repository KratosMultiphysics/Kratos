//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//


#ifndef KRATOS_MPM_PARTICLE_GENERATOR_UTILITY
#define KRATOS_MPM_PARTICLE_GENERATOR_UTILITY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "particle_mechanics_application_variables.h"


namespace Kratos
{
namespace MPMParticleGeneratorUtility
{

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    /**
     * @brief Function that return matrix of shape function value for 16 particles.
     * @details It is only possible to be used in 2D Triangular.
     */
    Matrix MP16ShapeFunctions();


    /**
     * @brief Function that return matrix of shape function value for 33 particles.
     * @details It is only possible to be used in 2D Triangular.
     */
    Matrix MP33ShapeFunctions();


    /**
     * @brief Construct material points or particles from given initial mesh
     * @details Generating particles using a designated shape functions
     */
    void GenerateMaterialPointElement(  ModelPart& rBackgroundGridModelPart,
                                        ModelPart& rInitialModelPart,
                                        ModelPart& rMPMModelPart,
                                        bool IsAxisSymmetry = false,
                                        bool IsMixedFormulation = false);


    /**
     * @brief Function to Initiate material point condition.
     * @details Generating particle condition using a designated shape functions
     */
    void GenerateMaterialPointCondition(    ModelPart& rBackgroundGridModelPart,
                                            ModelPart& rInitialModelPart,
                                            ModelPart& rMPMModelPart);

}; // end namespace MPMParticleGeneratorUtility
} // end namespace Kratos

#endif // KRATOS_MPM_PARTICLE_GENERATOR_UTILITY


