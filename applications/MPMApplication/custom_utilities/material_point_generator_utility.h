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


#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/quadrature_points_utility.h"
#include "mpm_application_variables.h"
#include "custom_utilities/mpm_math_utilities.h"


namespace Kratos::MaterialPointGeneratorUtility
{

    using IndexType =  std::size_t;

    using SizeType =  std::size_t;

    using GeometryType =  Geometry<Node>;

    using IntegrationMethod =  GeometryData::IntegrationMethod;

    /// Get integration weights of the geometry for the given integration method
    void GetIntegrationPointVolumes(
            const GeometryType& rGeom,
            const IntegrationMethod IntegrationMethod,
            Vector& rIntVolumes);

    /// Get integration method and shape function values for the given element
    void DetermineIntegrationMethodAndShapeFunctionValues(
            const GeometryType& rGeom,
            const SizeType MaterialPointsPerElement,
            IntegrationMethod& rIntegrationMethod,
            Matrix& rN);

    /// Get integration method and shape function values for the given condition
    void DetermineGeometryIntegrationMethod(
            const GeometryType& rGeom,
            const SizeType MaterialPointsPerCondition,
            IndexType& rNumPointsPerSpan);

    /**
     * @brief Construct material points from given initial mesh
     * @details Generating material points using a designated shape functions
     */
    template<SizeType TDimension>
    void GenerateMaterialPointElement(
            ModelPart& rBackgroundGridModelPart,
            ModelPart& rInitialModelPart,
            ModelPart& rMPMModelPart,
            bool IsMixedFormulation=false);

    /**
     * @brief Function to Initiate material point condition.
     * @details Generating material point condition using a designated shape functions
     */
    template<SizeType TDimension>
    void GenerateMaterialPointCondition(
            ModelPart& rBackgroundGridModelPart,
            ModelPart& rInitialModelPart,
            ModelPart& rMPMModelPart);

    /**
     * @brief Function to add dofs to elements for Lagrange multiplier.
     * @details non-conforming Perturbed Lagrangian or Lagrange multiplier conditions
     */
    void GenerateLagrangeNodes(ModelPart& rBackgroundGridModelPart);
    /**
     * @brief Function to Initiate material point condition.
     * @details Generating material point condition using a designated shape functions
     */
    void KRATOS_API(MPM_APPLICATION)
        GenerateMaterialPointCondition(
                ModelPart& rBackgroundGridModelPart,
                ModelPart& rInitialModelPart,
                ModelPart& rMPMModelPart);

} // end namespace Kratos::MaterialPointGeneratorUtility
