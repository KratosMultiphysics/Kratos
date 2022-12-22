//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"

namespace Kratos
{

/**
 * @class IntegrationUtilities
 * @ingroup KratosCore
 * @brief Utilities to integrate in different cases
 * @author Vicente Mataix Ferrandiz
 */
class IntegrationUtilities
{
public:
    /** 
     * @brief This method returns the integration order for the exact mass matrix evaluation
     * @param rGeometry The geometry considered
     * @tparam TPointType The Point type
     * @return The integration order
     */
    template<class TPointType>
    static GeometryData::IntegrationMethod GetIntegrationMethodForExactMassMatrixEvaluation(const Geometry<TPointType>& rGeometry)
    {
        GeometryData::IntegrationMethod integration_method = rGeometry.GetDefaultIntegrationMethod();
        if (integration_method == GeometryData::IntegrationMethod::GI_GAUSS_1)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        else if(integration_method == GeometryData::IntegrationMethod::GI_GAUSS_2)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_3;
        else if(integration_method == GeometryData::IntegrationMethod::GI_GAUSS_3)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;
        else if(integration_method == GeometryData::IntegrationMethod::GI_GAUSS_4)
            integration_method = GeometryData::IntegrationMethod::GI_GAUSS_5;
        return integration_method;
    }

    /** 
     * @brief This method calculates and returns the domain size of the geometry from any geometry in a generic manner
     * @param rGeometry The geometry considered
     * @param rIntegrationPoints The integration points considered
     * @tparam TPointType The Point type
     * @return double value contains volume.
     */
    template<class TPointType>
    static inline double ComputeDomainSize(
        const Geometry<TPointType>& rGeometry,
        const typename Geometry<TPointType>::IntegrationPointsArrayType& rIntegrationPoints
        )
    {
        const auto number_gp = rIntegrationPoints.size();
        Vector temp(number_gp);
        temp = rGeometry.DeterminantOfJacobian(temp);
        double domain_size = 0.0;
        for (unsigned int i = 0; i < number_gp; ++i) {
            domain_size += temp[i] * rIntegrationPoints[i].Weight();
        }
        return domain_size;
    }

    /** 
     * @brief This method calculates and returns the volume of the geometry from a 3D geometry
     * @param rGeometry The geometry considered
     * @tparam TPointType The Point type
     * @return double value contains volume.
     */
    template<class TPointType>
    static inline double ComputeArea2DGeometry(const Geometry<TPointType>& rGeometry)
    {
        const auto integration_method = rGeometry.GetDefaultIntegrationMethod();
        const auto& r_integration_points = rGeometry.IntegrationPoints( integration_method );
        double volume = 0.0;
        Matrix J(2, 2);
        for ( unsigned int i = 0; i < r_integration_points.size(); i++ ) {
            rGeometry.Jacobian( J, i, integration_method);
            volume += (J(0,0)*J(1,1)-J(0,1)*J(1,0)) * r_integration_points[i].Weight();
        }

        return volume;
    }

    /** 
     * @brief This method calculates and returns the volume of the geometry from a 3D geometry
     * @param rGeometry The geometry considered
     * @tparam TPointType The Point type
     * @return double value contains volume.
     */
    template<class TPointType>
    static inline double ComputeVolume3DGeometry(const Geometry<TPointType>& rGeometry)
    {
        const auto integration_method = rGeometry.GetDefaultIntegrationMethod();
        const auto& r_integration_points = rGeometry.IntegrationPoints( integration_method );
        double volume = 0.0;
        Matrix J(3, 3);
        for ( unsigned int i = 0; i < r_integration_points.size(); i++ ) {
            rGeometry.Jacobian( J, i, integration_method);
            volume += J(0,0)*(J(1,1)*J(2,2) - J(1,2)*J(2,1)) - J(0,1)*(J(1,0)*J(2,2) - J(1,2)*J(2,0)) + J(0,2)*(J(1,0)*J(2,1) - J(1,1)*J(2,0)) * r_integration_points[i].Weight();
        }

        return volume;
    }

};

}  // namespace Kratos.
