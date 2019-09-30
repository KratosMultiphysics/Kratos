//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License 
//                   Kratos default license: kratos/license.txt
//
//  Main author:    
//

#if !defined(KRATOS_QUADRATURE_POINTS_UTILITY_H_INCLUDED)
#define KRATOS_QUADRATURE_POINTS_UTILITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "geometries/quadrature_point_geometry.h"

namespace Kratos
{
    ///@name Kratos Classes
    ///@{

    /// A Class for the creation of integration points
    template<class TPointType>
    class CreateQuadraturePointsUtility
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef Geometry<TPointType> GeometryType;

        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        ///@}
        ///@name Operations
        ///@{

        static std::vector<typename GeometryType::Pointer> Create(
            typename GeometryType::Pointer pGeometry) {
            KRATOS_TRY;

            auto integration_points = pGeometry->IntegrationPoints();
            auto default_method = pGeometry->GetDefaultIntegrationMethod();
            auto r_N = pGeometry->ShapeFunctionsValues();

            std::vector<typename GeometryType::Pointer> geometry_pointer_vector(integration_points.size());

            for (IndexType i = 0; i < integration_points.size(); ++i)
            {
                Matrix N_i = ZeroMatrix(1, pGeometry->size());
                for (IndexType j = 0; j < pGeometry->size(); ++j)
                {
                    N_i(0, j) = r_N(0, j);
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    default_method,
                    integration_points[i],
                    N_i,
                    pGeometry->ShapeFunctionLocalGradient(i));

                if (pGeometry->WorkingSpaceDimension() == 1 && pGeometry->LocalSpaceDimension() == 1)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 1>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
                if(pGeometry->WorkingSpaceDimension() == 2 && pGeometry->LocalSpaceDimension() == 2)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 2>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
                else if (pGeometry->WorkingSpaceDimension() == 3 && pGeometry->LocalSpaceDimension() == 2)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 3, 2>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
                else if (pGeometry->WorkingSpaceDimension() == 3 && pGeometry->LocalSpaceDimension() == 3)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 3>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
            }
            return geometry_pointer_vector;

            KRATOS_CATCH("");
        }

        static std::vector<typename GeometryType::Pointer> Create(
            typename GeometryType::Pointer pGeometry,
            GeometryData::IntegrationMethod ThisIntegrationMethod) {
            KRATOS_TRY;

            auto integration_points = pGeometry->IntegrationPoints(ThisIntegrationMethod);
            auto r_N = pGeometry->ShapeFunctionsValues(ThisIntegrationMethod);

            std::vector<typename GeometryType::Pointer> geometry_pointer_vector(integration_points.size());

            for (IndexType i = 0; i < integration_points.size(); ++i)
            {
                Matrix N_i = ZeroMatrix(1, pGeometry->size());
                for (IndexType j = 0; j < pGeometry->size(); ++j)
                {
                    N_i(0, j) = r_N(0, j);
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    ThisIntegrationMethod,
                    integration_points[i],
                    N_i,
                    pGeometry->ShapeFunctionLocalGradient(i, ThisIntegrationMethod));

                if (pGeometry->WorkingSpaceDimension() == 1 && pGeometry->LocalSpaceDimension() == 1)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 1>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
                if (pGeometry->WorkingSpaceDimension() == 2 && pGeometry->LocalSpaceDimension() == 2)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 2>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
                else if (pGeometry->WorkingSpaceDimension() == 3 && pGeometry->LocalSpaceDimension() == 2)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 3, 2>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
                else if (pGeometry->WorkingSpaceDimension() == 3 && pGeometry->LocalSpaceDimension() == 3)
                    geometry_pointer_vector[i] = typename Geometry<TPointType>::Pointer(
                        Kratos::make_shared<
                        QuadraturePointGeometry<TPointType, 3>>(
                            pGeometry->Points(),
                            data_container,
                            pGeometry.get()));
            }
            return geometry_pointer_vector;

            KRATOS_CATCH("");
        }

        ///@}

    };
    ///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINTS_UTILITY_H_INCLUDED defined
