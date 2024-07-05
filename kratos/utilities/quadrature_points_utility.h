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
#include "geometries/quadrature_point_curve_on_surface_geometry.h"
#include "geometries/quadrature_point_surface_in_volume_geometry.h"
#include "geometries/quadrature_point_coupling_geometry_2d.h"

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
        typedef typename Geometry<TPointType>::Pointer GeometryPointerType;
        typedef typename QuadraturePointCurveOnSurfaceGeometry<TPointType>::Pointer QuadraturePointCurveOnSurfaceGeometryPointer;

        typedef std::size_t SizeType;
        typedef std::size_t IndexType;

        typedef PointerVector<TPointType> PointsArrayType;

        ///@}
        ///@name Operations
        ///@{

        static typename GeometryType::Pointer CreateFromCoordinates(
            typename GeometryType::Pointer pGeometry,
            const array_1d<double, 3>& rCoordinates,
            double integration_weight)
        {
            KRATOS_TRY;

            array_1d<double, 3> local_coordinates;
            pGeometry->PointLocalCoordinates(local_coordinates, rCoordinates);

            return CreateFromLocalCoordinates(*(pGeometry.get()), local_coordinates, integration_weight);

            KRATOS_CATCH("");
        }

        static typename GeometryType::Pointer CreateFromLocalCoordinates(
            GeometryType& rGeometry,
            const array_1d<double, 3>& rLocalCoordinates,
            double integration_weight)
        {
            KRATOS_TRY;

            IntegrationPoint<3> int_p(rLocalCoordinates, integration_weight);
            Vector N;
            rGeometry.ShapeFunctionsValues(N, rLocalCoordinates);
            Matrix N_matrix(1, N.size());
            for (IndexType i = 0; i < N.size(); ++i)
            {
                N_matrix(0, i) = N[i];
            }

            Matrix DN_De;
            rGeometry.ShapeFunctionsLocalGradients(DN_De, rLocalCoordinates);

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                rGeometry.GetDefaultIntegrationMethod(),
                int_p,
                N_matrix,
                DN_De);

            return CreateQuadraturePoint(
                rGeometry.WorkingSpaceDimension(),
                rGeometry.LocalSpaceDimension(),
                data_container,
                rGeometry.Points(),
                &rGeometry);

            KRATOS_CATCH("");
        }


        static GeometryPointerType CreateQuadraturePointCurveOnSurface(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            double LocalTangentU,
            double LocalTangentV,
            GeometryType* pGeometryParent)
        {
            return Kratos::make_shared<
                QuadraturePointCurveOnSurfaceGeometry<TPointType>>(
                    rPoints,
                    rShapeFunctionContainer,
                    LocalTangentU,
                    LocalTangentV,
                    pGeometryParent);
        }

        static GeometryPointerType CreateQuadraturePointCurveOnSurface(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            double LocalTangentU,
            double LocalTangentV)
        {
            return Kratos::make_shared<
                QuadraturePointCurveOnSurfaceGeometry<TPointType>>(
                    rPoints,
                    rShapeFunctionContainer,
                    LocalTangentU,
                    LocalTangentV);
        }

        static GeometryPointerType CreateQuadraturePointCouplingGeometry2D(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainerMaster,
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainerSlave,
            PointsArrayType rPointsMaster,
            PointsArrayType rPointsSlave,
            double LocalTangentMasterU,
            double LocalTangentMasterV,
            double LocalTangentSlaveU,
            double LocalTangentSlaveV,
            GeometryType* pGeometryParentMaster,
            GeometryType* pGeometryParentSlave)
        {
            QuadraturePointCurveOnSurfaceGeometryPointer quadraturePointMaster = Kratos::make_shared<QuadraturePointCurveOnSurfaceGeometry<TPointType>>(
                rPointsMaster,
                rShapeFunctionContainerMaster,
                LocalTangentMasterU,
                LocalTangentMasterV,
                pGeometryParentMaster);

            QuadraturePointCurveOnSurfaceGeometryPointer quadraturePointSlave = Kratos::make_shared<QuadraturePointCurveOnSurfaceGeometry<TPointType>>(
                rPointsSlave,
                rShapeFunctionContainerSlave,
                LocalTangentSlaveU,
                LocalTangentSlaveV,
                pGeometryParentSlave);

            // return Kratos::make_shared<QuadraturePointCouplingGeometry2D<TPointType>>(quadraturePointMaster, quadraturePointSlave,
            //                                                                     rPointsMaster,
            //                                                                     rShapeFunctionContainerMaster,
            //                                                                     pGeometryParentMaster);

            return Kratos::make_shared<QuadraturePointCouplingGeometry2D<TPointType>>(quadraturePointMaster,
                                                                                      quadraturePointSlave);
            

        }

        static GeometryPointerType CreateQuadraturePointSurfaceInVolume(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            Matrix LocalTangentMatrix,
            GeometryType* pGeometryParent)
        {
            return Kratos::make_shared<
                QuadraturePointSurfaceInVolumeGeometry<TPointType>>(
                    rPoints,
                    rShapeFunctionContainer,
                    LocalTangentMatrix,
                    pGeometryParent);
        }

        static GeometryPointerType CreateQuadraturePointSurfaceInVolume(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            Matrix LocalTangentMatrix,
            Vector Normal,
            GeometryType* pGeometryParent)
        {
            return Kratos::make_shared<
                QuadraturePointSurfaceInVolumeGeometry<TPointType>>(
                    rPoints,
                    rShapeFunctionContainer,
                    LocalTangentMatrix,
                    pGeometryParent,
                    Normal);
        }

        static GeometryPointerType CreateQuadraturePointSurfaceInVolume(
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            Matrix LocalTangentMatrix)
        {
            return Kratos::make_shared<
                QuadraturePointSurfaceInVolumeGeometry<TPointType>>(
                    rPoints,
                    rShapeFunctionContainer,
                    LocalTangentMatrix );
        }

        static GeometryPointerType CreateQuadraturePoint(
            SizeType WorkingSpaceDimension,
            SizeType LocalSpaceDimension,
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints,
            GeometryType* pGeometryParent)
        {
            if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                    QuadraturePointGeometry<TPointType, 1>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                    QuadraturePointGeometry<TPointType, 2, 1>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                    QuadraturePointGeometry<TPointType, 3, 1>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                    QuadraturePointGeometry<TPointType, 2>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                    QuadraturePointGeometry<TPointType, 3, 2>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
                return Kratos::make_shared<
                    QuadraturePointGeometry<TPointType, 3>>(
                        rPoints,
                        rShapeFunctionContainer,
                        pGeometryParent);
            else{
                KRATOS_ERROR << "Working/Local space dimension combinations are "
                    << "not provided for QuadraturePointGeometry. WorkingSpaceDimension: "
                    << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
                    <<  std::endl;
            }
        }

        static GeometryPointerType CreateQuadraturePoint(
            SizeType WorkingSpaceDimension,
            SizeType LocalSpaceDimension,
            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod>& rShapeFunctionContainer,
            PointsArrayType rPoints)
        {
            if (WorkingSpaceDimension == 1 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointGeometry<TPointType, 1>>(
                    rPoints,
                    rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointGeometry<TPointType, 2, 1>>(
                    rPoints,
                    rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 1)
                return Kratos::make_shared<
                QuadraturePointGeometry<TPointType, 3, 1>>(
                    rPoints,
                    rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 2 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointGeometry<TPointType, 2>>(
                    rPoints,
                    rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 2)
                return Kratos::make_shared<
                QuadraturePointGeometry<TPointType, 3, 2>>(
                    rPoints,
                    rShapeFunctionContainer);
            else if (WorkingSpaceDimension == 3 && LocalSpaceDimension == 3)
                return Kratos::make_shared<
                QuadraturePointGeometry<TPointType, 3>>(
                    rPoints,
                    rShapeFunctionContainer);
            else {
                KRATOS_ERROR << "Working/Local space dimension combinations are "
                    << "not provided for QuadraturePointGeometry. WorkingSpaceDimension: "
                    << WorkingSpaceDimension << ", LocalSpaceDimension: " << LocalSpaceDimension
                    << std::endl;
            }
        }

        static std::vector<GeometryPointerType> Create(
            GeometryPointerType pGeometry) {
            KRATOS_TRY;

            auto integration_points = pGeometry->IntegrationPoints();
            auto default_method = pGeometry->GetDefaultIntegrationMethod();
            auto r_N = pGeometry->ShapeFunctionsValues();

            std::vector<GeometryPointerType> geometry_pointer_vector(integration_points.size());

            for (IndexType i = 0; i < integration_points.size(); ++i)
            {
                Matrix N_i = ZeroMatrix(1, pGeometry->size());
                for (IndexType j = 0; j < pGeometry->size(); ++j)
                {
                    N_i(0, j) = r_N(i, j);
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    default_method,
                    integration_points[i],
                    N_i,
                    pGeometry->ShapeFunctionLocalGradient(i));

                geometry_pointer_vector[i] = CreateQuadraturePoint(
                    pGeometry->WorkingSpaceDimension(), pGeometry->LocalSpaceDimension(),
                    data_container, pGeometry->Points(), pGeometry.get());
            }
            return geometry_pointer_vector;

            KRATOS_CATCH("");
        }

        static std::vector<GeometryPointerType> Create(
            GeometryPointerType pGeometry,
            GeometryData::IntegrationMethod ThisIntegrationMethod) {
            KRATOS_TRY;

            auto integration_points = pGeometry->IntegrationPoints(ThisIntegrationMethod);
            auto r_N = pGeometry->ShapeFunctionsValues(ThisIntegrationMethod);

            std::vector<GeometryPointerType> geometry_pointer_vector(integration_points.size());

            for (IndexType i = 0; i < integration_points.size(); ++i)
            {
                Matrix N_i = ZeroMatrix(1, pGeometry->size());
                for (IndexType j = 0; j < pGeometry->size(); ++j)
                {
                    N_i(0, j) = r_N(i, j);
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    ThisIntegrationMethod,
                    integration_points[i],
                    N_i,
                    pGeometry->ShapeFunctionLocalGradient(i, ThisIntegrationMethod));

                geometry_pointer_vector[i] = CreateQuadraturePoint(
                    pGeometry->WorkingSpaceDimension(), pGeometry->LocalSpaceDimension(),
                    data_container, pGeometry->Points(), pGeometry.get());
            }
            return geometry_pointer_vector;

            KRATOS_CATCH("");
        }

        /// creates a quadrature point geometry on a provided location.
        static void Create(
            GeometryType& rGeometry,
            typename GeometryType::GeometriesArrayType& rResultGeometries,
            typename GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
            SizeType NumberOfShapeFunctionDerivatives)
        {
            KRATOS_ERROR_IF(NumberOfShapeFunctionDerivatives > 1)
                << "Create can only compute shape functions up to an derivative order of 1. "
                << "Demanded derivative order: " << NumberOfShapeFunctionDerivatives << std::endl;

            // Resize containers.
            if (rResultGeometries.size() != rIntegrationPoints.size())
                rResultGeometries.resize(rIntegrationPoints.size());

            auto default_method = rGeometry.GetDefaultIntegrationMethod();

            Vector N;
            Matrix DN_De;
            for (IndexType i = 0; i < rIntegrationPoints.size(); ++i)
            {
                rGeometry.ShapeFunctionsValues(N, rIntegrationPoints[i]);

                Matrix N_matrix = ZeroMatrix(1, N.size());
                if (NumberOfShapeFunctionDerivatives >= 0) {
                    for (IndexType j = 0; j < N.size(); ++j)
                    {
                        N_matrix(0, j) = N[j];
                    }
                }

                /// Get Shape Function Derivatives DN_De, ...
                if (NumberOfShapeFunctionDerivatives > 0) {
                    rGeometry.ShapeFunctionsLocalGradients(DN_De, rIntegrationPoints[i]);
                }

                GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                    default_method, rIntegrationPoints[i],
                    N_matrix, DN_De);

                rResultGeometries(i) = CreateQuadraturePointsUtility<TPointType>::CreateQuadraturePoint(
                    rGeometry.WorkingSpaceDimension(), rGeometry.LocalSpaceDimension(),
                    data_container, rGeometry);
            }
        }

        ///@}
        ///@name Update functions
        ///@{

        /* @brief This function updates the location of the respective
        *         QuadraturePointGeometry and resets the point vector and the parent.
         */
        static void UpdateFromLocalCoordinates(
            typename GeometryType::Pointer pGeometry,
            const array_1d<double, 3>& rLocalCoordinates,
            const double rIntegrationWeight,
            GeometryType& rParentGeometry)
        {
            pGeometry->SetGeometryParent(&rParentGeometry);
            pGeometry->Points() = rParentGeometry.Points();

            IntegrationPoint<3> int_p(rLocalCoordinates, rIntegrationWeight);

            Vector N;
            pGeometry->ShapeFunctionsValues(N, rLocalCoordinates);
            Matrix N_matrix(1, N.size());
            for (IndexType i = 0; i < N.size(); ++i) {
                N_matrix(0, i) = N[i];
            }

            Matrix DN_De;
            pGeometry->ShapeFunctionsLocalGradients(DN_De, rLocalCoordinates);

            GeometryShapeFunctionContainer<GeometryData::IntegrationMethod> data_container(
                pGeometry->GetDefaultIntegrationMethod(), int_p, N_matrix, DN_De);

            pGeometry->SetGeometryShapeFunctionContainer(data_container);
        }

    };
    ///@} // Kratos Classes
} // namespace Kratos.

#endif // KRATOS_QUADRATURE_POINTS_UTILITY_H_INCLUDED defined
