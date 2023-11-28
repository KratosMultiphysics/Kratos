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
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"

namespace Kratos
{
/**
 * @class GeometryUtils
 * @ingroup KratosCore
 * @brief This function provides basic routines for working with simplicial meshes
 * @details It is faster than using Geometry as it is more specialized
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) GeometryUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// The size type definition
    using SizeType = std::size_t;

    /// The index type definition
    using IndexType = std::size_t;

    /// Definition of the geometry
    using GeometryType = Geometry<Node>;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function returns a string equivalent for the geometry type
     * @param TypeOfGeometry The geometry type
     */
    static std::string GetGeometryName(const GeometryData::KratosGeometryType TypeOfGeometry);

    /**
     * @brief Returns the local coordinates of a given arbitrary point for a given linear tetrahedra
     * @details Based on https://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf. Section 9.1.6
     * @param rGeometry The geometry to be considered
     * @param rResult The vector containing the local coordinates of the point
     * @param rPoint The point in global coordinates
     * @return The vector containing the local coordinates of the point
     */
    template<class TGeometryType>
    static inline typename TGeometryType::CoordinatesArrayType& PointLocalCoordinatesPlanarFaceTetrahedra(
        const TGeometryType& rGeometry,
        typename TGeometryType::CoordinatesArrayType& rResult,
        const typename TGeometryType::CoordinatesArrayType& rPoint
        )
    {
        // Debug check that it is at least a tetrahedra
        KRATOS_DEBUG_ERROR_IF_NOT(rGeometry.GetGeometryFamily() == GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) << "Geometry should be a tetrahedra in order to use PointLocalCoordinatesPlanarFaceTetrahedra" << std::endl;

        // Compute RHS
        array_1d<double,4> X;
        X[0] = 1.0;
        X[1] = rPoint[0];
        X[2] = rPoint[1];
        X[3] = rPoint[2];

        // Auxiliary coordinates
        const auto& r_coordinates_0 = rGeometry[0].Coordinates();
        const auto& r_coordinates_1 = rGeometry[1].Coordinates();
        const auto& r_coordinates_2 = rGeometry[2].Coordinates();
        const auto& r_coordinates_3 = rGeometry[3].Coordinates();
        const double x1 = r_coordinates_0[0];
        const double y1 = r_coordinates_0[1];
        const double z1 = r_coordinates_0[2];
        const double x2 = r_coordinates_1[0];
        const double y2 = r_coordinates_1[1];
        const double z2 = r_coordinates_1[2];
        const double x3 = r_coordinates_2[0];
        const double y3 = r_coordinates_2[1];
        const double z3 = r_coordinates_2[2];
        const double x4 = r_coordinates_3[0];
        const double y4 = r_coordinates_3[1];
        const double z4 = r_coordinates_3[2];

        // Auxiliary diff
        const double x12 = x1 - x2;
        const double x13 = x1 - x3;
        const double x14 = x1 - x4;
        const double x21 = x2 - x1;
        const double x24 = x2 - x4;
        const double x31 = x3 - x1;
        const double x32 = x3 - x2;
        const double x34 = x3 - x4;
        const double x42 = x4 - x2;
        const double x43 = x4 - x3;
        const double y12 = y1 - y2;
        const double y13 = y1 - y3;
        const double y14 = y1 - y4;
        const double y21 = y2 - y1;
        const double y24 = y2 - y4;
        const double y31 = y3 - y1;
        const double y32 = y3 - y2;
        const double y34 = y3 - y4;
        const double y42 = y4 - y2;
        const double y43 = y4 - y3;
        const double z12 = z1 - z2;
        const double z13 = z1 - z3;
        const double z14 = z1 - z4;
        const double z21 = z2 - z1;
        const double z24 = z2 - z4;
        const double z31 = z3 - z1;
        const double z32 = z3 - z2;
        const double z34 = z3 - z4;
        const double z42 = z4 - z2;
        const double z43 = z4 - z3;

        // Compute LHS
        BoundedMatrix<double, 4,4> invJ;
        const double aux_volume = 1.0/(6.0*rGeometry.Volume());
        invJ(0,0) = aux_volume * (x2*(y3*z4-y4*z3)+x3*(y4*z2-y2*z4)+x4*(y2*z3-y3*z2));
        invJ(1,0) = aux_volume * (x1*(y4*z3-y3*z4)+x3*(y1*z4-y4*z1)+x4*(y3*z1-y1*z3));
        invJ(2,0) = aux_volume * (x1*(y2*z4-y4*z2)+x2*(y4*z1-y1*z4)+x4*(y1*z2-y2*z1));
        invJ(3,0) = aux_volume * (x1*(y3*z2-y2*z3)+x2*(y1*z3-y3*z1)+x3*(y2*z1-y1*z2));
        invJ(0,1) = aux_volume * (y42*z32 - y32*z42);
        invJ(1,1) = aux_volume * (y31*z43 - y34*z13);
        invJ(2,1) = aux_volume * (y24*z14 - y14*z24);
        invJ(3,1) = aux_volume * (y13*z21 - y12*z31);
        invJ(0,2) = aux_volume * (x32*z42 - x42*z32);
        invJ(1,2) = aux_volume * (x43*z31 - x13*z34);
        invJ(2,2) = aux_volume * (x14*z24 - x24*z14);
        invJ(3,2) = aux_volume * (x21*z13 - x31*z12);
        invJ(0,3) = aux_volume * (x42*y32 - x32*y42);
        invJ(1,3) = aux_volume * (x31*y43 - x34*y13);
        invJ(2,3) = aux_volume * (x24*y14 - x14*y24);
        invJ(3,3) = aux_volume * (x13*y21 - x12*y31);

        // Compute result
        const array_1d<double,4> result = prod(invJ, X);

        // Resize if needed
        if (rResult.size() != 3) {
            rResult.resize(3,false);
        }

        // Copy result
        rResult[0] = result[1];
        rResult[1] = result[2];
        rResult[2] = result[3];

        return rResult;
    }

    /**
     * @brief This function is designed to compute the shape function derivatives, shape functions and volume in 3D
     * @param rGeometry it is the array of nodes. It is expected to be a tetrahedra
     * @param rDN_DX a stack matrix of size 4*3 to store the shape function's derivatives
     * @param rN an array_1d to store the shape functions at the barycenter
     * @param rVolume the volume of the element
     */
    static inline void CalculateGeometryData(
        const GeometryType& rGeometry,
        BoundedMatrix<double,4,3>& rDN_DX,
        array_1d<double,4>& rN,
        double& rVolume
        )
    {
        const double x10 = rGeometry[1].X() - rGeometry[0].X();
        const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
        const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

        const double x20 = rGeometry[2].X() - rGeometry[0].X();
        const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
        const double z20 = rGeometry[2].Z() - rGeometry[0].Z();

        const double x30 = rGeometry[3].X() - rGeometry[0].X();
        const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
        const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

        const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        rDN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        rDN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        rDN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        rDN_DX(1,0) = y20 * z30 - y30 * z20;
        rDN_DX(1,1) = z20 * x30 - x20 * z30;
        rDN_DX(1,2) = x20 * y30 - y20 * x30;
        rDN_DX(2,0) = -y10 * z30 + z10 * y30;
        rDN_DX(2,1) = x10 * z30 - z10 * x30;
        rDN_DX(2,2) = -x10 * y30 + y10 * x30;
        rDN_DX(3,0) = y10 * z20 - z10 * y20;
        rDN_DX(3,1) = -x10 * z20 + z10 * x20;
        rDN_DX(3,2) = x10 * y20 - y10 * x20;

        rDN_DX /= detJ;

        rN[0] = 0.25;
        rN[1] = 0.25;
        rN[2] = 0.25;
        rN[3] = 0.25;

        rVolume = detJ*0.1666666666666666666667;
    }

    /**
     * @brief This function computes the element's volume (with sign)
     * @param rGeometry it is the array of nodes. It expects a tetrahedra
     * @deprecated This method can be replaced by geometry function without loosing performance
     */
    KRATOS_DEPRECATED_MESSAGE("Please use the Volume() method from the geometry")
    static inline double CalculateVolume3D(const GeometryType& rGeometry)
    {
        const double x10 = rGeometry[1].X() - rGeometry[0].X();
        const double y10 = rGeometry[1].Y() - rGeometry[0].Y();
        const double z10 = rGeometry[1].Z() - rGeometry[0].Z();

        const double x20 = rGeometry[2].X() - rGeometry[0].X();
        const double y20 = rGeometry[2].Y() - rGeometry[0].Y();
        const double z20 = rGeometry[2].Z() - rGeometry[0].Z();

        const double x30 = rGeometry[3].X() - rGeometry[0].X();
        const double y30 = rGeometry[3].Y() - rGeometry[0].Y();
        const double z30 = rGeometry[3].Z() - rGeometry[0].Z();

        const double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return detJ*0.1666666666666666666667;
    }

    /**
     * @brief This function is designed to compute the shape function derivatives, shape functions and area of a triangle
     * @param rGeometry it is the array of nodes. It is expected to be a triangle
     * @param DN_DX a stack matrix of size 3*2 to store the shape function's derivatives
     * @param N an array_1d to store the shape functions at the barycenter
     * @param Area the volume of the element
     */
    static inline void CalculateGeometryData(
        const GeometryType& rGeometry,
        BoundedMatrix<double,3,2>& DN_DX,
        array_1d<double,3>& N,
        double& rArea
        )
    {
        const double x10 = rGeometry[1].X() - rGeometry[0].X();
        const double y10 = rGeometry[1].Y() - rGeometry[0].Y();

        const double x20 = rGeometry[2].X() - rGeometry[0].X();
        const double y20 = rGeometry[2].Y() - rGeometry[0].Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|    |x1-x0   x2-x0|
        //J=|               |=   |             |
        //  |dy/dxi  dy/deta|    |y1-y0   y2-y0|

        double detJ = x10 * y20-y10 * x20;

        DN_DX(0,0) = -y20 + y10;
        DN_DX(0,1) = x20 - x10;
        DN_DX(1,0) =  y20       ;
        DN_DX(1,1) = -x20     ;
        DN_DX(2,0) = -y10       ;
        DN_DX(2,1) = x10       ;

        DN_DX /= detJ;
        N[0] = 0.333333333333333;
        N[1] = 0.333333333333333;
        N[2] = 0.333333333333333;

        rArea = 0.5*detJ;
    }

    /**
     * @brief This function computes the element's volume (with sign)
     * @param rGeometry it is the array of nodes. It expects a triangle
     * @deprecated This method can be replaced by geometry function without loosing performance
     */
    KRATOS_DEPRECATED_MESSAGE("Please use the Area() method from the geometry")
    static inline double CalculateVolume2D(const GeometryType& rGeometry)
    {
        double x10 = rGeometry[1].X() - rGeometry[0].X();
        double y10 = rGeometry[1].Y() - rGeometry[0].Y();

        double x20 = rGeometry[2].X() - rGeometry[0].X();
        double y20 = rGeometry[2].Y() - rGeometry[0].Y();

        double detJ = x10 * y20-y10 * x20;
        return 0.5*detJ;
    }

    /**
     * @brief This function compute the maximum and minimum edge lenghts
     * @details It assumes that it is a triangle
     * @param rGeometry it is the array of nodes. It expects a triangle
     * @param hmin The minimum length
     * @param hmax The maximum length
     */
    static inline void SideLenghts2D(
        const GeometryType& rGeometry,
        double& hmin,
        double& hmax
        )
    {
        const double x10 = rGeometry[1].X() - rGeometry[0].X();
        const double y10 = rGeometry[1].Y() - rGeometry[0].Y();

        const double x20 = rGeometry[2].X() - rGeometry[0].X();
        const double y20 = rGeometry[2].Y() - rGeometry[0].Y();

        double l = std::pow(x20, 2) + std::pow(y20, 2);
        hmax = l;
        hmin = l;

        if(l>hmax) hmax = l;
        else if(l<hmin) hmin = l;

        l = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10);
        if(l>hmax) hmax = l;
        else if(l<hmin) hmin = l;

        hmax = std::sqrt(hmax);
        hmin = std::sqrt(hmin);
    }

    /**
     * @brief This function is designed to compute the shape function derivatives, shape functions and length
     * @param rGeometry it is the array of nodes. It is expected to be a line
     * @param rDN_DX a stack matrix of size 3*2 to store the shape function's derivatives
     * @param rN an array_1d to store the shape functions at the barycenter
     * @param rLength the volume of the element
     */
    static inline void CalculateGeometryData(
        const GeometryType& rGeometry,
        BoundedMatrix<double,2,1>& rDN_DX,
        array_1d<double,2>& rN,
        double& rLength
        )
    {
        const double lx = rGeometry[0].X() - rGeometry[1].X();
        const double ly = rGeometry[0].Y() - rGeometry[1].Y();
        const double detJ = 0.5 * std::sqrt(std::pow(lx, 2) + std::pow(ly, 2));

        rDN_DX(0,0) = -0.5;
        rDN_DX(1,0) = 0.5;
        rDN_DX /= detJ;

        rN[0] = 0.5;
        rN[1] = 0.5;

        rLength = 2.0 * detJ;
    }

    /**
     * @brief Calculate the exact distances to the interface TRIANGLE defined by a set of initial distances.
     * @param rGeometry The tetrahedra itself. Note: If the geometry is not a tetrahedra the result is undefined and may cause memory error.
     * @param rDistances The distances which define the isosurface as input and the same argument is used to give the calculated exact distance
     */
    template<std::size_t TSize>
    static void CalculateTetrahedraDistances(
        const GeometryType& rGeometry, array_1d<double, TSize>& rDistances)
    {
        // Calculating the intersection points
        array_1d<Point, 4> intersection_points;
        int number_of_intersection_points = CalculateTetrahedraIntersectionPoints(rGeometry, rDistances, intersection_points);

        if(number_of_intersection_points == 0) {
            KRATOS_WARNING("CalculateTetrahedraDistances") << "WARNING:: The intersection with interface hasn't found!" << std::endl << "The distances are: " << rDistances << std::endl;
        } else if(number_of_intersection_points == 1) {
            // There is one point with zero distance. The distance of the nodes are their distance to this point
            array_1d<double,3> temp;
            // Loop over nodes to calculate their distance to the zero distance node.
            for(unsigned int i_node = 0; i_node < rGeometry.size() ; ++i_node) {
                noalias(temp) = intersection_points[0] - rGeometry[i_node];
                rDistances[i_node] = norm_2(temp);
            }
        } else if(number_of_intersection_points == 2) {
            // Loop over nodes to calculate their distance to the zero distance line.
            for(unsigned int i_node = 0; i_node < rGeometry.size() ; ++i_node) {
                rDistances[i_node] = PointDistanceToLineSegment3D(intersection_points[0], intersection_points[1], rGeometry[i_node]);
            }
        } else if(number_of_intersection_points == 3) {
            // Loop over nodes to calculate their distance to the zero distance triangle.
            for(unsigned int i_node = 0; i_node < rGeometry.size() ; ++i_node) {
                rDistances[i_node] = PointDistanceToTriangle3D(intersection_points[0], intersection_points[1], intersection_points[2], rGeometry[i_node]);
//                 rDistances[i_node] = std::abs(rGeometry[i_node].Z()); // TODO: To be removed. Pooyan.
            }

        } else if(number_of_intersection_points == 4) {
            // Loop over nodes to calculate their distance to the each zero distance triangle.
            for(unsigned int i_node = 0; i_node < rGeometry.size() ; ++i_node) {
                // Here I'm taking in account the order of edges where I'm looking for intersection
                double d1 = PointDistanceToTriangle3D(intersection_points[0], intersection_points[1], intersection_points[3], rGeometry[i_node]);
                double d2 = PointDistanceToTriangle3D(intersection_points[0], intersection_points[3], intersection_points[2], rGeometry[i_node]);

                rDistances[i_node] = (d1 > d2) ? d2 : d1;
            }
        }
    }

    /**
     * @brief Calculate the exact distances to the interface SEGMENT defined by a set of initial distances.
     * @param rGeometry The Triangle itself. Note: If the geometry is not a triangle the result is undefined and may cause memory error.
     * @param rDistances The distances which define the isosurface as input and the same argument is used to give the calculated exact distance
     */
    template<std::size_t TSize>
    static void CalculateTriangleDistances(
        const GeometryType& rGeometry,
        array_1d<double, TSize>& rDistances
        )
    {
        // Calculating the intersection points
        array_1d<Point, 4> intersection_points;
        int number_of_intersection_points = CalculateTetrahedraIntersectionPoints(rGeometry, rDistances, intersection_points);

        if(number_of_intersection_points == 0) {
            KRATOS_WARNING("CalculateTriangleDistances") << "WARNING:: The intersection with interface hasn't found!" << std::endl << "The distances are: " << rDistances << std::endl;
        } else if(number_of_intersection_points == 1) {   // There is one point with zero distance. The distance of the nodes are their distance to this point
            array_1d<double,3> temp;
            // Loop over nodes to calculate their distance to the zero distance node.
            for(unsigned int i_node = 0; i_node < rGeometry.size() ; ++i_node) {
                noalias(temp) = intersection_points[0] - rGeometry[i_node];
                rDistances[i_node] = norm_2(temp);
            }
        } else if(number_of_intersection_points == 2) {
            // loop over nodes to calculate their distance to the zero distance line.
            for(unsigned int i_node = 0; i_node < rGeometry.size() ; ++i_node) {
                rDistances[i_node] = PointDistanceToLineSegment3D(intersection_points[0], intersection_points[1], rGeometry[i_node]);
            }
        } else {
            KRATOS_WARNING("CalculateTriangleDistances") << "WARNING:: This is a triangle with more than two intersections!" << std::endl << "Too many intersections: " << number_of_intersection_points << std::endl << "The distances are: " << rDistances << std::endl;
        }
    }

    /**
     * @brief Calculate the exact distances to the plane interface defined by a set of initial distances.
     * @details Same argument is used to give the calculated exact distances back
     * @param rThisGeometry Geometry can be either a triangle or a tetrahedra
     * @param rDistances The distances which define the isosurface as input.
     */
    template<std::size_t TSize>
    static void CalculateExactDistancesToPlane(
        const GeometryType& rThisGeometry,
        array_1d<double, TSize>& rDistances
        )
    {
        array_1d<Point, TSize> intersection_points;
        int number_of_intersection_points = CalculateTetrahedraIntersectionPoints(rThisGeometry, rDistances, intersection_points);

        if(number_of_intersection_points == 0) {
            KRATOS_WARNING("GeometryUtilities") << "Warning: The intersection with interface hasn't found! The distances are" << rDistances << std::endl;
        } else {
            BoundedMatrix<double,TSize,TSize-1> DN_DX;
            array_1d<double, TSize> N;
            double volume;
            GeometryUtils::CalculateGeometryData(rThisGeometry, DN_DX, N, volume);
            array_1d<double, TSize-1> distance_gradient = prod(trans(DN_DX), rDistances);
            double distance_gradient_norm = norm_2(distance_gradient);
            if (distance_gradient_norm < 1e-15) distance_gradient_norm = 1e-15; //avoid division by zero
            distance_gradient /= distance_gradient_norm;
            // We use the first intersection point as reference
            const auto &ref_point = intersection_points[0].Coordinates();
            for (unsigned int i = 0; i < TSize; i++) {
                double d = 0.0;
                const auto &i_coords = rThisGeometry[i].Coordinates();
                for (unsigned int Dim = 0; Dim < TSize - 1; Dim++) {
                    d += (i_coords[Dim] - ref_point[Dim]) * distance_gradient[Dim];
                }
                d = std::abs(d);
                rDistances[i] = std::min(std::abs(rDistances[i]), d);
            }
        }
    }

    /**
     * @brief This function calculates the coordinates of the intersecion points between edges of tetrahedra and a isosurface given by the distances in its corners
     * @param rGeometry The tetrahedra itself. Note: If the geometry is not a tetrahedra the result is undefined and may cause memory error.
     * @param rDistances The distances of the 4 nodes of the tetrahedra to the iso surface.
     * @param rIntersectionPoints The result intersection points
     * @return Number of intersection points.
     */
    template<std::size_t TSize1, std::size_t TSize2>
    static int CalculateTetrahedraIntersectionPoints(
        const GeometryType& rGeometry,
        array_1d<double, TSize1>& rDistances,
        array_1d<Point, TSize2>& rIntersectionPoints
        )
    {
        const double epsilon = 1e-15; //1.00e-9;

        int number_of_intersection_points = 0;
        for(unsigned int i = 0 ; i < TSize1 ; i++) {
            if(std::abs(rDistances[i]) < epsilon) {
                noalias(rIntersectionPoints[number_of_intersection_points].Coordinates()) = rGeometry[i].Coordinates();

                number_of_intersection_points++;
                continue;
            }
            for(unsigned int j = i + 1 ; j < TSize1 ; j++) {
                if(std::abs(rDistances[j]) < epsilon)
                    continue; // we will add it to the intersections by the i index to be unique

                if(rDistances[i] * rDistances[j] < 0.00) { // The interface passes through the edge
                    const double delta_d = std::abs(rDistances[i]) + std::abs(rDistances[j]);  // we know that both distances are greater than epsilon.

                    const double di = std::abs(rDistances[i]) / delta_d;
                    const double dj = std::abs(rDistances[j]) / delta_d;

                    noalias(rIntersectionPoints[number_of_intersection_points].Coordinates()) = dj * rGeometry[i].Coordinates();
                    noalias(rIntersectionPoints[number_of_intersection_points].Coordinates()) += di * rGeometry[j].Coordinates();

                    number_of_intersection_points++;
                }
            }
        }

        return number_of_intersection_points;
    }

    /**
     * @brief This function calculates the distance of a 3D point to a 3D line segment
     * @param rLinePoint1 First point of the line segment
     * @param rLinePoint2 End point of the line segment
     * @param rToPoint The point which distance is required
     * @return The distance between the point and the line
     */
    static double PointDistanceToLineSegment3D(
        const Point& rLinePoint1,
        const Point& rLinePoint2,
        const Point& rToPoint
        );

    /**
     * @brief This function calculates the distance of a 3D point to a 3D triangle
     * @details The implementation is done using following reference:
     *          http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
     * @param rTrianglePoint1 First point of triangle
     * @param rTrianglePoint2 Second point of triangle
     * @param rTrianglePoint3 Third point of triangle
     * @param rPoint The point which distance is required
     * @return The distance between the point and the triangle
     */
    static double PointDistanceToTriangle3D(
        const Point& rTrianglePoint1,
        const Point& rTrianglePoint2,
        const Point& rTrianglePoint3,
        const Point& rPoint
        );

    /**
     * @brief This function calculates the distance of a 3D point to a 3D quadratic triangle
     * @details The implementation is done by decomposing the quadratic triangle into 3 triangles and calling PointDistanceToTriangle3D
     * @param rTrianglePoint1 First point of triangle
     * @param rTrianglePoint2 Second point of triangle
     * @param rTrianglePoint3 Third point of triangle
     * @param rTrianglePoint4 Fourth point of triangle
     * @param rTrianglePoint5 Fifth point of triangle
     * @param rTrianglePoint6 Sixth point of triangle
     * @param rPoint The point which distance is required
     * @return The distance between the point and the triangle
     */
    static double PointDistanceToTriangle3D(
        const Point& rTrianglePoint1,
        const Point& rTrianglePoint2,
        const Point& rTrianglePoint3,
        const Point& rTrianglePoint4,
        const Point& rTrianglePoint5,
        const Point& rTrianglePoint6,
        const Point& rPoint
        );

    /**
     * @brief This function calculates the distance of a 3D point to a 3D quadrilateral
     * @details The implementation is done by decomposing the quadrilateral into 2 triangles and calling PointDistanceToTriangle3D
     * @param rQuadrilateralPoint1 First point of quadrilateral
     * @param rQuadrilateralPoint2 Second point of quadrilateral
     * @param rQuadrilateralPoint3 Third point of quadrilateral
     * @param rQuadrilateralPoint4 Third point of quadrilateral
     * @param rPoint The point which distance is required
     * @return The distance between the point and the quadrilateral
     */
    static double PointDistanceToQuadrilateral3D(
        const Point& rQuadrilateralPoint1,
        const Point& rQuadrilateralPoint2,
        const Point& rQuadrilateralPoint3,
        const Point& rQuadrilateralPoint4,
        const Point& rPoint
        );

    /**
     * @brief Calculate the gradients of shape functions.
     * @param rDN_De local gradient of shape functions.
     * @param rInvJ inverse of the element Jacobian.
     * @param rDN_DX gradient of shape functions.
     */
    template<class TMatrix1, class TMatrix2, class TMatrix3>
    static void ShapeFunctionsGradients(
        TMatrix1 const& rDN_De,
        TMatrix2 const& rInvJ,
        TMatrix3& rDN_DX
        )
    {
        if (rDN_DX.size1() != rDN_De.size1() || rDN_DX.size2() != rInvJ.size2())
            rDN_DX.resize(rDN_De.size1(), rInvJ.size2(), false);

        noalias(rDN_DX) = prod(rDN_De, rInvJ);
    }

    /**
     * @brief Calculate the deformation gradient.
     * @details See, e.g., P. Wriggers, Nonlinear Finite Element Methods, Springer, 2008.
     * @param rJ element Jacobian.
     * @param rInvJ0 inverse of the element Jacobian of the initial configuration.
     * @param rF deformation gradient.
     */
    template<class TMatrix1, class TMatrix2, class TMatrix3>
    static void DeformationGradient(
        TMatrix1 const& rJ,
        TMatrix2 const& rInvJ0,
        TMatrix3& rF
        )
    {
        if (rF.size1() != rJ.size1() || rF.size2() != rInvJ0.size2())
            rF.resize(rJ.size1(), rInvJ0.size2(), false);

        noalias(rF) = prod(rJ, rInvJ0);
    }

    /**
     * @brief Calculate the Jacobian on the initial configuration.
     * @param rGeom element geometry.
     * @param rCoords local coordinates of the current integration point.
     * @param rJ0 Jacobian on the initial configuration.
     */
    static void JacobianOnInitialConfiguration(
        GeometryType const& rGeom,
        GeometryType::CoordinatesArrayType const& rCoords,
        Matrix& rJ0
        )
    {
        Matrix delta_position(rGeom.PointsNumber(), rGeom.WorkingSpaceDimension());
        for (std::size_t i = 0; i < rGeom.PointsNumber(); ++i)
            for (std::size_t j = 0; j < rGeom.WorkingSpaceDimension(); ++j)
                delta_position(i, j) = rGeom[i].Coordinates()[j] -
                                       rGeom[i].GetInitialPosition().Coordinates()[j];
        rGeom.Jacobian(rJ0, rCoords, delta_position);
    }

    /**
     * @brief Calculate the Jacobian on the initial configuration.
     * @param rGeom element geometry.
     * @param rCoords local coordinates of the current integration point.
     * @param rJ0 Jacobian on the initial configuration.
     */
    template<class TMatrix>
    static void DirectJacobianOnCurrentConfiguration(
        GeometryType const& rGeometry,
        GeometryType::CoordinatesArrayType const& rCoords,
        TMatrix& rJ
        )
    {
        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType points_number = rGeometry.PointsNumber();

        Matrix shape_functions_gradients(points_number, local_space_dimension);
        rGeometry.ShapeFunctionsLocalGradients( shape_functions_gradients, rCoords );

        rJ.clear();
        for (IndexType i = 0; i < points_number; ++i ) {
            const array_1d<double, 3>& r_coordinates = rGeometry[i].Coordinates();
            for(IndexType j = 0; j< working_space_dimension; ++j) {
                const double value = r_coordinates[j];
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rJ(j,m) += value * shape_functions_gradients(i,m);
                }
            }
        }
    }

    /**
     * @brief Calculate the Jacobian on the initial configuration.
     * @param rGeom element geometry
     * @param rJ0 Jacobian on the initial configuration
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    template<class TMatrix>
    static void DirectJacobianOnInitialConfiguration(
        GeometryType const& rGeometry,
        TMatrix& rJ0,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        )
    {
        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType points_number = rGeometry.PointsNumber();

        const Matrix& rDN_De = rGeometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[PointNumber];

        rJ0.clear();
        for (IndexType i = 0; i < points_number; ++i ) {
            const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition().Coordinates();
            for(IndexType j = 0; j< working_space_dimension; ++j) {
                const double value = r_coordinates[j];
                for(IndexType m = 0; m < local_space_dimension; ++m) {
                    rJ0(j,m) += value * rDN_De(i,m);
                }
            }
        }
    }

    /**
     * @brief Evaluates variable value at gauss point
     *
     * This method evaluates variable value at gauss point given by gauss point shape function values.
     *
     * @tparam TDataType                        Data type
     * @param rOutput                           Output which holds evaluated value
     * @param rGeometry                         Geometry from which gauss point values are interpolated
     * @param rVariable                         Variable for value interpolation
     * @param rGaussPointShapeFunctionValues    Shape function values evaluated at gauss point
     * @param Step                              Step to be used in historical variable value interpolation
     */
    template <class TDataType>
    static void EvaluateHistoricalVariableValueAtGaussPoint(
        TDataType& rOutput,
        const GeometryType& rGeometry,
        const Variable<TDataType>& rVariable,
        const Vector& rGaussPointShapeFunctionValues,
        const int Step = 0);

    /**
     * @brief Evaluates gradient of scalar at gauss point
     *
     * This method evaluates gradient of a scalar variable for given shape function derivative values.
     * For 2D, it returns 3rd component as zero.
     *
     * @param rOutput                                   Output 3d variable containing gradients at gauss point
     * @param rGeometry                                 Geometry from which gauss point values are interpolated
     * @param rVariable                                 Variable for value interpolation
     * @param rGaussPointShapeFunctionDerivativeValues  Shape function derivatives evaluated at gauss point
     * @param Step                                      Step to be used in historical variable value interpolation
     */
    static void EvaluateHistoricalVariableGradientAtGaussPoint(
        array_1d<double, 3>& rOutput,
        const GeometryType& rGeometry,
        const Variable<double>& rVariable,
        const Matrix& rGaussPointShapeFunctionDerivativeValues,
        const int Step = 0);

    /**
     * @brief Evaluates gradient of a 3d vector at gauss point
     *
     * This method evaluates gradient of a vector variable for given shape function derivative values.
     * \[
     *      rOutput(i, j) = \frac{\partial u_i}{\partial x_j}
     * \]
     * Where $u_i$ is the component of 3D rVariable, and $x_j$ is the cartesian coordinate component.
     * In 2D, it returns 3rd row and 3rd column with zeros.
     *
     * @param rOutput                                   Output matrix containing gradients at gauss point
     * @param rGeometry                                 Geometry from which gauss point values are interpolated
     * @param rVariable                                 Variable for value interpolation
     * @param rGaussPointShapeFunctionDerivativeValues  Shape function derivatives evaluated at gauss point
     * @param Step                                      Step to be used in historical variable value interpolation
     */
    static void EvaluateHistoricalVariableGradientAtGaussPoint(
        BoundedMatrix<double, 3, 3>& rOutput,
        const GeometryType& rGeometry,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rGaussPointShapeFunctionDerivativeValues,
        const int Step = 0);

    /**
     * @brief This method gives the transform of shape functions second derivatives from isoparametric to natural coordinates evaluated in all integration points
     *  @param rResult the transform of the second derivative of the shape function on the integration point
     */
    static void ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
        DenseVector<DenseVector<Matrix>>& rResult,
        const GeometryType& rGeometry,
        const GeometryType::IntegrationMethod& rIntegrationMethod );


    /**
     * @brief This method gives the transform of shape functions second derivatives from isoparametric to natural coordinates evaluated on an integration point.
     * @details The method used can be found here. https://scicomp.stackexchange.com/questions/25196/implementing-higher-order-derivatives-for-finite-element
     *  @param rLocalIntegrationPointCoordinates the local coordinates of the integration point
     *  @param rResult the transform of the second derivative of the shape function on the integration point
     */
    static void ShapeFunctionsSecondDerivativesTransformOnIntegrationPoint(
        const Matrix& DN_DX,
        const GeometryType& rGeometry,
        const GeometryType::CoordinatesArrayType& rLocalIntegrationPointCoordinates,
        DenseVector<Matrix>& rResult);

    /**
     * @brief Checks if given point in global space coordinates is inside the geometry boundaries.
     * @details This function computes the local coordinates and checks then if this point lays within the boundaries after projecting the points.
     * @param rPointGlobalCoordinates the global coordinates of the external point.
     * @param rResult the local coordinates of the point.
     * @param Tolerance the tolerance to the boundary.
     * @return true if the point is inside, false otherwise
     */
    static bool ProjectedIsInside(
        const GeometryType& rGeometry,
        const GeometryType::CoordinatesArrayType& rPointGlobalCoordinates,
        GeometryType::CoordinatesArrayType& rResult,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        );

    /**
    * @brief Computes the distance between an point in global coordinates and the closest point of this geometry.
    * @param rGeometry the geometry to compute the distance to.
    * @param rPointGlobalCoordinates the point to which the closest point has to be found.
    * @param Tolerance accepted orthogonal error.
    * @return Distance to geometry.
    *         positive -> outside of to the geometry (for 2D and solids)
    *         0        -> on/ in the geometry.
    */
    template <class TGeometryType>
    static double CalculateDistanceFrom3DGeometry(
        const TGeometryType& rGeometry,
        const typename TGeometryType::CoordinatesArrayType& rPointGlobalCoordinates,
        const double Tolerance = std::numeric_limits<double>::epsilon()
        )
    {
        typename TGeometryType::CoordinatesArrayType aux_coordinates;
        if (rGeometry.IsInside(rPointGlobalCoordinates, aux_coordinates, Tolerance)) {
            return 0.0;
        }

        // Generate faces
        std::vector<double> distances(rGeometry.FacesNumber());
        unsigned int i = 0;
        for (auto& r_face : rGeometry.GenerateFaces()) {
            distances[i] = r_face.CalculateDistance(rPointGlobalCoordinates, Tolerance);
            ++i;
        }
        const auto min = std::min_element(distances.begin(), distances.end());
        return *min;
    }
};

}  // namespace Kratos.


