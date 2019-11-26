//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_GEOMETRY_UTILITIES_INCLUDED )
#define  KRATOS_GEOMETRY_UTILITIES_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"

namespace Kratos
{
/**
 * @class GeometryUtils
 * @ingroup KratosCore
 * @brief This function provides basic routines for working with simplicial meshes
 * @details It is faster than using Geometry as it is more specialized
 * @author Riccardo Rossi
 */
class GeometryUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// The size type definition
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// Definition of the node
    typedef Node<3> NodeType;

    /// Definition of the geometry
    typedef Geometry<NodeType> GeometryType;

    ///@}

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
        //J=|                |=    |              |
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
        Point const& rLinePoint1,
        Point const& rLinePoint2,
        Point const& rToPoint
        )
    {
        const double epsilon = 1e-15; //1.00e-9;

        const array_1d<double,3> v1 = rLinePoint2 - rLinePoint1;
        const array_1d<double,3> v2 = rLinePoint1 - rToPoint;
        array_1d<double,3> v3;

        const double square_distance = inner_prod(v1,v1);

        if(square_distance < epsilon) // near zero length line
            return norm_2(v2); // we return the distance to the first point of line

        const double t = - inner_prod(v1,v2) / square_distance;

        if(t < 0.0) { // it is before point 1
            // We return the distance to point 1
            noalias(v3) = rLinePoint1 - rToPoint;

            return norm_2(v3);
        }

        if(t > 1.00) { // it is after point 2
            // We return the distance to point 2
            noalias(v3) = rLinePoint2 - rToPoint;

            return norm_2(v3);
        }

        // The projection point is between point 1 and 2 of the line segment
        noalias(v3) = rLinePoint1 * (1.0 - t) + rLinePoint2 * t;

        return norm_2(v3 - rToPoint);

    }

    /**
     * @brief This function calculates the distance of a 3D point to a 3D triangle
     * @details The implementation is done using following reference:
     *          http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
     * @param TrianglePoint1 First point of triangle
     * @param TrianglePoint2 Second point of triangle
     * @param TrianglePoint3 Third point of triangle
     * @param ToPoint The point which distance is required
     * @return The distance between the point and the triangle
     */
    static double PointDistanceToTriangle3D(
        Point const& TrianglePoint1,
        Point const& TrianglePoint2,
        Point const& TrianglePoint3,
        Point const& ToPoint
    )
    {
        const array_1d<double, 3> e0 = TrianglePoint2 - TrianglePoint1;
        const array_1d<double, 3> e1 = TrianglePoint3 - TrianglePoint1;
        const array_1d<double, 3> dd = TrianglePoint1 - ToPoint;

        const double a = inner_prod(e0, e0);
        const double b = inner_prod(e0, e1);
        const double c = inner_prod(e1, e1);
        const double d = inner_prod(e0, dd);
        const double e = inner_prod(e1, dd);
        const double f = inner_prod(dd, dd);

        const double det = a*c-b*b;
        double s = b*e-c*d;
        double t = b*d-a*e;

        double square_distance = 0.00;

        if ( s + t <= det ) {
            if ( s < 0.0 ) {
                if ( t < 0.0 ) { // region 4
                    if (d < 0) {
                        t = 0;
                        if (-d >= a) {
                            s = 1;
                            square_distance = a + 2*d + f;
                        } else {
                            s = -d/a;
                            square_distance = d*s + f;
                        }
                    } else {
                        s = 0;
                        if (e >= 0) {
                            t = 0;
                            square_distance = f;
                        } else {
                            if (-e >= c) {
                                t = 1;
                                square_distance = c + 2*e + f;
                            } else {
                                t = -e/c;
                                square_distance = e*t + f;
                            }
                        }
                    }
                } else { // region 3
                    s = 0.0;
                    if(e >= 0.0) {
                        t = 0.0;
                        square_distance = f;
                    } else {
                        if (-e >= c) {
                            t = 1.00;
                            square_distance = c + 2*e +f;
                        } else {
                            t = -e/c;
                            square_distance = e*t + f;
                        }
                    }

                }
            } else if ( t < 0.00 ) { // region 5
                t = 0;
                if (d >= 0) {
                    s = 0;
                    square_distance = f;
                } else {
                    if (-d >= a) {
                        s = 1;
                        square_distance = a + 2.0 * d + f;
                    } else {
                        s = -d / a;
                        square_distance = d * s + f;
                    }
                }
            } else { // region 0
                double inv_det = 1.0 / det;
                s *= inv_det;
                t *= inv_det;
                square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
            }
        } else {
            if ( s < 0.00 ) {
                // Region 2
                const double temp0 = b + d;
                const double temp1 = c + e;
                if (temp1 > temp0)  { // Minimum on edge s+t=1
                    const double numer = temp1 - temp0;
                    const double denom = a - 2*b + c;
                    if(numer >= denom) {
                        s = 1.0;
                        t = 0.0;
                        square_distance = a + 2*d + f;
                    } else {
                        s = numer/denom;
                        t = 1.0 - s;
                        square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                } else { // Minimum on edge s=0
                    s = 0.0;
                    if(temp1 <= 0.0) {
                        t = 1;
                        square_distance = c + 2*e + f;
                    } else {
                        if(e >= 0.0) {
                            t = 0.0;
                            square_distance = f;
                        } else {
                            t = -e/c;
                            square_distance = e*t + f;
                        }
                    }
                }
            } else if ( t < 0.0 ) {
                // Region 6
                double temp0 = b + e;
                double temp1 = a + d;
                if (temp1 > temp0) {
                    double numer = temp1 - temp0;
                    double denom = a - 2*b + c;
                    if(numer >= denom) {
                        s = 0.0;
                        t = 1.0;
                        square_distance = c + 2*e + f;
                    } else {
                        t = numer/denom;
                        s = 1.0 - t;
                        square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                } else {
                    t = 0.0;
                    if(temp1 <= 0.0) {
                        s = 1;
                        square_distance = a + 2*d + f;
                    } else {
                        if(d >= 0.0) {
                            s = 0.0;
                            square_distance = f;
                        } else {
                            s = -d/a;
                            square_distance = d*s + f;
                        }
                    }
                }
            } else {
                // Region 1
                double numer = c + e - b - d;

                if (numer <= 0.0) {
                    s = 0.0;
                    t = 1.0;
                    square_distance = c + 2.0 * e + f;
                } else {
                    double denom = a - 2.0 * b + c;
                    if (numer >= denom) {
                        s = 1.0;
                        t = 0.0;
                        square_distance = a + 2.0 * d + f;
                    } else {
                        s = numer / denom;
                        t = 1.0 - s;
                        square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
                    }
                }
            }
        }

        if(square_distance < 0.0)
            return 0.0; // avoiding -0 case!!

        return std::sqrt(square_distance);
    }

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
    #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        KRATOS_WARNING_IF("ShapeFunctionsGradients", rDN_DX.size1() != rDN_De.size1() || rDN_DX.size2() != rInvJ.size2()) << "ShapeFunctionsGradients has detected an incorrect size of your DN_DX matrix. Please resize before compute" << std::endl;
    #else
        if (rDN_DX.size1() != rDN_De.size1() || rDN_DX.size2() != rInvJ.size2())
            rDN_DX.resize(rDN_De.size1(), rInvJ.size2(), false);
    #endif // KRATOS_USE_AMATRIX
        
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
    #ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        KRATOS_WARNING_IF("DeformationGradient", rF.size1() != rJ.size1() || rF.size2() != rInvJ0.size2()) << "DeformationGradient has detected an incorrect size of your F matrix. Please resize before compute" << std::endl;
    #else
        if (rF.size1() != rJ.size1() || rF.size2() != rInvJ0.size2())
            rF.resize(rJ.size1(), rInvJ0.size2(), false);
    #endif // KRATOS_USE_AMATRIX
        
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
};

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined


