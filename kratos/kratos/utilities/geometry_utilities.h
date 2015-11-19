/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
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
///this function provides basic routines for working with simplicial meshes.
///It is faster than using Geometry as it is more specialized
class GeometryUtils
{
public:

    /**this function is designed to compute the shape function derivatives, shape functions and volume in 3D
     * @param geom it is the array of nodes. It is expected to be a tetrahedra
     * @param a stack matrix of size 4*3 to store the shape function's derivatives
     * @param an array_1d to store the shape functions at the barycenter
     * @param the volume of the element
     */
    static inline void CalculateGeometryData(
        Element::GeometryType& geom,
        boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX,
        array_1d<double,4>& N,
        double& Volume)
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();
        double z10 = geom[1].Z() - geom[0].Z();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();
        double z20 = geom[2].Z() - geom[0].Z();

        double x30 = geom[3].X() - geom[0].X();
        double y30 = geom[3].Y() - geom[0].Y();
        double z30 = geom[3].Z() - geom[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;

        N[0] = 0.25;
        N[1] = 0.25;
        N[2] = 0.25;
        N[3] = 0.25;

        Volume = detJ*0.1666666666666666666667;
    }

    /**this function computes the element's volume (with sign)
     * @param geom it is the array of nodes. It expects a tetrahedra
     */
    static inline double CalculateVolume3D(
        Element::GeometryType& geom)
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();
        double z10 = geom[1].Z() - geom[0].Z();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();
        double z20 = geom[2].Z() - geom[0].Z();

        double x30 = geom[3].X() - geom[0].X();
        double y30 = geom[3].Y() - geom[0].Y();
        double z30 = geom[3].Z() - geom[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return  detJ*0.1666666666666666666667;
    }

    //********************************************************************************
    //********************************************************************************
    /**this function is designed to compute the shape function derivatives, shape functions and volume in 3D
     * @param geom it is the array of nodes. It is expected to be a triangle
     * @param a stack matrix of size 3*2 to store the shape function's derivatives
     * @param an array_1d to store the shape functions at the barycenter
     * @param the volume of the element
     */
    static inline void CalculateGeometryData(
        Element::GeometryType& geom,
        boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,
        array_1d<double,3>& N,
        double& Area)
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();

        //Jacobian is calculated:
        //  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
        //J=|				|=	|			  |
        //  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


        double detJ = x10 * y20-y10 * x20;

        DN_DX(0,0) = -y20 + y10;
        DN_DX(0,1) = x20 - x10;
        DN_DX(1,0) =  y20	   ;
        DN_DX(1,1) = -x20     ;
        DN_DX(2,0) = -y10	   ;
        DN_DX(2,1) = x10	   ;

        DN_DX /= detJ;
        N[0] = 0.333333333333333;
        N[1] = 0.333333333333333;
        N[2] = 0.333333333333333;

        Area = 0.5*detJ;
    }


    //********************************************************************************
    //********************************************************************************
    /**this function computes the element's volume (with sign)
     * @param geom it is the array of nodes. It expects a triangle
     */
    static inline double CalculateVolume2D(
        Element::GeometryType& geom)
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();

        double detJ = x10 * y20-y10 * x20;
        return 0.5*detJ;
    }

    //********************************************************************************
    //********************************************************************************
    /** this function compute the maximum and minimum edge lenghts */
    static inline void SideLenghts2D(
        Element::GeometryType& geom,
        double& hmin, double& hmax)
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();

        double l = x20*x20 + y20*y20;
        hmax = l;
        hmin = l;

        if(l>hmax) hmax = l;
        else if(l<hmin) hmin = l;

        l = (x20-x10)*(x20-x10) + (y20-y10)*(y20-y10);
        if(l>hmax) hmax = l;
        else if(l<hmin) hmin = l;

        hmax = sqrt(hmax);
        hmin = sqrt(hmin);
    }



    static inline void CalculateGeometryData(
        Element::GeometryType& geom,
        boost::numeric::ublas::bounded_matrix<double,2,1>& DN_DX,
        array_1d<double,2>& N,
        double& Area)
    {
        double x10 = fabs(geom[1].X() - geom[0].X());


        double detJ = x10;

        DN_DX(0,0) = 1.0/(geom[0].X() - geom[1].X());
        DN_DX(0,1) = 1.0/(geom[1].X() - geom[0].X());


        //DN_DX /= detJ;
        N[0] = 0.5;
        N[1] = 0.5;

        Area = detJ;
    }

    /**
     * Calculate the exact distances to an isosurface define by a set of initial
     * distances
     * @param ThisGeometryThe tetrahedra itself. Note: If the geometry is not a
     * tetrahedra the result is undefined and may cause memory error.
     * @param Distances The distances which define the isosurface as input and
     * the same argument is used to give the calculated exact distance
     */
    template<std::size_t TSize>
    static void CalculateTetrahedraDistances(Element::GeometryType& ThisGeometry, array_1d<double, TSize>& Distances)
    {
        // Calculating the intersection points
        array_1d<Point<3>, 4> intersection_points;
        int number_of_intersection_points = CalculateTetrahedraIntersectionPoints(ThisGeometry, Distances, intersection_points);

//        for(int i = 0 ; i < number_of_intersection_points ; i++)
//            KRATOS_WATCH(intersection_points[i]);


		if(number_of_intersection_points == 0)
		{
			std::cout << "Warning: The intersection with interface hasn't found!" << std::endl;
			std::cout << "Warning: The distances are: " << Distances << std::endl;
		}
		else if(number_of_intersection_points == 1)
		{ // There is one point with zero distance. The distance of the nodes are their distance to this point
//                    std::cout << "1 intersection point" << std::endl;
			array_1d<double,3> temp;
			// loop over nodes to calculate their distance to the zero distance node.
                        for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
                        {
				noalias(temp) = intersection_points[0] - ThisGeometry[i_node];
				Distances[i_node] = norm_2(temp);
//                        KRATOS_WATCH(ThisGeometry[i_node].Id());
//                        KRATOS_WATCH(intersection_points[0]);
//                        KRATOS_WATCH(ThisGeometry[i_node].Coordinates());
//                        KRATOS_WATCH(Distances[i_node]);

			}
//                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
//                       {
//                           Distances[i_node] = fabs(ThisGeometry[i_node].Z()); // To be removed. Pooyan.
//                       }

		}
		else if(number_of_intersection_points == 2)
		{

//                    std::cout << "2 intersection points" << std::endl;
			// loop over nodes to calculate their distance to the zero distance line.
                        for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
                        {
				Distances[i_node] = PointDistanceToLineSegment3D(intersection_points[0], intersection_points[1], ThisGeometry[i_node]);
//                        KRATOS_WATCH(intersection_points[0]);
//                        KRATOS_WATCH(intersection_points[1]);
//                        KRATOS_WATCH(ThisGeometry[i_node]);
//                        KRATOS_WATCH(Distances[i_node]);
			}
//                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
//                       {
//                           Distances[i_node] = ThisGeometry[i_node].Z(); // To be removed. Pooyan.
//                       }

		}
		else if(number_of_intersection_points == 3)
		{
//                    std::cout << "3 intersection points" << std::endl;
			// loop over nodes to calculate their distance to the zero distance triangle.
                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
                       {
				Distances[i_node] = PointDistanceToTriangle3D(intersection_points[0], intersection_points[1], intersection_points[2], ThisGeometry[i_node]);
//                           Distances[i_node] = fabs(ThisGeometry[i_node].Z()); // To be removed. Pooyan.
                       }

		}
		else if(number_of_intersection_points == 4)
		{
                    //std::cout << "4 intersection points" << std::endl;
//                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
//                       {
//                           Distances[i_node] = fabs(ThisGeometry[i_node].Z()); // To be removed. Pooyan.
//                        }

			// loop over nodes to calculate their distance to the each zero distance triangle.
                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
                       {   // here I'm taking in account the order of edges where I'm looking for intersection
                           double d1 = PointDistanceToTriangle3D(intersection_points[0], intersection_points[1], intersection_points[3], ThisGeometry[i_node]);
                           double d2 = PointDistanceToTriangle3D(intersection_points[0], intersection_points[3], intersection_points[2], ThisGeometry[i_node]);

//                           KRATOS_WATCH(d1);
//                           KRATOS_WATCH(d2);
//                           KRATOS_WATCH(Distances[i_node] );

 			   Distances[i_node] = (d1 > d2) ? d2 : d1;
                       }

		}



    }

    /**
     * Calculate the exact distances to an isosurface define by a set of initial
     * distances
     * @param ThisGeometryThe Triangle itself. Note: If the geometry is not a
     * tetrahedra the result is undefined and may cause memory error.
     * @param Distances The distances which define the isosurface as input and
     * the same argument is used to give the calculated exact distance
     */
    template<std::size_t TSize>
    static void CalculateTriangleDistances(Element::GeometryType& ThisGeometry, array_1d<double, TSize>& Distances)
    {
        // Calculating the intersection points
        array_1d<Point<3>, 4> intersection_points;
        int number_of_intersection_points = CalculateTetrahedraIntersectionPoints(ThisGeometry, Distances, intersection_points);

//        for(int i = 0 ; i < number_of_intersection_points ; i++)
//            KRATOS_WATCH(intersection_points[i]);


		if(number_of_intersection_points == 0)
		{
			std::cout << "Warning: The intersection with interface hasn't found!" << std::endl;
			std::cout << "Warning: The distances are: " << Distances << std::endl;
		}
		else if(number_of_intersection_points == 1)
		{ // There is one point with zero distance. The distance of the nodes are their distance to this point
//                    std::cout << "1 intersection point" << std::endl;
			array_1d<double,3> temp;
			// loop over nodes to calculate their distance to the zero distance node.
                        for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
                        {
				noalias(temp) = intersection_points[0] - ThisGeometry[i_node];
				Distances[i_node] = norm_2(temp);
//                        KRATOS_WATCH(ThisGeometry[i_node].Id());
//                        KRATOS_WATCH(intersection_points[0]);
//                        KRATOS_WATCH(ThisGeometry[i_node].Coordinates());
//                        KRATOS_WATCH(Distances[i_node]);

			}
//                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
//                       {
//                           Distances[i_node] = fabs(ThisGeometry[i_node].Z()); // To be removed. Pooyan.
//                       }

		}
		else if(number_of_intersection_points == 2)
		{

//                    std::cout << "2 intersection points" << std::endl;
			// loop over nodes to calculate their distance to the zero distance line.
                        for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
                        {
				Distances[i_node] = PointDistanceToLineSegment3D(intersection_points[0], intersection_points[1], ThisGeometry[i_node]);
//                        KRATOS_WATCH(intersection_points[0]);
//                        KRATOS_WATCH(intersection_points[1]);
//                        KRATOS_WATCH(ThisGeometry[i_node]);
//                        KRATOS_WATCH(Distances[i_node]);
			}
//                       for(unsigned int i_node = 0; i_node < ThisGeometry.size() ; i_node++)
//                       {
//                           Distances[i_node] = ThisGeometry[i_node].Z(); // To be removed. Pooyan.
//                       }

		}
		else
		{
			std::cout << "This is a triangle with more than two intersections!" << std::endl;	
			std::cout << "Warning: Too many intersections: " << number_of_intersection_points << std::endl;		
			std::cout << "Warning: The distances are: " << Distances << std::endl;
			
		}

    }


    /**
     * This function calculates the coordinates of the intersecion points
     * between edges of tetrahedra and a isosurface given by the distances in
     * its corners
     * @param ThisGeometry The tetrahedra itself. Note: If the geometry is not a
     * tetrahedra the result is undefined and may cause memory error.
     * @param Distances The distances of the 4 nodes of the tetrahedra to the
     * iso surface.
     * @param IntersectionPoints The result intersection points
     * @return Number of intersection points.
     */
    template<std::size_t TSize1, std::size_t TSize2>
    static int CalculateTetrahedraIntersectionPoints(Element::GeometryType& ThisGeometry, array_1d<double, TSize1>& Distances, array_1d<Point<3>, TSize2>& IntersectionPoints)
    {
        const double epsilon = 1e-15; //1.00e-9;

        int number_of_intersection_points = 0;
        for(unsigned int i = 0 ; i < TSize1 ; i++)
        {
            if(fabs(Distances[i]) < epsilon)
            {
                noalias(IntersectionPoints[number_of_intersection_points].Coordinates()) = ThisGeometry[i].Coordinates();

                number_of_intersection_points++;
                continue;
            }
            for(unsigned int j = i + 1 ; j < TSize1 ; j++)
            {
                if(fabs(Distances[j]) < epsilon)
                    continue; // we will add it to the intersections by the i index to be unique

                if(Distances[i] * Distances[j] < 0.00)  // The interface passes through the edge
                {

                    double delta_d = fabs(Distances[i]) + fabs(Distances[j]);  // we know that both distances are greater than epsilon.

                    double di = fabs(Distances[i]) / delta_d;
                    double dj = fabs(Distances[j]) / delta_d;

                    noalias(IntersectionPoints[number_of_intersection_points].Coordinates()) = dj * ThisGeometry[i].Coordinates();
                    noalias(IntersectionPoints[number_of_intersection_points].Coordinates()) += di * ThisGeometry[j].Coordinates();

                    number_of_intersection_points++;
                }
            }
        }

		return number_of_intersection_points;
	}


    /**
     * This function calculates the distance of a 3D point to a 3D line segment
     * @param LinePoint1 First point of the line segment
     * @param LinePoint2 End point of the line segment
     * @param ToPoint The point which distance is required
     * @return The distance between the point and the line
     */
    static double PointDistanceToLineSegment3D(Point<3> const& LinePoint1,
                                 Point<3> const& LinePoint2,
                                  Point<3> const& ToPoint)
    {
        const double epsilon = 1e-15; //1.00e-9;

		array_1d<double,3> v1 = LinePoint2 - LinePoint1;
		array_1d<double,3> v2 = LinePoint1 - ToPoint;
		array_1d<double,3> v3;

//                KRATOS_WATCH(LinePoint1);
//                KRATOS_WATCH(LinePoint2);
//                KRATOS_WATCH(ToPoint.Coordinates());
//
//

		double square_distance = inner_prod(v1,v1);

		if(square_distance < epsilon) // near zero length line
			return norm_2(v2); // we return the distance to the first point of line

		double t = - inner_prod(v1,v2) / square_distance;

//                KRATOS_WATCH(t);

		if(t < 0.00) // it is before point 1
		{ // we return the distance to point 1
			v3 = LinePoint1 - ToPoint;

			return norm_2(v3);
		}

		if(t > 1.00) // it is after point 2
		{ // we return the distance to point 2
			v3 = LinePoint2 - ToPoint;

			return norm_2(v3);
		}

		// The projection point is between point 1 and 2 of the line segment
		v3 = LinePoint1 * (1.00 - t) + LinePoint2 * t;
//
//                KRATOS_WATCH(v3);
//                KRATOS_WATCH(v3 - ToPoint);

		return norm_2(v3 - ToPoint);

    }

    /**
     * This function calculates the distance of a 3D point to a 3D triangle
     * @param TrianglePoint1 First point of triangle
     * @param TrianglePoint2 Second point of triangle
     * @param TrianglePoint3 Third point of triangle
     * @param ToPoint The point which distance is required
     * @return The distance between the point and the triangle
     */
    static double PointDistanceToTriangle3D(Point<3> const& TrianglePoint1,
                                 Point<3> const& TrianglePoint2,
                                 Point<3> const& TrianglePoint3,
                                 Point<3> const& ToPoint)
    {
		// The implementation is done using following reference:
		// http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf



		array_1d<double, 3> e0 = TrianglePoint2 - TrianglePoint1;
		array_1d<double, 3> e1 = TrianglePoint3 - TrianglePoint1;
		array_1d<double, 3> dd = TrianglePoint1 - ToPoint;

		double a = inner_prod(e0, e0);
		double b = inner_prod(e0, e1);
		double c = inner_prod(e1, e1);
		double d = inner_prod(e0, dd);
		double e = inner_prod(e1, dd);
		double f = inner_prod(dd, dd);

		double det = a*c-b*b;
		double s = b*e-c*d;
		double t = b*d-a*e;

		double square_distance = 0.00;

		if ( s + t <= det )
		{
			if ( s < 0.00 )
			{
				if ( t < 0.00 )
				{ // region 4
					if (d < 0)
					{
						t = 0;
						if (-d >= a)
						{
							s = 1;
							square_distance = a + 2*d + f;
						}
						else
						{
							s = -d/a;
							square_distance = d*s + f;
						}
					}
					else
					{
						s = 0;
						if (e >= 0)
						{
							t = 0;
							square_distance = f;
						}
						else
						{
							if (-e >= c)
							{
								t = 1;
								square_distance = c + 2*e + f;
							}
							else
							{
								t = -e/c;
								square_distance = e*t + f;
							}
						}
					}
				}
				else
				{ // region 3
					s = 0.00;
					if(e >= 0.00)
					{
						t = 0.00;
						square_distance = f;
					}
					else
					{
						if (-e >= c)
						{
							t = 1.00;
							square_distance = c + 2*e +f;
						}
						else
						{
							t = -e/c;
							square_distance = e*t + f;
						}
					}

				}
			}
			else if ( t < 0.00 )
                        { // region 5
                            t = 0;
                            if (d >= 0)
                            {
                                s = 0;
                                square_distance = f;
                            }
                            else
                            {
                                if (-d >= a)
                                {
                                    s = 1;
                                    square_distance = a + 2.00 * d + f;
                                }
                                else
                                {
                                    s = -d / a;
                                    square_distance = d * s + f;
                                }
                            }
                        }
                        else
			{ // region 0
				double inv_det = 1.00 / det;
				s *= inv_det;
				t *= inv_det;
				square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
			}
		}
		else
		{
			if ( s < 0.00 )
			{ // region 2
				double temp0 = b + d;
				double temp1 = c + e;
				if (temp1 > temp0)  // minimum on edge s+t=1
				{
					double numer = temp1 - temp0;
					double denom = a - 2*b + c;
					if(numer >= denom)
					{
						s = 1.00;
						t = 0.00;
						square_distance = a + 2*d + f;
					}
					else
					{
						s = numer/denom;
						t = 1.00-s;
						square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
					}
				}
				else          // minimum on edge s=0
				{
					s = 0.00;
					if(temp1 <= 0.00)
					{
						t = 1;
						square_distance = c + 2*e + f;
					}
					else
					{
						if(e >= 0.00)
						{
							t = 0.00;
							square_distance = f;
						}
						else
						{
							t = -e/c;
							square_distance = e*t + f;
						}
					}
				}
			}
			else if ( t < 0.00 )
			{ // region 6
				double temp0 = b + e;
				double temp1 = a + d;
				if (temp1 > temp0)
				{
					double numer = temp1 - temp0;
					double denom = a - 2*b + c;
					if(numer >= denom)
					{
						s = 0.00;
						t = 1.00;
						square_distance = c + 2*e + f;
					}
					else
					{
						t = numer/denom;
						s = 1.00-t;
						square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
					}
				}
				else
				{
					t = 0.00;
					if(temp1 <= 0.00)
					{
						s = 1;
						square_distance = a + 2*d + f;
					}
					else
					{
						if(d >= 0.00)
						{
							s = 0.00;
							square_distance = f;
						}
						else
						{
							s = -d/a;
							square_distance = d*s + f;
						}
					}
				}
			}
			else
			{ // region 1
				double numer = c + e - b - d;

				if (numer <= 0.00)
				{
					s = 0.00;
					t = 1.00;
					square_distance = c + 2.00 * e + f;
				}
				else
				{
					double denom = a - 2.00 * b + c;
					if (numer >= denom)
					{
						s = 1.00;
						t = 0.00;
						square_distance = a + 2.00 * d + f;
					}
					else
					{
						s = numer / denom;
						t = 1.00 - s;
						square_distance = s*(a*s + b*t + 2*d) + t*(b*s + c*t + 2*e) + f;
					}
				}
			}
		}

                if(square_distance < 0.00)
                    return 0.00; // avoiding -0 case!!

        return std::sqrt(square_distance);
    }



};

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_UTILITIES_INCLUDED  defined 


