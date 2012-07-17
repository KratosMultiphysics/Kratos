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

/*
 * File:   signed_distance_calculation_utils.h
 * Author: rrossi
 *
 * Created on September 3, 2009, 8:34 AM
 */

#ifndef _SIGNED_DISTANCE_CALCULATION_UTILS_H
#define	_SIGNED_DISTANCE_CALCULATION_UTILS_H
/* System includes */


/* External includes */


/* Project includes */
//#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/body_distance_calculation_utils.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.

  \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	  \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

		\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


			\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

			  \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

				\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

				  \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
template< unsigned int TDim>
class SignedDistanceCalculationUtils
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */


    /** Destructor.
    */

    /*@} */
    /**@name Operators
    */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    void CalculateDistances(
        ModelPart& r_model_part,
        Variable<double>& rDistanceVar,
        const double max_distance
    )
    {
        if(TDim == 2)
            CalculateDistances2D(r_model_part, rDistanceVar, max_distance);
        else if (TDim == 3)
            CalculateDistances3D(r_model_part, rDistanceVar, max_distance);

     //   r_model_part.GetCommunicator().SynchronizeCurrentDataToMin(rDistanceVar);
    }
    
    //***********************************************************************
    //***********************************************************************
    void CalculateDistances2D(
        ModelPart& r_model_part,
        Variable<double>& rDistanceVar,
        const double max_distance
    )
    {
        KRATOS_TRY

        double large_distance = 1e6;
        double tol = 1.0/large_distance;

        //copy rDistance var to GetValue database
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            it->GetValue(rDistanceVar) = it->FastGetSolutionStepValue(rDistanceVar);
            it->FastGetSolutionStepValue(rDistanceVar) = large_distance;
            it->GetValue(IS_VISITED) = 0;
        }

        boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
        array_1d<double,TDim+1> N, distances;
        array_1d<double,TDim> grad_d;
        array_1d<double,3> coord_on_0;
        array_1d<double,3> temp;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> free_surface_points;

        //fill the list of the first elements to be solved for the "distance"
        for(ElementsArrayType::iterator it =  r_model_part.ElementsBegin(); it !=r_model_part.ElementsEnd(); it++)
        {
            Element::GeometryType& geom = it->GetGeometry();

            unsigned int n_positive = 0;
            unsigned int n_negative = 0;
            for(unsigned int kk = 0; kk<TDim+1 ; kk++)
            {
                const double dist = geom[kk].GetValue(rDistanceVar);
                distances[kk] = dist;
                if(dist < 0)
                    n_negative += 1;
                else
                    n_positive += 1;
            }

            if(n_negative > 0 && n_positive > 0) //ELEMENT IS CROSSED BY THE INTERFACE!
            {
//                                    //determine the positions at which the edges are cut by the free surface
//                                    unsigned int counter = 0;
//                                    for (unsigned int i = 0; i < TDim + 1; i++) {
//                                        for (unsigned int j = i + 1; j < TDim + 1; j++) {
//                                            if (distances[i] * distances[j] <= 0.0) {
//                                                double fact = fabs(distances[i]) / (fabs(distances[i]) + fabs(distances[j]));
//                                                //std::cout << i << " " << j << std::endl;
//                                                for (unsigned int k = 0; k < TDim; k++)
//                                                    free_surface_points(counter, k) = geom[i].Coordinates()[k] + (geom[j].Coordinates()[k] - geom[i].Coordinates()[k]) * fact;
//
//                                                counter++;
//                                                //                                KRATOS_WATCH(counter);
//                                            }
//
//                                            if (counter >= TDim)
//                                                break;
//                                        }
//                                        if (counter >= TDim)
//                                            break;
//                                    }
//
//                                    //generate list of virtual points on the free surface (for the element divided)
//                                    //and compute the distance of all of the nodes in the element from such list
//                                    if(TDim == 2)
//                                        KRATOS_ERROR(std::logic_error,"2d not yet implemented","")
//                                    else
//                                    {
//                                        array_1d<double,3> N;
//                                        //internal points
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.2, 0.2, N, free_surface_points, rDistanceVar );
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.6, 0.2, N,free_surface_points, rDistanceVar );
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.2, 0.6, N,free_surface_points, rDistanceVar );
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.33333333333333333333, 0.33333333333333333333, N, free_surface_points, rDistanceVar );
//
//                                        //points on faces
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.5, 0.0, N,free_surface_points, rDistanceVar );
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.0, 0.5, N,free_surface_points, rDistanceVar );
//                                        ComputeDistancesFromVirtualPoint3D(geom,0.5, 0.5, N,free_surface_points, rDistanceVar );
//
//                                        //points on edges
////                                        ComputeDistancesFromVirtualPoint3D(geom,1.0, 0.0, N,free_surface_points, rDistanceVar );
////                                        ComputeDistancesFromVirtualPoint3D(geom,0.0, 1.0, N,free_surface_points, rDistanceVar );
////                                        ComputeDistancesFromVirtualPoint3D(geom,0.0, 0.0, N,free_surface_points, rDistanceVar );
//
//                                        for(unsigned int kk = 0; kk<TDim+1 ; kk++)
//                                            geom[kk].GetValue(IS_VISITED) = 1;
//                                    }




//
                double Volume;
                GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

                //compute the gradient of the distance and normalize it
                noalias(grad_d) = prod(trans(DN_DX),distances);
                double norm = norm_2(grad_d);
                grad_d /= norm;

                //find one division point on one edge
                for(unsigned int i = 1; i<TDim+1; i++)
                {
                    if(distances[0]*distances[i]<=0) //if the edge is divided
                    {
                        double delta_d = fabs(distances[i]) + fabs(distances[0]);

                        if(delta_d>1e-20)
                        {
                            double Ni = fabs(distances[0]) / delta_d;
                            double N0 = fabs(distances[i]) / delta_d;

                            noalias(coord_on_0) = N0 * geom[0].Coordinates();
                            noalias(coord_on_0) += Ni * geom[i].Coordinates();
                        }
                        else
                            noalias(coord_on_0) = geom[0].Coordinates();

                        break;

                    }
                }

                //now calculate the distance of all the nodes from the elemental free surface
                for(unsigned int i = 0; i<TDim+1; i++)
                {
                    noalias(temp) = geom[i].Coordinates();
                    noalias(temp) -= coord_on_0 ;

                    double real_distance = 0.0;
                    for(unsigned int k=0; k<TDim; k++)
                        real_distance += temp[k]*grad_d[k];
                    real_distance = fabs(real_distance);

                    double& dist_i = geom[i].FastGetSolutionStepValue(rDistanceVar);
                    if(real_distance < dist_i)
                        dist_i = real_distance;

                    //select the nodes for the computation of the distance
                    geom[i].GetValue(IS_VISITED) = 1;
                }

            }



        }

        //loop on all nodes to treat correctly the case of the surface coinciding with a node
        //copy rDistance var to GetValue database
        array_1d<double,3> aux;
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(fabs(it->GetValue(rDistanceVar)) < tol)
            {
                it->FastGetSolutionStepValue(rDistanceVar) = 0.0;
                it->GetValue(IS_VISITED) = 1;

                //const array_1d<double,3>& center_coords = it->Coordinates();

                //now loop all of its neighbours and calculate the distance value
//                for (WeakPointerVector< Node<3> >::iterator in = it->GetValue(NEIGHBOUR_NODES).begin();
//                        in != it->GetValue(NEIGHBOUR_NODES).end(); in++)
//                {
//                    const array_1d<double,3>& coords = in->Coordinates();
//                    noalias(aux) = coords;
//                    noalias(aux) -= center_coords;
//                    double dist = norm_2(aux);
//
//                    if(in->FastGetSolutionStepValue(rDistanceVar) > dist)
//                        in->FastGetSolutionStepValue(rDistanceVar) = dist;
//
//                    in->GetValue(IS_VISITED)=1;
//                }

            }
        }

        //compute the distance using the element based approach
        BodyDistanceCalculationUtils util;
        util.CalculateDistances<TDim>(r_model_part.Elements(),rDistanceVar,max_distance);

        //finally change the sign to the distance as needed
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(it->GetValue(rDistanceVar) < 0)
                it->FastGetSolutionStepValue(rDistanceVar) = -it->FastGetSolutionStepValue(rDistanceVar);
        }

        //check if something is wrong
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(fabs(it->FastGetSolutionStepValue(rDistanceVar)) == large_distance)
            {
                KRATOS_WATCH("error in the calculation of the distance for node")
                KRATOS_WATCH(it->Id());

            }
        }




        KRATOS_CATCH("")

    }

    void CalculateDistances3D(
        ModelPart& r_model_part,
        Variable<double>& rDistanceVar,
        const double max_distance
    )
    {
        KRATOS_TRY

        const double large_distance = 1e9;

        std::vector<array_1d<double,TDim + 1> > distances;
        std::vector<Element::Pointer> interface_elements;
        
         //copy rDistance var to GetValue database
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            it->GetValue(rDistanceVar) = it->FastGetSolutionStepValue(rDistanceVar);
            it->GetValue(IS_VISITED) = 0;
        }

        // Make a list of all elements crossed by interface and get their nodal distances
        for(ElementsArrayType::iterator i_element =  r_model_part.ElementsBegin(); i_element !=r_model_part.ElementsEnd(); i_element++)
        {
            // Get a reference to the geometry
            Element::GeometryType& element_geometry = i_element->GetGeometry();
            array_1d<double,TDim + 1> element_distance;

            bool positive = false;
            bool negative = false;
            
            // loop over nodes to see if they have positive or negative distance.
            for(unsigned int i_node = 0; i_node < element_geometry.size() ; i_node++)
            {
                const double distance = element_geometry[i_node].GetSolutionStepValue(rDistanceVar);
                element_distance[i_node] = distance;
                negative |= (distance < 0);
                positive |= (distance >= 0);
            }

            if(negative && positive) //ELEMENT IS CROSSED BY THE INTERFACE!
            {
                interface_elements.push_back(*(i_element.base()));
                distances.push_back(element_distance);
            }
        }
        
        // Reset the distance to a large value maintaining the sign
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            double& distance = it->FastGetSolutionStepValue(rDistanceVar);

            if(distance >= 0.00)
                distance = large_distance;
            else
                distance = -large_distance;
        }

        // Now calculating the new distance for interface elements
        int size = static_cast<int>(interface_elements.size());
        for(int i = 0 ; i < size ; i++)
        {
            // Get a reference to the geometry
            Element::GeometryType& element_geometry = interface_elements[i]->GetGeometry();

            CalculateTetrahedraDistances(element_geometry, distances[i]);

            // loop over nodes and apply the new distances.
            for(unsigned int i_node = 0; i_node < element_geometry.size() ; i_node++)
            {
                double& distance = element_geometry[i_node].GetSolutionStepValue(rDistanceVar);
                double new_distance = distances[i][i_node];

//                if(fabs(fabs(distance) - fabs(new_distance)) > 1.00e-9)
//                {
//                    KRATOS_WATCH(interface_elements[i]->Id());
//                    KRATOS_WATCH(element_geometry[i_node]);
//                    KRATOS_WATCH(distance);
//                    KRATOS_WATCH(new_distance);
//                }



                if(fabs(distance) > fabs(new_distance))
                {
                    if(distance >= 0)
                        distance = new_distance;
                    else
                        distance = new_distance;
                }
                element_geometry[i_node].GetValue(IS_VISITED) = 1;
            }
       }
        //compute the distance using the element based approach
        BodyDistanceCalculationUtils util;
        util.CalculateDistances<TDim>(r_model_part.Elements(),rDistanceVar,max_distance);

        //finally change the sign to the distance as needed
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(it->GetValue(rDistanceVar) < 0)
                it->FastGetSolutionStepValue(rDistanceVar) = -it->FastGetSolutionStepValue(rDistanceVar);
        }

        //check if something is wrong
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(fabs(it->FastGetSolutionStepValue(rDistanceVar)) == large_distance)
            {
                KRATOS_WATCH("error in the calculation of the distance for node")
                KRATOS_WATCH(it->Id());

            }
        }

        KRATOS_CATCH("")

    }

    void CalculateTetrahedraDistances(Element::GeometryType& ThisGeometry, array_1d<double, TDim + 1>& Distances/*, array_1d<double,4>& NewDistances*/)
    {
        // Calculating the intersection points
        array_1d<Point<3>, 4> intersection_points;
        int number_of_intersection_points = CalculateTetrahedraIntersectionPoints(ThisGeometry, Distances, intersection_points);
        
//        for(int i = 0 ; i < number_of_intersection_points ; i++)
//            KRATOS_WATCH(intersection_points[i]);


		if(number_of_intersection_points == 0)
		{
			std::cout << "Warning: The intersection with interface hasn't found!" << std::endl;
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

    int CalculateTetrahedraIntersectionPoints(Element::GeometryType& ThisGeometry, array_1d<double, TDim + 1>& Distances, array_1d<Point<3>, 4>& IntersectionPoints)
    {
        const double epsilon = 1e-15; //1.00e-9;

        int number_of_intersection_points = 0;
        for(int i = 0 ; i < 4 ; i++)
        {
            if(fabs(Distances[i]) < epsilon)
            {
                noalias(IntersectionPoints[number_of_intersection_points].Coordinates()) = ThisGeometry[i].Coordinates();
                
                number_of_intersection_points++;
                continue;
            }
            for(int j = i + 1 ; j < 4 ; j++)
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

	double PointDistanceToLineSegment3D(Point<3> const& LinePoint1, 
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

    double PointDistanceToTriangle3D(Point<3> const& TrianglePoint1, 
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


    //**********************************************************************************+
    //**********************************************************************************+
    double FindMaximumEdgeSize(ModelPart& r_model_part)
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& rNodes = r_model_part.Nodes();

        double h_max = 0.0;

        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            double xc = in->X();
            double yc = in->Y();
            double zc = in->Z();

            double h = 0.0;
            for (WeakPointerVector< Node < 3 > >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                    i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
            {
                double x = i->X();
                double y = i->Y();
                double z = i->Z();
                double l = (x - xc)*(x - xc);
                l += (y - yc)*(y - yc);
                l += (z - zc)*(z - zc);

                if (l > h) h = l;
            }
            h = sqrt(h);

            if(h > h_max) h_max = h;

        }

        return h_max;

        KRATOS_CATCH("");
    }

    /*@} */
    /**@name Acces */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */
//                 void ComputeDistancesFromVirtualPoint3D(
//                                 Geometry< Node<3> >& geom,
//                                 double xi, double eta,
//                                 array_1d<double, 3 > N,
//                                 const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> free_surface_points,
//                                 Variable<double>& rDistanceVar )
//                 {
//                     Point<3> aux(0.0, 0.0, 0.0);
//                     array_1d<double, 3 >& coords = aux.Coordinates();
//
//                     //setting the coordinates for the virtual points to the expected value
//                     N[0] = 1.0 - xi - eta;
//                     N[1] = xi;
//                     N[2] = eta;
//                     for (unsigned int i = 0; i < TDim; i++)
//                     {
//                         for (unsigned int j = 0; j < TDim; j++)
//                             coords[j] += N[i] * free_surface_points(i,j);
//                     }
//
//                     //now compute and verify the distance
//                     for (unsigned int i = 0; i < TDim+1; i++)
//                     {
// 		      const array_1d<double,3>& node_coords = geom[i].Coordinates();
//                         double dist = 0.0;
//                         for (unsigned int j = 0; j < TDim; j++)
//                             dist += pow( coords[j] - node_coords[j] , 2);
//
//                         dist = sqrt(dist);
//
//                         double& node_dist = geom[i].FastGetSolutionStepValue(rDistanceVar);
//                         if(dist < node_dist)
//                             node_dist = dist;
//                     }
//                 }


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Acces */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */



    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/


#endif	/* _SIGNED_DISTANCE_CALCULATION_UTILS_H */

