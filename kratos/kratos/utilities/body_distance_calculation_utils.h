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


/* *********************************************************   
 *
 *   Last Modified by:    $Author: antonia $
 *   Date:                $Date: 2008-05-13 14:12:19 $
 *   Revision:            $Revision: 1.5 $
 *
 * ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_BODY_DISTANCE_CALCULATION_UTILS )
#define  KRATOS_BODY_DISTANCE_CALCULATION_UTILS


/* System includes */


/* External includes */


/* Project includes */
//#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"


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
    class BodyDistanceCalculationUtils
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

        //***********************************************************************
        //***********************************************************************
        //this function calculates the "area normal" (vector oriented as the normal
        //with a dimension proportional to the area. This is done basing on the volume discretization.

        template< unsigned int TDim>
        void CalculateDistances(
                ElementsArrayType& rElements,
                Variable<double>& rDistanceVar,
                const double max_distance)
        {
            KRATOS_TRY


            std::cout << "dimension in distance computation " << TDim << std::endl;

            //defining work arrays
            PointerVector< Element > elements_to_solve;
            elements_to_solve.reserve(rElements.size());

            PointerVector< Node < 3 > > active_nodes;
            active_nodes.reserve(rElements.size());

            std::vector<double> node_distance_values;
            node_distance_values.reserve(rElements.size());

            PointerVector< Node < 3 > > failed_nodes;



            //fill the list of the first elements to be solved for the "distance"
            for (ElementsArrayType::iterator it = rElements.begin(); it != rElements.end(); it++)
            {
                //reset is visited flag for the elements
                it->GetValue(IS_VISITED) = 0.0;

                unsigned int visited_nodes = 0;
                Element::GeometryType& geom = it->GetGeometry();

                for (unsigned int kk = 0; kk < TDim + 1; kk++)
                {
                    if (geom[kk].GetValue(IS_VISITED) == 1)
                        visited_nodes += 1;
                    else
                        geom[kk].FastGetSolutionStepValue(rDistanceVar) = max_distance;
                }

                //save in the elements to solve if just one node is missing to be determined
                if (visited_nodes == TDim)
                {
                    elements_to_solve.push_back(*(it.base()));
                    it->GetValue(IS_VISITED) = 1;
                }

            }
            //            KRATOS_WATCH(elements_to_solve.size());

            //this is the "total" solution loop
            boost::numeric::ublas::bounded_matrix<double, TDim + 1, TDim> DN_DX;
            array_1d<double, TDim + 1 > N;

            array_1d<double, TDim> d;

            while (elements_to_solve.size() != 0)
            {
//                KRATOS_WATCH(elements_to_solve.size());
                //compute all of the candidate elements
                for (unsigned int current_position = 0; current_position != elements_to_solve.size(); current_position++)
                {
                    PointerVector< Element >::iterator current_element = (elements_to_solve.begin() + current_position);
                    Geometry< Node < 3 > >& geom = current_element->GetGeometry();

                      unsigned int unknown_node_index = 0;



                    double Volume;
                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

                    //compute discriminant
                    noalias(d) = ZeroVector(TDim);
                    for (unsigned int iii = 0; iii < TDim + 1; iii++)
                    {
                        const double distance = geom[iii].FastGetSolutionStepValue(rDistanceVar);
                        double node_is_known = geom[iii].GetValue(IS_VISITED);
                        if (node_is_known == 1) //identyfing the unknown node
                            for (unsigned int jjj = 0; jjj < TDim; jjj++)
                                d[jjj] += DN_DX(iii, jjj) * distance;
                        else
                            unknown_node_index = iii;
                    }

                    //finalizing computation of discriminant
                    double c = -1.0;
                    double a = 0.0;
                    double b = 0.0;
                    for (unsigned int jjj = 0; jjj < TDim; jjj++)
                    {
                        a += DN_DX(unknown_node_index, jjj) * DN_DX(unknown_node_index, jjj);
                        b += d[jjj] * DN_DX(unknown_node_index, jjj);
                        c += d[jjj] * d[jjj];
                    }
                    b *= 2.0;

                    double discriminant = b * b - 4.0 * a*c;

                    if (discriminant < 0.0) //element distance computation fails - we may have to tread specially this node
                    {
                        failed_nodes.push_back(Node < 3 > ::Pointer(geom(unknown_node_index)));
                    } else
                    {
                        //(accurate) computation of the distance
                        //requires the solution of a*x^2+b*x+c=0
                        double q, root1, root2, distance;
                        if (a != 0.0)
                        {
                            if (b > 0) q = -0.5 * (b + sqrt(discriminant));
                            else q = -0.5 * (b - sqrt(discriminant));
                            root1 = q / a;
                            root2 = c / q;
                            if (root1 > root2) distance = root1;
                            else distance = root2;
                        } else //in this case we have a linear equation
                        {
                            distance = -c / b;
                        }

                        //saving the distance
                        active_nodes.push_back(Node < 3 > ::Pointer(geom(unknown_node_index))); //save the pointer to the node to update
                        node_distance_values.push_back(distance); //save a value
                        //                                                    geom[unknown_node_index].FastGetSolutionStepValue(rDistanceVar) = distance;
                        //                                                    geom[unknown_node_index].GetValue(IS_VISITED) = 1.0;
                    }
                }
                //                }

                //now loop over all of the active nodes, and assign the minimum value of distance
                for (unsigned int k = 0; k < active_nodes.size(); k++)
                {
//                                        std::cout << " " << active_nodes[k].Id() << " " << active_nodes[k].FastGetSolutionStepValue(rDistanceVar) <<" " << node_distance_values[k];
                    if (active_nodes(k)->FastGetSolutionStepValue(rDistanceVar) > node_distance_values[k])
                        active_nodes(k)->FastGetSolutionStepValue(rDistanceVar) = node_distance_values[k];
//                                        std::cout << " " << active_nodes[k].FastGetSolutionStepValue(rDistanceVar) << std::endl;
                    active_nodes(k)->GetValue(IS_VISITED) = 1;
                }


  //CHAPUZA TEST
//                 for (PointerVector< Node < 3 > >::iterator it = failed_nodes.begin(); it != failed_nodes.end(); it++)
//                 {
//                     if (it->GetValue(IS_VISITED) != 1) //it was not possible to calculate the distance
//                     {
//                         for (WeakPointerVector< Element >::iterator ie = it->GetValue(NEIGHBOUR_ELEMENTS).begin();
//                                 ie != it->GetValue(NEIGHBOUR_ELEMENTS).end(); ie++)
//                             ie->GetValue(IS_VISITED)=1;
//                     }
//                 }

//                KRATOS_WATCH(active_nodes.size());
//                KRATOS_WATCH(node_distance_values.size());
// 
//                unsigned int k=0;
//                for (PointerVector< Node < 3 > >::iterator it = active_nodes.begin(); it != active_nodes.end(); it++)
//                {
//                    if (it->FastGetSolutionStepValue(rDistanceVar) > node_distance_values[k])
//                        it->FastGetSolutionStepValue(rDistanceVar) = node_distance_values[k];
//
//                    it->GetValue(IS_VISITED) = 1.0;
//                    k++;
//                }



                elements_to_solve.clear();
          

                //now loop over all of the active nodes, and assign the minimum value of distance
                for (PointerVector< Node < 3 > >::iterator it = active_nodes.begin(); it != active_nodes.end(); it++)
                {
                    if(it->FastGetSolutionStepValue(rDistanceVar) < max_distance)
                    {
                        //loop over neighbour elements and add them to the todo list
                        for (WeakPointerVector< Element >::iterator ie = it->GetValue(NEIGHBOUR_ELEMENTS).begin();
                                ie != it->GetValue(NEIGHBOUR_ELEMENTS).end(); ie++)
                        {
                            unsigned int visited_nodes = 0;
                            Element::GeometryType& geom = ie->GetGeometry();

                            for (unsigned int kk = 0; kk < TDim + 1; kk++)
                            {
                                if (geom[kk].GetValue(IS_VISITED) == 1)
                                    visited_nodes += 1;
                            }

                            if (visited_nodes == TDim)
                                if( ie->GetValue(IS_VISITED) != 1 ) //it is to be used for the next step (an was not added before by another element)
                                {
                                    ie->GetValue(IS_VISITED) = 1;
                                    elements_to_solve.push_back(Element::Pointer(*(ie.base())));
                                }

    //                        if(visited_nodes == TDim+1)
    //                            std::cout << "should not add the element" << ie->Id() << std::endl;

                        }
                    }
                }



                //erase working arrays
                active_nodes.clear();
                node_distance_values.clear();
 

            }



            //approximate computation of distance on failed nodes
            unsigned int confirmed_failures = 0;
            for (PointerVector< Node < 3 > >::iterator it = failed_nodes.begin(); it != failed_nodes.end(); it++)
            {
                if (it->GetValue(IS_VISITED) != 1 ) //it was not possible to calculate the distance
                {
                    confirmed_failures++;
                    double davg = 0.0;
                    double counter = 0.0;
                    for (WeakPointerVector< Node < 3 > >::iterator in = it->GetValue(NEIGHBOUR_NODES).begin(); in != it->GetValue(NEIGHBOUR_NODES).end(); in++)
                    {
                        if (in->GetValue(IS_VISITED) == 1)
                        {
                            davg += in->FastGetSolutionStepValue(rDistanceVar);
                            counter += 1.0;
                        }
                    }
//CHAPUZA WITH MAX DISTANCE
                    if (counter != 0.0)
                    {
                        double estimated_distance = davg / counter;
                        if(estimated_distance < max_distance)
                            it->FastGetSolutionStepValue(rDistanceVar) = estimated_distance;
			else
				it->FastGetSolutionStepValue(rDistanceVar) = max_distance;	
                    }
                    else
                    {
                        KRATOS_WATCH("distance computation failed for node:");
                        KRATOS_WATCH(it->Id());
                        // 						KRATOS_WATCH("distance is set to zero:");
                        //it->FastGetSolutionStepValue(rDistanceVar) = 0.0;
                        KRATOS_ERROR(std::logic_error, "no neighbour nodes was succesfully computed ... impossible to recover", "");
                    }
                }
            }

            //set all of the elements as not visited




            KRATOS_CATCH("")

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
        //this function adds the Contribution of one of the geometries
        //to the corresponding nodes

        /*@} */
        /**@name Private Operators*/

        /*@{ */


        class CompareElementDistance : public std::binary_function<const Element::Pointer, const Element::Pointer, bool >
        {
        public:

            CompareElementDistance(const Variable<double>& rDistanceVar) : mrDistanceVar(rDistanceVar) { }

            bool operator()(const Element::Pointer a,
                    const Element::Pointer b) const
            {
                return a->GetValue(mrDistanceVar) < b->GetValue(mrDistanceVar);
            }

        private:
            const Variable<double>& mrDistanceVar;

        };
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

        //BodyDistanceCalculationUtils(void);

        //BodyDistanceCalculationUtils(BodyDistanceCalculationUtils& rSource);


        /*@} */

    }; /* Class ClassName */

    /*@} */

    /**@name Type Definitions */
    /*@{ */


    /*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_BODY_DISTANCE_CALCULATION_UTILS defined */

