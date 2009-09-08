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

		//***********************************************************************
		//***********************************************************************
		//this function calculates the "area normal" (vector oriented as the normal
		//with a dimension proportional to the area. This is done basing on the volume discretization.
		void CalculateDistances(
			ModelPart& r_model_part,
			Variable<double>& rDistanceVar)
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
                                    double Volume;
                                    GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

                                    //compute the gradient of the distance and normalize it
                                    noalias(grad_d) = prod(trans(DN_DX),distances);
                                    double norm = norm_2(grad_d);
                                    grad_d /= norm;

                                    //find one division point on one edge
                                    for(unsigned int i = 1; i<TDim+1; i++)
                                    {
                                        if(distances[0]*distances[i]<0) //if the edge is divided
                                         {
                                                double delta_d = fabs(distances[i]) + fabs(distances[0]);

                                                double Ni = fabs(distances[0]) / delta_d;
                                                double N0 = fabs(distances[i]) / delta_d;

                                                noalias(coord_on_0) = N0 * geom[0].Coordinates();
                                                noalias(coord_on_0) += Ni * geom[i].Coordinates();

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

//                        //loop on all nodes to treat correctly the case of the surface coinciding with a node
//                        //copy rDistance var to GetValue database
//			for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
//			{
//                            if(fabs(it->GetValue(rDistanceVar)) < tol)
//                            {
//                                it->FastGetSolutionStepValue(rDistanceVar) = 0.0;
//                                it->GetValue(IS_VISITED) = 1;
//                            }
//                        }
                        
                        //compute the distance using the element based approach
                        BodyDistanceCalculationUtils util;
                        util.CalculateDistances<TDim>(r_model_part.Elements(),rDistanceVar,true);

                        //finally change the sign to the distance as needed
			for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
			{
                            if(it->GetValue(rDistanceVar) < 0)
                                it->FastGetSolutionStepValue(rDistanceVar) = -it->FastGetSolutionStepValue(rDistanceVar);
                        }




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

		/*@} */
		/**@name Private Operators*/
		/*@{ */


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

