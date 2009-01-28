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
			bool reorder )
        {
			KRATOS_TRY
					
					
			std::cout << "dimension in distance computation " << TDim << std::endl;
					
					
			//defining work arrays
			PointerVector< Element > elements_to_solve;
			PointerVector< Node<3> > failed_nodes;
//			elements_to_solve.reserve( int(rElements.size() ) );
								
			//fill the list of the first elements to be solved for the "distance"
			for(ElementsArrayType::iterator it =  rElements.begin(); it !=rElements.end(); it++)
			{
				unsigned int visited_nodes = 0;
				Element::GeometryType& geom = it->GetGeometry();
				
				for(unsigned int kk = 0; kk<TDim+1 ; kk++)
				{
					if( geom[kk].GetValue(IS_VISITED) == 1)
						visited_nodes += 1;
				}
				
				//saving the number of nodes already calculated for this element
				it->GetValue(IS_VISITED) = visited_nodes;
				
				//save in the elements to solve if just one node is missing to be determined
				if( visited_nodes == TDim )
					elements_to_solve.push_back( *(it.base()) );
				
				//resetting elemental distance
				it->GetValue(rDistanceVar) = 0.0;
				
			}
	KRATOS_WATCH(elements_to_solve.size());
			//this is the "total" solution loop
			boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX; 
			array_1d<double,TDim+1> N; 
			
			
			//TO DO: improve!!
			unsigned int reorder_frequency = elements_to_solve.size(); //100;
//ANTONIA change begin
			//unsigned int next_reorder = reorder_frequency*2;
			unsigned int next_reorder = reorder_frequency;
//ANTONIA change end

			array_1d<double,TDim> d;
			unsigned int current_position = 0;
			
// PointerVector<Element> temporaneo;
			//element_to_solve.size() is increasing at each loop that's why current_position can be  < elements_to_solve.size()  and  >= next_reorder (which is 2*elements_to_solve.size()INITIAL!!!)
			while(current_position < elements_to_solve.size() )
			{
				
				if(current_position >= next_reorder && reorder==true) //we reorder the elements to be computed depending on their distance
				{
					//saving the maximum distance
					for(PointerVector< Element >::iterator ie = elements_to_solve.begin()+current_position; ie != elements_to_solve.end(); ie++)
					{
						Geometry< Node<3> >& geom = ie->GetGeometry();
						double max_dist = 0.0;
						for(unsigned int iii=0; iii<TDim+1; iii++)
						{
							if( geom[iii].GetValue(IS_VISITED) == 1)
							{
								double node_dist = geom[iii].FastGetSolutionStepValue(rDistanceVar);
								if(node_dist > max_dist) max_dist = node_dist;
							}							
						}
						ie->GetValue(rDistanceVar) = max_dist;
/*KRATOS_WATCH(max_dist);
temporaneo.push_back( Element::Pointer( *(ie.base()) ) );*/
						
					}
					
					
						
					//sorting the elements depending on their max dist
					CompareElementDistance comparor(rDistanceVar);
					std::sort(elements_to_solve.ptr_begin()+current_position,    
							elements_to_solve.ptr_end(), comparor    );
/*					std::sort(elements_to_solve.ptr_begin()+current_position,    
							elements_to_solve.ptr_end(), CompareElementDistance()    );*/
							
					reorder_frequency = (unsigned int)( (elements_to_solve.size() - current_position)*0.95);
							 
					next_reorder = current_position + reorder_frequency;
				}
				
				PointerVector< Element >::iterator current_element = (elements_to_solve.begin()+current_position);
			
				//identify if the element is to be solved
				if(current_element->GetValue(IS_VISITED) == TDim)
				{
					unsigned int unknown_node_index = 0;
					
					Geometry< Node<3> >& geom = current_element->GetGeometry();
					
					double Volume;
					GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
					
					//compute discriminant
					noalias(d) = ZeroVector(TDim);
					for(unsigned int iii=0; iii<TDim+1; iii++)
					{
						double distance = geom[iii].FastGetSolutionStepValue(rDistanceVar);
						
						double node_is_known = geom[iii].GetValue(IS_VISITED);
						
						if( node_is_known == 1) //identyfing the unknown node
						{
							for(unsigned int jjj=0; jjj<TDim; jjj++)
								d[jjj] += DN_DX(iii,jjj) * distance;
						}
						else
						{
							unknown_node_index = iii;
						}	
					}
					
					//finalizing computation of discriminant
					double c = -1.0;
					double a = 0.0;
					double b = 0.0;
					for(unsigned int jjj=0; jjj<TDim; jjj++)
					{
						a += DN_DX(unknown_node_index,jjj) * DN_DX(unknown_node_index,jjj);
						b += d[jjj] * DN_DX(unknown_node_index,jjj);
						c += d[jjj]*d[jjj];						
					}
					b *= 2.0;
					
					double discriminant = b*b - 4.0*a*c;
					
					if(discriminant < 0.0) //element distance computation fails - we may have to tread specially this node
					{
						failed_nodes.push_back(Node<3>::Pointer(geom(unknown_node_index)));
					   //failed_elements.push_back( Element::WeakPointer(*(current_element.base())));	
					}
					else
					{
						//(accurate) computation of the distance
						//requires the solution of a*x^2+b*x+c=0
						double q,root1,root2,distance;
						if( a != 0.0)
						{							
							if (b > 0) q = -0.5 * (b + sqrt(discriminant));
							else       q = -0.5 * (b - sqrt(discriminant));
							root1 = q / a;
							root2 = c / q;
							if(root1 > root2) distance = root1;
							else distance = root2;
						}
						else //in this case we have a linear equation
						{
							distance = -c/b;
						}
						
						//saving the distance
						geom[unknown_node_index].FastGetSolutionStepValue(rDistanceVar) = distance;
						geom[unknown_node_index].GetValue(IS_VISITED) = 1.0;
						
						//loop over neighbour elements and add them to the todo list
						for(WeakPointerVector< Element >::iterator ie = geom[unknown_node_index].GetValue(NEIGHBOUR_ELEMENTS).begin(); ie != geom[unknown_node_index].GetValue(NEIGHBOUR_ELEMENTS).end(); ie++)
						{
							unsigned int times_visited = (unsigned int) ie->GetValue(IS_VISITED);
							
							//increase visitation counter
							if(times_visited != TDim+1)
							{
								times_visited++;
								ie->GetValue(IS_VISITED) = times_visited;
								
							}
							
							if(times_visited == TDim) //it is to be used for the next step
								elements_to_solve.push_back( Element::Pointer( *(ie.base()) ) );
								
								
							
						}
					}
				}
				
				
				//update the counter
				current_position++;
			}
			
			unsigned int confirmed_failures = 0;			
			for(PointerVector< Node<3> >::iterator it=failed_nodes.begin(); it!=failed_nodes.end(); it++)
			{
				if(it->GetValue(IS_VISITED) != 1) //it was not possible to calculate the distance
				{
					confirmed_failures++;
					double davg = 0;
					double counter = 0;
					for(WeakPointerVector< Node<3> >::iterator in = it->GetValue(NEIGHBOUR_NODES).begin(); in != it->GetValue(NEIGHBOUR_NODES).end(); in++)
					{
						if(in->GetValue(IS_VISITED) == 1)
						{
							davg += in->FastGetSolutionStepValue(rDistanceVar);
							counter += 1.0;
						}
					}
					it->FastGetSolutionStepValue(rDistanceVar) = davg/counter;
					if(counter == 0)
					{
						KRATOS_WATCH("distance computation failed for node:");
						KRATOS_WATCH(it->Id());
// 						KRATOS_WATCH("distance is set to zero:");
						//it->FastGetSolutionStepValue(rDistanceVar) = 0.0;
						KRATOS_ERROR(std::logic_error,"no neighbour nodes was succesfully computed ... impossible to recover","");
					}
				}
			}
//			std::cout << "confirmed failures in distance computation = " << confirmed_failures << std::endl;




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

		
		class CompareElementDistance : public std::binary_function<const Element::Pointer , const Element::Pointer, bool >{
			public:
				CompareElementDistance(const Variable<double>& rDistanceVar):mrDistanceVar(rDistanceVar){}
				
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
	
}  /* namespace Kratos.*/

#endif /* KRATOS_BODY_DISTANCE_CALCULATION_UTILS defined */

