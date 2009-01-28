/*
==============================================================================
KratosPFEMApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
*   Last Modified by:    $Author: pooyan $
*   Date:                $Date: 2008-11-13 12:12:16 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_NORMAL_TO_WALL_CALCULATION_UTILS )
#define  KRATOS_NORMAL_TO_WALL_CALCULATION_UTILS


/* System includes */


/* External includes */


/* Project includes */
#include "utilities/math_utils.h"


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
	class NormalToWallCalculationUtils
    {
    public:
		/**@name Type Definitions */       
		/*@{ */
		typedef ModelPart::NodesContainerType NodesArrayType; 
		typedef ModelPart::ConditionsContainerType ConditionsArrayType;
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
		//assumes the normals are already calculated
        void CalculateNormalToWall(
			ConditionsArrayType& rConditions,
			int dimension)
        {
			KRATOS_TRY

			//resetting the normals
			array_1d<double,3> zero;
			noalias(zero) = ZeroVector(3);
				
			for(ConditionsArrayType::iterator it =  rConditions.begin(); it !=rConditions.end(); it++)
			{
				Element::GeometryType& rNodes = it->GetGeometry();
				for(unsigned int in = 0; in<rNodes.size(); in++)
					noalias((rNodes[in]).FastGetSolutionStepValue(NORMAL_TO_WALL)) = zero;
			}

			//calculating the normals to the walls by adding up the contributions of elements which are exclusively on the wall
			for(ConditionsArrayType::iterator it =  rConditions.begin(); 
							it !=rConditions.end(); it++)
			{
				Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

				double coeff = 1.00/double(pGeometry.size());
				const array_1d<double,3>& An = it->GetValue(NORMAL);

				//count the nodes on the boundary which are "of wall"
				unsigned int structural_nodes = 0;
				for(unsigned int i = 0; i<pGeometry.size(); i++)
					if(pGeometry[i].FastGetSolutionStepValue(IS_STRUCTURE) == 1)
						structural_nodes += 1;

				if(structural_nodes == pGeometry.size())
				{
					for(unsigned int i = 0; i<pGeometry.size(); i++)
					{
						noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL_TO_WALL)) += coeff * An;
					}
				}
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
		//this function adds the Contribution of one of the geometries 
		//to the corresponding nodes

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
		
        //NormalToWallCalculationUtils(void);
		
        //NormalToWallCalculationUtils(NormalToWallCalculationUtils& rSource);
		
		
		/*@} */   
		
    }; /* Class ClassName */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_NORMAL_TO_WALL_CALCULATION_UTILS  defined */

