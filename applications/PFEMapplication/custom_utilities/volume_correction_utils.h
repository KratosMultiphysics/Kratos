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
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2007-03-06 10:30:31 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_VOLUME_CORRECTION_UTILS )
#define  KRATOS_VOLUME_CORRECTION_UTILS


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
	class VolumeCorrectionUtils
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
        void CorrectVolume(double old_volume, double new_volume, ModelPart::NodesContainerType& nodes)
        {
		KRATOS_TRY

		double delta_vol = new_volume - old_volume;

		//calculate the total area of the free surface by using the NORMAL
		double area_free_surf = 0.0;

		for(ModelPart::NodesContainerType::iterator it =  nodes.begin(); it!=nodes.end(); it++)
		{
			if(it->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1 && 
				(it->FastGetSolutionStepValue(NEIGHBOUR_NODES)).size() != 0 )
			{
				const array_1d<double,3>& An = it->FastGetSolutionStepValue(NORMAL);
				area_free_surf += norm_2(An);
			}
		}

		double h_node = -delta_vol / area_free_surf;
KRATOS_WATCH(delta_vol);
KRATOS_WATCH(area_free_surf);
KRATOS_WATCH(h_node);
		for(ModelPart::NodesContainerType::iterator it =  nodes.begin(); it!=nodes.end(); it++)
		{
			if(it->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1 && 
				(it->FastGetSolutionStepValue(NEIGHBOUR_NODES)).size() != 0 )
			{
				const array_1d<double,3>& An = it->FastGetSolutionStepValue(NORMAL);
				double Aloc = norm_2(An);

				if(Aloc > 1.0e-15)
				{
					array_1d<double,3>& disp = it->FastGetSolutionStepValue(DISPLACEMENT);
	
					double temp = h_node/Aloc;
					noalias(disp) += temp * An;
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
		
        //VolumeCorrectionUtils(void);
		
        //VolumeCorrectionUtils(VolumeCorrectionUtils& rSource);
		
		
		/*@} */   
		
    }; /* Class ClassName */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_VOLUME_CORRECTION_UTILS  defined */

