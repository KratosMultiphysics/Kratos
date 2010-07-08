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
*   Last Modified by:    $Author: antonia $
*   Date:                $Date: 2009-10-19 12:12:16 $
*   Revision:            $Revision: 0.1 $
*
* ***********************************************************/

#if !defined(KRATOS_DRAG_UTILS )
#define  KRATOS_DRAG_UTILS


/* System includes */
#include <string>
#include <iostream> 
#include <algorithm> 

/* External includes */


/* Project includes */
#include "utilities/math_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
#include "includes/node.h"
// #include "includes/element.h"

#include "PFEM_application.h"


//Database includes
// #include "spatial_containers/spatial_containers.h"

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
	
	class DragUtils
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

		
		
		
		void CalculateFluidDrag(
			ModelPart& rFluid_ModelPart,
			Variable< array_1d<double,3> >& rSeepageDragVar,
			const double& mu 
			)
 		{

			KRATOS_TRY

			for( ModelPart::NodesContainerType::iterator inode = rFluid_ModelPart.NodesBegin();
				   inode != rFluid_ModelPart.NodesEnd();
				   inode++)	
			{
			      if(inode->FastGetSolutionStepValue(DISTANCE) < 0.0)
			      {
				  const double& porosity = inode->FastGetSolutionStepValue(POROSITY);
				  const double& diameter = inode->FastGetSolutionStepValue(DIAMETER);
				  //Nodal velocity
				  const array_1d<double,3>& vel = inode->FastGetSolutionStepValue(VELOCITY);
				  double vel_norm = norm_2(vel);
				  array_1d<double,3>& darcy_term = inode->FastGetSolutionStepValue(rSeepageDragVar);
// 				  KRATOS_WATCH(inode->Id())
// 				  KRATOS_WATCH(vel)
// 				  KRATOS_WATCH(vel_norm)
				  
				  double lin_coeff = 150 * mu * (1-porosity)*(1-porosity)/(porosity*porosity*diameter*diameter);
				  double non_lin_coeff = 1.75 * (1-porosity)/(porosity*porosity*diameter);
				  noalias(darcy_term) = (lin_coeff + non_lin_coeff * vel_norm) * vel;
// 				  KRATOS_WATCH(darcy_term);
			      }
			}
			
			
			
			KRATOS_CATCH("")

		}


		void AddDrag(
			ModelPart& rStructural_ModelPart,
			const Variable< array_1d<double,3> >& rSeepageDragVar,
			Variable< array_1d<double,3> >& rBodyForce,
			const array_1d<double,3>& rGravityAcc
			)
 		{

			KRATOS_TRY
			for( ModelPart::NodesContainerType::iterator inode = rStructural_ModelPart.NodesBegin();
				   inode != rStructural_ModelPart.NodesEnd();
				   inode++)	
			{
			  if(inode->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
			  {
			  
			    array_1d<double,3>& bodyforce_drag = inode->FastGetSolutionStepValue(rBodyForce);
			    bodyforce_drag  = inode->FastGetSolutionStepValue(rSeepageDragVar);
			    bodyforce_drag  +=  rGravityAcc;
// 			    KRATOS_WATCH(bodyforce_drag)
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
		///@name Static Member rVariables 
		///@{ 


		///@} 
		///@name Member rVariables 
		///@{ 
			

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
		
        //DragUtils(void);
		
        //DragUtils(DragUtils& rSource);
		
		
		/*@} */   
		
    }; /* Class ClassName */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_DRAG_UTILS  defined */

