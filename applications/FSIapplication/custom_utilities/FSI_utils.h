
/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2009-01-13 15:39:55 $
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_FSI_UTILS )
#define  KRATOS_FSI_UTILS


/* System includes */


/* External includes */


/* Project includes */
#include "utilities/math_utils.h"
#include "fsi_application.h"
//#include "custom_conditions/coupling_face2D.h"


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
	class FSIUtils
    {
    public:
		/**@name Type Definitions */       
		/*@{ */
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
 /*       void GenerateCouplingElements(ModelPart& structural_model_part, int dim)
        {
			KRATOS_TRY

			if(dim == 2)
			{
				for(ModelPart::ConditionsContainerType::iterator it = structural_model_part.ConditionsBegin();
																it != structural_model_part.ConditionsEnd();
																it ++ )
				{
					Geometry< Node<3> >::Pointer pgeom =it->pGetGeometry();
					Properties::Pointer props = structural_model_part.GetMesh().pGetProperties(1);

					unsigned int id = (structural_model_part.Elements().end() - 1)->Id() + 1;
					Element::Pointer newel = Element::Pointer( new CouplingFace2D(id, pgeom,props) );

					(structural_model_part.Elements()).push_back(newel);
				}
			}

			KRATOS_CATCH("")			
		}
*/
				

		//***********************************************************************
		//***********************************************************************
 		bool CheckPressureConvergence(ModelPart::NodesContainerType& fluid_interface_nodes, double toll, double abs_toll)
		{
			KRATOS_TRY

			double dp_norm = 0.00;
			double p_norm = 0.00;
			double dp=0.0;
			double p=0.0;

			//verifies that the pressure variation is sufficiently small
			for(ModelPart::NodesContainerType::iterator it = fluid_interface_nodes.begin(); it != fluid_interface_nodes.end(); it++)
			{
				p = it->FastGetSolutionStepValue(PRESSURE);
				dp = p - (it->FastGetSolutionStepValue(PRESSURE_OLD_IT));
				dp_norm += dp*dp;
				p_norm += p*p;
			}

			if( p_norm < 1e-20)
				p_norm = 1.00;

			double ratio = sqrt(dp_norm/p_norm);
			std::cout << "dp_norm = " << dp_norm << std::endl;
			std::cout << "FSI convergence ratio = " << ratio << std::endl;

                        if( dp < abs_toll*abs_toll )
                            return true;
                        else
                        {
                            if(ratio < toll)
                                    return true;
                            else
                                    return false;
                        }
			KRATOS_CATCH("")
		}
		
		//***********************************************************************
		//***********************************************************************
		void StructuralPressurePrediction( Variable<double>& rVariable, ModelPart::NodesContainerType& structure_interface_nodes, const int structure_buffer_size, const int order )
		{
			KRATOS_TRY
					
			if( order > structure_buffer_size-1)
				KRATOS_ERROR(std::logic_error,"using a force prediction order higher than the buffer size ... increase the buffer size or reduce the prediction order","");
					
			if(order == 1)
			{
				std::cout << "first order force prediction" << std::endl;
				for(ModelPart::NodesContainerType::iterator it = structure_interface_nodes.begin(); it != structure_interface_nodes.end(); it++)
				{
					it->FastGetSolutionStepValue(rVariable) = it->FastGetSolutionStepValue(rVariable,1);
				}
			}
			else if (order == 2)
			{
				std::cout << "second order force prediction" << std::endl;
				for(ModelPart::NodesContainerType::iterator it = structure_interface_nodes.begin(); it != structure_interface_nodes.end(); it++)
				{
					it->FastGetSolutionStepValue(rVariable) = 2.0*it->FastGetSolutionStepValue(rVariable,1) - it->FastGetSolutionStepValue(rVariable,2);
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
		
        //FSIUtils(void);
		
        //FSIUtils(FSIUtils& rSource);
		
		
		/*@} */   
		
    }; /* Class ClassName */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_FSI_UTILS  defined */

