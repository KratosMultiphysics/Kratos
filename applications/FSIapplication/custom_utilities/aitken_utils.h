
/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2007-08-17 11:48:30 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_AITKEN_UTILS )
#define  KRATOS_AITKEN_UTILS


/* System includes */


/* External includes */


/* Project includes */
#include "utilities/math_utils.h"
#include "fsi_application.h"

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
	class AitkenUtils
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
	//save Variables needed for Aitken in the Non Historical Database
	double ComputeAitkenFactor(ModelPart::NodesContainerType& rInterfaceNodes, const double old_mu)
	{
		KRATOS_TRY
				
/*double rel_it = 0.0;
double rel = 0.0;*/
		double num = 0.0;
		double denom = 0.0;
		array_1d<double,3> dd_i;
		array_1d<double,3> dd_i1;
		array_1d<double,3> difference;
		for(ModelPart::NodesContainerType::iterator it = rInterfaceNodes.begin(); it!=rInterfaceNodes.end(); it++)
		{
			noalias(dd_i) = it->GetValue(RELAXED_DISP);
			noalias(dd_i) -= it->GetValue(DISPLACEMENT);
			
			noalias(dd_i1) = it->FastGetSolutionStepValue(RELAXED_DISP);
			noalias(dd_i1) -= it->FastGetSolutionStepValue(DISPLACEMENT);

			noalias(difference) = dd_i;
			noalias(difference) -= dd_i1;
			
			num += inner_prod(difference,dd_i1);
			denom += inner_prod(difference,difference);
			
/*rel_it += inner_prod(it->GetValue(RELAXED_DISP),it->GetValue(RELAXED_DISP));
rel += inner_prod(it->FastGetSolutionStepValue(RELAXED_DISP),it->FastGetSolutionStepValue(RELAXED_DISP));*/
		}
/*KRATOS_WATCH(rel);		
KRATOS_WATCH(rel_it);*/
		double new_mu = old_mu + (old_mu - 1.0)*num/denom;
		
//		KRATOS_WATCH(num);
// 		KRATOS_WATCH(denom);
		
		return new_mu;
		KRATOS_CATCH("");
	}
// 	double ComputeAitkenFactor(ModelPart::NodesContainerType& rInterfaceNodes, const double old_mu)
// 	{
// 		KRATOS_TRY
// 				
// 		double num = 0.0;
// 		double denom = 0.0;
// 		array_1d<double,3> dd_i;
// 		array_1d<double,3> dd_i1;
// 		array_1d<double,3> difference;
// 		array_1d<double,3> aaa;
// 		for(ModelPart::NodesContainerType::iterator it = rInterfaceNodes.begin(); it!=rInterfaceNodes.end(); it++)
// 		{
// 			noalias(dd_i) = it->GetValue(RELAXED_DISP);
// 			noalias(dd_i) -= it->GetValue(DISPLACEMENT);
// 			
// 			noalias(dd_i1) = it->FastGetSolutionStepValue(RELAXED_DISP);
// 			noalias(dd_i1) -= it->FastGetSolutionStepValue(DISPLACEMENT);
// 			
// 			noalias(aaa) = it->FastGetSolutionStepValue(RELAXED_DISP);
// 			noalias(aaa) -= it->GetValue(RELAXED_DISP);
// 			
// 			noalias(difference) = dd_i;
// 			noalias(difference) -= dd_i1;
// 			
// 			num += inner_prod(difference,aaa);
// 			denom += inner_prod(difference,difference);
// 		}
// 		
// 		double omega = num/denom;
// 		
// 		KRATOS_WATCH("aaaa");
// 		KRATOS_WATCH(num);
// 		KRATOS_WATCH(denom);
// 		
// 		return omega;
// 		KRATOS_CATCH("");
// 	}		
	//***********************************************************************
	//***********************************************************************
	//save Variables needed for Aitken in the Non Historical Database
	void ComputeRelaxedDisplacement(ModelPart::NodesContainerType& rInterfaceNodes, const double omega)
	{
		KRATOS_TRY
				
/*double rel_before = 0.0;
for(ModelPart::NodesContainerType::iterator it = rInterfaceNodes.begin(); it!=rInterfaceNodes.end(); it++)
{
	array_1d<double,3>& relaxed = it->FastGetSolutionStepValue(RELAXED_DISP);
	rel_before += inner_prod(relaxed,relaxed);
}*/
		

		for(ModelPart::NodesContainerType::iterator it = rInterfaceNodes.begin(); it!=rInterfaceNodes.end(); it++)
		{
			array_1d<double,3>& relaxed = it->FastGetSolutionStepValue(RELAXED_DISP);
//	KRATOS_WATCH(relaxed);		
			noalias(relaxed) = omega * it->FastGetSolutionStepValue(DISPLACEMENT);
			noalias(relaxed) +=  (1.0 - omega) * it->GetValue(RELAXED_DISP);
//			KRATOS_WATCH(relaxed);
//			KRATOS_WATCH(" ")
		}
		
/*		double rel_after = 0.0;
		double rel_saved = 0.0;
		for(ModelPart::NodesContainerType::iterator it = rInterfaceNodes.begin(); it!=rInterfaceNodes.end(); it++)
		{
			array_1d<double,3>& relaxed = it->FastGetSolutionStepValue(RELAXED_DISP);
			array_1d<double,3>& rel_s = it->GetValue(RELAXED_DISP);
			rel_after += inner_prod(relaxed,relaxed);
			rel_saved += inner_prod(rel_s,rel_s);
		}
std::cout << "before= "<< rel_before << " after = " << rel_after << " saved=" << rel_saved << std::endl;*/
		
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

#endif /* KRATOS_AITKEN_UTILS  defined */

