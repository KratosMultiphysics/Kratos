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
*   Last Modified by:    $Author: pooyan $
*   Date:                $Date: 2008-11-13 12:12:17 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_NORMAL_CALCULATION_UTILS )
#define  KRATOS_NORMAL_CALCULATION_UTILS


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
	class NormalCalculationUtils
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
		//this function calculates the "area normal" (vector oriented as the normal 
		//with a dimension proportional to the area. This is done on the base of the 
		//Conditions provided
		//which should be intended as the surface elements of the area of interest.
        void CalculateOnSimplex(
			ConditionsArrayType& rConditions,
			int dimension)
        {
			KRATOS_TRY

			//resetting the normals
			array_1d<double,3> zero = Vector(3);
			noalias(zero) = ZeroVector(3);
				
			for(ConditionsArrayType::iterator it =  rConditions.begin(); 
											it !=rConditions.end(); it++)
			{
				Element::GeometryType& rNodes = it->GetGeometry();
				for(unsigned int in = 0; in<rNodes.size(); in++)
					noalias((rNodes[in]).GetSolutionStepValue(NORMAL)) = zero;
			}


			//calculating the normals and storing on the conditions
			array_1d<double,3> An;
			if(dimension == 2)
			{
				for(ConditionsArrayType::iterator it =  rConditions.begin(); 
											it !=rConditions.end(); it++)
				{
					CalculateNormal2D(it,An);
				}
			}
			else if(dimension == 3)
			{
				array_1d<double,3> v1;
				array_1d<double,3> v2;
				for(ConditionsArrayType::iterator it =  rConditions.begin(); 
											it !=rConditions.end(); it++)
				{
					//calculate the normal on the given condition
					CalculateNormal3D(it,An,v1,v2);
				}
			}

			//adding the normals to the nodes
			for(ConditionsArrayType::iterator it =  rConditions.begin(); 
							it !=rConditions.end(); it++)
			{
				Geometry<Node<3> >& pGeometry = (it)->GetGeometry();
				double coeff = 1.00/pGeometry.size();
				const array_1d<double,3>& An = it->GetValue(NORMAL);

				for(unsigned int i = 0; i<pGeometry.size(); i++)
				{
					noalias(pGeometry[i].FastGetSolutionStepValue(NORMAL)) += coeff * An;
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
		static void CalculateNormal2D(ConditionsArrayType::iterator it, array_1d<double,3>& An)
		{
				Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

				An[0] =    pGeometry[1].Y() - pGeometry[0].Y();
				An[1] = - (pGeometry[1].X() - pGeometry[0].X());
				An[2] =    0.00;

 				array_1d<double,3>& normal = (it)->GetValue(NORMAL);
 				noalias(normal) = An;

// 				(it)->SetValue(NORMAL,An);
		}	

		static void CalculateNormal3D(ConditionsArrayType::iterator it, array_1d<double,3>& An,
			array_1d<double,3>& v1,array_1d<double,3>& v2 )
		{
				Geometry<Node<3> >& pGeometry = (it)->GetGeometry();

				v1[0] = pGeometry[1].X() - pGeometry[0].X();
				v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
				v1[2] = pGeometry[1].Z() - pGeometry[0].Z();
						
				v2[0] = pGeometry[2].X() - pGeometry[0].X();
				v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
				v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

				MathUtils<double>::CrossProduct(An,v1,v2);
				An *= 0.5;

 				array_1d<double,3>& normal = (it)->GetValue(NORMAL);
 				noalias(normal) = An;
// 				noalias((it)->GetValue(NORMAL)) = An;
		}		
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
		
        //NormalCalculationUtils(void);
		
        //NormalCalculationUtils(NormalCalculationUtils& rSource);
		
		
		/*@} */   
		
    }; /* Class ClassName */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_NORMAL_CALCULATION_UTILS  defined */

