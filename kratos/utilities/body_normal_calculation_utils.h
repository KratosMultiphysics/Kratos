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
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2007-10-25 10:14:31 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/
#include "includes/model_part.h"

#if !defined(KRATOS_BODY_NORMAL_CALCULATION_UTILS )
#define  KRATOS_BODY_NORMAL_CALCULATION_UTILS


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
	class BodyNormalCalculationUtils
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
        void CalculateBodyNormals(
			ModelPart& r_model_part,
			int dimension)
        {
			KRATOS_TRY

			ModelPart::ElementsContainerType& rElements = r_model_part.Elements();

			//resetting the normals - only for the nodes on which we will do the calculate
			array_1d<double,3> zero = ZeroVector(3);
				
			for(ModelPart::ElementsContainerType::iterator it =  rElements.begin(); 
											it !=rElements.end(); it++)
			{
				Element::GeometryType& rNodes = it->GetGeometry();
				for(unsigned int in = 0; in<rNodes.size(); in++)
					noalias((rNodes[in]).GetSolutionStepValue(NORMAL)) = zero;
			}


			//calculating the normals and storing on the conditions
			array_1d<double,3> An;
			if(dimension == 2)
			{
				boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
				array_1d<double,3> N;
				double Volume;
				for(ModelPart::ElementsContainerType::iterator it =  rElements.begin(); it !=rElements.end(); it++)
				{				
					Element::GeometryType& geom = it->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

					for(unsigned int i = 0; i<geom.size(); i++)
					{
						array_1d<double,3>& normal = geom[i].FastGetSolutionStepValue(NORMAL);
						for(unsigned int j=0; j<2; j++)
						{	
							normal[j] += Volume*DN_DX(i,j);
						}
					}
				}
			}
			else if(dimension == 3)
			{
				boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
				array_1d<double,4> N;
				double Volume;
				for(ModelPart::ElementsContainerType::iterator it =  rElements.begin(); it !=rElements.end(); it++)
				{
					Element::GeometryType& geom = it->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

					for(unsigned int i = 0; i<geom.size(); i++)
					{
						array_1d<double,3>& normal = geom[i].FastGetSolutionStepValue(NORMAL);
						for(unsigned int j=0; j<3; j++)
						{	
							normal[j] += Volume*DN_DX(i,j);
						}
					}
				}
			}
			
			r_model_part.GetCommunicator().AssembleCurrentData(NORMAL);


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
		
        //BodyNormalCalculationUtils(void);
		
        //BodyNormalCalculationUtils(BodyNormalCalculationUtils& rSource);
		
		
		/*@} */   
		
    }; /* Class ClassName */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_BODY_NORMAL_CALCULATION_UTILS defined */

