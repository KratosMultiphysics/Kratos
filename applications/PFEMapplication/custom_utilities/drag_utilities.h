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
	
	///@name Kratos Globals 
	///@{
	
	
	///@} 
	///@name Type Definitions       
	///@{
	
	
	///@}
	
	
	///@name  Enum's */       
	///@{
	
	
	///@}
	///@name  Functions */       
	///@{
	
	
	///@}
	///@name Kratos Classes */
	///@{
	
	/// This class allows the calculation of the seepage drag that the fluid transmits to the solid porous matrix.
	  /** @author  Antonia Larese De Tetto <antoldt@cimne.upc.edu>
	  * 
	  * This class includes some usefull utilities for the calculation of the drag that water
	  * makes on the solid porous matrix in a coupled problem.
	  * This utilities are used in the pfem_nonnewtonian_coupled gid problem type 
	  * where an edgebased model is used for the calculation of evolution of the fluid @see edgebased_levelset.h
	  * and it is coupled with a pfem non newtonian visco rigid model is used to avaluate the structural response  @see nonewtonian_asgs_2d.h, nonewtonian_asgs_3d.h
	  */
	
	class DragUtils
         {
          public:
		///@name Type Definitions      
		///@{ 
		typedef ModelPart::NodesContainerType NodesArrayType; 
		typedef ModelPart::ElementsContainerType ElementsArrayType;
		///@} 
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

		
		
		/// Provides the global indices for each one of this element's local rows
		/**
		* Seepage drag is evaluated in the fluid edgebased model part. Ready to be projected to the structural pfem mesh.
		* 
		* @param rFluid_ModelPart: The fluid model part where Seepage drag has been evaluated for the current time step @see nonewtonian_asgs_2d.h, nonewtonian_asgs_3d.h
		* @param rSeepageDragVar: The seepage drag vector that is evaluated in the fluid model part
		* @param mu: Fluid viscosity
		* @param fluid_density: Fluid density
		* @param solid_density: Solid density

		*/
		void CalculateFluidDrag(
			ModelPart& rFluid_ModelPart,
			Variable< array_1d<double,3> >& rSeepageDragVar,
			const double & fluid_nu ,
			const double & fluid_density,
			const double & solid_density
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
				  
				  //dimensionally accelerations. A, B not multiplied by porosity like in the fluid momentum equation
				  double lin_coeff = 150 * fluid_nu * (1-porosity)*(1-porosity)/(porosity*porosity*porosity*diameter*diameter);
				  double non_lin_coeff = 1.75 * (1-porosity)/(porosity*porosity*porosity*diameter);
				  noalias(darcy_term) = (lin_coeff + non_lin_coeff * vel_norm) * vel * fluid_density / solid_density;

			      }
			}

			KRATOS_CATCH("")

		}

		/// Adds the drag to the body force vector of the structural model
		/**
		* After projecting the seepage drag on the structural pfem mesh
		* body force vector of the structural model part is constructed with the gravity  and the seepage_drag contributions 
		* 
		* @param rStructural_ModelPart: The structural model part
		* @param rSeepageDragVar: the seepage drag nodal vector
		* @param rBodyForce: the body force vector that will be reset and calculated as the sum of the gravity contribution and the seepage drag contribution	
		* @param rGravityAcc: the gravity acceleration vector
		*/
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
// 			  if(inode->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
// 			  {
			    const double& porosity = inode->FastGetSolutionStepValue(POROSITY);
			    const array_1d<double,3>& press_grad = inode->FastGetSolutionStepValue(PRESS_PROJ);
			    const double& solid_density = inode->FastGetSolutionStepValue(DENSITY);

  			    array_1d<double,3>& bodyforce_drag = inode->FastGetSolutionStepValue(rBodyForce);
			    bodyforce_drag  = inode->FastGetSolutionStepValue(rSeepageDragVar);

			    ////subtracting buoyancy force
			    if(inode->FastGetSolutionStepValue(WATER_PRESSURE) > 0.0) //TO IDENTIFY FLUID NODES...we don't have any information over distance
			    {
				noalias(bodyforce_drag) -= (1-porosity) * press_grad / solid_density;
// 			    	array_1d<double,3> temp  = (1-porosity) * press_grad / solid_density;
// 			    	KRATOS_WATCH(porosity);
// 			    	KRATOS_WATCH(press_grad);
// 			    	KRATOS_WATCH(solid_density);
			    }
			    
			    bodyforce_drag  +=  rGravityAcc;
// 			    KRATOS_WATCH(bodyforce_drag)
// 			  } 
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
			

		///@} 
		///@name Private Operators*/
		///@{ 
			

		///@} 
		///@name Private Operations*/
		///@{ 
			

		///@} 
		///@name Private  Acces */
		///@{ 
			

		///@}     
		///@name Private Inquiry */
		///@{ 
			

		///@}   
		///@name Un accessible methods */
		///@{ 
			

		///@}   
		
    }; /* Class ClassName */
	
	///@} 
	
	/**@name Type Definitions */       
	///@{ 
	
	
	///@}
	
}  /* namespace Kratos.*/

#endif /* KRATOS_DRAG_UTILS  defined */

