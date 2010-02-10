/*
==============================================================================
KratosIncompressibleFluidApplication 
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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(ASSIGN_NO_SLIP_CONDITION )
#define  ASSIGN_NO_SLIP_CONDITION



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point.h"
#include "incompressible_fluid_application.h"
#include "custom_conditions/no_slip_condition_2d.h"
//#include "includes/variables.h"



namespace Kratos
{
 	

	class AssignNoSlipCondition
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		//
		//template<unsigned int TDim>
		void AssignNoSlipCondition2D(ModelPart& ThisModelPart)
		{			
			KRATOS_TRY;


			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			int id = ThisModelPart.Conditions().size();
                        //const char* ConditionName = "NoSlipCondition2D";

			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			 {
				if( it->FastGetSolutionStepValue(IS_BOUNDARY) == 1.0)
				   {	
	      KRATOS_WATCH(">>>>>>>>>>>>>>>>>>>>> NOSLIPBOUNDARY «««««««««««««««««««««");	
				       Condition::NodesArrayType temp;	
				       temp.reserve(1);
				       temp.push_back(*(it.base()));	
		
									
			Condition::Pointer p_cond = (KratosComponents<Condition>::Get("NoSlipCondition2D")).Create(id, temp, properties);	
				
					//Condition::Pointer p_cond = Condition::Pointer(new NoSlipCondition2D(id, temp, properties) );
					(ThisModelPart.Conditions()).push_back(p_cond);

					id++;
				  }

			}
			ThisModelPart.Conditions().Sort();

			

			KRATOS_CATCH("")
		}



	private:


	};

}  // namespace Kratos.

#endif // ASSIGN_NO_SLIP_CONDITION  defined 


