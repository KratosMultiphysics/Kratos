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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_LAGRANGIAN_UTILITIES_INCLUDED )
#define  KRATOS_LAGRANGIAN_UTILITIES_INCLUDED



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
#include "geometries/tetrahedra_3d_4.h"
#include "PFEM_application.h"

namespace Kratos
{
	class LagrangianUtils
	{
	public:


		//**********************************************************************************************
		//**********************************************************************************************
		void ExplicitLagrangianPrediction(ModelPart& rModelPart)
		{
			KRATOS_TRY;

			array_1d<double, 3> delta_disp;
			double dt =  rModelPart.GetProcessInfo()[DELTA_TIME];
			ModelPart::NodesContainerType rNodes = rModelPart.Nodes();

			for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
			{
				if(i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
				{
					noalias(delta_disp) = i->FastGetSolutionStepValue(VELOCITY);
					delta_disp *= dt;
			
					noalias(i->FastGetSolutionStepValue(DISPLACEMENT)) += delta_disp;
				}
			}
			
			KRATOS_CATCH("");   
		}

		//**********************************************************************************************
		//**********************************************************************************************
		// TAKE CARE! this is an eulerian step => the mesh velocity should be previously set to zero
		void ImplicitLagrangianPrediction(ModelPart& rModelPart)
		{
			KRATOS_TRY;

			array_1d<double, 3> delta_disp;
			double dt =  rModelPart.GetProcessInfo()[DELTA_TIME];
			ModelPart::NodesContainerType rNodes = rModelPart.Nodes();

			for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
			{
				if(i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
				{				noalias(delta_disp) = i->FastGetSolutionStepValue(FRACT_VEL);
				delta_disp *= dt;

				noalias(i->FastGetSolutionStepValue(DISPLACEMENT)) += delta_disp;
				}
			}
			
			KRATOS_CATCH("");   
		}

		//**********************************************************************************************
		//**********************************************************************************************
		//calculates lagrangian displacements corrsponding to the first step of the fractional step solution
		//assumes the FractionalStepVelocity at the old iteration to be saved in VAUX
		void CalculateStep1DisplacementCorrection(ModelPart& rModelPart)
		{
			KRATOS_TRY;

			array_1d<double, 3> delta_disp;
			double dt =  rModelPart.GetProcessInfo()[DELTA_TIME];
			ModelPart::NodesContainerType rNodes = rModelPart.Nodes();

			for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
			{
				if(i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
				{
				noalias(delta_disp) = i->FastGetSolutionStepValue(FRACT_VEL);
				noalias(delta_disp) -= i->GetValue(VAUX);
				delta_disp *= dt;

				noalias(i->FastGetSolutionStepValue(DISPLACEMENT)) += delta_disp;
				}
			}
			
			KRATOS_CATCH("");   
		}
	
		//**********************************************************************************************
		//**********************************************************************************************

		//calculates lagrangian displacements corrsponding to the first step of the fractional step solution
		//assumes the FractionalStepVelocity at the old iteration to be saved in VAUX
		void CalculateFinalDisplacementCorrection(ModelPart& rModelPart)
		{
			KRATOS_TRY;

			array_1d<double, 3> delta_disp;
			double dt =  rModelPart.GetProcessInfo()[DELTA_TIME];
			ModelPart::NodesContainerType rNodes = rModelPart.Nodes();

			for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
			{
				if(i->FastGetSolutionStepValue(IS_STRUCTURE) == 0)
				{
					noalias(delta_disp) = i->FastGetSolutionStepValue(VELOCITY);
					noalias(delta_disp) -= i->GetValue(VAUX);
					delta_disp *= dt;

					noalias(i->FastGetSolutionStepValue(DISPLACEMENT)) += delta_disp;
				}
			}
			
			KRATOS_CATCH("");   
		}
	private:

	};

}  // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_UTILITIES_INCLUDED  defined 


