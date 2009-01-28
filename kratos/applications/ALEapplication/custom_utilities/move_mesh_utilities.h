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
//   Date:                $Date: 2008-02-21 09:36:17 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_MESHMOVING_UTILITIES_H_INCLUDED )
#define  KRATOS_MESHMOVING_UTILITIES_H_INCLUDED



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
#include "ale_application.h"
#include "ale_application.h"

namespace Kratos
{

	class MoveMeshUtilities
	{
		public:

		MoveMeshUtilities(){};
		~MoveMeshUtilities(){};
		
		void BDF_MoveMesh(unsigned int time_order, ModelPart& rModelPart)
		{
			KRATOS_TRY
			//calculating time integration coefficients
			ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
			double dt = CurrentProcessInfo[DELTA_TIME];
		
			Vector BDFcoeffs(time_order+1);
		
			if(time_order == 2)
			{
				if(rModelPart.GetBufferSize() < 3)
					KRATOS_ERROR(std::logic_error,"insufficient buffer size for BDF2","")
	
				BDFcoeffs[0] =	1.5 / dt;	//coefficient for step n+1
				BDFcoeffs[1] =	-2.0 / dt;//coefficient for step n
				BDFcoeffs[2] =	0.5 / dt;//coefficient for step n-1
			}
			else
			{
				BDFcoeffs[0] =	1.0 / dt;	//coefficient for step n+1
				BDFcoeffs[1] =	-1.0 / dt;//coefficient for step n
			}
			
		//update nodal coordinates
			array_1d<double,3> mesh_vel;
			for(ModelPart::NodesContainerType::iterator i = rModelPart.NodesBegin();  
						 i!=rModelPart.NodesEnd(); i++) 
			{
				const array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
				i->X() = i->X0() + disp[0];
				i->Y() = i->Y0() + disp[1];
				i->Z() = i->Z0() + disp[2];
				
				//calculating the mesh velocity
				noalias(mesh_vel) = BDFcoeffs[0] * disp;
				for(unsigned int step=1; step<time_order+1; step++)
					noalias(mesh_vel) += BDFcoeffs[step]*i->FastGetSolutionStepValue(DISPLACEMENT,step);
					
				//saving the mesh velocity
				noalias(i->FastGetSolutionStepValue(MESH_VELOCITY)) = mesh_vel;
			}
			
			
			KRATOS_CATCH("");
		}
	

 



    private:
	    
	};



}  // namespace Kratos.

#endif // KRATOS_MESHMOVING_UTILITIES_H_INCLUDED  defined 


