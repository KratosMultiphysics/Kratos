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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_FLUID_INTERPOLATOR_INCLUDED )
#define  KRATOS_FLUID_INTERPOLATOR_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "custom_utilities/custom_interpolator.h"


namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	class FluidInterpolator : public CustomInterpolator
	{
	public:

		KRATOS_CLASS_POINTER_DEFINITION(FluidInterpolator);

		FluidInterpolator(ModelPart& model_part)
			: CustomInterpolator(model_part)
		{}

		~FluidInterpolator()
		{}


		virtual void InterpolateNode(Node<3>::Pointer new_node, 
			const Geometry< Node<3> >& geom,
			const Vector& N)
		{
			KRATOS_TRY
				
			//add dofs
			new_node->pAddDof(PRESSURE);
			new_node->pAddDof(VELOCITY_X);
			new_node->pAddDof(VELOCITY_Y);
			new_node->pAddDof(VELOCITY_Z);
			new_node->pAddDof(FRACT_VEL_X);
			new_node->pAddDof(FRACT_VEL_Y);
			new_node->pAddDof(FRACT_VEL_Z);

			//interpolate scalar variables
			std::vector< const Variable<double>* > double_vars;
			double_vars.push_back(&PRESSURE);
			double_vars.push_back(&NODAL_H);
			AdaptivityUtils::Interpolate<double>(new_node,geom,N,double_vars,mr_model_part.GetBufferSize());

			//interpolate vector variables
			std::vector<const Variable<array_1d<double,3> >* > vector_vars;
			vector_vars.push_back(&VELOCITY);
			AdaptivityUtils::Interpolate< array_1d<double,3> >(new_node,geom,N,vector_vars,mr_model_part.GetBufferSize());

			KRATOS_CATCH("")
		}


	protected:


	private:



	};

}  // namespace Kratos.

#endif // KRATOS_FLUID_INTERPOLATOR_INCLUDED  defined 


