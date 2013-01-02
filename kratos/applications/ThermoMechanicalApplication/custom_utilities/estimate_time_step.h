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


#if !defined(KRATOS_ESTIMATE_TIME_STEP )
#define  KRATOS_ESTIMATE_TIME_STEP



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
#include "thermo_mechanical_application.h"
// #include "custom_conditions/environment_contact.h"
//#include "includes/variables.h"



namespace Kratos
{
 	

	class EstimateTimeStep
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		//
		//template<unsigned int TDim>
		double ComputeDt(ModelPart& ThisModelPart, const double dt_min, const double dt_max  )
		{			
		  KRATOS_TRY

		  array_1d<double, 3 > dx, dv;
		  double deltatime = dt_max;
		  double dvside, lside;

		  for (ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
			  i != ThisModelPart.ElementsEnd(); i++)
		  {
		      //calculating shape functions values
		      Geometry< Node < 3 > >& geom = i->GetGeometry();

		      for (unsigned int i1 = 0; i1 < geom.size() - 1; i1++)
		      {
			  for (unsigned int i2 = i1 + 1; i2 < geom.size(); i2++)
			  {
			      dx[0] = geom[i2].X() - geom[i1].X();
			      dx[1] = geom[i2].Y() - geom[i1].Y();
			      dx[2] = geom[i2].Z() - geom[i1].Z();

			      lside = inner_prod(dx, dx);

			      noalias(dv) = geom[i2].FastGetSolutionStepValue(VELOCITY);
			      noalias(dv) -= geom[i1].FastGetSolutionStepValue(VELOCITY);

			      dvside = inner_prod(dx, dv);

			      double dt;
			      if (dvside < 0.0) //otherwise the side is "getting bigger" so the time step has not to be diminished
			      {
				  dt = fabs(lside / dvside);
				  if (dt < deltatime) deltatime = dt;
			      }

			  }
		      }
		  }

		  if (deltatime < dt_min)
		  {
		      std::cout << "ATTENTION dt_min is being used" << std::endl;
		      deltatime = dt_min;
		  }

		  return deltatime;

		  KRATOS_CATCH("")	
		}



	private:


	};

}  // namespace Kratos.

#endif // ASSIGN_NO_SLIP_CONDITION  defined 


