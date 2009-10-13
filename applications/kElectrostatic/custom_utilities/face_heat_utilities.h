/*
==============================================================================
KratosR1ConvectionDiffusionApplication 
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


#if !defined(KRATOS_FACE_HEAT_UTILITIES_INCLUDED )
#define  KRATOS_FACE_HEAT_UTILITIES_INCLUDED



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
#include "kElectrostatic.h"

namespace Kratos
{
	class FaceHeatUtilities
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		void ApplyFaceHeat(
			ModelPart::ConditionsContainerType& conditions , 
			double face_heat_source
			)
		{
			KRATOS_TRY
			for(ModelPart::ConditionsContainerType::iterator iii = conditions.begin(); iii != conditions.end(); iii++)
			{
				Geometry< Node<3> >& geom = iii->GetGeometry();

				for(unsigned int k = 0; k<geom.size(); k++)
					geom[k].FastGetSolutionStepValue(FACE_HEAT_FLUX) = face_heat_source;
				
			}
			std::cout << "Conditions are generated" << std::endl;
			
			KRATOS_CATCH("");		
		}



	private:

	};

}  // namespace Kratos.

#endif // KRATOS_FACE_HEAT_UTILITIES_INCLUDED  defined 


