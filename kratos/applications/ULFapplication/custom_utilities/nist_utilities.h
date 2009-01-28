/*
==============================================================================
KratosULFApplication 
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
//   Date:                $Date: 2007-03-07 15:44:25 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_NIST_UTILITIES_INCLUDED )
#define  KRATOS_NIST_UTILITIES_INCLUDED



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
#include "ULF_application.h"

namespace Kratos
{
	class NistUtils
	{
	public:


		//**********************************************************************************************
		//**********************************************************************************************
		void GenerateModelPart(
			ModelPart& OriginModelPart , 
			ModelPart& DestinationModelPart,
			Element const& rReferenceElement, 
			Condition const& rReferenceBoundaryCondition
			)
		{
			KRATOS_TRY;

			//assigning the nodes to the new model part
			DestinationModelPart.Nodes().clear();
			DestinationModelPart.Nodes() = OriginModelPart.Nodes();
			
			//generating the elements
			int id = 1;
			Properties::Pointer properties = OriginModelPart.GetMesh().pGetProperties(1);			
			for(ModelPart::ElementsContainerType::iterator iii = OriginModelPart.ElementsBegin(); iii != OriginModelPart.ElementsEnd(); iii++)
			{
				Geometry< Node<3> >& geom = iii->GetGeometry();
				Element::Pointer p_element = rReferenceElement.Create(id, geom ,properties);
				DestinationModelPart.Elements().push_back(p_element);
				id = id + 1;
			}
			std::cout << "Elements are generated" << std::endl;

			//generating the conditions
			id = 1;
			for(ModelPart::ConditionsContainerType::iterator iii = OriginModelPart.ConditionsBegin(); iii != OriginModelPart.ConditionsEnd(); iii++)
			{
				Geometry< Node<3> >& geom = iii->GetGeometry();
				double nfree_surf = 0;
				for(unsigned int k = 0; k<geom.size(); k++)
					nfree_surf += geom[k].FastGetSolutionStepValue(IS_FREE_SURFACE);

				if(nfree_surf > 1)
				{
					Condition::Pointer p_condition = rReferenceBoundaryCondition.Create(id, geom,properties);
						DestinationModelPart.Conditions().push_back(p_condition);
						id = id + 1;
				}
			}
			std::cout << "Conditions are generated" << std::endl;
			
			KRATOS_CATCH("");   
		}


		//**********************************************************************************************
		//**********************************************************************************************
		void ApplyInitialTemperature(
			ModelPart& ThisModelPart , 
			double wall_temperature
			)
		{
			KRATOS_TRY;
			for(ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); 
				in!=ThisModelPart.NodesEnd(); in++)
			{
				if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 1 &&
				  (in->GetValue(NEIGHBOUR_ELEMENTS)).size() == 0 )
				{
					in->FastGetSolutionStepValue(TEMPERATURE) = wall_temperature;
				}
			}

			KRATOS_CATCH("");   
		}

		//**********************************************************************************************
		//**********************************************************************************************
		double FindFluidLevel(
			ModelPart::NodesContainerType nodes
			)
		{
			KRATOS_TRY;

			double level = 0.00;
			for(ModelPart::NodesContainerType::iterator in = nodes.begin(); 
				in!=nodes.end(); in++)
			{
				if(  (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0 )
				{
					if( in->Y() > level) level = in->Y();
				}
			}
			return level;

			KRATOS_CATCH("");   
		}

	private:

	};

}  // namespace Kratos.

#endif // KRATOS_NIST_UTILITIES_INCLUDED  defined 


