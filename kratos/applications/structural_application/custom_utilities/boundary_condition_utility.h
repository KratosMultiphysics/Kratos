/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2007-09-20 10:41:02 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_BOUNDARY_CONDITION_UTILITY_INCLUDED )
#define  KRATOS_BOUNDARY_CONDITION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/variables.h"


namespace Kratos
{
	/**
	 * Deactivation and reactivation of elements
	 * This process handles deactivation and reactivation
	 * of elements and associated conditions.
	 * In order to use this process, a variable calles
	 * ACTIVATION_LEVEL of type int has to be defined.
	 * The following values can be assigned to ACTIVATION_LEVEL:
	 * ACTIVATION_LEVEL == 0: element is always active
	 * ACTIVATION_LEVEL > 0: element will be deactivated on demand
	 * ACTIVATION_LEVEL < 0: element is initially deactivated and
	 *                       will be reactivated on demand
	 */
	class BoundaryConditionUtility
	{
		public:
			
			typedef PointerVectorSet<Element> ElementsArrayType;
			typedef PointerVectorSet<Node<3> > NodesArrayType;
			typedef PointerVectorSet<Condition> ConditionsArrayType;
			typedef Geometry<Node<3> > GeometryType;
			/**
			 * class pointer definition
			 */
			KRATOS_CLASS_POINTER_DEFINITION( BoundaryConditionUtility );
			
			/**
			 * Constructor.
			 * The constructor takes the current model_part as argument.
			 * Please note that reactivation of elements does only work
			 * as long as the process that deactivated the elements before
			 * is living.
			 */
			BoundaryConditionUtility()
			{
			}
			
			/**
			 * Destructor.
			 */
			virtual ~BoundaryConditionUtility()
			{
			}
			
			/**
			 * Initializes all elements before performing any calculation.
			 * This is done even for those elements that are deactivated
			 * in the beginning of the simulation
			 */
			void Initialize( ModelPart& model_part )
			{
				//initializing elements
				for( ModelPart::NodeIterator it=model_part.NodesBegin();
								 it !=model_part.NodesEnd(); ++it )
				{
					it->AddDof(DISPLACEMENT_X);
					it->AddDof(DISPLACEMENT_Y);
					it->AddDof(DISPLACEMENT_Z);
				}
			}
			
			/**
			 * applies boundary conditions to a tunnelling simulation model
			 */
			void ApplyBoundaryConditions( ModelPart& model_part, 
										  double x_inf, double x_sup,
										  double y_inf, double y_sup,
										  double z_inf, double z_sup
										)
			{
				std::cout << "Applying boundary conditions..." << std::endl;
				for( ModelPart::NodeIterator it=model_part.NodesBegin();
								 it !=model_part.NodesEnd(); ++it )
				{
					if( it->HasDofFor(DISPLACEMENT) )
						std::cout << "has dof for disp." << std::endl;
					if( it->HasDofFor(DISPLACEMENT_X) )
						std::cout << "has dof for X" << std::endl;
					if( it->HasDofFor(DISPLACEMENT_Y) )
						std::cout << "has dof for Y" << std::endl;
					if( it->HasDofFor(DISPLACEMENT_Z) )
						std::cout << "has dof for Z" << std::endl;
					if( fabs(it->GetInitialPosition()[0]-x_inf) < 1.0e-2 )
					{
// 						if( (*it)->HasDofFor(DISPLACEMENT_X) )
// 						{
							std::cout << "fixing x" << std::endl;
							it->Fix(DISPLACEMENT_X);
// 						}
// 						else
// 							std::cout << "DOF x does not exist" << std::endl;
					}
					if( fabs(it->GetInitialPosition()[0]-x_sup) < 1.0e-2 )
					{
// 						if( (*it)->HasDofFor(DISPLACEMENT_X) )
// 						{
							std::cout << "fixing x" << std::endl;
							it->Fix(DISPLACEMENT_X);
// 						}
// 						else
// 							std::cout << "DOF x does not exist" << std::endl;
					}
					if( fabs(it->GetInitialPosition()[1]-y_inf) < 1.0e-2 )
					{
// 						if( (*it)->HasDofFor(DISPLACEMENT_Y) )
// 						{
							std::cout << "fixing y" << std::endl;
							it->Fix(DISPLACEMENT_Y);
// 						}
// 						else
// 							std::cout << "DOF y does not exist" << std::endl;
					}
					if( fabs(it->GetInitialPosition()[1]-y_sup) < 1.0e-2 )
					{
// 						if( (*it)->HasDofFor(DISPLACEMENT_Y) )
// 						{
							std::cout << "fixing y" << std::endl;
							it->Fix(DISPLACEMENT_Y);
// 						}
// 						else
// 							std::cout << "DOF y does not exist" << std::endl;
					}
					if( fabs(it->GetInitialPosition()[2]-z_inf) < 1.0e-2 )
					{
// 						if( (*it)->HasDofFor(DISPLACEMENT_Z) )
// 						{
							std::cout << "fixing z" << std::endl;
							it->Fix(DISPLACEMENT_Z);
// 						}
// 						else
// 							std::cout << "DOF z does not exist" << std::endl;
					}
				}
			}
			
			
		
		private:
			
			/**
			 * Assignment operator
			 */
 			//BoundaryConditionUtility& operator=(BoundaryConditionUtility const& rOther);
			
			/**
			 * Copy constructor
			 */
			//BoundaryConditionUtility(BoundaryConditionUtility const& rOther);
	
	};//class BoundaryConditionUtility
}  // namespace Kratos.

#endif //KRATOS_BOUNDARY_CONDITION_UTILITY_INCLUDED  defined 
