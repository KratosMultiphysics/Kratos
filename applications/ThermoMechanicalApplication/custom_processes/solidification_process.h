//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $
//
//  this process save structural elements in a separate list

#if !defined(KRATOS_SOLIDIFICATION_PROCESS )
#define  KRATOS_SOLIDIFICATION_PROCESS



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
#include "thermo_mechanical_application.h"
namespace Kratos
{

	///@name Kratos Globals
	///@{

	///@}
	///@name Type Definitions
	///@{


	///@}
	///@name  Enum's
	///@{

	///@}
	///@name  Functions
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
		Update the PRESSURE_FORCE on the nodes


	*/

	class SolidificationProcess
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(SolidificationProcess);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
// 		DuplicateInterfaceNodesCreateConditionsProcess()
// 		{
// 		}
	       SolidificationProcess(ModelPart& ThisModelPart,const double solidificarion_temperature )
		:Process(), mr_model_part(ThisModelPart), mr_solid_temp(solidificarion_temperature)
		{
		}

		/// Destructor.
		virtual ~SolidificationProcess()
		{
		}


		///@}
		///@name Operators
		///@{

		void operator()()
		{
		  Execute();
		}

		///@}
		///@name Operations
		///@{

	   void Execute() override
		 {
	            KRATOS_WATCH("INSIDE SolidificationProcess");

		   //loop over nodes to assign zero fixe velocity base on solidification temperature
		    for(ModelPart::NodesContainerType::iterator ind = mr_model_part.NodesBegin();
			    ind!=mr_model_part.NodesEnd(); ind++)
		    {
		      const double TT = ind->FastGetSolutionStepValue(TEMPERATURE);
		      const double dist = ind->FastGetSolutionStepValue(DISTANCE);
		      if( TT <= mr_solid_temp && dist<0.0 ){
			ind->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
			ind->Fix(VELOCITY_X);
			ind->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
			ind->Fix(VELOCITY_Y);
			ind->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
			ind->Fix(VELOCITY_Z);
		      }
		    }

		}//end of execute




		private:
			  ModelPart& mr_model_part;
			  const double mr_solid_temp;

		//functions



	};//end of class


}//end of namespace Kratos

#endif // KRATOS_DUPLICATE_INTERFACE_NODES_CREATE_CONDITIONS_PROCESS  defined


