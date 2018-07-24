//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pavel Ryzhakov and Alex Jarauta

#if !defined(ASSIGN_POINT_NEUMANN3D_CONDITION )
#define  ASSIGN_POINT_NEUMANN3D_CONDITION



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"

#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
//#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point.h"
#include "ULF_application.h"
#include "custom_conditions/Point_Neumann3D.h"
#include "custom_conditions/Point_Neumann2D.h"
//#include "custom_conditions/Point_Neumann_Monolithic2D.h"
//#include "includes/variables.h"



namespace Kratos
{
 	

	class AssignPointNeumannConditions
	{
	public:

		//**********************************************************************************************
		//**********************************************************************************************
		//
		//template<unsigned int TDim>
		void AssignPointNeumannConditionsDisp(ModelPart& ThisModelPart)
		{			
			KRATOS_TRY;
			//KRATOS_WATCH("Inside of AssignPointNeumannConditions UTILITY")
			const int TDim=ThisModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			int id = ThisModelPart.Conditions().size();
			//KRATOS_WATCH(ThisModelPart.Conditions().size())
			Condition::NodesArrayType temp;	
			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			 {
				if( it->FastGetSolutionStepValue(FLAG_VARIABLE) == 1.0)
				   {	
	      			       temp.reserve(1);
				       temp.push_back(*(it.base()));	
		
					if (TDim==3)
						{
						Condition::Pointer p_cond = (KratosComponents<Condition>::Get("PointNeumann3D")).Create(id, temp, properties);	
						(ThisModelPart.Conditions()).push_back(p_cond);
						}
					else if (TDim==2)
						{
						Condition::Pointer p_cond = (KratosComponents<Condition>::Get("PointNeumann2D")).Create(id, temp, properties);	
						(ThisModelPart.Conditions()).push_back(p_cond);
						
						}
					
					id++;
					temp.clear();
					
				  }

			}

			//ThisModelPart.Conditions().Sort();


			KRATOS_CATCH("")
		}

		void AssignPointNeumannConditionsDispAxisym(ModelPart& ThisModelPart)
		{
			KRATOS_TRY;
			KRATOS_WATCH("Inside of AssignPointNeumannConditions UTILITY")
			const int TDim=ThisModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			int id = ThisModelPart.Conditions().size();
			//KRATOS_WATCH(ThisModelPart.Conditions().size())
			Condition::NodesArrayType temp;	
			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			 {
			 	
			       	temp.reserve(1);
			       	temp.push_back(*(it.base()));	

				if( it->FastGetSolutionStepValue(FLAG_VARIABLE) == 1.0)
				   {	
	      			       	Condition::Pointer p_cond = (KratosComponents<Condition>::Get("PointNeumannAxisym")).Create(id, temp, properties);	
					(ThisModelPart.Conditions()).push_back(p_cond);
					
					
					id++;
					temp.clear();
					
				  }

			}
			//KRATOS_WATCH(ThisModelPart.Conditions().size())
			//KRATOS_WATCH("BEFORE SORTING")
			//ThisModelPart.Conditions().Sort();
			//KRATOS_WATCH(ThisModelPart.Conditions().size())
			

			KRATOS_CATCH("")

		}

		void AssignPointNeumannConditionsVel(ModelPart& ThisModelPart)
		{			
			KRATOS_TRY;
			//KRATOS_WATCH("Inside of AssignPointNeumannConditions UTILITY")
			const int TDim=ThisModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			int id = ThisModelPart.Conditions().size();
			//KRATOS_WATCH(ThisModelPart.Conditions().size())
			Condition::NodesArrayType temp;	
			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			 {
				if( it->FastGetSolutionStepValue(FLAG_VARIABLE) == 1.0)
				   {	
	      			       temp.reserve(1);
				       temp.push_back(*(it.base()));	
		
					if (TDim==3)
						{
						Condition::Pointer p_cond = (KratosComponents<Condition>::Get("PointNeumann3D_vel")).Create(id, temp, properties);	
						(ThisModelPart.Conditions()).push_back(p_cond);
						}
					else if (TDim==2)
						{
						Condition::Pointer p_cond = (KratosComponents<Condition>::Get("PointNeumann2D_vel")).Create(id, temp, properties);	
						(ThisModelPart.Conditions()).push_back(p_cond);
						
						}
					
					id++;
					temp.clear();
					
				  }

			}

			//ThisModelPart.Conditions().Sort();


			KRATOS_CATCH("")
		}
		 //THIS ONE IS FOR THE VELOCITY-BASED FORMUALTION
		/* IT WORKS, BUT IS COMMENTED
		void AssignPointNeumannConditionsMonolithic2D(ModelPart& ThisModelPart)
		{			
			KRATOS_TRY;
			KRATOS_WATCH("INSide of AssignPointNeumannConditions UTILITY Point_Neumann_Monolithic2D")
			const int TDim=ThisModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			int id = ThisModelPart.Conditions().size();
                        //const char* ConditionName = "NoSlipCondition2D";
			KRATOS_WATCH(ThisModelPart.Conditions().size())
			for(ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin(); 
				it!=ThisModelPart.NodesEnd(); it++)
			 {
				if( it->FastGetSolutionStepValue(FLAG_VARIABLE) == 1.0)
				   {	
	      			       //KRATOS_WATCH(">>>>>>>>>>>>>>>>>>>>> ASSIGNING THE NEUMANN POINT CONDITIONS -> EXTERNAL PRESSURE WILL BE APPLIED <<<<<<<<<<<<<<<<<<<<<");
				       Condition::NodesArrayType temp;	
				       temp.reserve(1);
				       temp.push_back(*(it.base()));	
		
					if (TDim==3)
						{
						KRATOS_ERROR(std::logic_error,"Method not implemented for 3D","");
						}
					else if (TDim==2)
						{
						Condition::Pointer p_cond = (KratosComponents<Condition>::Get("PointNeumannMonolithic2D")).Create(id, temp, properties);	
						(ThisModelPart.Conditions()).push_back(p_cond);
						}
					
					id++;
					
				  }

			}
KRATOS_WATCH(ThisModelPart.Conditions().size())
			//KRATOS_WATCH("BEFORE SORTING")
			ThisModelPart.Conditions().Sort();
KRATOS_WATCH(ThisModelPart.Conditions().size())
			

			KRATOS_CATCH("")
		}

		*/
	private:


	};

}  // namespace Kratos.

#endif // ASSIGN_POINT_NEUMANN3D_CONDITION  defined 
