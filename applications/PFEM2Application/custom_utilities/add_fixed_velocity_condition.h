

#if !defined(KRATOS_ADD_FIXED_VELOCITY_CONDITION_UTILITY_INCLUDED )
#define KRATOS_ADD_FIXED_VELOCITY_CONDITION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "custom_conditions/fixed_velocity_2d.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"


namespace Kratos
{
	//this class is to be modified by the user to customize the interpolation process
	//template< unsigned int TDim>
	class AddFixedVelocityCondition2D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddFixedVelocityCondition2D);

		AddFixedVelocityCondition2D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 2d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~AddFixedVelocityCondition2D()
		{}

		
		void AddThem() 
		{
			KRATOS_TRY
			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number += lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if ((inode->IsFixed(VELOCITY_X)) || (inode->IsFixed(VELOCITY_Y)) || (inode->IsFixed(VELOCITY_Z)) )
				{
					Condition const& rReferenceCondition = KratosComponents<Condition>::Get("FixedVelocity2D");         //condition type
					Point2D<Node > geometry(Node::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node::Pointer( *inode.base() ));
					Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
					Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding conditions on fixed velocity boundaries" << condition_number << std::endl;
			
			
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:
		ModelPart& mr_model_part;

	};
	
	
	
	//for 3D
	class AddFixedVelocityCondition3D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddFixedVelocityCondition3D);

		AddFixedVelocityCondition3D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 3d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~AddFixedVelocityCondition3D()
		{}

		
		void AddThem() 
		{
			KRATOS_TRY
			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number += lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if ((inode->IsFixed(VELOCITY_X)) || (inode->IsFixed(VELOCITY_Y)) || (inode->IsFixed(VELOCITY_Z)) )
				{
					Condition const& rReferenceCondition = KratosComponents<Condition>::Get("FixedVelocity3D");         //condition type
					Point3D<Node > geometry(Node::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node::Pointer( *inode.base() ));
					Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
					Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding conditions on fixed velocity boundaries" << condition_number << std::endl;
			
			
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:
		ModelPart& mr_model_part;

	};
	
	
	
	
	
	class AddFixedPressureCondition2D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddFixedPressureCondition2D);

		AddFixedPressureCondition2D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Pressure 2d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~AddFixedPressureCondition2D()
		{}

		
		void AddThem() 
		{
			KRATOS_TRY
			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number += lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if (inode->IsFixed(PRESSURE) && ( (inode->IsFixed(VELOCITY_X))==false || (inode->IsFixed(VELOCITY_Y))  ) )
				{
					Condition const& rReferenceCondition = KratosComponents<Condition>::Get("FixedPressure2D");         //condition type
					Point2D<Node > geometry(Node::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node::Pointer( *inode.base() ));
					Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
					Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding conditions on fixed pressure boundaries" << condition_number << std::endl;
			
			
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:
		ModelPart& mr_model_part;

	};
	
	
	
	//for 3D
	class AddFixedPressureCondition3D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddFixedPressureCondition3D);

		AddFixedPressureCondition3D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Pressure 3d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~AddFixedPressureCondition3D()
		{}

		
		void AddThem() 
		{
			KRATOS_TRY
			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number += lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if (inode->IsFixed(PRESSURE) && ( (inode->IsFixed(VELOCITY_X))==false || (inode->IsFixed(VELOCITY_Y)) || (inode->IsFixed(VELOCITY_Z)) ) )
				{
					Condition const& rReferenceCondition = KratosComponents<Condition>::Get("FixedPressure3D");         //condition type
					Point3D<Node > geometry(Node::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node::Pointer( *inode.base() ));
					Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
					Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding conditions on fixed pressure boundaries" << condition_number << std::endl;
			
			
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:
		ModelPart& mr_model_part;

	};
	
	
	class AddWaterFixedVelocityCondition2D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddWaterFixedVelocityCondition2D);

		AddWaterFixedVelocityCondition2D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Water Fixed Velocity 2d Utility" << std::endl;
			KRATOS_CATCH("")	
		}
		

		~AddWaterFixedVelocityCondition2D()
		{}

		
		void AddThem() 
		{
			KRATOS_TRY
			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number += lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if ((inode->IsFixed(WATER_VELOCITY_X)) || (inode->IsFixed(WATER_VELOCITY_Y)) || (inode->IsFixed(WATER_VELOCITY_Z)) )
				{
					Condition const& rReferenceCondition = KratosComponents<Condition>::Get("WaterFixedVelocity2D");         //condition type
					Point2D<Node > geometry(Node::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node::Pointer( *inode.base() ));
					Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
					Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding conditions on fixed water velocity boundaries" << condition_number << std::endl;
			
			
			KRATOS_CATCH("")
		} 
		
		
	protected:


	private:
		ModelPart& mr_model_part;

	};
	

}  // namespace Kratos.

#endif // KRATOS_ADD_FIXED_VELOCITY_CONDITION_UTILITY_INCLUDED  defined 


