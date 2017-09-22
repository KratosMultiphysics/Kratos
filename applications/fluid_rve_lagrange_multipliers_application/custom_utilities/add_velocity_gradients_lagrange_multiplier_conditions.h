

#if !defined(KRATOS_ADD_VELOCITY_GRADIENTS_LAGRANGE_MULTIPLIER_CONDITION_UTILITY_INCLUDED )
#define KRATOS_ADD_VELOCITY_GRADIENTS_LAGRANGE_MULTIPLIER_CONDITION_UTILITY_INCLUDED

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "custom_conditions/lagrange_multiplier_velocity_gradients_2d.h"
#include "fluid_rve_lagrange_multipliers_application.h"
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
	class AddVelocityGradientsLagrangeMultiplierConditions2D
	{
	public:
	
		KRATOS_CLASS_POINTER_DEFINITION(AddVelocityGradientsLagrangeMultiplierConditions2D);

		AddVelocityGradientsLagrangeMultiplierConditions2D(ModelPart& model_part)
			: mr_model_part(model_part) 
		{
			KRATOS_TRY	
			//std::cout << "Hello, I am the constructor of the Fixed Velocity 2d Utility" << std::endl;
			Check();		
			KRATOS_CATCH("")	
		}
		

		~AddVelocityGradientsLagrangeMultiplierConditions2D()
		{}

		
		void AddThem(const double x_larger_than ,const double x_smaller_than, const double y_larger_than, const double y_smaller_than  ) 
		{
			KRATOS_TRY
			
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number = lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); 
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				//Condition const& rReferenceCondition = KratosComponents<Condition>::Get("VelocityGradientsLagrangeMultiplierCondition2D");         //condition type
				//Point2D<Node<3> > geometry(Node<3>::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node<3>::Pointer( *inode.base() ));
				//Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
				//Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
				
				//we only add conditions if all the nodes Ys are higher than "y_larger_than":
				Geometry<Node<3> >& geom = ielem->GetGeometry();
				if( geom[0].X()>=x_larger_than &&  geom[1].X()>=x_larger_than && geom[2].X()>=x_larger_than	 &&
					geom[0].X()<=x_smaller_than &&  geom[1].X()<=x_smaller_than && geom[2].X()<=x_smaller_than &&
					geom[0].Y()>=y_larger_than &&  geom[1].Y()>=y_larger_than && geom[2].Y()>=y_larger_than  &&
					geom[0].Y()<=y_smaller_than &&  geom[1].Y()<=y_smaller_than && geom[2].Y()<=y_smaller_than )
				{
					Condition::Pointer p_condition = Condition::Pointer(new VelocityGradientsLagrangeMultiplierCondition2D(
						  condition_number,
						  (*ielem).pGetGeometry(),
						  (*ielem).pGetProperties() ) );
					p_condition->GetValue(NEIGHBOUR_NODES).push_back(Node<3>::WeakPointer( *(mr_model_part.NodesBegin()).base() ));     
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding lagrange multipliers for the velocity gradients (ELEMENT BASED)" << condition_number << std::endl;
			


			
			KRATOS_CATCH("")
		} 
		
		
		void AddThem2(const double x_larger_than ,const double x_smaller_than, const double y_larger_than, const double y_smaller_than , const double input_y_division ) 
		{
			KRATOS_TRY
			
			const double y_division = input_y_division-1.0e-4;
			unsigned int condition_number=1;
			if(mr_model_part.Conditions().size()!=0)
			{
				ModelPart::ConditionsContainerType::iterator lastcondition = (mr_model_part.ConditionsEnd()-1);
				condition_number = lastcondition->Id();
			}
			
			//Condition const& rReferenceCondition = FixedVelocity2D();         //condition type
			//getting data for the given geometry
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); 
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				//Condition const& rReferenceCondition = KratosComponents<Condition>::Get("VelocityGradientsLagrangeMultiplierCondition2D");         //condition type
				//Point2D<Node<3> > geometry(Node<3>::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node<3>::Pointer( *inode.base() ));
				//Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
				//Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
				
				//we only add conditions if all the nodes Ys are higher than "y_larger_than":
				Geometry<Node<3> >& geom = ielem->GetGeometry();
				if( geom[0].X()>=x_larger_than &&  geom[1].X()>=x_larger_than && geom[2].X()>=x_larger_than	 &&
					geom[0].X()<=x_smaller_than &&  geom[1].X()<=x_smaller_than && geom[2].X()<=x_smaller_than &&
					geom[0].Y()>=y_larger_than &&  geom[1].Y()>=y_larger_than && geom[2].Y()>=y_larger_than  &&
					geom[0].Y()<=y_smaller_than &&  geom[1].Y()<=y_smaller_than && geom[2].Y()<=y_smaller_than )
				{
					int multiplier_node = 1;
					if(geom[0].Y()>y_division && geom[1].Y()>y_division && geom[2].Y()>y_division)
						multiplier_node = 0;
					
					Condition::Pointer p_condition = Condition::Pointer(new VelocityGradientsLagrangeMultiplierCondition2D(
						  condition_number,
						  (*ielem).pGetGeometry(),
						  (*ielem).pGetProperties() ) );
					p_condition->GetValue(NEIGHBOUR_NODES).push_back(Node<3>::WeakPointer( *(mr_model_part.NodesBegin()+multiplier_node).base()  ));     
					mr_model_part.Conditions().push_back(p_condition);
					++condition_number;
				}
			}
			std::cout << "Finished adding lagrange multipliers for the velocity gradients (ELEMENT BASED)" << condition_number << std::endl;
			


			
			KRATOS_CATCH("")
		} 
		
		
		
	protected:

		void Check()
		{
			//if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NORMAL) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data","");
			//if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS variable on solution step data","");
		}

	private:
		ModelPart& mr_model_part;

	};
	
	
	
	

}  // namespace Kratos.

#endif // KRATOS_MEAN_VELOCITY_GRADIENTS_LAGRANGE_MULTIPLIER_CONDITION_UTILITY_INCLUDED


