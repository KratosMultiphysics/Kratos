//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Pavel Ryzhakov


#if !defined(KRATOS_ADD_MEAN_VEL_CONDITION_PROCESS_INCLUDED )
#define  KRATOS_ADD_MEAN_VEL_CONDITION_PROCESS_INCLUDED



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
//#include "custom_elements/updated_lagrangian_fluid.h"
//#include "custom_elements/updated_lagrangian_fluid3D.h"
//#include "custom_elements/updated_lagrangian_fluid_inc.h"
//#include "custom_elements/updated_lagrangian_fluid3D_inc.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "includes/ublas_interface.h"
#include "includes/variables.h" 
#include "includes/model_part.h"
#include "includes/condition.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"

#include "ULF_application.h"
#include "custom_conditions/lagrange_multiplier_mean_velocity_2D.h"
#include "custom_conditions/lagrange_multiplier_mean_velocity_3D.h"



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

class AddMeanVelocityConditionProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PushStructureProcess
    KRATOS_CLASS_POINTER_DEFINITION(AddMeanVelocityConditionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AddMeanVelocityConditionProcess(ModelPart& model_part): mr_model_part(model_part) 
    //ModelPart& fluid_model_part, ModelPart& structure_model_part, ModelPart& combined_model_part)
    //: mr_fluid_model_part(fluid_model_part), mr_structure_model_part(structure_model_part), mr_combined_model_part(combined_model_part)
    {
	Check();
	//KRATOS_WATCH(" INSIDE ADD WALL NODES CONSTRUCTOR") 
    }

    /// Destructor.
    ~AddMeanVelocityConditionProcess() //override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    //	void operator()()
    //	{
    //		MergeParts();
    //	}


    ///@}
    ///@name Operations
    ///@{

 
    //works on a cube. Note that we dont prescribe periodics on the corner nodes!
    void AddMeanVelCond3D(const double x_larger_than ,const double x_smaller_than, const double y_larger_than, const double y_smaller_than, const double z_larger_than, const double z_smaller_than  ) 
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
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if( inode->X()>=x_larger_than &&
					inode->X()<=x_smaller_than &&
					inode->Y()>=y_larger_than && 
					inode->Y()<=y_smaller_than &&
					inode->Z()>=z_larger_than && 
					inode->Z()<=z_smaller_than )
				{
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("MeanVelocityLagrangeMultiplierCondition3D");         //condition type
				Point3D<Node<3> > geometry(Node<3>::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node<3>::Pointer( *inode.base() ));
				Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
				Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
				p_condition->GetValue(NEIGHBOUR_NODES).push_back(Node<3>::WeakPointer( *(mr_model_part.NodesBegin()).base() ));
				mr_model_part.Conditions().push_back(p_condition);
				++condition_number;
				}
			}
			std::cout << "Finished adding Lagrange multiplier for mean velocity to all nodes!" << condition_number << std::endl;
			


			//now setting the nodal areas:
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				inode->GetSolutionStepValue(NODAL_AREA)=0.0;
			}
			//here we consider that the nodal_areas of all the elements are approximately equal (domain_size/number of elements)
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); 
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				Geometry< Node<3> >& geom = ielem->GetGeometry();
				const double area = geom.DomainSize(); 
				const double nodal_area = area/double(geom.size());
				for(unsigned int i=0; i<geom.size(); i++)
				{
					geom[i].FastGetSolutionStepValue(NODAL_AREA)+=nodal_area;
				}
				
			}

			
			KRATOS_CATCH("")
		} 

	//////////////////////////////////////////////////////////////////////////////////////
	void AddMeanVelCond2D(const double x_larger_than ,const double x_smaller_than, const double y_larger_than, const double y_smaller_than  ) 
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
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				if( inode->X()>=x_larger_than &&
					inode->X()<=x_smaller_than &&
					inode->Y()>=y_larger_than && 
					inode->Y()<=y_smaller_than )
				{
				Condition const& rReferenceCondition = KratosComponents<Condition>::Get("MeanVelocityLagrangeMultiplierCondition2D");         //condition type
				Point2D<Node<3> > geometry(Node<3>::Pointer( *inode.base() ));//mr_model_part.Nodes().(inode->Id()));//Node<3>::Pointer( *inode.base() ));
				Properties::Pointer properties = mr_model_part.GetMesh().pGetProperties(0); 		//this will allow us later to turn this layer on/off in GID
				Condition::Pointer p_condition = rReferenceCondition.Create(condition_number, geometry, properties); 
				p_condition->GetValue(NEIGHBOUR_NODES).push_back(Node<3>::WeakPointer( *(mr_model_part.NodesBegin()).base() ));
				mr_model_part.Conditions().push_back(p_condition);
				++condition_number;
				}
			}
			std::cout << "Finished adding lagrange multiplier to all nodes!" << condition_number << std::endl;
			


			//now setting the nodal areas:
			for(ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); 
				inode!=mr_model_part.NodesEnd(); inode++)
			{
				inode->GetSolutionStepValue(NODAL_AREA)=0.0;
			}
			for(ModelPart::ElementsContainerType::iterator ielem = mr_model_part.ElementsBegin(); 
				ielem!=mr_model_part.ElementsEnd(); ielem++)
			{
				Geometry< Node<3> >& geom = ielem->GetGeometry();
				const double area = geom.DomainSize(); 
				const double nodal_area = area/double(geom.size());
				for(unsigned int i=0; i<geom.size(); i++)
				{
					geom[i].FastGetSolutionStepValue(NODAL_AREA)+=nodal_area;
				}
				
			}

			
			KRATOS_CATCH("")
		} 

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "AddMeanVelocityConditionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AddMeanVelocityConditionProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
	void Check()
		{
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NORMAL) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data","");
			//if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NEIGHBOUR_NODES) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NEIGHBOUR_NODES variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data","");
			if(mr_model_part.NodesBegin()->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing LAGRANGE_MULTIPLIER_VELOCITY variable 				on solution step data","");
		}
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
		ModelPart& mr_model_part;
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    //ModelPart& mr_fluid_model_part;
    //ModelPart& mr_structure_model_part;
    //ModelPart& mr_combined_model_part;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//		AddMeanVelocityConditionProcess& operator=(AddMeanVelocityConditionProcess const& rOther);

    /// Copy constructor.
//		AddMeanVelocityConditionProcess(AddMeanVelocityConditionProcess const& rOther);


    ///@}

}; // Class AddMeanVelocityConditionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AddMeanVelocityConditionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AddMeanVelocityConditionProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ADD_MEAN_VEL_CONDITION_PROCESS_INCLUDED  defined 


