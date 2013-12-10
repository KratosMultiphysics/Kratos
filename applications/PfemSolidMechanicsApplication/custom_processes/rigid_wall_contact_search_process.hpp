//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_RIGID_CONTACT_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_RIGID_WALL_CONTACT_SEARCH_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "includes/model_part.h"
#include "custom_conditions/point_rigid_contact_penalty_2D_condition.hpp"


namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for all processes in Kratos.
/** The process is the base class for all processes and defines a simple interface for them.
    Execute method is used to execute the Process algorithms. While the parameters of this method
  can be very different from one Process to other there is no way to create enough overridden
  versions of it. For this reason this method takes no argument and all Process parameters must
  be passed at construction time. The reason is that each constructor can take different set of
  argument without any dependency to other processes or the base Process class.
*/
class RigidWallContactSearchProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RigidWallContactSearchProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    typedef Point2D<ModelPart::NodeType>       Point2DType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidWallContactSearchProcess(ModelPart& rModelPart): mrModelPart(rModelPart) {}


    RigidWallContactSearchProcess(SpatialBoundingBox::Pointer pRigidWall, ModelPart& rModelPart) 
      : mrModelPart(rModelPart)
    {
      mpRigidWall = pRigidWall;
    } 

    /// Destructor.
    virtual ~RigidWallContactSearchProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the Process algorithms.
    virtual void Execute() {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {

      KRATOS_TRY

      ProcessInfo& CurrentProcessInfo= mrModelPart.GetProcessInfo();	  
      double Time      = CurrentProcessInfo[TIME];   

      mpRigidWall->Center() = mpRigidWall->OriginalCenter() +  mpRigidWall->Velocity() * Time;

      if (Time == 0)
	KRATOS_ERROR( std::logic_error, "detected time = 0 in the Solution Scheme ... check if the time step is created correctly for the current model part", "" )

      ModelPart::NodesContainerType& NodesArray = mrModelPart.Nodes();
      
      //Create Rigid Contact Conditions
      int MeshId = 0;
      int id = mrModelPart.Conditions(MeshId).back().Id() + 1;

      mConditionsNumber =  mrModelPart.Conditions(MeshId).size();
      

      for ( ModelPart::NodesContainerType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  if((*nd)->FastGetSolutionStepValue(RIGID_WALL)==true){

	    //(*nd)->Set(STRUCTURE);
	    (*nd)->Set(RIGID);

	    //set new coordinates (MOVE MESH)
	    (*nd)->X() = (*nd)->X0() + mpRigidWall->Velocity()[0] * Time;
            (*nd)->Y() = (*nd)->Y0() + mpRigidWall->Velocity()[1] * Time;
            (*nd)->Z() = (*nd)->Z0() + mpRigidWall->Velocity()[2] * Time;

	    //std::cout<<" node "<<(*nd)->Id()<<" Position ("<<(*nd)->X()<<", "<<(*nd)->Y()<<" "<<(*nd)->Z()<<") "<<std::endl;
	  }

	  //set point rigid wall condition : usually in non rigid_wall points

	  Vector Point(3);
	  Point[0] = (*nd)->X();
	  Point[1] = (*nd)->Y();
	  Point[2] = (*nd)->Z();

	  if( (*nd)->Is(BOUNDARY) && mpRigidWall->IsInside(Point,Time) ){
	    
	    int number_properties = mrModelPart.NumberOfProperties();

	    PropertiesType::Pointer p_properties = mrModelPart.pGetProperties(number_properties-1);

	    GeometryType::Pointer p_geometry = GeometryType::Pointer(new Point2DType( (*nd) ));

	    ConditionType::Pointer p_cond = ModelPart::ConditionType::Pointer(new PointRigidContactPenalty2DCondition(id, p_geometry, p_properties, mpRigidWall) ); 
	    //pcond->SetValue(mpRigidWall); the boundingbox of the rigid wall must be passed to the condition

	    mrModelPart.Conditions(MeshId).push_back(p_cond);

	    id +=1;
	  }
	  
	}
      
      std::cout<<" [ rigid contacts : "<<mrModelPart.Conditions(MeshId).size() - mConditionsNumber<<" ]"<<std::endl;

      KRATOS_CATCH( "" )
	
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
      KRATOS_TRY	  

      // To write correct displacements
      ProcessInfo& CurrentProcessInfo= mrModelPart.GetProcessInfo();
      double Time      = CurrentProcessInfo[TIME];
      //double DeltaTime = CurrentProcessInfo[DELTA_TIME];

      ModelPart::NodesContainerType& NodesArray = mrModelPart.Nodes();

      for ( ModelPart::NodesContainerType::ptr_iterator nd = NodesArray.ptr_begin(); nd != NodesArray.ptr_end(); ++nd)
	{
	  if((*nd)->Is(RIGID)){
	    array_1d<double, 3 > &CurrentDisplacement  = (*nd)->FastGetSolutionStepValue(DISPLACEMENT);	    
	    CurrentDisplacement[0] = mpRigidWall->Velocity()[0] * Time;
	    CurrentDisplacement[1] = mpRigidWall->Velocity()[1] * Time;
	    CurrentDisplacement[2] = mpRigidWall->Velocity()[2] * Time;
	  }
	}


      //Clean Rigid Contact Conditions
      ModelPart::ConditionsContainerType NonRigidContactConditions;

      int MeshId = 0;
      unsigned int id=0;
	    
      //std::cout<<" [ NUMBER OF CONDITIONS before rigid contact update: "<<mrModelPart.Conditions(MeshId).size()<<" ]"<<std::endl;
    
      for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(MeshId); ic!= mrModelPart.ConditionsEnd(MeshId); ic++)
	{
	  if( id == mConditionsNumber )
	    break;
	  
	  NonRigidContactConditions.push_back(*(ic.base()));  

	  id +=1;
	}
    
      mrModelPart.Conditions(MeshId).swap( NonRigidContactConditions );

      //std::cout<<" [ NUMBER OF CONDITIONS after  rigid contact update: "<<mrModelPart.Conditions(MeshId).size()<<" ]"<<std::endl;

      //calculate elemental contribution
      KRATOS_CATCH( "" )      
    }


    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
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
    virtual std::string Info() const
    {
        return "RigidWallContactSearchProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RigidWallContactSearchProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Static Member Variables
    ///@{
    ModelPart&  mrModelPart;

    SpatialBoundingBox::Pointer mpRigidWall;

    unsigned int mConditionsNumber;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RigidWallContactSearchProcess& operator=(RigidWallContactSearchProcess const& rOther);

    /// Copy constructor.
    //Process(Process const& rOther);


    ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  RigidWallContactSearchProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RigidWallContactSearchProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_RIGID_CONTACT_SEARCH_PROCESS_H_INCLUDED  defined 


