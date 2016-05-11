//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_REFINE_MESH_BOUNDARY_PROCESS_H_INCLUDED )
#define  KRATOS_REFINE_MESH_BOUNDARY_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//Data:     DOMAIN_LABEL(nodes)(set)
//StepData: NODAL_H, NORMAL, CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked) BOUNDARY, TO_SPLIT
//          (set)     BOUNDARY(nodes), TO_ERASE(conditions), NEW_ENTITY(conditions,nodes)(set), TO_SPLIT(conditions)->locally
//          (modified)  
//          (reset)   TO_SPLIT
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Boundary Process
/** The process labels the boundary conditions (TO_SPLIT)
    Dependencies: RemoveMeshNodesProcess.Execute()  is needed as a previous step
    
    Determines if new conditions must be inserted in boundary.
    If boundary must to be kept (CONSTRAINED), 
    New conditions will be rebuild (splitting the old ones and inserting new nodes)
    Old conditions will be erased at the end.
    
*/

class RefineMeshBoundaryProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RefineMeshBoundaryProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RefineMeshBoundaryProcess(ModelPart& rModelPart,
			      ModelerUtilities::MeshingParameters& rRemeshingParameters,
			      ModelPart::IndexType MeshId,
			      int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {     
      mMeshId = MeshId;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~RefineMeshBoundaryProcess() {}


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
    virtual void Execute() 
    {
      
      KRATOS_TRY

      if( this->mEchoLevel > 0 ){
        std::cout<<" [ REFINE BOUNDARY : "<<std::endl;
	//std::cout<<"   Nodes and Conditions : "<<mrModelPart.Nodes(mMeshId).size()<<", "<<mrModelPart.Conditions(mMeshId).size()<<std::endl;
      }
      
      mrRemesh.Info->InsertedConditions    = mrModelPart.NumberOfConditions(mMeshId);
      mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes(mMeshId);


     //if the insert switches are activated, we check if the boundaries got too coarse
      if( (mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_INSERT_NODES) || mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_ADD_NODES)) && mrRemesh.Refine->RefiningOptions.Is(ModelerUtilities::REFINE_BOUNDARY) )
     {

        std::vector<Point<3> > list_of_points;
        std::vector<Condition::Pointer> list_of_conditions;

	unsigned int conditions_size = mrModelPart.Conditions(mMeshId).size();
        list_of_points.reserve(conditions_size);
        list_of_conditions.reserve(conditions_size);


        //std::cout<<"   List of Conditions Reserved Size: "<<conditions_size<<std::endl;


        //*********************************************************************************
        // REFINE BOUNDARY CONDITIONS START
        //*********************************************************************************
	RefineBoundary(list_of_points,list_of_conditions);
	//*********************************************************************************
	// REFINE BOUNDARY CONDITIONS END
	//*********************************************************************************



        //*********************************************************************************
        //                   DOFS AND NEW CONDITIONS REBUILD START                       //
        //*********************************************************************************
	BuildNewConditions(list_of_points,list_of_conditions);
        //*********************************************************************************
        //                   DOFS AND NEW CONDITIONS REBUILD END                         //
        //*********************************************************************************




        //*********************************************************************************
        //                   CLEAN CONDITIONS AND FLAGS START                            //
	//*********************************************************************************
	CleanConditionsAndFlags();
        //*********************************************************************************
        //                   CLEAN CONDITIONS AND FLAGS END                              //
        //*********************************************************************************




     } // REFINE END;



     mrRemesh.Info->InsertedConditions    = mrModelPart.NumberOfConditions(mMeshId)-mrRemesh.Info->InsertedConditions;
     mrRemesh.Info->InsertedBoundaryNodes = mrModelPart.NumberOfNodes(mMeshId)-mrRemesh.Info->InsertedBoundaryNodes;

     if( this->mEchoLevel > 0 ){
        std::cout<<"   [ CONDITIONS ( inserted : "<<mrRemesh.Info->InsertedConditions<<" ) ]"<<std::endl;
        std::cout<<"   [ NODES      ( inserted : "<<mrRemesh.Info->InsertedBoundaryNodes<<" ) ]"<<std::endl;
        std::cout<<"   REFINE BOUNDARY ]; "<<std::endl;
     }

     KRATOS_CATCH(" ")
    
    }
    

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
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
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
        return "RefineMeshBoundaryProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RefineMeshBoundaryProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{

    ///@}


protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
 
    ModelerUtilities::MeshingParameters& mrRemesh;

    ModelerUtilities mModelerUtilities;
  
    ModelPart::IndexType mMeshId; 

    int mEchoLevel;

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    void RefineBoundary(std::vector<Point<3> >& list_of_points, std::vector<Condition::Pointer>& list_of_conditions){};
  
    void BuildNewConditions( std::vector<Point<3> >& list_of_points, std::vector<Condition::Pointer>& list_of_conditions){};
  
    void CleanConditionsAndFlags(){};
  
    /// Assignment operator.
    RefineMeshBoundaryProcess& operator=(RefineMeshBoundaryProcess const& rOther);


    /// this function is a private function


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
                                  RefineMeshBoundaryProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const RefineMeshBoundaryProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_REFINE_MESH_BOUNDARY_PROCESS_H_INCLUDED  defined 


