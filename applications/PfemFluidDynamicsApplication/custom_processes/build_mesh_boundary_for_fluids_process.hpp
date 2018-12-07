//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:          AFranci $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:     September 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_BUILD_MESH_BOUNDARY_FOR_FLUIDS_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_MESH_BOUNDARY_FOR_FLUIDS_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_processes/build_model_part_boundary_process.hpp"

///VARIABLES used:
//Data:     MASTER_ELEMENTS(set), MASTER_NODES(set)
//StepData:
//Flags:    (checked) TO_ERASE, TO_REFINE, CONTACT, NEW_ENTITY
//          (set)     BOUNDARY(set),  [TO_REFINE(nodes), TO_ERASE(condition)]->locally to not preserve condition
//          (modified)
//          (reset)
// (set):=(set in this process)

namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef  ModelPart::NodesContainerType NodesContainerType;
  typedef  ModelPart::ElementsContainerType ElementsContainerType;
  typedef  ModelPart::ConditionsContainerType ConditionsContainerType;

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
   */
  class BuildMeshBoundaryForFluidsProcess
    : public BuildModelPartBoundaryProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildMeshBoundaryForFluidsProcess
    KRATOS_CLASS_POINTER_DEFINITION( BuildMeshBoundaryForFluidsProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuildMeshBoundaryForFluidsProcess(ModelPart& rModelPart,
			     MesherUtilities::MeshingParameters& rRemeshingParameters,
			     int EchoLevel)
      : BuildModelPartBoundaryProcess(rModelPart, rModelPart.Name(), EchoLevel),
	mrRemesh(rRemeshingParameters)
    {

    }

    /// Destructor.
    virtual ~BuildMeshBoundaryForFluidsProcess()
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
      KRATOS_TRY

      bool success=false;

      boost::timer auxiliary;

      if( mEchoLevel > 0 )
	std::cout<<" [ Build Boundary on ModelPart ["<<mrModelPart.Name()<<"] ]"<<std::endl;

      success=UniqueSkinSearch(mrModelPart);

      if(!success)
	{
          std::cout<<"  ERROR:  BOUNDARY BUILD FAILED ModelPart : ["<<mrModelPart<<"] "<<std::endl;
        }
      else
	{
          if( mEchoLevel >= 1 )
	    std::cout<<" [ Search performed in Time = "<<auxiliary.elapsed()<<" ]"<<std::endl;
            //PrintSkin(mrModelPart);
        }


      KRATOS_CATCH(" ")
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
      return "BuildMeshBoundaryForFluidsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "BuildMeshBoundaryForFluidsProcess";
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
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{



    //**************************************************************************
    //**************************************************************************



     
    bool UniqueSkinSearch( ModelPart& rModelPart ) 
    { 
 
      KRATOS_TRY 
 
  if( mEchoLevel > 0 ){ 
    std::cout<<" [ Initial Conditions : "<<rModelPart.Conditions().size()<<std::endl; 
  } 
 
      if( !rModelPart.Elements().size() || (rModelPart.Is(ACTIVE)) ){ 
  if( mEchoLevel > 0 ){ 
    std::cout<<" [ Final Conditions   : "<<rModelPart.Conditions().size()<<std::endl; 
  } 
  return true; 
      } 
       
      //reset the boundary flag in all nodes and check if a remesh process has been performed 
      bool any_node_to_erase = false; 
      for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); in++) 
  { 
    in->Reset(BOUNDARY); 
    in->Reset(FREE_SURFACE); 
 
    if( any_node_to_erase == false ) 
      if( in->Is(TO_ERASE) ) 
        any_node_to_erase = true; 
       
  } 
      SetBoundaryAndFreeSurface(rModelPart);

      return true; 
 
      KRATOS_CATCH( "" ) 
	} 
     



      
    // bool SetBoundaryAndFreeSurface( ModelPart& rModelPart, ModelPart::ConditionsContainerType& rTemporaryConditions, std::vector<int>& rPreservedConditions, unsigned int& rConditionId )
    bool SetBoundaryAndFreeSurface( ModelPart& rModelPart)
    {

      KRATOS_TRY

	
      //properties to be used in the generation
      int number_properties = rModelPart.GetParentModelPart()->NumberOfProperties();
      Properties::Pointer properties = rModelPart.GetParentModelPart()->pGetProperties(number_properties-1);

      ModelPart::ElementsContainerType::iterator elements_begin  = mrModelPart.ElementsBegin();
      ModelPart::ElementsContainerType::iterator elements_end    = mrModelPart.ElementsEnd();

      for(ModelPart::ElementsContainerType::iterator ie = elements_begin; ie != elements_end ; ie++)
	{

	  Geometry< Node<3> >& rElementGeometry = ie->GetGeometry();

	  if( rElementGeometry.FacesNumber() >= 3 ){ //3 or 4

	    /*each face is opposite to the corresponding node number so in 2D triangle
	      0 ----- 1 2
	      1 ----- 2 0
	      2 ----- 0 1
	    */

	    /*each face is opposite to the corresponding node number so in 3D tetrahedron
	      0 ----- 1 2 3
	      1 ----- 2 0 3
	      2 ----- 0 1 3
	      3 ----- 0 2 1
	    */

	    //finding boundaries and creating the "skin"
	    //
	    //********************************************************************

	    boost::numeric::ublas::matrix<unsigned int> lpofa; //connectivities of points defining faces
	    boost::numeric::ublas::vector<unsigned int> lnofa; //number of points defining faces

	    WeakPointerVector<Element >& rE = ie->GetValue(NEIGHBOUR_ELEMENTS);

	    //get matrix nodes in faces
	    rElementGeometry.NodesInFaces(lpofa);
	    rElementGeometry.NumberNodesInFaces(lnofa);

	    //loop on neighbour elements of an element
	    unsigned int iface=0;
	    for(WeakPointerVector< Element >::iterator ne = rE.begin(); ne!=rE.end(); ne++)
	      {

		unsigned int NumberNodesInFace = lnofa[iface];

		if (ne->Id() == ie->Id())
		  {

		    //if no neighbour is present => the face is free surface
		    bool freeSurfaceFace=false;
		    for(unsigned int j=1; j<=NumberNodesInFace; j++)
		      {
			rElementGeometry[lpofa(j,iface)].Set(BOUNDARY);
			if(rElementGeometry[lpofa(j,iface)].IsNot(RIGID)){
			  freeSurfaceFace=true;
			}
		      }
		    if(freeSurfaceFace==true){
		      for(unsigned int j=1; j<=NumberNodesInFace; j++)
			{
			  rElementGeometry[lpofa(j,iface)].Set(FREE_SURFACE);
			}
		    }


		  } //end face condition

		iface+=1;

	      } //end loop neighbours
	  }
	}

      return true;

      KRATOS_CATCH( "" )
    }




    //**************************************************************************
    //**************************************************************************

    bool CheckAcceptedCondition(ModelPart& rModelPart, Condition& rCondition) override
    {
      KRATOS_TRY

      bool node_not_preserved = false;
      bool condition_not_preserved = false;

      Geometry< Node<3> >& rConditionGeometry = rCondition.GetGeometry();

      for(unsigned int j=0; j<rConditionGeometry.size(); j++)
	{
	  if( rConditionGeometry[j].Is(TO_ERASE) || rConditionGeometry[j].Is(TO_REFINE) )
	    node_not_preserved = true;

	  if( rConditionGeometry[j].Is(ISOLATED) || rConditionGeometry[j].IsNot(BOUNDARY) )
	    condition_not_preserved = true;
	}

      if( rCondition.Is(TO_ERASE) )
	condition_not_preserved = true;

      if( rCondition.Is(BOUNDARY) ) //flag for composite condition
	condition_not_preserved = true;

      if(node_not_preserved == true || condition_not_preserved == true)
	return false;
      else
	return true;

      KRATOS_CATCH( "" )
    }


    //**************************************************************************
    //**************************************************************************


    void AddConditionToModelPart(ModelPart& rModelPart, Condition::Pointer pCondition) override
    {
      KRATOS_TRY

      //rModelPart.AddCondition(pCondition); //if a Local Id corresponds to a Global Id not added
      rModelPart.Conditions().push_back(pCondition);

      KRATOS_CATCH( "" )
    }



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
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    MesherUtilities::MeshingParameters& mrRemesh;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    //**************************************************************************
    //**************************************************************************


    bool FindNodeInCondition(Geometry< Node<3> >& rConditionGeometry,Geometry< Node<3> >& rElementGeometry , boost::numeric::ublas::matrix<unsigned int>& lpofa, boost::numeric::ublas::vector<unsigned int>& lnofa, unsigned int& iface)
    {
      KRATOS_TRY

      // not equivalent geometry sizes for boundary conditions:
      if( rConditionGeometry.size() != lnofa[iface] )
	return false;

      // line boundary condition:
      if( lnofa[iface] == 2 )
	{
	  if( rConditionGeometry[0].Id() == rElementGeometry[lpofa(1,iface)].Id()  ||
	      rConditionGeometry[1].Id() == rElementGeometry[lpofa(2,iface)].Id()  ||
	      rConditionGeometry[0].Id() == rElementGeometry[lpofa(2,iface)].Id()  ||
	      rConditionGeometry[1].Id() == rElementGeometry[lpofa(1,iface)].Id()  )
	  {
	    return true;
	  }
	  else
	  {
	    return false;
	  }

	}

      //3D faces:
      if(  lnofa[iface] == 3 )
	{
	  if( rConditionGeometry[0].Id() == rElementGeometry[lpofa(1,iface)].Id() ||
	      rConditionGeometry[1].Id() == rElementGeometry[lpofa(2,iface)].Id() ||
	      rConditionGeometry[2].Id() == rElementGeometry[lpofa(3,iface)].Id() ||
	      rConditionGeometry[0].Id() == rElementGeometry[lpofa(3,iface)].Id() ||
	      rConditionGeometry[1].Id() == rElementGeometry[lpofa(1,iface)].Id() ||
	      rConditionGeometry[2].Id() == rElementGeometry[lpofa(2,iface)].Id() ||
	      rConditionGeometry[0].Id() == rElementGeometry[lpofa(2,iface)].Id() ||
	      rConditionGeometry[1].Id() == rElementGeometry[lpofa(3,iface)].Id() ||
	      rConditionGeometry[2].Id() == rElementGeometry[lpofa(1,iface)].Id()  )
	  {
	    return true;
	  }
	  else
	  {
	    return false;
	  }

	}

      if(  lnofa[iface] > 3 )
	{
	  KRATOS_THROW_ERROR( std::logic_error, "Wrong Condition Number of Face Nodes",*this );
	}

      return false;

      KRATOS_CATCH(" ")
    }

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
    BuildMeshBoundaryForFluidsProcess& operator=(BuildMeshBoundaryForFluidsProcess const& rOther);

    /// Copy constructor.
    //BuildMeshBoundaryForFluidsProcess(BuildMeshBoundaryForFluidsProcess const& rOther);


    ///@}

  }; // Class BuildMeshBoundaryForFluidsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    BuildMeshBoundaryForFluidsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const BuildMeshBoundaryForFluidsProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.



#endif // KRATOS_BUILD_MESH_BOUNDARY_FOR_FLUIDS_PROCESS_H_INCLUDED  defined
