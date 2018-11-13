//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_BUILD_MESH_BOUNDARY_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_MESH_BOUNDARY_PROCESS_H_INCLUDED


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
  class BuildMeshBoundaryProcess
    : public BuildModelPartBoundaryProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildMeshBoundaryProcess
    KRATOS_CLASS_POINTER_DEFINITION( BuildMeshBoundaryProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuildMeshBoundaryProcess(ModelPart& rModelPart,
			     ModelerUtilities::MeshingParameters& rRemeshingParameters,
			     int EchoLevel)
      : BuildModelPartBoundaryProcess(rModelPart, rModelPart.Name(), EchoLevel),
	mrRemesh(rRemeshingParameters)
    { 

    }

    /// Destructor.
    virtual ~BuildMeshBoundaryProcess()
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

    virtual void Execute()
    {
		KRATOS_TRY

		bool success=false;

		boost::timer auxiliary;

		if( mEchoLevel > 0 )
			std::cout<<" [ Build Boundary on ModelPart ["<<mrModelPart.Name()<<"] ]"<<std::endl;

		success=this->UniqueSkinSearch(mrModelPart);

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
    virtual std::string Info() const
    {
      return "BuildMeshBoundaryProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "BuildMeshBoundaryProcess";
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

    bool BuildCompositeConditions( ModelPart& rModelPart, ModelPart::ConditionsContainerType& rTemporaryConditions, std::vector<int>& rPreservedConditions, unsigned int& rConditionId )
    {

      KRATOS_TRY

      //master conditions must be deleted and set them again in the build
      this->ClearMasterEntities(rModelPart, rTemporaryConditions);
      
      //properties to be used in the generation
      int number_properties = rModelPart.GetParentModelPart()->NumberOfProperties();
      Properties::Pointer properties = rModelPart.GetParentModelPart()->pGetProperties(number_properties-1);

      ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
      
      ModelPart::ElementsContainerType::iterator elements_begin  = mrModelPart.ElementsBegin();
      ModelPart::ElementsContainerType::iterator elements_end    = mrModelPart.ElementsEnd();
    
      rConditionId=0;
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
	    
	    //Get the standard ReferenceCondition
	    const Condition & rReferenceCondition = mrRemesh.GetReferenceCondition();
	    
	    //loop on neighbour elements of an element
	    unsigned int iface=0;
	    for(WeakPointerVector< Element >::iterator ne = rE.begin(); ne!=rE.end(); ne++)
	      {

		unsigned int NumberNodesInFace = lnofa[iface];
		
		if (ne->Id() == ie->Id())
		  {
		    
		    //if no neighbour is present => the face is free surface
		    for(unsigned int j=1; j<=NumberNodesInFace; j++)
		      {
			rElementGeometry[lpofa(j,iface)].Set(BOUNDARY);
		      }
	
	
		    //Get the correct ReferenceCondition
		    Condition::Pointer pBoundaryCondition;
		    bool condition_found = false;
		    bool point_condition = false; 

		    bool inserted = false;
		    for(ModelPart::ConditionsContainerType::iterator ic = rTemporaryConditions.begin(); ic!= rTemporaryConditions.end(); ic++)
		      {
			Geometry< Node<3> >& rConditionGeometry = ic->GetGeometry();

			if( ic->IsNot(TO_ERASE) ){

			  if( ic->IsNot(CONTACT) ){
			  
			    if(ic->Is(NEW_ENTITY)){
			      inserted = false;
			    }
			    else{
			      // remeshing rebuild
			      for( unsigned int i=0; i<rConditionGeometry.size(); i++ )
				{
				  if( rConditionGeometry[i].Is(TO_ERASE)){
				    inserted = true;
				    break;
				  }
				}
			      // remeshing rebuild
			    }

			    if( !inserted ){
			      
			      if( rPreservedConditions[ic->Id()-1] == 0 ){

				ModelerUtilities ModelerUtils;
				condition_found = ModelerUtils.FindCondition(rConditionGeometry,rElementGeometry,lpofa,lnofa,iface);
				
				if( condition_found ){
				
				  pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
				  rPreservedConditions[ic->Id()-1] += 1; //add each time is used

				  if( rConditionGeometry.PointsNumber() == 1 )
				    point_condition = true;

				  //break;
				}
			      }
			    
			    }
			    else{
			    
			      if( rPreservedConditions[ic->Id()-1] < 2 ){
				
				condition_found = this->FindNodeInCondition(rConditionGeometry,rElementGeometry,lpofa,lnofa,iface);
				
				if( condition_found ){
				
				  pBoundaryCondition = (*(ic.base())); //accessing boost::shared_ptr  get() to obtain the raw pointer
				  rPreservedConditions[ic->Id()-1] += 1; //add each time is used

				  if( rConditionGeometry.PointsNumber() == 1 )
				    point_condition = true;
				  
				  //break;
				}	

				// std::cout<<" INSERTED COND "<<ic->Id()<<std::endl;
			      }			
			  
			    }
			    			    			    			    
			  }
			  else{
			    
			    rPreservedConditions[ic->Id()-1] += 1;  //will not be restored
			    //std::cout<<" Condition Contact "<<ic->Id()<<std::endl;

			  }
			
			}
			  

			if(condition_found==true){
			  // std::cout<<" Condition Found:  "<<ic->Id()<<" ("<<ic->GetGeometry()[0].Id()<<", "<<ic->GetGeometry()[1].Id()<<") == ("<<rGeom[lpofa(1,i)].Id()<<" "<<rGeom[lpofa(2,i)].Id()<<") ->  Used : "<<rPreservedConditions[ic->Id()-1]<<" times "<<std::endl;
			  break;
			}
		      }


		    // Set new conditions:  start
		    if( !point_condition ){

		      //1.- create geometry: points array and geometry type
		      
		      Condition::NodesArrayType        FaceNodes;
		      Condition::GeometryType::Pointer ConditionVertices;
		      
		      FaceNodes.reserve(NumberNodesInFace);

		      for(unsigned int j=1; j<=NumberNodesInFace; j++)
			{
			  FaceNodes.push_back(rElementGeometry(lpofa(j,iface)));
			}

		      rConditionId +=1;

		      //Create a condition
		      Condition::Pointer p_cond;
		      if(condition_found){
			
			p_cond = pBoundaryCondition->Clone(rConditionId, FaceNodes);
		      
			//p_cond->Data() = pBoundaryCondition->Data();

			//std::cout<<" _IDa_ "<<p_cond->Id()<<" MASTER ELEMENT "<<ie->Id()<<" MASTER NODE "<<rElementGeometry[lpofa(0,iface)].Id()<<" or "<<rElementGeometry[lpofa(NumberNodesInFace,iface)].Id()<<std::endl;
			
			WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
			MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
			p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

			//p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rElementGeometry(lpofa(0,i)) ) );			
			WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
			MasterNodes.push_back( Node<3>::WeakPointer( rElementGeometry(lpofa(0,iface)) ) );
			p_cond->SetValue(MASTER_NODES,MasterNodes);

		      }
		      else{
		  
			if( mEchoLevel > 1 ){
			  std::cout<<"   NOT FOUND CONDITION :: CREATED-> ["<<rConditionId<<"] (";
			  std::cout<<FaceNodes[0].Id();
			  for(unsigned int f=1; f<FaceNodes.size(); f++)
			    std::cout<<", "<<FaceNodes[f].Id();
			    
			  std::cout<<")"<<std::endl;				
			}

			// something not implemented in geometry or condition PrintData
			//std::cout<<" ReferenceCondition "<<rReferenceCondition<<std::endl;

			p_cond = rReferenceCondition.Create(rConditionId, FaceNodes, properties);
		      
			//if a condition is created new nodes must be labeled TO_REFINE
			for(unsigned int j=0; j<FaceNodes.size(); j++)
			  {
			    FaceNodes[j].Set(TO_REFINE);
			  }

			MeshDataTransferUtilities TransferUtilities;

			TransferUtilities.InitializeBoundaryData(p_cond, *(mrRemesh.Transfer), rCurrentProcessInfo);
			
			//std::cout<<" _IDb_ "<<p_cond->Id()<<" MASTER ELEMENT "<<ie->Id()<<" MASTER NODE "<<rElementGeometry[lpofa(0,iface)].Id()<<" or "<<rElementGeometry[lpofa(NumberNodesInFace,iface)].Id()<<std::endl;
						
			WeakPointerVector< Element >& MasterElements = p_cond->GetValue(MASTER_ELEMENTS);
			MasterElements.push_back( Element::WeakPointer( *(ie.base()) ) );
			p_cond->SetValue(MASTER_ELEMENTS,MasterElements);

			//p_cond->GetValue(MASTER_NODES).push_back( Node<3>::WeakPointer( rElementGeometry(lpofa(0,i)) ) );			
			WeakPointerVector< Node<3> >& MasterNodes = p_cond->GetValue(MASTER_NODES);
			MasterNodes.push_back( Node<3>::WeakPointer( rElementGeometry(lpofa(0,iface)) ) );
			p_cond->SetValue(MASTER_NODES,MasterNodes);

		      }

		      mrModelPart.Conditions().push_back(p_cond);
		      // Set new conditions: end

		    } //end no point condition
		    
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

    bool CheckAcceptedCondition(ModelPart& rModelPart, Condition& rCondition)
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


    void AddConditionToModelPart(ModelPart& rModelPart, Condition::Pointer pCondition)
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

    ModelerUtilities::MeshingParameters& mrRemesh;

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
    BuildMeshBoundaryProcess& operator=(BuildMeshBoundaryProcess const& rOther);

    /// Copy constructor.
    //BuildMeshBoundaryProcess(BuildMeshBoundaryProcess const& rOther);


    ///@}

  }; // Class BuildMeshBoundaryProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    BuildMeshBoundaryProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const BuildMeshBoundaryProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.



#endif // KRATOS_BUILD_MESH_BOUNDARY_PROCESS_H_INCLUDED  defined 
