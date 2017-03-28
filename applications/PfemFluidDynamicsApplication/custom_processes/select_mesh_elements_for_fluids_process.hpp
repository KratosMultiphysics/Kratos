//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined( KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_PROCESS_H_INCLUDED )
#define KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"
#include "geometries/tetrahedra_3d_4.h"

///VARIABLES used:
//Data:     
//StepData: NODAL_H, CONTACT_FORCE
//Flags:    (checked) TO_ERASE, BOUNDARY, NEW_ENTITY
//          (set)     
//          (modified)  
//          (reset)   
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process labels the elements to be refined in the mesher
    it applies a size constraint to elements that must be refined.
    
*/
class SelectMeshElementsForFluidsProcess
  : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( SelectMeshElementsForFluidsProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SelectMeshElementsForFluidsProcess(ModelPart& rModelPart,
			      ModelerUtilities::MeshingParameters& rRemeshingParameters,
			      int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
    
      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~SelectMeshElementsForFluidsProcess() {}


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

	if( mEchoLevel > 1 ){
	  std::cout<<" [ SELECT MESH ELEMENTS in PfemFluid: ("<<mrRemesh.OutMesh.GetNumberOfElements()<<") "<<std::endl;
	  std::cout<<"MODEL PART InNumberOfElements "<<mrRemesh.InMesh.GetNumberOfElements()<<std::endl;
	  std::cout<<"MODEL PART InNumberOfPoints "<<mrRemesh.InMesh.GetNumberOfPoints()<<std::endl;
	  std::cout<<"MODEL PART OutNumberOfElements "<<mrRemesh.OutMesh.GetNumberOfElements()<<std::endl;
	  std::cout<<"MODEL PART OutNumberOfPoints "<<mrRemesh.OutMesh.GetNumberOfPoints()<<std::endl;
	}
      int& OutNumberOfElements = mrRemesh.OutMesh.GetNumberOfElements();
      mrRemesh.PreservedElements.clear();
      mrRemesh.PreservedElements.resize(OutNumberOfElements);
      std::fill( mrRemesh.PreservedElements.begin(), mrRemesh.PreservedElements.end(), 0 );
      mrRemesh.MeshElementsSelectedFlag = true;

      
      mrRemesh.Info->NumberOfElements=0;
    
      const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      bool firstMesh=false;
      if(currentTime<2*timeInterval){
	firstMesh=true;
      }
  
      bool box_side_element = false;
      bool wrong_added_node = false;

      int number_of_slivers = 0;
      int sliverNodes = 0;

      if(mrRemesh.ExecutionOptions.IsNot(ModelerUtilities::SELECT_TESSELLATION_ELEMENTS))
	{
	  for(int el=0; el<OutNumberOfElements; el++)
	    {
	      mrRemesh.PreservedElements[el]=1;
	      mrRemesh.Info->NumberOfElements+=1;
	    }
	}
      else
	{
	  if( mEchoLevel > 0 )
	    std::cout<<"   Start Element Selection "<<OutNumberOfElements<<std::endl;

	  ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
	  const unsigned int nds = element_begin->GetGeometry().size();
	  const unsigned int dimension = element_begin->GetGeometry().WorkingSpaceDimension();

	  int* OutElementList = mrRemesh.OutMesh.GetElementList();
	 
	  ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes(mMeshId);

	  int el = 0;
	  int number = 0;
	  // std::cout<<" num nodes "<<rNodes.size()<<std::endl;
	  // std::cout<<" NodalPreIdsSize "<<mrRemesh.NodalPreIds.size()<<std::endl;

	  //#pragma omp parallel for reduction(+:number) private(el)
	  for(el=0; el<OutNumberOfElements; el++)
	    {
	      Geometry<Node<3> > vertices;

	      unsigned int  numfreesurf =0;
	      unsigned int  numboundary =0;	      
	      unsigned int  numrigid =0;	      
	      unsigned int  numinternalsolid =0;	      
	      unsigned int  numsolid =0;	      
	      unsigned int  numinlet =0;	      
	      unsigned int  numisolated =0;	      
	      // unsigned int  numinsertednodes =0;	      


	      box_side_element = false;
	      for(unsigned int pn=0; pn<nds; pn++)
		{
		  // std::cout<<" pn "<<pn<<std::endl;
		  //set vertices
		  if(mrRemesh.NodalPreIds[OutElementList[el*nds+pn]]<0){
		    if(mrRemesh.Options.IsNot(ModelerUtilities::CONTACT_SEARCH))
		      std::cout<<" ERROR: something is wrong: nodal id < 0 "<<std::endl;
		    box_side_element = true;
		    break;
		  }
		  
		  if(OutElementList[el*nds+pn]<=0)
		    std::cout<<" ERROR: something is wrong: nodal id < 0 "<<el<<std::endl;


		  if( (unsigned int)OutElementList[el*nds+pn] > mrRemesh.NodalPreIds.size() ){
		    wrong_added_node = true;
		    std::cout<<" ERROR: something is wrong: node out of bounds "<<std::endl;
		    break;
		  }
		
		  //vertices.push_back( *((rNodes).find( OutElementList[el*nds+pn] ).base() ) );
		  vertices.push_back(rNodes(OutElementList[el*nds+pn]));

		  //check flags on nodes

		  if(vertices.back().Is(ISOLATED)){
		    numisolated++;
		  }

		  if(vertices.back().Is(BOUNDARY)){
		    numboundary++;
		    // std::cout<<" BOUNDARY COORDINATES: "<<vertices.back().Coordinates()<<std::endl;
		  }
		  if(vertices.back().Is(RIGID) || vertices.back().Is(SOLID)){
		    numrigid++;
		    // std::cout<<" rigid COORDINATES: "<<vertices.back().Coordinates()<<std::endl;
		  }
		  if(vertices.back().Is(SOLID) && vertices.back().IsNot(BOUNDARY)){
		    numinternalsolid++;
		    // std::cout<<" internal solid COORDINATES: "<<vertices.back().Coordinates()<<std::endl;
		  }
		  if(vertices.back().Is(SOLID)){
		    numsolid++;
		    // std::cout<<"solid COORDINATES: "<<vertices.back().Coordinates()<<std::endl;
		  }
		  if(vertices.back().IsNot(RIGID) && vertices.back().Is(BOUNDARY)){
		    numfreesurf++;
		    // std::cout<<"FREE SURFACE COORDINATES: "<<vertices.back().Coordinates()<<std::endl;
		  }

		  if(vertices.back().Is(INLET)){
		    // vertices.back().Reset(INLET);
		    numinlet++;
		  }
		}
	      
	      
	      if(box_side_element || wrong_added_node){
		std::cout<<" ,,,,,,,,,,,,,,,,,,,,,,,,,,,,, Box_Side_Element "<<std::endl;
		continue;
	      }
	      

	      double Alpha =  mrRemesh.AlphaParameter; //*nds;

	      // if((numinlet+numrigid)>=nds || (numinlet+numsolid)>=nds){
	      // 	Alpha*=0;
	      // }
	      // if(numinlet==3){
	      // 	Alpha*=1.5;
	      // } else if(numinlet==2){
	      // 	Alpha*=1.15;
	      // }
	 
	      // if(numfreesurf==nds || (numisolated+numfreesurf)==nds){
	      // 	Alpha*=0.75;
	      // }else if(numfreesurf==2 || (numisolated+numfreesurf)==2){
	      // 	Alpha*=0.95;
	      // }
	      // if(numrigid==0 && numfreesurf==0 && numisolated==0){
	      // 	Alpha*=1.5;
	      // }


	      // if(numisolated==0 && numfreesurf==0 && firstMesh==false){
	      // if(numisolated==0 && numinlet>0 && (numfreesurf>0 && (numisolated+numfreesurf+numrigid)<nds) && firstMesh==false){
	      // if((numfreesurf>0 && (numisolated+numfreesurf+numrigid)<nds) && firstMesh==false){
	      // 	if(dimension==2){
	      // 	  if(numrigid==0 && numfreesurf==0){
	      // 	    // Alpha*=2.0;
	      // 	    Alpha*=1.5;
	      // 	  }else{
	      // 	    //Alpha*=1.75;
	      // 	    Alpha*=1.35;
	      // 	  }

	      // 	}
	      // 	if(dimension==3){
	      // 	  if(numrigid==0){
	      // 	    // Alpha*=1.15;
	      // 	    Alpha*=1.5;
	      // 	    // std::cout<<"RIGID 0";
	      // 	  }else{
	      // 	    // Alpha*=1.25;
	      // 	    Alpha*=1.5;
	      // 	    // std::cout<<"Alpha*=1.5 "<<Alpha;
	      // 	  }
	      // 	}
	      // }

	      if(dimension==2){
		if(numfreesurf==nds || (numisolated+numfreesurf)==nds){
		  Alpha*=0.7;
		}else if((numrigid+numisolated+numfreesurf)==nds){
		  Alpha*=0.9;
		}else if(numfreesurf==2 || (numisolated+numfreesurf)==2){
		  Alpha*=0.95;
		}
		if(numrigid==0 && numfreesurf==0 && numisolated==0){
		  Alpha*=1.75;
		}else if(numfreesurf==0 && numisolated==0){
		  Alpha*=1.25;
		}
	      }else  if(dimension==3){
		if(numfreesurf==nds || (numisolated+numfreesurf)==nds){
		  Alpha*=0.85;
		}else if((numrigid+numisolated+numfreesurf)==nds){
		  Alpha*=0.95;
		}// else if(numfreesurf==3 || (numisolated+numfreesurf)==3){
		//   Alpha*=0.975;
		// }
		if(numrigid==0 && numfreesurf==0 && numisolated==0){
		  Alpha*=1.75;
		}else if(numfreesurf==0 && numisolated==0){
		  Alpha*=1.25;
		}
	      }
	      if(firstMesh==true){
		Alpha*=1.15;
	      }
	      
    
	      bool accepted=false;
	      
	      ModelerUtilities ModelerUtils;
	      
	      if(mrRemesh.Options.Is(ModelerUtilities::CONTACT_SEARCH))
		{
		  accepted=ModelerUtils.ShrankAlphaShape(Alpha,vertices,mrRemesh.OffsetFactor,dimension);
		}
	      else
		{

		  double MeanMeshSize=mrRemesh.Refine->CriticalRadius;
		  accepted=ModelerUtils.AlphaShape(Alpha,vertices,dimension,MeanMeshSize);

		}


	      //3.1.-
	      bool self_contact = false;
	      if(mrRemesh.Options.Is(ModelerUtilities::CONTACT_SEARCH))
		self_contact = ModelerUtils.CheckSubdomain(vertices);
	    	    
	      //4.- to control that the element is inside of the domain boundaries
	      if(accepted)
		{
		  if(mrRemesh.Options.Is(ModelerUtilities::CONTACT_SEARCH))
		    {
		      accepted=ModelerUtils.CheckOuterCentre(vertices,mrRemesh.OffsetFactor, self_contact);
		    }
		}

	      if(numrigid==nds){
		accepted=false;
	      }
	      //5.- to control that the element has a good shape
	      if(accepted)
		{
		  if(nds==4){
		    Geometry<Node<3> >* tetrahedron = new Tetrahedra3D4<Node<3> > (vertices);
		    double Volume = tetrahedron->Volume();
		    // int sliver=0;
		    // accepted=ModelerUtils.CheckGeometryShape(*tetrahedron,sliver);
		    double CriticalVolume=0.01*mrRemesh.Refine->MeanVolume;
		    if(Volume<CriticalVolume){
		      // accepted=ModelerUtils.CheckGeometryShape(*tetrahedron,sliver);
		      // std::cout<<"SLIVER! VOLUME IS "<<Volume<<" VS CRITICAL VOLUME "<<CriticalVolume<<std::endl;
		      // if( sliver){
		      for( unsigned int n=0; n<nds; n++)
		    	{
		    	  vertices[n].Set(INTERFACE);
		    	  sliverNodes++;
       		    	}
		      // accepted = false;
		      number_of_slivers++;
		    }

		    delete tetrahedron;
		  }
		}

	      if(accepted)
		{
		  //std::cout<<" Element ACCEPTED after cheking Center "<<number<<std::endl;
		  number+=1;
		  mrRemesh.PreservedElements[el] = number;
		}
	      // else{
	      
	      //   std::cout<<" Element DID NOT pass INNER/OUTER check "<<std::endl;
	      // }


	    }

	  mrRemesh.Info->NumberOfElements=number;

	}


      std::cout<<" sliver nodes= "<<sliverNodes<<std::endl;
      // mrRemesh.Info->RemovedNodes+=sliverNodes;
      std::cout<<"  fluid Number of Preserved Elements "<<mrRemesh.Info->NumberOfElements<<" (slivers detected: "<<number_of_slivers<<") "<<std::endl;
      std::cout<<"  TOTAL removed nodes "<<mrRemesh.Info->RemovedNodes<<std::endl;

      if(mrRemesh.ExecutionOptions.IsNot(ModelerUtilities::KEEP_ISOLATED_NODES)){


	ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin(mMeshId);	  
	const unsigned int nds = (*element_begin).GetGeometry().size();
      
	int* OutElementList = mrRemesh.OutMesh.GetElementList();
      
	ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes(mMeshId);

	//check engaged nodes
	for(int el=0; el<OutNumberOfElements; el++)
	  {
	    if( mrRemesh.PreservedElements[el] ){
	      for(unsigned int pn=0; pn<nds; pn++)
		{
		  //set vertices
		  rNodes[OutElementList[el*nds+pn]].Set(BLOCKED);

		}

	    }
	    
	  }

	int count_release = 0;
	for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; i_node++)
	  {

	    if( i_node->IsNot(BLOCKED)  ){
	      if(!(i_node->Is(FREE_SURFACE) || i_node->Is(RIGID))){
		i_node->Set(TO_ERASE);
		if( mEchoLevel > 0 )
		  std::cout<<" NODE "<<i_node->Id()<<" RELEASE "<<std::endl;
		if( i_node->Is(BOUNDARY) )
		  std::cout<<" ERROR: node "<<i_node->Id()<<" IS BOUNDARY RELEASE "<<std::endl;
	      }
	    }
	    if( i_node->Is(TO_ERASE)  ){
	      count_release++;
	    }

	    i_node->Reset(BLOCKED);
	  }
	  
	if( mEchoLevel > 0 )
	  std::cout<<"   fluid NUMBER OF RELEASED NODES "<<count_release<<std::endl;

      }
      else{
	
	ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes(mMeshId);
	int count_sliver = 0;

	for(ModelPart::NodesContainerType::iterator i_node = rNodes.begin() ; i_node != rNodes.end() ; i_node++)
	  { 
	    i_node->Reset(BLOCKED);
	      
	    if( i_node->Is(INTERFACE)){
	      count_sliver++;
	    }

	  }
	// std::cout<<"count_sliver++ "<<count_sliver<<std::endl;

      }


      mrRemesh.InputInitializedFlag = false;
      // mModelerUtilities.SetNodes(mrModelPart,mrRemesh,mMeshId);
      mModelerUtilities.SetNodes(mrModelPart,mrRemesh);
      mrRemesh.InputInitializedFlag = true;

      if( mEchoLevel > 0 ){
	std::cout<<"   Generated_Elements :"<<OutNumberOfElements<<std::endl;
	std::cout<<"   Passed_AlphaShape  :"<<mrRemesh.Info->NumberOfElements<<std::endl;
	std::cout<<"   SELECT MESH ELEMENTS ]; "<<std::endl;
      }

      KRATOS_CATCH( "" )

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
        return "SelectMeshElementsForFluidsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SelectMeshElementsForFluidsProcess";
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
    ModelPart& mrModelPart;
 
    ModelerUtilities::MeshingParameters& mrRemesh;

    ModelerUtilities mModelerUtilities;  

    ModelPart::IndexType mMeshId; 

    int mEchoLevel;

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
    SelectMeshElementsForFluidsProcess& operator=(SelectMeshElementsForFluidsProcess const& rOther);


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
                                  SelectMeshElementsForFluidsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SelectMeshElementsForFluidsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SELECT_MESH_ELEMENTS_FOR_FLUIDS_PROCESS_H_INCLUDED defined 

