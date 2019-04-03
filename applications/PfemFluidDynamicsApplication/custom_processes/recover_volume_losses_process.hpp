//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//


#if !defined(KRATOS_RECOVER_VOLUME_LOSSES_PROCESS_H_INCLUDED )
#define  KRATOS_RECOVER_VOLUME_LOSSES_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes 
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:      
//StepData: CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked) 
//          (set)     
//          (modified)  
//          (reset)   
//(set):=(set in this process)

namespace Kratos
{

  ///@name Kratos Classes
  ///@{

  /// Refine Mesh Elements Process 2D and 3D
  /** The process labels the nodes to be refined (TO_REFINE)
      if the ThresholdVariable  is larger than a ReferenceThreshold
  */

  class RecoverVolumeLossesProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION( RecoverVolumeLossesProcess );

    typedef ModelPart::NodeType                   NodeType;
    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RecoverVolumeLossesProcess(ModelPart& rModelPart,
			       MesherUtilities::MeshingParameters& rRemeshingParameters,
			       int EchoLevel) 
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
    {
      std::cout<<" generate_new_nodes_process_for_homogeneous_mesh CONSTRUCTOR "<<std::endl;

      mEchoLevel = EchoLevel;
    }


    /// Destructor.
    virtual ~RecoverVolumeLossesProcess() {}


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
    void Execute() override
    {

      KRATOS_TRY

	if( mEchoLevel > 0 )
	  std::cout<<" [ RECOVER VOLUME LOSSES: "<<std::endl;


      if( mrModelPart.Name() != mrRemesh.SubModelPartName )
	std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

      MesherUtilities MesherUtils;
      const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];
      // double TotalVolumeLoss=0;
      double initialVolume=0;
      double currentVolume=0;
      double volumeLoss=0;
	
      if(currentTime<=2*timeInterval){
	initialVolume=MesherUtils.ComputeModelPartVolume(mrModelPart);
	if( mEchoLevel > 0 )
	  std::cout<<"setting                                InitialVolume "<<initialVolume<<std::endl;
	mrRemesh.Info->SetInitialMeshVolume(initialVolume);
      }

      if(currentTime>2*timeInterval){

	initialVolume=mrRemesh.Info->InitialMeshVolume;
	currentVolume=MesherUtils.ComputeModelPartVolume(mrModelPart);
	volumeLoss=initialVolume-currentVolume;
	if( mEchoLevel > 0 ){
	  std::cout<<" InitialVolume "<<initialVolume<<" currentVolume "<<currentVolume<<"-->  volumeLoss "<<volumeLoss<<std::endl;
	}
	double freeSurfaceLength=0;

	/////////////////////////        compute the free-surface length         /////////////////////////
	ModelPart::ElementsContainerType::iterator element_begin = mrModelPart.ElementsBegin();	  
	//ModelPart::NodesContainerType::iterator nodes_begin = mrModelPart.NodesBegin();
	const unsigned int nds = element_begin->GetGeometry().size();
	for(ModelPart::ElementsContainerType::const_iterator ie = element_begin; ie != mrModelPart.ElementsEnd(); ie++)
	  {
	    unsigned int freesurfaceNodes=0;
	    double freeSurfaceElementalSize=0;
	    for(unsigned int pn=0; pn<nds; pn++)
	      {
		// if(ie->GetGeometry()[pn].Is(FREE_SURFACE)){
		if(ie->GetGeometry()[pn].Is(BOUNDARY) && ie->GetGeometry()[pn].IsNot(RIGID) && ie->GetGeometry()[pn].Is(FLUID)){
		  freesurfaceNodes++;
		}
	      }
	    if(dimension==2 && freesurfaceNodes==2){
	      array_1d<double,3> CoorDifference(3,0.0);
	      if(ie->GetGeometry()[0].Is(BOUNDARY) && ie->GetGeometry()[0].IsNot(RIGID) &&
		 ie->GetGeometry()[1].Is(BOUNDARY) && ie->GetGeometry()[1].IsNot(RIGID))
		noalias(CoorDifference) = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[0].Coordinates();
		// CoorDifference = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[0].Coordinates();
	      if(ie->GetGeometry()[0].Is(BOUNDARY) && ie->GetGeometry()[0].IsNot(RIGID) &&
		 ie->GetGeometry()[2].Is(BOUNDARY) && ie->GetGeometry()[2].IsNot(RIGID))
		noalias(CoorDifference) = ie->GetGeometry()[2].Coordinates() - ie->GetGeometry()[0].Coordinates();
		// CoorDifference = ie->GetGeometry()[2].Coordinates() - ie->GetGeometry()[0].Coordinates();
	      if(ie->GetGeometry()[1].Is(BOUNDARY) && ie->GetGeometry()[1].IsNot(RIGID) &&
		 ie->GetGeometry()[2].Is(BOUNDARY) && ie->GetGeometry()[2].IsNot(RIGID))
		noalias(CoorDifference) = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[2].Coordinates();
		// CoorDifference = ie->GetGeometry()[1].Coordinates() - ie->GetGeometry()[2].Coordinates();

	      double SquaredLength = CoorDifference[0]*CoorDifference[0] + CoorDifference[1]*CoorDifference[1];
	      freeSurfaceElementalSize=sqrt(SquaredLength);
	      // std::cout<<"length "<<freeSurfaceElementalSize<<std::endl;

	    }
	    if(dimension==3){
	      if(ie->GetGeometry()[0].Is(BOUNDARY) && ie->GetGeometry()[0].IsNot(RIGID) &&
		 ie->GetGeometry()[1].Is(BOUNDARY) && ie->GetGeometry()[1].IsNot(RIGID) &&
		 ie->GetGeometry()[2].Is(BOUNDARY) && ie->GetGeometry()[2].IsNot(RIGID)){
		freeSurfaceElementalSize=Compute3dArea(ie->GetGeometry()[0].Coordinates(),
						       ie->GetGeometry()[1].Coordinates(),
						       ie->GetGeometry()[2].Coordinates());
	      }
	      if(ie->GetGeometry()[0].Is(BOUNDARY) && ie->GetGeometry()[0].IsNot(RIGID) &&
		 ie->GetGeometry()[1].Is(BOUNDARY) && ie->GetGeometry()[1].IsNot(RIGID) &&
		 ie->GetGeometry()[3].Is(BOUNDARY) && ie->GetGeometry()[3].IsNot(RIGID)){
		freeSurfaceElementalSize=Compute3dArea(ie->GetGeometry()[0].Coordinates(),
						       ie->GetGeometry()[1].Coordinates(),
						       ie->GetGeometry()[3].Coordinates());
	      }
	      if(ie->GetGeometry()[0].Is(BOUNDARY) && ie->GetGeometry()[0].IsNot(RIGID) &&
		 ie->GetGeometry()[2].Is(BOUNDARY) && ie->GetGeometry()[2].IsNot(RIGID) &&
		 ie->GetGeometry()[3].Is(BOUNDARY) && ie->GetGeometry()[3].IsNot(RIGID)){
		freeSurfaceElementalSize=Compute3dArea(ie->GetGeometry()[0].Coordinates(),
						       ie->GetGeometry()[2].Coordinates(),
						       ie->GetGeometry()[3].Coordinates());
	      }
	      if(ie->GetGeometry()[1].Is(BOUNDARY) && ie->GetGeometry()[1].IsNot(RIGID) &&
		 ie->GetGeometry()[2].Is(BOUNDARY) && ie->GetGeometry()[2].IsNot(RIGID) &&
		 ie->GetGeometry()[3].Is(BOUNDARY) && ie->GetGeometry()[3].IsNot(RIGID)){
		freeSurfaceElementalSize=Compute3dArea(ie->GetGeometry()[1].Coordinates(),
						       ie->GetGeometry()[2].Coordinates(),
						       ie->GetGeometry()[3].Coordinates());
	      }
	    }
	    freeSurfaceLength+=freeSurfaceElementalSize;
	  }

	if(freeSurfaceLength==0){
	  freeSurfaceLength=1.0;
	}
	double offset=volumeLoss/freeSurfaceLength;
	// if( mEchoLevel > 0 )
	std::cout<<"freeSurface length "<<freeSurfaceLength<<"  offset "<<offset<<std::endl;

	for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; i_node++)
	  {
	    if(i_node->Is(BOUNDARY) && i_node->IsNot(RIGID)){
	      array_1d<double, 3>  Normal(3,0.0);
	      Normal    =i_node->FastGetSolutionStepValue(NORMAL);
	      i_node->X()+=Normal[0]*offset;
	      i_node->Y()+=Normal[1]*offset;
	      if(dimension==3){
		i_node->Z()+=Normal[2]*offset;
	      }
	    }
	  }
      }

 
      if( mEchoLevel > 0 )
	std::cout<<" RECOVER VOLUME LOSSES ]; "<<std::endl;

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
      return "RecoverVolumeLossesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "RecoverVolumeLossesProcess";
    }

    /// Print object's data.s
    void PrintData(std::ostream& rOStream) const override
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
 
    MesherUtilities::MeshingParameters& mrRemesh;

    MesherUtilities mMesherUtilities;  

    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    double Compute3dArea(array_1d<double,3> PointA, array_1d<double,3> PointB, array_1d<double,3> PointC){
      double a = MathUtils<double>::Norm3(PointA - PointB);
      double b = MathUtils<double>::Norm3(PointB - PointC);
      double c = MathUtils<double>::Norm3(PointC - PointA);
      double s = (a+b+c) / 2.0;
      double Area=std::sqrt(s*(s-a)*(s-b)*(s-c));
      return Area;
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
    RecoverVolumeLossesProcess& operator=(RecoverVolumeLossesProcess const& rOther);


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
				    RecoverVolumeLossesProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const RecoverVolumeLossesProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_RECOVER_VOLUME_LOSSES_PROCESS_H_INCLUDED  defined 


