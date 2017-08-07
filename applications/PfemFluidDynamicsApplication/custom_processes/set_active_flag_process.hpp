//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:         August 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_ACTIVE_FLAG_PROCESS_H_INCLUDED )
#define  KRATOS_SET_ACTIVE_FLAG_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes


#include "spatial_containers/spatial_containers.h"

#include "custom_processes/set_active_flag_process.hpp"
#include "custom_utilities/modeler_utilities.hpp"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"

///VARIABLES used:
//Data:     
//StepData: 
//Flags:    (checked) 
//          (set)     
//          (modified)  
//          (reset)   


namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef  ModelPart::NodesContainerType                      NodesContainerType;
  typedef  ModelPart::ElementsContainerType                ElementsContainerType;
  typedef  ModelPart::MeshType::GeometryType::PointsArrayType    PointsArrayType;

 
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
  class SetActiveFlagProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetActiveFlagProcess
    KRATOS_CLASS_POINTER_DEFINITION( SetActiveFlagProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetActiveFlagProcess(ModelPart& rModelPart,
			 bool unactivePeakElements,
			 bool unactiveSliverElements,
			 int EchoLevel)
      : mrModelPart(rModelPart)
    {
      mUnactivePeakElements = unactivePeakElements;
      mUnactiveSliverElements = unactiveSliverElements;
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~SetActiveFlagProcess()
    {
    }

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
	// std::cout<<" Execute() in SetActiveFlagProcess"<<std::endl;
#pragma omp parallel
	{
	  double tolerance=0.0000000001;
	  ModelPart::ElementIterator ElemBegin;
	  ModelPart::ElementIterator ElemEnd;
	  OpenMPUtils::PartitionedIterators(mrModelPart.Elements(),ElemBegin,ElemEnd);
	  for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
	    {
	      const unsigned int dimension = (itElem)->GetGeometry().WorkingSpaceDimension();
	      bool peakElementsEliminationCriteria=false;
	      bool sliverEliminationCriteria=false;

	      if(mUnactiveSliverElements == true){
		ModelerUtilities ModelerUtils;
		double ElementalVolume =  0;
		if(dimension==2){
		  ElementalVolume =  (itElem)->GetGeometry().Area();
		}else if(dimension==3){
		  ElementalVolume =  (itElem)->GetGeometry().Volume();
		}else{
		  ElementalVolume = 0;
		}
		double ModelPartVolume=ModelerUtils.ComputeModelPartVolume(mrModelPart);
		double CriticalVolume=0.01*ModelPartVolume/double(mrModelPart.Elements().size());
		if(ElementalVolume<CriticalVolume && ElementalVolume>0){
		  sliverEliminationCriteria=true;
		  // std::cout<<"RESET ACTIVE FOR THIS SLIVER! \t";
		  // std::cout<<"its volume is "<<ElementalVolume<<" vs CriticalVolume "<<CriticalVolume<<std::endl;
		}
	      }

	      if(mUnactivePeakElements == true && sliverEliminationCriteria==false){
		ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
		double timeInterval = rCurrentProcessInfo[DELTA_TIME];
		double scalarProduct=1.0;
		bool doNotErase=false;
		for(unsigned int i=0; i<itElem->GetGeometry().size(); i++)
		  {
		    if(itElem->GetGeometry()[i].IsNot(RIGID) && itElem->GetGeometry()[i].IsNot(FREE_SURFACE)){
		      peakElementsEliminationCriteria=false;
		      doNotErase=true;
		      break;
		    }else if(itElem->GetGeometry()[i].Is(RIGID) && itElem->GetGeometry()[i].Is(FREE_SURFACE) && doNotErase==false){
		      peakElementsEliminationCriteria=true;
		      const array_1d<double,3> &wallVelocity = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
		      double normWallVelocity=norm_2(wallVelocity);
		      if(normWallVelocity==0){// up to now this is for fixed walls only
			for(unsigned int j=0; j<itElem->GetGeometry().size(); j++)
			  {

			    if(itElem->GetGeometry()[j].IsNot(RIGID) && itElem->GetGeometry()[j].Is(FREE_SURFACE)){
			      Point<3> freeSurfaceToRigidNodeVector=itElem->GetGeometry()[i].Coordinates()-itElem->GetGeometry()[j].Coordinates();
			      const array_1d<double,3> &freeSurfaceVelocity = itElem->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);

			      double freeSurfaceToRigidNodeDistance=sqrt(freeSurfaceToRigidNodeVector[0]*freeSurfaceToRigidNodeVector[0]+
									 freeSurfaceToRigidNodeVector[1]*freeSurfaceToRigidNodeVector[1]+
									 freeSurfaceToRigidNodeVector[2]*freeSurfaceToRigidNodeVector[2]);
			      double displacementFreeSurface=timeInterval*(sqrt(freeSurfaceVelocity[0]*freeSurfaceVelocity[0]+
										freeSurfaceVelocity[1]*freeSurfaceVelocity[1]+
										freeSurfaceVelocity[2]*freeSurfaceVelocity[2]));
			      if(dimension==2){
				scalarProduct=freeSurfaceToRigidNodeVector[0]*freeSurfaceVelocity[0]+freeSurfaceToRigidNodeVector[1]*freeSurfaceVelocity[1];
			      }else if(dimension==3){
				scalarProduct=freeSurfaceToRigidNodeVector[0]*freeSurfaceVelocity[0]+freeSurfaceToRigidNodeVector[1]*freeSurfaceVelocity[1]+freeSurfaceToRigidNodeVector[2]*freeSurfaceVelocity[2];
			      }
			      if(scalarProduct>tolerance && displacementFreeSurface>(0.01*freeSurfaceToRigidNodeDistance)){
				// if(scalarProduct>tolerance){
				peakElementsEliminationCriteria=false;
				doNotErase=true;
				break;
			      }else{
				// I will not unactive the element if the free-surface node is sorrounded by rigd nodes only
				WeakPointerVector<Node<3> >& rN = itElem->GetGeometry()[j].GetValue(NEIGHBOUR_NODES);
				unsigned int rigidNodes=0;
				unsigned int freeSurfaceNodes=0;
				for(unsigned int i = 0; i < rN.size(); i++)
				  {
				    if(rN[i].Is(RIGID))
				      rigidNodes += 1;
				    if(rN[i].Is(FREE_SURFACE) && rN[i].IsNot(RIGID))
				      freeSurfaceNodes += 1;
				  }
				if(dimension==2){
				  if(rigidNodes==rN.size()){
				    peakElementsEliminationCriteria=false;
				    doNotErase=true;
				    break;
				  }
				}else if(dimension==3){
				  if(rigidNodes==rN.size() || freeSurfaceNodes==1 || (scalarProduct>tolerance && freeSurfaceNodes<4) ){
				    peakElementsEliminationCriteria=false;
				    doNotErase=true;
				    break;
				  }
				}

			      }
			    }
			  }
		      }
		    }
		    // if(peakElementsEliminationCriteria==false)
		    //   break;
		    // if(scalarProduct<tolerance && goAhead==true){
		    //   peakElementsEliminationCriteria=true;
		    //   // std::cout<<"RESET ACTIVE FOR THIS PEAK ELEMENT! \t";
		    //   // std::cout<<"scalarProduct= "<<scalarProduct<<std::endl;
		    // }
		  }  
	      }

	      if(peakElementsEliminationCriteria==true ||  sliverEliminationCriteria==true){
		(itElem)->Set(ACTIVE,false);
	      }else{
		(itElem)->Set(ACTIVE,true);
	      }

	    }

	}


      KRATOS_CATCH(" ")
	};

    ///@}
    ///@name Operators
    ///@{ 

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
      return "SetActiveFlagProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "SetActiveFlagProcess";
    }

    virtual void ExecuteInitialize()
    {
    }

    virtual void ExecuteFinalize()
    {
      KRATOS_TRY
#pragma omp parallel
	{
	  ModelPart::ElementIterator ElemBegin;
	  ModelPart::ElementIterator ElemEnd;
	  OpenMPUtils::PartitionedIterators(mrModelPart.Elements(),ElemBegin,ElemEnd);
	  for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
	    {
	      if((itElem)->IsNot(ACTIVE)){
		for(unsigned int i=0; i<itElem->GetGeometry().size(); i++)
		  {
		    if(itElem->GetGeometry()[i].Is(RIGID) && itElem->GetGeometry()[i].Is(FREE_SURFACE)){
		      WeakPointerVector<Element >& neighb_elems = itElem->GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS);
		      bool doNotSetNullPressure=false;
		      for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
			{
			  if((ne)->Is(ACTIVE)){
			    doNotSetNullPressure=true;
			    break;
			  }
			}
		      if(doNotSetNullPressure==false)
			itElem->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0;
		    }
		  }
	      }
	      // (itElem)->Set(ACTIVE,true);
	    }


	}
      KRATOS_CATCH(" ")    
	}

  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

 

    //*******************************************************************************************
    //*******************************************************************************************


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
    ModelPart& mrModelPart;

    int mEchoLevel;
    bool mUnactivePeakElements;
    bool mUnactiveSliverElements;


    ///@}
    ///@name Member Variables
    ///@{


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
    SetActiveFlagProcess& operator=(SetActiveFlagProcess const& rOther);

    /// Copy constructor.
    //SetActiveFlagProcess(SetActiveFlagProcess const& rOther);


    ///@}

  }; // Class SetActiveFlagProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    SetActiveFlagProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const SetActiveFlagProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_SET_ACTIVE_FLAG_PROCESS_H_INCLUDED  defined 

