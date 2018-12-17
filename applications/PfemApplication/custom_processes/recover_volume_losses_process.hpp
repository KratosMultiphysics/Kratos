//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:        August 2018 $
//   Revision:            $Revision:            0.0 $
//
//


#if !defined(KRATOS_RECOVER_VOLUME_LOSSES_PROCESS_H_INCLUDED )
#define  KRATOS_RECOVER_VOLUME_LOSSES_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "processes/process.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Move free surface to restore volume losses
/** The process moves free surface nodes artificially in order to recover the loss of volume
 */

class RecoverVolumeLossesProcess : public Process
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
                             int EchoLevel)
      : Process(Flags()), mrModelPart(rModelPart), mEchoLevel(EchoLevel)
  {
    mTotalVolume = 0;
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

  /// Execute method is used to execute the AssignPropertiesToNodesProcess algorithms.
  void Execute()  override
  {
  }

  /// this function is designed for being called at the beginning of the computations
  /// right after reading the model and the groups
  void ExecuteInitialize() override
  {
  }

  /// this function is designed for being execute once before the solution loop but after all of the
  /// solvers where built
  void ExecuteBeforeSolutionLoop() override
  {
    mTotalVolume = this->ComputeVolume(mrModelPart);
  }


  /// this function will be executed at every time step BEFORE performing the solve phase
  void ExecuteInitializeSolutionStep() override
  {
    KRATOS_TRY

    //recover volume losses due meshing
    mTotalVolume = this->RecoverVolume(mrModelPart);

    KRATOS_CATCH("")
  }

  /// this function will be executed at every time step AFTER performing the solve phase
  void ExecuteFinalizeSolutionStep() override
  {

    //recover volume losses due computation
    mTotalVolume = this->RecoverVolume(mrModelPart);

  }


  /// this function will be executed at every time step BEFORE  writing the output
  void ExecuteBeforeOutputStep() override
  {
  }


  /// this function will be executed at every time step AFTER writing the output
  void ExecuteAfterOutputStep() override
  {
  }


  /// this function is designed for being called at the end of the computations
  /// right after reading the model and the groups
  void ExecuteFinalize() override
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

 protected:

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

  ///@name Static Member Variables
  ///@{
  ///@}
  ///@name Member Variables
  ///@{

  ModelPart& mrModelPart;

  double mTotalVolume;

  int mEchoLevel;

  ///@}
  ///@name Private Operators
  ///@{
  ///@}
  ///@name Private Operations
  ///@{

  //*******************************************************************************************
  //*******************************************************************************************

  double RecoverVolume(ModelPart& rModelPart)
  {

    KRATOS_TRY

    double Tolerance = 1e-6;
    unsigned int NumberOfIterations = 5;
    unsigned int iteration = -1;

    double CurrentVolume = this->ComputeVolume(rModelPart);

    double Error = fabs(mTotalVolume-CurrentVolume);

    this->SetFreeSurfaceElements(rModelPart);
    double FreeSurfaceVolume = this->ComputeFreeSurfaceVolume(rModelPart);
    double FreeSurfaceArea   = this->ComputeFreeSurfaceArea(rModelPart);

    double VolumeIncrement = mTotalVolume-CurrentVolume;
    //initial prediction of the offset
    double Offset = VolumeIncrement/FreeSurfaceArea;

    FreeSurfaceVolume += VolumeIncrement;

    while( ++iteration < NumberOfIterations && Error > Tolerance )
    {

      this->MoveFreeSurface(rModelPart,Offset);

      double CurrentFreeSurfaceVolume = this->ComputeFreeSurfaceVolume(rModelPart);

      VolumeIncrement = (FreeSurfaceVolume - CurrentFreeSurfaceVolume);

      Offset = (VolumeIncrement / FreeSurfaceArea );

      Error = fabs(VolumeIncrement);

      //std::cout<<" Iteration: "<<iteration<<" Error in Volume "<<Error<<std::endl;
    }


    this->ResetFreeSurfaceElements(rModelPart);

    CurrentVolume = this->ComputeVolume(rModelPart);

    //std::cout<<" Recover Volume Losses perfomed in "<<iteration<<" iterations : Error "<<Error<<" CurrentVolume "<<CurrentVolume<<std::endl;

    return CurrentVolume;

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void MoveFreeSurface(ModelPart& rModelPart, const double& rOffset)
  {
    KRATOS_TRY

    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();
    int NumberOfNodes = rModelPart.NumberOfNodes();

    #pragma omp parallel for
    for (int i=0; i<NumberOfNodes; ++i)
    {
      ModelPart::NodesContainerType::iterator i_node = it_begin + i;
      if(i_node->Is(FREE_SURFACE) && (i_node->IsNot(RIGID) && i_node->IsNot(SOLID)) ){
        const array_1d<double,3>& rNormal = i_node->FastGetSolutionStepValue(NORMAL);
        i_node->Coordinates() += rOffset * rNormal;
        i_node->FastGetSolutionStepValue(DISPLACEMENT) += rOffset * rNormal;
        i_node->FastGetSolutionStepValue(DISPLACEMENT,1) += rOffset * rNormal;
      }
    }

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void SetFreeSurfaceElements(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
    {
      for(unsigned int i=0; i<i_elem->GetGeometry().size(); ++i)
      {
        if(i_elem->GetGeometry()[i].Is(FREE_SURFACE)){
          i_elem->Set(FREE_SURFACE,true);
          break;
        }
      }
    }

    KRATOS_CATCH(" ")
  }


  //*******************************************************************************************
  //*******************************************************************************************

  void ResetFreeSurfaceElements(ModelPart& rModelPart)
  {
    KRATOS_TRY

    for(ModelPart::ElementsContainerType::iterator i_elem = rModelPart.ElementsBegin() ; i_elem != rModelPart.ElementsEnd() ; ++i_elem)
    {
      i_elem->Set(FREE_SURFACE,false);
    }

    KRATOS_CATCH(" ")
  }

  //*******************************************************************************************
  //*******************************************************************************************

  double ComputeVolume(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double ModelPartVolume = 0;

    if( dimension == 2 ){

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:ModelPartVolume)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        if( i_elem->GetGeometry().Dimension() == 2 && i_elem->Is(FLUID) )
          ModelPartVolume += i_elem->GetGeometry().Area();
      }
    }
    else{ //dimension == 3

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:ModelPartVolume)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        if( i_elem->GetGeometry().Dimension() == 3 && i_elem->Is(FLUID) )
          ModelPartVolume += i_elem->GetGeometry().Volume();
      }
    }

    return ModelPartVolume;

    KRATOS_CATCH(" ")

  }

  //*******************************************************************************************
  //*******************************************************************************************

  double ComputeFreeSurfaceVolume(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double ModelPartVolume = 0;

    if( dimension == 2 ){

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:ModelPartVolume)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        if( i_elem->GetGeometry().Dimension() == 2 && i_elem->Is(FREE_SURFACE) && i_elem->Is(FLUID) )
          ModelPartVolume += i_elem->GetGeometry().Area();
      }
    }
    else{ //dimension == 3

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:ModelPartVolume)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        if( i_elem->GetGeometry().Dimension() == 3  && i_elem->Is(FREE_SURFACE) && i_elem->Is(FLUID) )
          ModelPartVolume += i_elem->GetGeometry().Volume();
      }
    }

    return ModelPartVolume;

    KRATOS_CATCH(" ")

  }


  //*******************************************************************************************
  //*******************************************************************************************

  double ComputeFreeSurfaceArea(ModelPart& rModelPart)
  {
    KRATOS_TRY

    const unsigned int dimension = rModelPart.GetProcessInfo()[SPACE_DIMENSION];
    double FreeSurfaceArea = 0;

    if( dimension == 2 ){

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for reduction(+:FreeSurfaceArea)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;
        GeometryType& rGeometry = i_elem->GetGeometry();
        if( rGeometry.Dimension() == 2 && i_elem->Is(FREE_SURFACE) && i_elem->Is(FLUID) ){
          for(unsigned int j=0; j<rGeometry.size()-1; ++j)
          {
            if(rGeometry[j].Is(FREE_SURFACE)){
              for(unsigned int k=j+1; k<rGeometry.size(); ++k)
              {
                if(rGeometry[k].Is(FREE_SURFACE)){
                  FreeSurfaceArea += norm_2( rGeometry[k].Coordinates() - rGeometry[j].Coordinates() );
                }
              }

            }
          }
        }
      }
    }
    else{ //dimension == 3

      DenseMatrix<unsigned int> lpofa; //connectivities of points defining faces
      DenseVector<unsigned int> lnofa; //number of points defining faces

      ModelPart::ElementsContainerType::iterator it_begin = rModelPart.ElementsBegin();
      int NumberOfElements = rModelPart.NumberOfElements();

      #pragma omp parallel for private(lpofa,lnofa) reduction(+:FreeSurfaceArea)
      for (int i=0; i<NumberOfElements; ++i)
      {
        ModelPart::ElementsContainerType::iterator i_elem = it_begin + i;

        GeometryType& rGeometry = i_elem->GetGeometry();

        if( rGeometry.Dimension() == 3 && i_elem->Is(FREE_SURFACE) && i_elem->Is(FLUID) ){

          rGeometry.NodesInFaces(lpofa);
          rGeometry.NumberNodesInFaces(lnofa);

          for(unsigned int iface=0; iface<rGeometry.FacesNumber(); ++iface)
          {
            unsigned int free_surface = 0;
            for(unsigned int j=1; j<=lnofa[iface]; ++j)
              if(rGeometry[j].Is(FREE_SURFACE))
                ++free_surface;

            if(free_surface==lnofa[iface])
              FreeSurfaceArea+=Compute3DArea(rGeometry[lpofa(1,iface)].Coordinates(),
                                             rGeometry[lpofa(2,iface)].Coordinates(),
                                             rGeometry[lpofa(3,iface)].Coordinates());
          }
        }
      }
    }

    return FreeSurfaceArea;

    KRATOS_CATCH(" ")

  }

  //*******************************************************************************************
  //*******************************************************************************************

  double Compute3DArea(array_1d<double,3> PointA, array_1d<double,3> PointB, array_1d<double,3> PointC){
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
