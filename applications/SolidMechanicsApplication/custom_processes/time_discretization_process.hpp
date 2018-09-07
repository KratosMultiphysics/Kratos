//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:           September 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TIME_DISCRETIZATION_PROCESS_H_INCLUDED )
#define  KRATOS_TIME_DISCRETIZATION_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "processes/process.h"

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
class TimeDiscretizationProcessTimeDiscretizationProcess
    : public Process
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of TimeDiscretizationProcess
  KRATOS_CLASS_POINTER_DEFINITION( TimeDiscretizationProcess );

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  TimeDiscretizationProcess(ModelPart& rModelPart)
      : mrModelPart(rModelPart)
  {
  }

  /// Default constructor.
  TimeDiscretizationProcess(ModelPart& rModelPart,
                            Parameters rParameters) : Process(Flags()) , mrModelPart(rModelPart)
  {
    KRATOS_TRY

    Parameters default_parameters( R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "start_time": 0,
                "end_time": 1,
                "time_step": 1,
                "prediction_level": -1, 
                "increase_factor": 2,
                "decrease_factor": 2,
                "steps_update_delay": 4
            }  )" );


    // Validate against defaults -- this ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mTime.Initialize(rParameters["time_step"].GetDouble(), rParameters["start_time"].GetDouble(), rParameters["end_time"].GetDouble());

    mTime.SetFactors(rParameters["increase_factor"].GetDouble(), rParameters["decrease_factor"].GetDouble(), rParameters["steps_update_delay"]);

    mTime.PredictionLevel = rParameters["prediction_level"].GetInt();
    
  }

  /// Destructor.
  virtual ~TimeDiscretizationProcess()
  {
  }

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

    ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
   
    if(!rCurrentProcessInfo[CONVERGENCE_ACHIEVED])
    {
      this->ReduceTimeStep();
    }
    else{
      this->UpdateTimeStep();
    }

    rCurrentProcessInfo[TIME] += mTime.CurrentStep;
    rCurrentProcessInfo[STEP] += 1;

    // it sets TIME and DELTA_TIME internally and Clones the SolutionStepData and ProcessInfo buffers
    mrModelPart.CloneTimeStep(rCurrentProcessInfo[TIME]);

    // update total time
    mTime.Total = rCurrentProcessInfo[TIME];
    
    KRATOS_CATCH("");
  }


  // Deprecated methods:
  
  // void PreExecute()
  // {

  //   KRATOS_TRY
        
  //   const double initialTimeInterval = rCurrentProcessInfo[INITIAL_DELTA_TIME];
  //   const double currentTimeInterval = rCurrentProcessInfo[CURRENT_DELTA_TIME];

  //   double updatedTime               = rCurrentProcessInfo[TIME];
  //   double updatedTimeInterval       = rCurrentProcessInfo[DELTA_TIME];

  //   double deltaTimeToNewMilestone   = initialTimeInterval;
  //   double minimumTimeInterval       = initialTimeInterval*0.0001;

  //   rCurrentProcessInfo.SetValue(PREVIOUS_DELTA_TIME,currentTimeInterval);
  //   rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,false);

  //   bool milestoneTimeReached=true;
  //   bool increaseTimeInterval=true;
  //   bool timeIntervalReduced=false;

  //   double tolerance=0.0001;
  //   updatedTime -= initialTimeInterval;
  //   unsigned int previousMilestoneStep=updatedTime/initialTimeInterval;
  //   deltaTimeToNewMilestone=initialTimeInterval*(previousMilestoneStep+1)-updatedTime;

  //   updatedTimeInterval =currentTimeInterval;

  //   bool badVelocityConvergence=rCurrentProcessInfo[BAD_VELOCITY_CONVERGENCE];
  //   bool badPressureConvergence=rCurrentProcessInfo[BAD_PRESSURE_CONVERGENCE];

  //   if(updatedTimeInterval<2.0*minimumTimeInterval && mEchoLevel > 0 && mrModelPart.GetCommunicator().MyPID() == 0){
  //     std::cout<<"ATTENTION! time step much smaller than initial time step, I'll not reduce it"<<std::endl;
  //   }

  //   if((badPressureConvergence==true || badVelocityConvergence==true) && updatedTimeInterval>(2.0*minimumTimeInterval)){
  //     updatedTimeInterval *=0.5;
  //     /* std::cout<<"reducing time step (bad convergence at the previous step)"<<updatedTimeInterval<<std::endl; */
  //     rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,true);
  //     timeIntervalReduced=true;
  //   }

  //   if(deltaTimeToNewMilestone<(1.0+tolerance)*updatedTimeInterval && deltaTimeToNewMilestone>initialTimeInterval*tolerance){
  //     rCurrentProcessInfo.SetValue(DELTA_TIME,deltaTimeToNewMilestone);
  //     if(deltaTimeToNewMilestone<0.75*updatedTimeInterval){
  //       timeIntervalReduced=true;
  //       rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,true);
  //     }
  //     updatedTimeInterval =deltaTimeToNewMilestone;
  //     milestoneTimeReached=true;
  //   }else{
  //     milestoneTimeReached=false;
  //     rCurrentProcessInfo.SetValue(DELTA_TIME,updatedTimeInterval);
  //   }

  //   if(timeIntervalReduced==false){

  //     if(updatedTimeInterval>(2.0*minimumTimeInterval)){

  //       const unsigned int dimension =  mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
  //       if(dimension==2){
  //         CheckNodalConditionForTimeStepReduction(updatedTimeInterval,increaseTimeInterval,timeIntervalReduced);
  //         if(timeIntervalReduced==false){
  //           CheckElementalConditionForTimeStepReduction(increaseTimeInterval);
  //         }
  //       }
  //     }

  //     // if(increaseTimeInterval==true && initialTimeInterval>(1.0+tolerance)*updatedTimeInterval && badPressureConvergence==false && badVelocityConvergence==false ){
  //     if(increaseTimeInterval==true && initialTimeInterval>(1.0+tolerance)*updatedTimeInterval && badVelocityConvergence==false ){
  //       IncreaseTimeInterval(updatedTimeInterval,deltaTimeToNewMilestone,tolerance,increaseTimeInterval);
  //     }
  //     else{
  //       increaseTimeInterval=false;
  //     }

  //   }

  //   double newTimeInterval = rCurrentProcessInfo[DELTA_TIME];
  //   double milestoneGap=fabs(newTimeInterval-deltaTimeToNewMilestone);

  //   if(milestoneGap<0.49*newTimeInterval && milestoneTimeReached==false){
  //     /* std::cout<<"the milestone is very close, I add "<<milestoneGap<<" to "<<newTimeInterval<<std::endl;*/
  //     newTimeInterval+=milestoneGap;
  //     rCurrentProcessInfo.SetValue(DELTA_TIME,newTimeInterval);
  //     milestoneTimeReached=true;
  //   }

  //   updatedTime+=newTimeInterval;
  //   rCurrentProcessInfo.SetValue(TIME,updatedTime);
  //   rCurrentProcessInfo.SetValue(CURRENT_DELTA_TIME,newTimeInterval);

  //   /* if(newTimeInterval<(0.49*currentTimeInterval)){ */
  //   /* std::cout<<"ATTENTION! new time step is more than 2 times smaller than the previous one"<<std::endl; */
  //   /* } */

  //   if(increaseTimeInterval==false && milestoneTimeReached==true && fabs(newTimeInterval-initialTimeInterval)>tolerance && !(deltaTimeToNewMilestone>newTimeInterval*(1.0+tolerance))){
  //     rCurrentProcessInfo.SetValue(CURRENT_DELTA_TIME,currentTimeInterval);
  //   }


  //   if (newTimeInterval<initialTimeInterval){
  //     std::cout<<"current time "<<updatedTime<<" time step: new  "<<newTimeInterval<<" previous "<<currentTimeInterval<<" initial  "<<initialTimeInterval<<"\n"<<std::endl;
  //   }


  //   KRATOS_CATCH("");

  // };


//   void CheckNodalConditionForTimeStepReduction(double updatedTimeInterval,
//                                                bool &increaseTimeInterval,
//                                                bool &timeIntervalReduced)
//   {

//     ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

// #pragma omp parallel
//     {
//       ModelPart::NodeIterator NodeBegin;
//       ModelPart::NodeIterator NodeEnd;
//       OpenMPUtils::PartitionedIterators(mrModelPart.Nodes(),NodeBegin,NodeEnd);
//       for (ModelPart::NodeIterator itNode = NodeBegin; itNode != NodeEnd; ++itNode)
//       {
//         if(itNode->IsNot(TO_ERASE) && itNode->IsNot(ISOLATED) && itNode->IsNot(SOLID)){
//           const array_1d<double,3> &Vel = itNode->FastGetSolutionStepValue(VELOCITY);
//           double NormVelNode=0;
//           for (unsigned int d = 0; d < 3; ++d){
//             NormVelNode+=Vel[d] * Vel[d];
//           }
//           double motionInStep=sqrt(NormVelNode)*updatedTimeInterval;
//           double unsafetyFactor=0;
//           WeakPointerVector< Node < 3 > >& neighb_nodes = itNode->GetValue(NEIGHBOUR_NODES);
//           for (WeakPointerVector< Node <3> >::iterator nn = neighb_nodes.begin();nn != neighb_nodes.end(); ++nn)
//           {
//             array_1d<double,3>  CoorNeighDifference=itNode->Coordinates()-nn->Coordinates();
//             double squaredDistance=0;
//             for (unsigned int d = 0; d < 3; ++d){
//               squaredDistance+=CoorNeighDifference[d]*CoorNeighDifference[d];
//             }
//             double nodeDistance=sqrt(squaredDistance);
//             double tempUnsafetyFactor=motionInStep/nodeDistance;
//             if(tempUnsafetyFactor>unsafetyFactor){
//               unsafetyFactor=tempUnsafetyFactor;
//             }
//           }

//           if(unsafetyFactor>0.35){
//             increaseTimeInterval=false;
//             if(unsafetyFactor>1.0){
//               double temporaryTimeInterval = rCurrentProcessInfo[DELTA_TIME];
//               double reducedTimeInterval=0.5*updatedTimeInterval;
//               if(reducedTimeInterval<temporaryTimeInterval){
//                 rCurrentProcessInfo.SetValue(DELTA_TIME,reducedTimeInterval);
//                 /* std::cout<<"reducing time step (nodal criterion)"<<reducedTimeInterval<<std::endl; */
//                 rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,true);
//                 timeIntervalReduced=true;
//                 break;
//               }
//             }
//           }
//         }
//       }
//     }

//   }


//   void CheckElementalConditionForTimeStepReduction(bool &increaseTimeInterval)
//   {

//     ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

// #pragma omp parallel
//     {
//       ModelPart::ElementIterator ElemBegin;
//       ModelPart::ElementIterator ElemEnd;
//       OpenMPUtils::PartitionedIterators(mrModelPart.Elements(),ElemBegin,ElemEnd);
//       for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
//       {
//         double temporaryTimeInterval=rCurrentProcessInfo[DELTA_TIME];
//         double currentElementalArea =  0;
//         const unsigned int dimension = (itElem)->GetGeometry().WorkingSpaceDimension();
//         if(dimension==2){
//           currentElementalArea =  (itElem)->GetGeometry().Area();
//           Geometry<Node<3> >  updatedElementCoordinates;
//           bool solidElement=false;
//           for(unsigned int i=0; i<itElem->GetGeometry().size(); ++i)
//           {
//             if(itElem->GetGeometry()[i].Is(SOLID) || itElem->GetGeometry()[i].Is(TO_ERASE) || itElem->IsNot(ACTIVE)){
//               solidElement=true;
//             }

//             const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
//             Point updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval;
//             updatedElementCoordinates.push_back(Kratos::make_shared<Node<3> >(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z()));
//           }

//           double newArea=0;
//           if(itElem->GetGeometry().size()==3){
//             Triangle2D3<Node<3> > myGeometry(updatedElementCoordinates);
//             newArea=myGeometry.Area();
//           }else if(itElem->GetGeometry().size()==6){
//             Triangle2D6<Node<3> > myGeometry(updatedElementCoordinates);
//             newArea=myGeometry.Area();
//           }else{
//             std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
//           }

//           if(solidElement==true){
//             newArea=currentElementalArea;
//           }

//           if(newArea<0.001*currentElementalArea && currentElementalArea>0){
//             double reducedTimeInterval=0.5*temporaryTimeInterval;

//             if(reducedTimeInterval<temporaryTimeInterval){
//               rCurrentProcessInfo.SetValue(DELTA_TIME,reducedTimeInterval);
//               /* std::cout<<"reducing time step (elemental inversion)"<<reducedTimeInterval<<std::endl; */
//               rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,true);
//               increaseTimeInterval=false;
//               break;
//             }
//           }else{
//             Geometry<Node<3> >  updatedEnlargedElementCoordinates;

//             for(unsigned int i=0; i<itElem->GetGeometry().size(); ++i)
//             {
//               const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
//               Point updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval*2.5;
//               updatedEnlargedElementCoordinates.push_back(Kratos::make_shared<Node<3> >(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z()));

//             }

//             if(itElem->GetGeometry().size()==3){
//               Triangle2D3<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
//               newArea=myGeometry.Area();
//             }else if(itElem->GetGeometry().size()==6){
//               Triangle2D6<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
//               newArea=myGeometry.Area();
//             }else{
//               std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
//             }

//             if(newArea<0.001*currentElementalArea && currentElementalArea>0){
//               increaseTimeInterval=false;
//               /* std::cout<<"I'll not reduce the time step but I'll not allow to increase it"<<std::endl; */
//             }

//           }
//         }
//         else if(dimension==3){
//           double currentElementalVolume =  (itElem)->GetGeometry().Volume();
//           Geometry<Node<3> >  updatedElementCoordinates;
//           bool solidElement=false;
//           for(unsigned int i=0; i<itElem->GetGeometry().size(); ++i)
//           {
//             if(itElem->GetGeometry()[i].Is(SOLID) || itElem->IsNot(ACTIVE)){
//               solidElement=true;
//             }
//             const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
//             Point updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval;
//             updatedElementCoordinates.push_back(Kratos::make_shared<Node<3> >(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z()));
//           }

//           double newVolume=0;
//           if(itElem->GetGeometry().size()==4){
//             Tetrahedra3D4<Node<3> > myGeometry(updatedElementCoordinates);
//             newVolume=myGeometry.Volume();
//           }else if(itElem->GetGeometry().size()==10){
//             Tetrahedra3D10<Node<3> > myGeometry(updatedElementCoordinates);
//             newVolume=myGeometry.Volume();
//           }else{
//             std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
//           }

//           if(solidElement==true){
//             newVolume=currentElementalVolume;
//           }

//           if(newVolume<0.001*currentElementalVolume && currentElementalVolume>0){
//             double reducedTimeInterval=0.5*temporaryTimeInterval;

//             if(reducedTimeInterval<temporaryTimeInterval){
//               rCurrentProcessInfo.SetValue(DELTA_TIME,reducedTimeInterval);
//               /* std::cout<<"reducing time step (elemental inversion)"<<reducedTimeInterval<<std::endl; */
//               rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,true);
//               increaseTimeInterval=false;
//               break;
//             }
//           }else{
//             Geometry<Node<3> >  updatedEnlargedElementCoordinates;

//             for(unsigned int i=0; i<itElem->GetGeometry().size(); ++i)
//             {
//               const array_1d<double,3> &Vel = itElem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
//               Point updatedNodalCoordinates=itElem->GetGeometry()[i].Coordinates()+Vel*temporaryTimeInterval*2.5;
//               updatedEnlargedElementCoordinates.push_back(Kratos::make_shared<Node<3> >(i,updatedNodalCoordinates.X(),updatedNodalCoordinates.Y(),updatedNodalCoordinates.Z()));
//             }

//             if(itElem->GetGeometry().size()==4){
//               Tetrahedra3D4<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
//               newVolume=myGeometry.Volume();
//             }else if(itElem->GetGeometry().size()==10){
//               Tetrahedra3D10<Node<3> > myGeometry(updatedEnlargedElementCoordinates);
//               newVolume=myGeometry.Volume();
//             }else{
//               std::cout<<"GEOMETRY NOT DEFINED"<<std::endl;
//             }

//             if(newVolume<0.001*currentElementalVolume && currentElementalVolume>0){
//               increaseTimeInterval=false;
//               /* std::cout<<"I'll not reduce the time step but I'll not allow to increase it"<<std::endl; */
//             }


//           }



//         }

//       }

//     }
//   }

  // void IncreaseTimeInterval(double updatedTimeInterval,
  //                           double deltaTimeToNewMilestone,
  //                           double tolerance,
  //                           bool &increaseTimeInterval)
  // {
  //   ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();

  //   double increasedTimeInterval = 2.0 * updatedTimeInterval;

  //   if(increasedTimeInterval<deltaTimeToNewMilestone*(1.0+tolerance))
  //   {
  //     rCurrentProcessInfo.SetValue(DELTA_TIME,increasedTimeInterval);
  //     rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED,true);
  //   }
  //   else{
  //     increaseTimeInterval=false;
  //   }

  // }

  ///@}
  ///@name Operators
  ///@{
  ///@}
  ///@name Access
  ///@{
  ///@}
  ///@name Inquiry
  ///@
  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const override
  {
    return "TimeDiscretizationProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "TimeDiscretizationProcess";
  }


 protected:
  ///@name Protected static Member Variables
  ///@{
  ///@}
  ///@name Protected member Variables
  ///@{

  struct TimeParameters {

    int PredictionLevel;
    
    double InitialStep;
    double PreviousStep;
    double CurrentStep;

    double PredictedStep;
    
    double Total;

    double Start;
    double End;
    double MileStone;

    double IncreaseFactor;
    double DecreaseFactor;

    int UpdateDelay;
    int DelayCounter;

    void Initialize(const double& rTimeStep, const double& rStartTime, const double& rEndTime)
    {
      Predict = false;
      PredictionLevel = 0;
      InitialStep = rTimeStep;
      PreviousStep = rTimeStep;
      CurrentStep = rTimeStep;
      Start = rStartTime;
      End = rEndTime;
      MileStone = rTimeStep+rTimeStep;
    }

    void SetFactors(const double& rIncreaseFactor, const double& rDecreaseFactor, const double& rDelay)
    {
      IncreaseFactor = rIncreaseFactor;
      DecreaseFactor = rDecreaseFactor;
      UpdateDelay = rDelay;
    }

    bool Increase()
    {
      if( !ActiveDelay() ){
        PreviousStep = CurrentStep;
        CurrentStep  = PreviousStep + IncreaseFactor;
        if( InitialStep <= CurrentStep )
          CurrentStep == InitialStep;
        if( Total+CurrentStep >= MileStone ){
          CurrentStep == MileStone-Total;
          MileStone += InitialStep;
        }
      }
      return CheckSameTime(PreviousStep, CurrentStep);      
    }

    bool Decrease()
    {
      if( !ActiveDelay() ){
        PreviousStep = CurrentStep;
        CurrentStep  = PreviousStep - DecreaseFactor;
        if( CurrentStep <= 1e-2*InitialStep )
          CurrentStep = 1e-2*InitialStep;
      }
      return CheckSameTime(PreviousStep, CurrentStep); 
    }

    bool PredictActive()
    {
      if( PredictionLevel >= 0)
      {
        if( PredictionLevel == 0 && Total > Start )
          return false;
        else
          return true;
      }
    }

    bool Update()
    {
      if( PredictedStep > CurrentStep )
        return Increase();
      else
        return Decrease();
    }
    
    bool ActiveDelay()
    {
      if( DelayCounter == UpdateDelay )
      {
        DelayCounter = 0;
        return false;
      }
      else{
        ++DelayCounter;
        return true;
      }
    }

    bool CheckSameTime(const double& rTime, const double& rNewTime)
    {
      double tolerance = Initial * 1e-5;
      if( rNewTime > rTime-tolerance && rNewTime < rTime+tolerance )
        return true;
      else
        return false;      
    }

  }

  ///@}
  ///@name Protected Operators
  ///@{
  ///@}
  ///@name Protected Operations
  ///@{

  void ReduceTimeStep()
  {
    KRATOS_TRY   
    
    rCurrentProcessInfo[TIME] -= mTime.CurrentStep;
    rCurrentProcessInfo[STEP] -= 1;
    
    mrModelPart.ReduceTimeStep(rCurrentProcessInfo[TIME]);
    
    rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED, mTime.Decrease());
   
    KRATOS_CATCH("");
  }

  void UpdateTimeStep()
  {
    KRATOS_TRY
        
    this->PredictTimeStep();

    rCurrentProcessInfo.SetValue(DELTA_TIME_CHANGED, mTime.Update());
   
    KRATOS_CATCH("");
  }
  
  void PredictTimeStep()
  {
    KRATOS_TRY

    if( mTime.PredictActive() ){
      this->PredictTimeStep(mTime.PredictedStep);
    }
    
    KRATOS_CATCH("");
  }
  
  void PredictTimeStep(double& rTimeStep)
  {
    KRATOS_TRY

    KRATOS_CATCH("");
  }

  void CheckCriticalElement()
  {
    KRATOS_TRY
    
    KRATOS_CATCH("");
  }

    
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
  ///@name Private Static Member Variables
  ///@{

  ModelPart& mrModelPart;

  TimeParameters mTime;

  ///@}
  ///@name Private Member Variables
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
  TimeDiscretizationProcess& operator=(TimeDiscretizationProcess const& rOther);

  /// Copy constructor.
  //TimeDiscretizationProcess(TimeDiscretizationProcess const& rOther);


  ///@}

}; // Class TimeDiscretizationProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TimeDiscretizationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TimeDiscretizationProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TIME_DISCRETIZATION_PROCESS_H_INCLUDED  defined
