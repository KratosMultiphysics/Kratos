//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:          July 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_MANAGE_ISOLATED_NODES_PROCESS_H_INCLUDED)
#define  KRATOS_MANAGE_ISOLATED_NODES_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_bounding/spatial_bounding_box.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Process for managing the time integration of variables for isolated nodes
class ManageIsolatedNodesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ManageIsolatedNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(ManageIsolatedNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    ManageIsolatedNodesProcess(ModelPart& rModelPart) : Process(Flags()) , mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~ManageIsolatedNodesProcess() {}


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


    /// Execute method is used to execute the ManageIsolatedNodesProcess algorithms.
    void Execute() override
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
      KRATOS_TRY

      double Radius = 0.0;
      //BOUNDARY flag must be set in model part nodes
      mBoundingBox = SpatialBoundingBox(mrModelPart,Radius);

      KRATOS_CATCH("")
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
      KRATOS_TRY

      const int nnodes = mrModelPart.Nodes().size();

      if (nnodes != 0)
      {
        ModelPart::NodesContainerType::iterator it_begin =
            mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
          ModelPart::NodesContainerType::iterator it = it_begin + i;

          if( it->Is(ISOLATED) || (it->Is(RIGID) && (it->IsNot(SOLID) && it->IsNot(FLUID))) ){

            if(it->SolutionStepsDataHas(PRESSURE)){
              it->FastGetSolutionStepValue(PRESSURE,0) = 0.0;
              it->FastGetSolutionStepValue(PRESSURE,1) = 0.0;
            }
            if(it->SolutionStepsDataHas(PRESSURE_VELOCITY)){
              it->FastGetSolutionStepValue(PRESSURE_VELOCITY,0) = 0.0;
              it->FastGetSolutionStepValue(PRESSURE_VELOCITY,1) = 0.0;
            }
            if(it->SolutionStepsDataHas(PRESSURE_ACCELERATION)){
              it->FastGetSolutionStepValue(PRESSURE_ACCELERATION,0) = 0.0;
              it->FastGetSolutionStepValue(PRESSURE_ACCELERATION,1) = 0.0;
            }

          }
        }
      }

      KRATOS_CATCH("")
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
      KRATOS_TRY

      const int nnodes = mrModelPart.Nodes().size();
      const double TimeStep = mrModelPart.GetProcessInfo()[DELTA_TIME];
      if (nnodes != 0)
      {
        this->SetSemiIsolatedNodes(mrModelPart);

        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
          ModelPart::NodesContainerType::iterator it = it_begin + i;

          if( it->Is(ISOLATED) ){
            if(it->SolutionStepsDataHas(VOLUME_ACCELERATION)){
              const array_1d<double,3>& VolumeAcceleration = it->FastGetSolutionStepValue(VOLUME_ACCELERATION);
              noalias(it->FastGetSolutionStepValue(ACCELERATION)) = VolumeAcceleration;
            }
            // std::cout<<"PRE POSITION "<<it->Coordinates()<<"  ACCELERATION "<<it->FastGetSolutionStepValue(ACCELERATION)<<std::endl;
            // std::cout<<"DISPLACEMENT "<<it->FastGetSolutionStepValue(DISPLACEMENT)<<"  VELOCITY "<<it->FastGetSolutionStepValue(VELOCITY)<<std::endl;

            noalias(it->FastGetSolutionStepValue(VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY,1) + it->FastGetSolutionStepValue(ACCELERATION)*TimeStep;
            noalias(it->Coordinates()) -= it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(it->FastGetSolutionStepValue(DISPLACEMENT)) = it->FastGetSolutionStepValue(DISPLACEMENT,1) + it->FastGetSolutionStepValue(VELOCITY)*TimeStep + 0.5*it->FastGetSolutionStepValue(ACCELERATION)*TimeStep*TimeStep;
            noalias(it->Coordinates()) += it->FastGetSolutionStepValue(DISPLACEMENT);

            // std::cout<<"POST POSITION "<<it->Coordinates()<<"  ACCELERATION "<<it->FastGetSolutionStepValue(ACCELERATION)<<std::endl;
            // std::cout<<"DISPLACEMENT "<<it->FastGetSolutionStepValue(DISPLACEMENT)<<"  VELOCITY "<<it->FastGetSolutionStepValue(VELOCITY)<<std::endl;

            //std::cout<<" ISOLATED Node ["<<it->Id()<<"] Displacement"<<it->FastGetSolutionStepValue(DISPLACEMENT)<<std::endl;

            if( !mBoundingBox.IsInside( it->Coordinates() ) ){
              it->Set(TO_ERASE);
              std::cout<<" ISOLATED to erase "<<std::endl;
            }
          }
          else if( it->Is(VISITED) ){

            if(it->SolutionStepsDataHas(VOLUME_ACCELERATION)){
              const array_1d<double,3>& VolumeAcceleration = it->FastGetSolutionStepValue(VOLUME_ACCELERATION);
              noalias(it->FastGetSolutionStepValue(ACCELERATION)) = VolumeAcceleration;
            }

            noalias(it->FastGetSolutionStepValue(VELOCITY)) = 0.5 * (it->FastGetSolutionStepValue(VELOCITY) + it->FastGetSolutionStepValue(VELOCITY,1) + it->FastGetSolutionStepValue(ACCELERATION)*TimeStep);
            noalias(it->Coordinates()) -= it->FastGetSolutionStepValue(DISPLACEMENT);
            noalias(it->FastGetSolutionStepValue(DISPLACEMENT)) = 0.5 * (it->FastGetSolutionStepValue(DISPLACEMENT) + it->FastGetSolutionStepValue(DISPLACEMENT,1) + it->FastGetSolutionStepValue(VELOCITY)*TimeStep + 0.5*it->FastGetSolutionStepValue(ACCELERATION)*TimeStep*TimeStep);
            noalias(it->Coordinates()) += it->FastGetSolutionStepValue(DISPLACEMENT);

          }
        }

        this->ResetSemiIsolatedNodes(mrModelPart);

      }

      KRATOS_CATCH("")
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
        return "ManageIsolatedNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ManageIsolatedNodesProcess";
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

    /// Copy constructor.
    ManageIsolatedNodesProcess(ManageIsolatedNodesProcess const& rOther);

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

    SpatialBoundingBox mBoundingBox;

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{

    void SetSemiIsolatedNodes(ModelPart& rModelPart)
    {

     const int nnodes = mrModelPart.Nodes().size();

      if (nnodes != 0)
      {
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
          ModelPart::NodesContainerType::iterator it = it_begin + i;

          if( it->Is(FREE_SURFACE) ){

            WeakPointerVector<Node<3> >& rN = it->GetValue(NEIGHBOUR_NODES);
            unsigned int NumberOfNeighbours = rN.size();
            unsigned int rigid = 0;
            for(unsigned int j = 0; j < NumberOfNeighbours; ++j)
	    {
              if(rN[j].Is(RIGID))
                ++rigid;
	    }
            if( rigid == NumberOfNeighbours )
              it->Set(VISITED,true);
          }
        }
      }


    }

    void ResetSemiIsolatedNodes(ModelPart& rModelPart)
    {

      const int nnodes = mrModelPart.Nodes().size();

      if (nnodes != 0)
      {
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
          ModelPart::NodesContainerType::iterator it = it_begin + i;

          if(it->Is(VISITED))
            it->Set(VISITED,false);
        }
      }

    }


    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    ManageIsolatedNodesProcess& operator=(ManageIsolatedNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ManageIsolatedNodesProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ManageIsolatedNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ManageIsolatedNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MANAGE_ISOLATED_NODES_PROCESS_H_INCLUDED  defined
