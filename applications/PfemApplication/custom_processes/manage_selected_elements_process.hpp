//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:          July 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_MANAGE_SELECTED_ELEMENTS_PROCESS_H_INCLUDED)
#define  KRATOS_MANAGE_SELECTED_ELEMENTS_PROCESS_H_INCLUDED


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

/// Process for managing the time integration of variables for selected elements
/**
   Selected elements are SLIVERS and not COMPUTED ELEMENTS
   Setting geometrical SLIVERS to NOT ACTIVE before computation
   Integrating variables for Not Computed elements and nodes
 */
class ManageSelectedElementsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef ModelPart::ElementType::GeometryType       GeometryType;

    /// Pointer definition of ManageSelectedElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ManageSelectedElementsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    ManageSelectedElementsProcess(ModelPart& rModelPart) : Process(Flags()) , mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~ManageSelectedElementsProcess() {}


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


    /// Execute method is used to execute the ManageSelectedElementsProcess algorithms.
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

      const int nelements = mrModelPart.NumberOfElements();

      if (nelements != 0)
      {
        ModelPart::ElementsContainerType::iterator it_begin = mrModelPart.ElementsBegin();

        #pragma omp parallel for
        for (int i = 0; i < nelements; ++i)
        {
          ModelPart::ElementsContainerType::iterator it = it_begin + i;

          if( it->Is(FLUID) ){

            GeometryType& rGeometry = it->GetGeometry();
            const unsigned int number_of_nodes = rGeometry.size();
            unsigned int selected_nodes = 0;
            for(unsigned int j=0; j<number_of_nodes; ++j)
            {
              if(rGeometry[j].Is(SELECTED))
                ++selected_nodes;
            }

            if(selected_nodes == number_of_nodes){
              it->Set(SELECTED,true);
              it->Set(ACTIVE,false);
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

      const int nelements = mrModelPart.NumberOfElements();

      if (nelements != 0)
      {
        ModelPart::ElementsContainerType::iterator it_begin = mrModelPart.ElementsBegin();

        //(set on nodes, not in parallel)
        for (int i = 0; i < nelements; ++i)
        {
          ModelPart::ElementsContainerType::iterator it = it_begin + i;

          if( it->Is(FLUID) && it->Is(SELECTED) ){

            GeometryType& rGeometry = it->GetGeometry();
            const unsigned int number_of_nodes = rGeometry.size();

            for(unsigned int j=0; j<number_of_nodes; ++j)
            {
              rGeometry[j].Set(SELECTED,true);
            }
            it->Set(SELECTED,false);
            it->Set(ACTIVE,true);
          }
        }
      }

      const int nnodes = mrModelPart.NumberOfNodes();

      if (nnodes != 0)
      {
        ModelPart::NodesContainerType::iterator it_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for (int i = 0; i < nnodes; ++i)
        {
          ModelPart::NodesContainerType::iterator it = it_begin + i;

          if( it->Is(SELECTED) ){

            if( norm_2(it->FastGetSolutionStepValue(VELOCITY)) > 1.5 * norm_2(it->FastGetSolutionStepValue(VELOCITY,1)) ||
                norm_2(it->FastGetSolutionStepValue(ACCELERATION)) > 1.5 * norm_2(it->FastGetSolutionStepValue(ACCELERATION,1)) ){

              // std::cout<<"PRE POSITION "<<it->Coordinates()<<"  ACCELERATION "<<it->FastGetSolutionStepValue(ACCELERATION)<<std::endl;
              // std::cout<<"DISPLACEMENT "<<it->FastGetSolutionStepValue(DISPLACEMENT)<<"  VELOCITY "<<it->FastGetSolutionStepValue(VELOCITY)<<std::endl;

              noalias(it->FastGetSolutionStepValue(VELOCITY)) = it->FastGetSolutionStepValue(VELOCITY,1);
              noalias(it->FastGetSolutionStepValue(ACCELERATION)) = it->FastGetSolutionStepValue(ACCELERATION,1);
              noalias(it->Coordinates()) -= it->FastGetSolutionStepValue(DISPLACEMENT);
              noalias(it->FastGetSolutionStepValue(DISPLACEMENT)) = it->FastGetSolutionStepValue(DISPLACEMENT,1);
              noalias(it->Coordinates()) += it->FastGetSolutionStepValue(DISPLACEMENT);

              // std::cout<<"POST POSITION "<<it->Coordinates()<<"  ACCELERATION "<<it->FastGetSolutionStepValue(ACCELERATION)<<std::endl;
              // std::cout<<"DISPLACEMENT "<<it->FastGetSolutionStepValue(DISPLACEMENT)<<"  VELOCITY "<<it->FastGetSolutionStepValue(VELOCITY)<<std::endl;

              //std::cout<<" SELECTED Node ["<<it->Id()<<"] Displacement"<<it->FastGetSolutionStepValue(DISPLACEMENT)<<std::endl;

            }

            //check if free surface nodes are outside the bounding box
            if( it->Is(FREE_SURFACE) ){
              if( !mBoundingBox.IsInside( it->Coordinates() ) ){
                it->Set(TO_ERASE,true);
                std::cout<<" SELECTED to erase "<<std::endl;
              }
            }

            it->Set(SELECTED,false);
          }


        }
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
        return "ManageSelectedElementsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ManageSelectedElementsProcess";
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
    ManageSelectedElementsProcess(ManageSelectedElementsProcess const& rOther);

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
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    ManageSelectedElementsProcess& operator=(ManageSelectedElementsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class ManageSelectedElementsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ManageSelectedElementsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ManageSelectedElementsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MANAGE_SELECTED_ELEMENTS_PROCESS_H_INCLUDED  defined
