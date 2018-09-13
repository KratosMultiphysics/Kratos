//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:          July 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_ASSIGN_PROPERTIES_TO_NODES_PROCESS_H_INCLUDED)
#define  KRATOS_ASSIGN_PROPERTIES_TO_NODES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

#include "pfem_application_variables.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for fixing scalar variable Dof or array_1d component Dof processes in Kratos.
/** This function fix the variable dof belonging to all of the nodes in a given mesh
*/
class AssignPropertiesToNodesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// NodeType
    typedef Node<3> NodeType;

    typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;
    typedef typename PropertiesContainerType::Pointer   PropertiesContainerPointerType;

    /// Pointer definition of AssignPropertiesToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignPropertiesToNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    AssignPropertiesToNodesProcess(ModelPart& model_part,
                                   Parameters rParameters
                                   ) : Process() , mrModelPart(model_part)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
            {
                "fluid_mixture": false,
                "solid_mixture": false
            }  )" );


        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mFluidMixture = rParameters["fluid_mixture"].GetBool();
        mSolidMixture = rParameters["solid_mixture"].GetBool();

        mpProperties  = mrModelPart.pProperties();

        KRATOS_CATCH("");
    }


    /// Destructor.
    virtual ~AssignPropertiesToNodesProcess() {}


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
      KRATOS_TRY

      this->AssignPropertiesToNodes();

      KRATOS_CATCH("")
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
      KRATOS_TRY

      this->AssignMaterialPercentageToNodes();

      KRATOS_CATCH("")
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
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
        return "AssignPropertiesToNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "AssignPropertiesToNodesProcess";
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
    AssignPropertiesToNodesProcess(AssignPropertiesToNodesProcess const& rOther);

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

    bool mFluidMixture;

    bool mSolidMixture;

    PropertiesContainerPointerType mpProperties;

    ///@}
    ///@name Private Operators
    ///@{


    void AssignPropertiesToNodes()
    {
        const int nnodes = mrModelPart.GetMesh().Nodes().size();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

            //#pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
		it->SetValue(PROPERTIES_VECTOR,mpProperties);
            }
        }
    }

    void AssignMaterialPercentageToNodes()
    {
        const int nnodes = mrModelPart.GetMesh().Nodes().size();

        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

            //#pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                Vector MaterialPercentage;
                this->CalculateMaterialPercentage(*it, MaterialPercentage);
		it->SetValue(MATERIAL_PERCENTAGE, MaterialPercentage);
            }
        }
    }


    void CalculateMaterialPercentage(NodeType& rNode, Vector& MaterialPercentage)
    {
      KRATOS_TRY

      unsigned int size = mpProperties->size();
      MaterialPercentage.resize(size,false);
      noalias(MaterialPercentage) = ZeroVector(size);

      double counter = 0;
      if( rNode.Is(FLUID) && mFluidMixture ){

        WeakPointerVector<Element >& rE = rNode.GetValue(NEIGHBOUR_ELEMENTS);

        for(unsigned int i = 0; i < rE.size(); i++)
        {
          if(rE[i].Is(FLUID)){
            unsigned int id = rE[i].GetProperties().Id();
            if( id < size ){
              MaterialPercentage[id] += 1;
              ++counter;
            }
          }
        }
      }
      else if( rNode.Is(SOLID) && mSolidMixture ){

        WeakPointerVector<Element >& rE = rNode.GetValue(NEIGHBOUR_ELEMENTS);
        for(unsigned int i = 0; i < rE.size(); i++)
        {
          if(rE[i].Is(SOLID)){
            unsigned int id = rE[i].GetProperties().Id();
            if( id < size ){
              MaterialPercentage[id] += 1;
              ++counter;
            }
          }
        }

      }

      double divider = 1.0;
      if( counter != 0 )
        divider = 1.0/counter;

      for(unsigned int i=0; i<size; ++i)
        MaterialPercentage[i] *= divider;

      KRATOS_CATCH("")
    }

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    AssignPropertiesToNodesProcess& operator=(AssignPropertiesToNodesProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignPropertiesToNodesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  AssignPropertiesToNodesProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AssignPropertiesToNodesProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_PROPERTIES_TO_NODES_PROCESS_H_INCLUDED  defined
