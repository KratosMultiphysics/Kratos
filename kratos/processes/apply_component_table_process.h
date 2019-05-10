//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Alejandro Cornejo
//

#if !defined(KRATOS_APPLY_COMPONENT_TABLE_PROCESS_H_INCLUDED)
#define  KRATOS_APPLY_COMPONENT_TABLE_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/table.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

// Variables types
typedef VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>  ComponentVariableType;
typedef Variable<double>                                                  DoubleVariableType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @ingroup KratosCore
 * @class ApplyComponentTableProcess
 * @brief This class applies a value from a table to a BC o load
 * @details This function applies a table value to a component
 * @author Ignasi de Pouplana
 * @author Alejandro Cornejo
 */
template<class TVariableType>
class KRATOS_API(KRATOS_CORE) ApplyComponentTableProcess
    : public Process
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyConstantScalarValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyComponentTableProcess);

    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the search of neighbours is performed
     * @param ThisParameters The parameters of configuration
     */
    ApplyComponentTableProcess(ModelPart& rModelPart, Parameters rParameters);

    /// Destructor
    virtual ~ApplyComponentTableProcess() override {}

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

    /**
     * @brief This method executes the algorithm that looks for neighbour nodes and elements in a  mesh of prismatic elements
     */
    void Execute() override;

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override;

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
        return "ApplyComponentTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyComponentTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPart;             /// The main model part
    std::string mVariableName;          /// The variable to be fixed with the value of the table
    std::string mTimeVariableName;      /// The variable that represents "time"
    bool mIsFixed;                      /// Boolean var to know wether the DoF is fixed
    double mInitialValue;               /// Initial value of the variable
    TableType::Pointer mpTable;         /// Table containing all the information

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Auxiliar implementation of ExecuteInitialize to reduce code duplication
     */
    void ImplementationExecuteInitialize(const TVariableType& rVariable);

    /**
     * @brief Auxiliar implementation of ExecuteInitializeSolutionStep to reduce code duplication
     */
    void ImplementationExecuteInitializeSolutionStep(const TVariableType& rVariable);

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    /// Copy constructor.
    ApplyComponentTableProcess(ApplyComponentTableProcess const& rOther);

    ///@}

private:

    /// Assignment operator.
    ApplyComponentTableProcess& operator=(ApplyComponentTableProcess const& rOther);

}; // Class ApplyComponentTableProcess

/// input stream function
template<class TVariableType>
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyComponentTableProcess<TVariableType>& rThis);

/// output stream function
template<class TVariableType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyComponentTableProcess<TVariableType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_COMPONENT_TABLE_PROCESS defined */
