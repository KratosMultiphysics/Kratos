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
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

// Variables types
typedef VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>  ComponentVariableType;
typedef Variable<double>                                                  DoubleVariableType;

/// The base class for all processes in Kratos.
/** This function applies a table value to a component
*/
template<class TVariableType>
class KRATOS_API(KRATOS_CORE) ApplyComponentTableProcess : public Process
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyConstantScalarValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyComponentTableProcess);
    
    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;
    
    /**
     * @brief Default constructor.
     * @param rModelPart The model part where the search of neighbours is performed
     * @param ThisParameters The parameters of configuration
     */
	ApplyComponentTableProcess(ModelPart& rModelPart, Parameters rParameters);
   
    /// Destructor
    virtual ~ApplyComponentTableProcess() override {}
    
    /// This function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

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

    /// Copy constructor.
    ApplyComponentTableProcess(ApplyComponentTableProcess const& rOther);

    /// Member Variables
    ModelPart& mrModelPart;  // The main model part
    std::string mVariableName;  // The variable to be fixed with the value of the table
    std::string mTimeVariableName;  // The variable that represents "time"
    bool mIsFixed;  // Boolean var to know wether the DoF is fixed
    double mInitialValue;  // Initial value of the variable
    TableType::Pointer mpTable;  // Table containing all the information
    
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
