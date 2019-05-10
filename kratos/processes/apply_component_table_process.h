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
//  Collaborator:    Alejandro Cornejo
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

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for all processes in Kratos.
/** This function applies a table value to a component
*/
class KRATOS_API(KRATOS_CORE) ApplyComponentTableProcess : public Process
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyConstantScalarValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyComponentTableProcess);
    
    /// Defining a table with double argument and result type as table type.
    typedef Table<double,double> TableType;
    
    /// Constructor
	ApplyComponentTableProcess(ModelPart& rModelPart, Parameters rParameters);
   
    /// Destructor
    virtual ~ApplyComponentTableProcess() override {}
    
    /// this function is designed for being called at the beginning of the computations
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
    ModelPart& mrModelPart;
    std::string mVariableName;
    bool mIsFixed;
    double mInitialValue;
    TableType::Pointer mpTable;
    
private:

    /// Assignment operator.
    ApplyComponentTableProcess& operator=(ApplyComponentTableProcess const& rOther);
    
}; // Class ApplyComponentTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyComponentTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyComponentTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_COMPONENT_TABLE_PROCESS defined */
