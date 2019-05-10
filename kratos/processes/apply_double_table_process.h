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

#if !defined(KRATOS_APPLY_DOUBLE_TABLE_PROCESS_H_INCLUDED)
#define  KRATOS_APPLY_DOUBLE_TABLE_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/apply_component_table_process.h"

namespace Kratos
{

/// The base class for all processes in Kratos.
/** This function applies a value according to a table defined in the mdpa
*/
class KRATOS_API(KRATOS_CORE) ApplyDoubleTableProcess : public ApplyComponentTableProcess
{
    
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyConstantScalarValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyDoubleTableProcess);

    /// Constructor
    ApplyDoubleTableProcess(ModelPart& rModelPart, Parameters Parameters);

    
    /// Destructor
    virtual ~ApplyDoubleTableProcess() override {}
    
    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ApplyDoubleTableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyDoubleTableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

protected:

    /// Copy constructor.
    ApplyDoubleTableProcess(ApplyDoubleTableProcess const& rOther);
private:

    /// Assignment operator.
    ApplyDoubleTableProcess& operator=(ApplyDoubleTableProcess const& rOther);
  
}; // Class ApplyDoubleTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyDoubleTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyDoubleTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif /* KRATOS_APPLY_DOUBLE_TABLE_PROCESS defined */
