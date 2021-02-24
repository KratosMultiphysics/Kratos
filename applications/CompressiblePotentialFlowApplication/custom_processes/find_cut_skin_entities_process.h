//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, Inigo Lopez based on R.Rossi and V.Mataix work
//
//

#if !defined(KRATOS_FIND_CUT_SKIN_ENTITIES_PROCESS_INCLUDED)
#define  KRATOS_FIND_CUT_SKIN_ENTITIES_PROCESS_INCLUDED

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{


class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) FindCutSkinEntitiesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindCutSkinEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindCutSkinEntitiesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    FindCutSkinEntitiesProcess(
        ModelPart& rModelPart,
        ModelPart& rSectionModelPart,
        const array_1d<double,3>& rVersor,
        const array_1d<double,3>& rOrigin);

    /// Destructor.
    ~FindCutSkinEntitiesProcess() override
    {
    }


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

    /**
     * Execute method is used to execute the Process algorithms.
     * In this process the gradient of a scalar variable will be computed
     */
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindCutSkinEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindCutSkinEntitiesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; // The main model part
    ModelPart& mrSectionModelPart; // The main model part
    const array_1d<double,3>& mrVersor; // The plane normal
    const array_1d<double,3>& mrOrigin; // A point of the plane

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindCutSkinEntitiesProcess& operator=(FindCutSkinEntitiesProcess const& rOther);

    /// Copy constructor.
    //FindCutSkinEntitiesProcess(FindCutSkinEntitiesProcess const& rOther);


    ///@}

}; // Class FindCutSkinEntitiesProcess

///@}
}  // namespace Kratos.

#endif // KRATOS_COMPUTE_NODAL_VALUE_PROCESS_INCLUDED  defined


