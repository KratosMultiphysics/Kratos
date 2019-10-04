//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Marc Nunez,
//

#ifndef KRATOS_DEFINE_EMBEDDED_WAKE_3D_PROCESS_H
#define KRATOS_DEFINE_EMBEDDED_WAKE_3D_PROCESS_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

	class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) DefineEmbeddedWake3DProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(DefineEmbeddedWake3DProcess);

    // Constructor for DefineEmbeddedWake3DProcess Process
    DefineEmbeddedWake3DProcess(ModelPart& rModelPart, ModelPart& rWakeModelPart);

    /// Destructor.
    ~DefineEmbeddedWake3DProcess() = default;

    /// Assignment operator.
    DefineEmbeddedWake3DProcess& operator=(DefineEmbeddedWake3DProcess const& rOther) = delete;

    /// Copy constructor.
    DefineEmbeddedWake3DProcess(DefineEmbeddedWake3DProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    void Execute() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DefineEmbeddedWake3DProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DefineEmbeddedWake3DProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->PrintInfo(rOStream);
    }

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    ModelPart& mrWakeModelPart;
    GlobalPointersVector<Element> mKuttaWakeElementCandidates;

    void ComputeDistanceToWake();

    void MarkWakeElements();

    void ComputeTrailingEdgeNode();

    void MarkKuttaWakeElements();

}; // Class Process
} // namespace Kratos


#endif // KRATOS_DEFINE_EMBEDDED_WAKE_PROCESS_H
