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

#include <unordered_set>
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{

	class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) DefineEmbeddedWakeProcess3D : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(DefineEmbeddedWakeProcess3D);

    // Constructor for DefineEmbeddedWakeProcess3D Process
    DefineEmbeddedWakeProcess3D(ModelPart& rModelPart, ModelPart& rWakeModelPart);

    /// Destructor.
    ~DefineEmbeddedWakeProcess3D() = default;

    /// Assignment operator.
    DefineEmbeddedWakeProcess3D& operator=(DefineEmbeddedWakeProcess3D const& rOther) = delete;

    /// Copy constructor.
    DefineEmbeddedWakeProcess3D(DefineEmbeddedWakeProcess3D const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    void ExecuteInitialize() override;

    void Execute() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DefineEmbeddedWakeProcess3D";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DefineEmbeddedWakeProcess3D";
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

    void ComputeDistanceToWake();

    void MarkWakeElements();

    void ComputeTrailingEdgeNode();

    bool TouchesWake(Element& rElem, std::unordered_set<std::size_t> visited_elements);

    // std::vector<IndexType> GetTrailingEdgeNodeList();

}; // Class Process
} // namespace Kratos


#endif // KRATOS_DEFINE_EMBEDDED_WAKE_PROCESS_H
