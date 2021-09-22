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

#ifndef KRATOS_DEFINE_EMBEDDED_WAKE_PROCESS_H
#define KRATOS_DEFINE_EMBEDDED_WAKE_PROCESS_H

#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include <unordered_set>

namespace Kratos
{

	class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) DefineEmbeddedWakeProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(DefineEmbeddedWakeProcess);

    typedef Node<3> NodeType;
    typedef Geometry<NodeType>::Pointer GeomPointerType;

    // Constructor for DefineEmbeddedWakeProcess Process
    DefineEmbeddedWakeProcess(ModelPart& rModelPart, ModelPart& rWakeModelPart);

    /// Destructor.
    ~DefineEmbeddedWakeProcess() = default;

    /// Assignment operator.
    DefineEmbeddedWakeProcess& operator=(DefineEmbeddedWakeProcess const& rOther) = delete;

    /// Copy constructor.
    DefineEmbeddedWakeProcess(DefineEmbeddedWakeProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    void Execute() override;

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "DefineEmbeddedWakeProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DefineEmbeddedWakeProcess";
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
    NodeType mrTrailingEdgeNode;
    array_1d<double,3> mWakeOrigin;

    void ComputeDistanceToWake();

    void MarkWakeElements();

    void ComputeTrailingEdgeNode();

    void MarkKuttaWakeElements();

    void RedefineWake();

    NodeType* FindNode(double restarted_search_maximum);

    bool CheckIfValid(NodeType* trailing_edge_candidate);

    bool TouchesWake(Element& rElem, std::unordered_set<std::size_t> visited_elements);

    ModifiedShapeFunctions::Pointer pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const;

}; // Class Process
} // namespace Kratos


#endif // KRATOS_DEFINE_EMBEDDED_WAKE_PROCESS_H
