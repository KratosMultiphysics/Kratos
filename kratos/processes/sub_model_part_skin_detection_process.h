//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_SUB_MODEL_PART_SKIN_DETECTION_PROCESS_H_INCLUDED)
#define KRATOS_SUB_MODEL_PART_SKIN_DETECTION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "skin_detection_process.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Create a SubModelPart covering a part of the outside skin of the computation domain where a condition is met.
/** For example, create the outer skin for the part of the domain belonging to a given SubModelPart.
*/
template<SizeType TDim>
class KRATOS_API(KRATOS_CORE) SubModelPartSkinDetectionProcess: public SkinDetectionProcess<TDim>
{

KRATOS_DEFINE_LOCAL_FLAG( NODE_SELECTED );

// Internal class used to select which faces to create.
class FaceSelector
{
public:
KRATOS_CLASS_POINTER_DEFINITION(FaceSelector);
virtual void Prepare(ModelPart& rMainModelPart) const = 0;
virtual bool IsSelected(const Geometry<Node<3>>::PointsArrayType&) const = 0;
};

// Select faces where all nodes belong to given SubModelPart.
class SelectIfAllNodesOnSubModelPart: public FaceSelector
{
std::string mName;
public:
SelectIfAllNodesOnSubModelPart(const std::string& rName): mName(rName) {}

void Prepare(ModelPart& rMainModelPart) const override
{
    ModelPart& r_model_part = rMainModelPart.GetSubModelPart(mName);
    auto node_begin = r_model_part.NodesBegin();
    const int num_nodes = r_model_part.NumberOfNodes();

    #pragma omp parallel for
    for (int k = 0; k < num_nodes; k++)
    {
        (node_begin+k)->Set(SubModelPartSkinDetectionProcess::NODE_SELECTED);
    }
}

bool IsSelected(const Geometry<Node<3>>::PointsArrayType& rNodes) const override
{
    bool select = true;
    for (auto i_node = rNodes.begin(); i_node != rNodes.end(); ++i_node)
    {
        select &= i_node->Is(SubModelPartSkinDetectionProcess::NODE_SELECTED);
    }
    return select;
}

};

public:
///@name Type Definitions
///@{

/// Pointer definition of SubModelPartSkinDetectionProcess
KRATOS_CLASS_POINTER_DEFINITION(SubModelPartSkinDetectionProcess);

using typename SkinDetectionProcess<TDim>::HashMapVectorIntType;
using typename SkinDetectionProcess<TDim>::HashMapVectorIntIdsType;
using typename SkinDetectionProcess<TDim>::VectorIndexType;

using ConditionCheckType = bool(const Geometry<Node<3>>::PointsArrayType&);

///@}
///@name Life Cycle
///@{

/// Constructor
SubModelPartSkinDetectionProcess(ModelPart& rModelPart, Parameters Settings);

/// Deleted default constructor.
SubModelPartSkinDetectionProcess() = delete;

/// Deleted copy constructor.
SubModelPartSkinDetectionProcess(SubModelPartSkinDetectionProcess const &rOther) = delete;

/// Destructor.
~SubModelPartSkinDetectionProcess() override = default;

///@}
///@name Operators
///@{

/// Deleted sssignment operator.
SubModelPartSkinDetectionProcess &operator=(SubModelPartSkinDetectionProcess const &rOther) = delete;

///@}
///@name Operations
///@{

void Execute() override;

///@}
///@name Input and output
///@{

std::string Info() const override
{
    return "SkinDetectionProcess";
}

/// Print information about this object.
void PrintInfo(std::ostream& rOStream) const override
{
    rOStream << "SkinDetectionProcess";
}

/// Print object's data.
void PrintData(std::ostream& rOStream) const override
{
}

///@}

protected:

///@name Protected Operations
///@{

void CreateConditions(
    ModelPart& rMainModelPart,
    ModelPart& rSkinModelPart,
    HashMapVectorIntType& rInverseFaceMap,
    HashMapVectorIntIdsType& rPropertiesFaceMap,
    std::unordered_set<IndexType>& rNodesInTheSkin,
    const std::string& rConditionName) const override;

Parameters GetDefaultSettings() const override;

///@}

private:

///@name Member Variables
///@{

typename FaceSelector::Pointer mpFaceSelector;

///@}
///@name Private Operations
///@{

static bool FaceIsNeeded(const Geometry<Node<3>>::PointsArrayType&)
{
    return true;
}

///@}

}; // Class SubModelPartSkinDetectionProcess

///@}
///@name Input and output
///@{

/// input stream function
template<SizeType TDim>
inline std::istream &operator>>(std::istream &rIStream,
                                SubModelPartSkinDetectionProcess<TDim> &rThis);

/// output stream function
template<SizeType TDim>
inline std::ostream &operator<<(std::ostream &rOStream,
                                const SubModelPartSkinDetectionProcess<TDim> &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@}

} // namespace Kratos.

#endif // KRATOS_SUB_MODEL_PART_SKIN_DETECTION_PROCESS_H_INCLUDED  defined
