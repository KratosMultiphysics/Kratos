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

/**
 * @class SubModelPartSkinDetectionProcess
 * @brief Create a SubModelPart covering a part of the outside skin of the computation domain where a condition is met.
 * @details For example, create the outer skin for the part of the domain belonging to a given SubModelPart.
 */
template<SizeType TDim>
class KRATOS_API(KRATOS_CORE) SubModelPartSkinDetectionProcess: public SkinDetectionProcess<TDim>
{

KRATOS_DEFINE_LOCAL_FLAG( NODE_SELECTED );

/**
 * @class FaceSelector
 * @brief Internal class used to select which faces to create.
 */
class FaceSelector
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FaceSelector);
    virtual ~FaceSelector() = default;
    virtual void Prepare(ModelPart& rMainModelPart) const = 0;
    virtual bool IsSelected(const Geometry<Node>::PointsArrayType&) const = 0;
};

/**
 * @class SelectIfOneNodeNotOnSubModelPart
 * @brief Select faces where all nodes belong to given SubModelPart.
 * @see FaceSelector
 */
class SelectIfAllNodesOnSubModelPart: public FaceSelector
{
    std::string mName;

public:
    SelectIfAllNodesOnSubModelPart(const std::string& rName): mName(rName) {}

    void Prepare(ModelPart& rMainModelPart) const override;

    bool IsSelected(const Geometry<Node>::PointsArrayType& rNodes) const override;
};

/**
 * @class SelectIfOneNodeNotOnSubModelPart
 * @brief Select faces where almost one node does not belong to the given SubModelParts.
 * @see FaceSelector
 * @author Miguel Maso
 */
class SelectIfOneNodeNotOnSubModelPart: public FaceSelector
{
    std::vector<std::string> mNames;

public:
    SelectIfOneNodeNotOnSubModelPart(const std::vector<std::string>& rNames): mNames(rNames) {}

    void Prepare(ModelPart& rMainModelPart) const override;

    bool IsSelected(const Geometry<Node>::PointsArrayType& rNodes) const override;
};

public:
///@name Type Definitions
///@{

/// Pointer definition of SubModelPartSkinDetectionProcess
KRATOS_CLASS_POINTER_DEFINITION(SubModelPartSkinDetectionProcess);

using typename SkinDetectionProcess<TDim>::HashMapVectorIntType;
using typename SkinDetectionProcess<TDim>::HashMapVectorIntIdsType;
using typename SkinDetectionProcess<TDim>::VectorIndexType;

using ConditionCheckType = bool(const Geometry<Node>::PointsArrayType&);

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

/**
 * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
 */
const Parameters GetDefaultParameters() const override;

///@}

private:

///@name Member Variables
///@{

typename FaceSelector::Pointer mpFaceSelector;

///@}
///@name Private Operations
///@{

static bool FaceIsNeeded(const Geometry<Node>::PointsArrayType&)
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
