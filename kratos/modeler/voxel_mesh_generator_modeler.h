//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part_io.h"
#include "modeler/modeler.h"
#include "internals/cartesian_mesh_colors.h"
#include "internals/skin_intersection.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) VoxelMeshGeneratorModeler
    : public Modeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(VoxelMeshGeneratorModeler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    VoxelMeshGeneratorModeler() ;

    /// Constructor.
    VoxelMeshGeneratorModeler(
        Model& rModel,
        Parameters ModelerParameters = Parameters());

    /// Destructor.
    ~VoxelMeshGeneratorModeler() override = default;

    /// Creates the Modeler Pointer
        /// Creates the Modeler Pointer
    Modeler::Pointer Create(
        Model& rModel, const Parameters ModelParameters) const override
    {
        return Kratos::make_shared<VoxelMeshGeneratorModeler>(rModel, ModelParameters);
    }

    ///@}
    ///@name Stages
    ///@{

    void SetupModelPart() override ;

    ModelPart& ReadModelParts();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Access
    ///@{

    /// Returns the key planes in given direction
    const std::vector<double>& GetKeyPlanes(std::size_t Direction) const {
        return mKeyPlanes[Direction];
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "VoxelMeshGeneratorModeler";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream & rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream & rOStream) const override
    {
    }

    const Parameters GetDefaultParameters() const override ;
    const Parameters GetDefaultKeyPlaneGeneratorParameters() const;

    ///@}
private:
    // Data that I store in the internal data structure for each node of the cartesian grid
    class CartesianNodalData
    {
        Node::Pointer mpNode = nullptr;

    public:
        CartesianNodalData() : mpNode(nullptr) {}
        Node::Pointer pGetNode() const { return mpNode; }
        void pSetNode(Node::Pointer pNode) { mpNode = pNode; }
        bool IsCreated() const { return mpNode != nullptr; }
    };

    // A cartesian data storage for the nodes of the cartesian mesh.
    // Is used as internal storage
    class CartesianData
    {
        std::vector<CartesianNodalData> mNodalData;
        array_1d<std::size_t, 3> mNumberOfDivisions;

    public:
        CartesianData() = default;

        void SetNumberOfDivisions(std::size_t XDivisions, std::size_t YDivisions, std::size_t ZDivisions)
        {
            mNumberOfDivisions[0] = XDivisions;
            mNumberOfDivisions[1] = YDivisions;
            mNumberOfDivisions[2] = ZDivisions;

            KRATOS_ERROR_IF((mNumberOfDivisions[0] < 1) || (mNumberOfDivisions[1] < 1) || (mNumberOfDivisions[2] < 1))
                << "The number of division cannot be zero" << std::endl;

            std::size_t size = mNumberOfDivisions[0] * mNumberOfDivisions[1] * mNumberOfDivisions[2];
            mNodalData.assign(size, CartesianNodalData());
        }

        std::size_t GetNodeIndex(std::size_t I, std::size_t J, std::size_t K) const
        {
            const std::size_t index = I + J * mNumberOfDivisions[0] + K * mNumberOfDivisions[1] * mNumberOfDivisions[0];
            return index;
        }

        CartesianNodalData &GetNodalData(std::size_t I, std::size_t J, std::size_t K)
        {
            const std::size_t index = GetNodeIndex(I,J,K);
            return mNodalData[index];
        }

        const array_1d<std::size_t, 3> &GetNumberOfDivisions()
        {
            return mNumberOfDivisions;
        }

        bool IsInitialized() const
        {
            return mNumberOfDivisions[0] && mNumberOfDivisions[1] && mNumberOfDivisions[2];
        }
    };

    struct ContactContainer
    {
        std::unordered_map<int,std::vector<std::size_t>> mNodeMap;
        std::unordered_map<int,std::vector<std::size_t>> mConditionsMap;
        std::unordered_map<int,ModelPart*> mModelPartMap;
        std::vector<int> mColorVector;
    };

    ///@name Private members
    ///@{

    std::size_t mStartNodeId = 0;
    std::size_t mStartElementId = 0;
    std::size_t mStartConditionId = 0;

    Model* mpModel = nullptr;
    ModelPart* mpInputModelPart = nullptr;
    array_1d<std::vector<double>,3> mKeyPlanes;
    std::map<std::size_t,ModelPart&> mModelPartsColors;
    Internals::CartesianMeshColors mColors;
    CartesianData mMeshingData;
    CartesianData mQuadraticNodeData{};
    std::map<int, Internals::SkinIntersection> mSkinIntersections;

    ///@}
    ///@name Private Operations
    ///@{

    /// This initializes de internal cartesian mesh data structure to be used for coloring
    void PreparingTheInternalDataStructure();

    /// Goes over the coloring array and apply them to the internal data structure
    void ApplyColoring(Parameters ColoringParameters, int OutsideColor);

    void SetStartIds(ModelPart& rTheVolumeModelPart);

    ModelPart& CreateAndGetModelPart(std::string const& FullName);

    void GenerateEntities(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters);

    Node::Pointer GenerateOrRetrieveNode(ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, const std::size_t I, const std::size_t J, const std::size_t K);

    Node::Pointer GenerateOrRetrieveQuadraticNode(ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, const std::size_t I, const std::size_t J, const std::size_t K);

    Node::Pointer GenerateNode(ModelPart& rTheVolumeModelPart, const Point& rCoordinates);

    void ApplyOperations(Parameters ThisParameters);

    void GenerateIntersectedSkinElements(
        ModelPart& rIntersectionModelPart,
        double Color,
        std::size_t PropertiesId,
        const Internals::SkinIntersection& rSkinIntersection);

    Internals::SkinIntersection& GetIntersections(int Color);

    const Internals::SkinIntersection& GetIntersections(int Color) const;

    bool IntersectionsGenerated(int Color) const;

    ///@}
    ///@name Friends
    ///@{

    friend class VoxelMesherColoring;

    friend class VoxelMesherKeyPlaneGeneration;

    friend class VoxelMesherEntityGeneration;

    friend class VoxelMesherOperation;

    ///@}
}; // Class VoxelMeshGeneratorModeler

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    VoxelMeshGeneratorModeler& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const VoxelMeshGeneratorModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.
