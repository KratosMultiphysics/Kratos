//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//
#if !defined(KRATOS_VOXEL_MESH_GENERATOR_MODELER_H_INCLUDED )
#define  KRATOS_VOXEL_MESH_GENERATOR_MODELER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part_io.h"
#include "modeler/modeler.h"
#include "internals/cartesian_mesh_colors.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(MESHING_APPLICATION) VoxelMeshGeneratorModeler
    : public Modeler
{
    
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Modeler
    KRATOS_CLASS_POINTER_DEFINITION(VoxelMeshGeneratorModeler);

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

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
    virtual ~VoxelMeshGeneratorModeler() = default;

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

    void SetupGeometryModel() override ;

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

    ///@}

protected:
    // Data that I store in the internal data structure for each node of the cartesian grid
    class CartesianNodalData
    {
        Node<3>::Pointer mpNode = nullptr;

    public:
        CartesianNodalData() : mpNode(nullptr) {}
        Node<3>::Pointer pGetNode() const { return mpNode; }
        void pSetNode(Node<3>::Pointer pNode) { mpNode = pNode; }
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
    };

    ///@name Private members
    ///@{

    std::size_t mStartNodeId;
    std::size_t mStartElementId;
    std::size_t mStartConditionId; 
    ModelPart* mpInputModelPart = nullptr;
    Model* mpModel = nullptr;   
    array_1d<std::vector<double>,3> mKeyPlanes;
    std::map<std::size_t,ModelPart&> mModelPartsColors;
    Internals::CartesianMeshColors mColors;
    CartesianData mMeshingData;

    ///@}
    ///@name Private Operations
    ///@{

    /// Creates the key planes positions in x,y,z directions
    /** It takes the min_point and max_point of the bounding box and 
     *  divides each direction into uniform divisions with length as 
     *  close as possible to the given size.
    */
    void VoxelMeshKeyPlaneGeneratorBySize(Parameters KeyPlaneGeneratorParameters);

    /// This initializes de internal cartesian mesh data structure to be used for coloring
    void PraparingTheInternalDataStructure();

    /// Goes over the coloring array and apply them to the internal data structure
    void ApplyColoring(Parameters ColoringParameters);

    /// Applyes the coloring to the cells which their center is inside the given skin model part
    void ApplyColorToCellsWithInsideCenter(ModelPart const& TheSkinModelPart, Parameters TheParameters, double OutsideColor);

    /// Applyes the coloring to the cells which are in touch with skin model part
    void ApplyColorToCellsInTouch(ModelPart const& TheSkinModelPart, Parameters TheParameters, double OutsideColor);

    void ApplyColorToCellsInTouchWithGeometry(Element::GeometryType& TheGeometry, double InsideColor, double OutsideColor);

	void ApplyColorToCellsFaces(ModelPart const& TheSkinModelPart, Parameters parameters, double OutsideColor);

	void ApplyColorToOuterFacesOfCellsWithColor(Parameters parameters, double OutsideColor);

    void ApplyColorIfOuterFace(double InterfaceColor, double CellColor, array_1d<std::size_t, 3> CellIndices);

    void ApplyColorToCellsWithCenterInSphereArroundNodes(Parameters parameters, double OutsideColor);

    void ApplyColorToConnectedCellsInTouch(ModelPart const& TheSkinModelPart, Parameters parameters);

    void ApplyColorToConnectedCellsInTouchWithGeometry(Element::GeometryType& TheGeometry, double InsideColor, double CellColor);

    void ColorConnectedCellsToThisCell(std::size_t I, std::size_t J, std::size_t K, double InsideColor, double CellColor);
   

    virtual void GenerateElementsWithCellColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters);

    void GenerateTetrahedralElementsWithCellColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters);

    void GenerateConditionsWithFaceColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters);

    void GenerateTriangularConditionsWithFaceColor(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters);

    void CheckGeometryType(const GeometryData::KratosGeometryType &geometry_type);

    void SetStartIds(ModelPart& rTheVolumeModelPart);

    ModelPart& CreateAndGetModelPart(std::string const& FullName);

    void GenerateEntities(ModelPart& rTheVolumeModelPart, Parameters EntityGeneratorParameters);

    Node<3>::Pointer GenerateOrRetriveNode(ModelPart& rTheVolumeModelPart, ModelPart::NodesContainerType& rThisNodes, std::size_t I, std::size_t J, std::size_t K);

    bool NodeBelongsToInsideElementsAndConditions(std::size_t I, std::size_t J, std::size_t K);

    static void CreateTetrahedraInCell(ModelPart& rTheVolumeModelPart, Element::NodesArrayType & rCellNodes, const std::size_t StartId, Properties::Pointer & pProperties );

 
    ///@}
    ///@name Serializer
    ///@{

    friend class Serializer;

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

#endif //KRATOS_VOXEL_MESH_GENERATOR_MODELER_H_INCLUDED defined
