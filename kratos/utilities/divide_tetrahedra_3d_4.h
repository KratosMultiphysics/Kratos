//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_DIVIDE_TETRAHEDRA_3D_4_UTILS)
#define KRATOS_DIVIDE_TETRAHEDRA_3D_4_UTILS

// System includes
#include <bitset>

// External includes

// Project includes

#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "utilities/divide_geometry.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{
template<class TPointType>
class KRATOS_API(KRATOS_CORE) DivideTetrahedra3D4 : public DivideGeometry<TPointType>
{
public:
    ///@name Type Definitions
    ///@{

    typedef DivideGeometry<TPointType>                                  BaseType;
    typedef typename BaseType::GeometryType                                      GeometryType;
    typedef typename BaseType::IndexedPointType                                  IndexedPointType;
    typedef typename BaseType::IndexedPointPointerType                           IndexedPointPointerType;
    typedef typename BaseType::IndexedPointGeometryType                          IndexedPointGeometryType;
    typedef typename BaseType::IndexedPointGeometryPointerType                   IndexedPointGeometryPointerType;
    typedef typename BaseType::IndexedPointGeometryType::GeometriesArrayType     IndexedGeometriesArrayType;
    typedef Triangle3D3 < IndexedPointType >                            IndexedPointTriangleType;
    typedef Tetrahedra3D4 < IndexedPointType >                          IndexedPointTetrahedraType;

    typedef Line3D2 < IndexedPointType >                                IndexedPointLineType; // Needed to construct the contact line

    /// Pointer definition of DivideTetrahedra3D4
    KRATOS_CLASS_POINTER_DEFINITION(DivideTetrahedra3D4);

    const std::vector<int> mEdgeNodeI = {0, 0, 0, 1, 1, 2};
    const std::vector<int> mEdgeNodeJ = {1, 2, 3, 2, 3, 3};
    std::vector<int> mSplitEdges = {0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1};

    // std::vector < unsigned int > mContactInterface;     // Zero or One, gives the interface that contacts the solid
    // std::vector < unsigned int > mContactEdge;          // Zero, One, or Two, gives the contact edge of the contact interface 
    // std::vector < IndexedPointGeometryPointerType > mContactLine; // Object to store the contact line(s) (intersection of the interface with solid).
    // //std::vector < unsigned int > mContactLineNodeIds;   // Object to store the contact line(s)' pair of node Ids in order (e.g. 1,2 , 1,3)
    // std::vector < unsigned int > mContactFace;          // Object to store the face (local) number corresponding to mContactLine.
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DivideTetrahedra3D4(const GeometryType& rInputGeometry, const Vector& rNodalDistances);

    //
    // DivideTetrahedra3D4(const GeometryType& rInputGeometry, const Vector& rNodalDistances, const std::vector<int> rStructureNodes);
    DivideTetrahedra3D4(const GeometryType& rInputGeometry, const Vector& rNodalDistances, const Vector& rStructureNodes);

    /// Destructor
    ~DivideTetrahedra3D4();

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Returns the member vector containing the edges node I ids.
     */
    const std::vector<int>& GetEdgeIdsI() const override;

    /**
     * Returns the member vector containing the edges node J ids.
     */
    const std::vector<int>& GetEdgeIdsJ() const override;

    /**
     * Returns the member vector containing the split edges member vector.
     */
    std::vector<int>& GetSplitEdges() override;

    /**
     * Divides the input geometry according to the provided distance data.
     */
    void GenerateDivision() override;

    /**
     * Generates a list containing the intersection interface geometries for either the positive or the negative element subdivisions.
     */
    void GenerateIntersectionsSkin() override;

    /**
     * Generates a list containing the exterior (boundary) faces geometries for either the positive or the negative element subdivisions.
     * @param rExteriorFacesVector Vector containing the generated exterior subfaces geometries
     * @param rExteriorFacesParentSubdivisionsIdsVector Vector containing the ids of the parent subdivision of each subface
     * @param rSubdivisionsContainer Positive or negative parent geometry subdivisions container
     */
    void GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer) override;

    /**
     * Given a father face id, generates a list containing the exterior (boundary)
     * faces geometries belonging to either the positive or negative side of that that father face.
     * @param rExteriorFacesVector Vector containing the generated exterior subfaces geometries
     * @param rExteriorFacesParentSubdivisionsIdsVector Vector containing the ids of the parent subdivision of each subface
     * @param rSubdivisionsContainer Positive or negative parent geometry subdivisions container
     * @param FatherFaceId Father face in where the positive exterior faces are to be obtained
     */
    void GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer,
        const unsigned int FatherFaceId) override;

    /**
     * Given two edge numbers, the common face (if available) is given
     * @param edgeIdI Id of the first edge
     * @param edgeIdJ Id of the second edge
     */
    int FindCommonFace(const int edgeIdI, const int edgeIdJ);

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const Vector& m_structure_node_id;
    
    std::bitset<4> mNodeIsCut{0x0};

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    bool NodeIsInterface(int NodeKey) const;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DivideTetrahedra3D4& operator=(DivideTetrahedra3D4 const& rOther);

    /// Copy constructor.
    //DivideTetrahedra3D4(DivideTetrahedra3D4 const& rOther)
        //: DivideGeometry<TPointType>(rOther.GetInputGeometry(), rOther.GetNodalDistances()) {};
    DivideTetrahedra3D4(DivideTetrahedra3D4 const& rOther)
        : DivideGeometry<TPointType>(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
          m_structure_node_id(rOther.m_structure_node_id) {};

    ///@}

};// class DivideTetrahedra3D4

}
#endif /* KRATOS_DIVIDE_TETRAHEDRA_3D_4_UTILS defined */
