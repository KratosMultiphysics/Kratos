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

#if !defined(KRATOS_DIVIDE_TRIANGLE_2D_3_UTILS)
#define KRATOS_DIVIDE_TRIANGLE_2D_3_UTILS

// System includes
#include <bitset>

// External includes

// Project includes

#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
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
class KRATOS_API(KRATOS_CORE) DivideTriangle2D3 : public DivideGeometry<TPointType>
{
public:
    ///@name Type Definitions
    ///@{

    typedef DivideGeometry<TPointType>                                           BaseType;
    typedef typename BaseType::GeometryType                                      GeometryType;
    typedef typename BaseType::IndexedPointType                                  IndexedPointType;
    typedef typename BaseType::IndexedPointPointerType                           IndexedPointPointerType;
    typedef typename BaseType::IndexedPointGeometryType                          IndexedPointGeometryType;
    typedef typename BaseType::IndexedPointGeometryType::GeometriesArrayType     IndexedGeometriesArrayType;
    typedef typename BaseType::IndexedPointGeometryPointerType                   IndexedPointGeometryPointerType;
    typedef Line2D2 < IndexedPointType >                                         IndexedPointLineType;
    typedef Triangle2D3 < IndexedPointType >                                     IndexedPointTriangleType;

    /// Pointer definition of DivideTriangle2D3
    KRATOS_CLASS_POINTER_DEFINITION(DivideTriangle2D3);

    const std::vector<int> mEdgeNodeI = {0, 1, 2};
    const std::vector<int> mEdgeNodeJ = {1, 2, 0};
    std::vector<int> mSplitEdges = {0, 1, 2, -1, -1, -1};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DivideTriangle2D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~DivideTriangle2D3();

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

    ///@}

    virtual IndexedPointGeometryPointerType GenerateAuxiliaryPartitionTriangle(
        const int I0,
        const int I1,
        const int I2);

    virtual IndexedPointGeometryPointerType GenerateIntersectionLine(
        const int I0,
        const int I1);

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::bitset<3> mNodeIsCut{0x0}; // If the cut passes through a node, store this information here

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
    DivideTriangle2D3& operator=(DivideTriangle2D3 const& rOther);

    /// Copy constructor.
    DivideTriangle2D3(DivideTriangle2D3 const& rOther)
        : DivideGeometry<TPointType>(rOther.GetInputGeometry(), rOther.GetNodalDistances()) {};

    ///@}

};// class DivideTriangle2D3

}
#endif /* KRATOS_DIVIDE_TRIANGLE_2D_3_UTILS defined */
