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

#if !defined(KRATOS_DIVIDE_GEOMETRY_UTILS)
#define KRATOS_DIVIDE_GEOMETRY_UTILS

// System includes

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

class KRATOS_API(KRATOS_CORE) DivideTriangle2D3 : public DivideGeometry
{
public:
    ///@name Type Definitions
    ///@{

    typedef DivideGeometry                                                    BaseType;
    typedef typename BaseType::GeometryType                               GeometryType;
    typedef typename BaseType::IndexedPointType                       IndexedPointType;
    typedef Line2D2 < IndexedPointType >                          IndexedPointLineType;
    typedef Triangle2D3 < IndexedPointType >                  IndexedPointTriangleType;

    /// Pointer definition of DivideTriangle2D3
    KRATOS_CLASS_POINTER_DEFINITION(DivideTriangle2D3);

    int mEdgeNodeI[3] = {0, 1, 2};
    int mEdgeNodeJ[3] = {1, 2, 0};
    int mSplitEdges[6] = {0, 1, 2, -1, -1, -1};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DivideTriangle2D3(GeometryType& rInputGeometry, Vector& rNodalDistances);

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
     * Divides the input geometry according to the provided distance data.
     * @return rAuxPoints: Reference to the pointer vector set containing the original nodes plus the intersection ones.
     * @return rPositiveSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the positive subdivisions.
     * @return rNegativeSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the negative subdivisions.
     */
    bool GenerateDivision(IndexedPointsContainerType& rAuxPoints,
                          std::vector < IndexedPointGeometryPointerType >& rPositiveSubdivisions,
                          std::vector < IndexedPointGeometryPointerType >& rNegativeSubdivisions) override;

    /**
     * Generates a list containing the intersection interface geometries for either the positive or the negative element subdivisions.
     * @return rInterfacesVector: Reference to a std::vector containing the generated interface geometries.
     * @param rSubdivisionsVector: std::vector of subdivisions point based geometries to its intersection geometries.
     */
    void GenerateIntersectionsSkin(std::vector < IndexedPointGeometryPointerType >& rInterfacesVector,
                                   IndexedPointsContainerType& rAuxPoints,
                                   const std::vector < IndexedPointGeometryPointerType >& rSubdivisionsVector) override;

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

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

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // DivideTriangle2D3& operator=(DivideTriangle2D3 const& rOther);
    //
    // /// Copy constructor.
    // DivideTriangle2D3(DivideTriangle2D3 const& rOther)
    //     : DivideGeometry(rOther.mrInputGeometry, rOther.mrNodalDistances) {};

    ///@}

};// class DivideTriangle2D3

}
#endif /* KRATOS_DIVIDE_GEOMETRY_UTILS defined */
