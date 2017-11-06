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

class KRATOS_API(KRATOS_CORE) DivideTetrahedra3D4 : public DivideGeometry
{
public:
    ///@name Type Definitions
    ///@{

    typedef DivideGeometry                                              BaseType;
    typedef BaseType::GeometryType                                      GeometryType;
    typedef BaseType::IndexedPointType                                  IndexedPointType;
    typedef BaseType::IndexedPointGeometryType::GeometriesArrayType     GeometriesArrayType;
    typedef Triangle3D3 < IndexedPointType >                            IndexedPointTriangleType;
    typedef Tetrahedra3D4 < IndexedPointType >                          IndexedPointTetrahedraType;

    /// Pointer definition of DivideTetrahedra3D4
    KRATOS_CLASS_POINTER_DEFINITION(DivideTetrahedra3D4);

    const std::vector<int> mEdgeNodeI = {0, 0, 0, 1, 1, 2};
    const std::vector<int> mEdgeNodeJ = {1, 2, 3, 2, 3, 3};
    std::vector<int> mSplitEdges = {0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DivideTetrahedra3D4(const GeometryType& rInputGeometry, const Vector& rNodalDistances);

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
     * Divides the input geometry according to the provided distance data.
     */
    void GenerateDivision() override;

    /**
     * Generates a list containing the intersection interface geometries for either the positive or the negative element subdivisions.
     */
    void GenerateIntersectionsSkin() override;

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
    // DivideTetrahedra3D4& operator=(DivideTetrahedra3D4 const& rOther);
    //
    // /// Copy constructor.
    // DivideTetrahedra3D4(DivideTetrahedra3D4 const& rOther)
    //     : DivideGeometry(rOther.mrInputGeometry, rOther.mrNodalDistances) {};

    ///@}

};// class DivideTetrahedra3D4

}
#endif /* KRATOS_DIVIDE_TETRAHEDRA_3D_4_UTILS defined */
