//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//                   Kratos default license: kratos/license.txt
//
//   Main authors:   Pablo Becker
//

#if !defined(KRATOS_DIVIDE_TRIANGLE_3D_3_UTILS)
#define KRATOS_DIVIDE_TRIANGLE_3D_3_UTILS

// System includes

// External includes

// Project includes

#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/divide_geometry.h"
#include "utilities/divide_triangle_2d_3.h"

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
class KRATOS_API(KRATOS_CORE) DivideTriangle3D3 : public DivideTriangle2D3<TPointType>
{
public:
    ///@name Type Definitions
    ///@{
    typedef DivideGeometry<TPointType>                                  BaseType;
    typedef typename BaseType::GeometryType                             GeometryType;
    typedef typename BaseType::IndexedPointType                         IndexedPointType;
    typedef Line3D2 < IndexedPointType >                                IndexedPointLineType;
    typedef Triangle3D3 < IndexedPointType >                            IndexedPointTriangleType;
    typedef typename BaseType::IndexedPointGeometryPointerType          IndexedPointGeometryPointerType;
    /// Pointer definition of DivideTriangle2D3
    KRATOS_CLASS_POINTER_DEFINITION(DivideTriangle3D3);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DivideTriangle3D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~DivideTriangle3D3();

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{
    IndexedPointGeometryPointerType GenerateAuxiliaryPartitionTriangle(
        const int I0,
        const int I1,
        const int I2) override;

    IndexedPointGeometryPointerType GenerateIntersectionLine(
        const int I0,
        const int I1) override;
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
    DivideTriangle3D3& operator=(DivideTriangle3D3 const& rOther);

    /// Copy constructor.
    DivideTriangle3D3(DivideTriangle3D3 const& rOther)
        : DivideTriangle2D3<TPointType>(rOther.GetInputGeometry(), rOther.GetNodalDistances()) {};

    ///@}

};// class DivideTriangle3D3

}
#endif /* KRATOS_DIVIDE_TRIANGLE_3D_3_UTILS defined */
