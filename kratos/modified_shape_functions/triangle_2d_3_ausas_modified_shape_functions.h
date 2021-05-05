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

#if !defined(KRATOS_TRIANGLE_2D_3_AUSAS_MODIFIED_SHAPE_FUNCTIONS)
#define KRATOS_TRIANGLE_2D_3_AUSAS_MODIFIED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/ausas_modified_shape_functions.h"

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

class KRATOS_API(KRATOS_CORE) Triangle2D3AusasModifiedShapeFunctions : public AusasModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Triangle2D3AusasModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(Triangle2D3AusasModifiedShapeFunctions);

    // General type definitions
    typedef AusasModifiedShapeFunctions                        BaseType;
    typedef BaseType::GeometryType                             GeometryType;
    typedef BaseType::GeometryPointerType                      GeometryPointerType;
    typedef BaseType::IntegrationMethodType                    IntegrationMethodType;
    typedef BaseType::ShapeFunctionsGradientsType              ShapeFunctionsGradientsType;

    typedef BaseType::IndexedPointGeometryType                 IndexedPointGeometryType;
    typedef BaseType::IndexedPointGeometryPointerType          IndexedPointGeometryPointerType;

    typedef BaseType::IntegrationPointType                     IntegrationPointType;
    typedef BaseType::IntegrationPointsArrayType               IntegrationPointsArrayType;
    typedef BaseType::IntegrationPointsContainerType           IntegrationPointsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    Triangle2D3AusasModifiedShapeFunctions(const GeometryPointerType rpInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~Triangle2D3AusasModifiedShapeFunctions();

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
    * Returns the member pointer to the splitting utility.
    */
    const DivideGeometry::Pointer pGetSplittingUtil() const override;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveSideInterfaceAreaNormal: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveSideInterfaceAreaNormals(
        std::vector<Vector> &rPositiveSideInterfaceAreaNormal,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeSideInterfaceAreaNormal: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeSideInterfaceAreaNormals(
        std::vector<Vector> &rNegativeSideInterfaceAreaNormal,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveExteriorFaceAreaNormals(
        std::vector<Vector> &rPositiveExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeExteriorFaceAreaNormals(
        std::vector<Vector> &rNegativeExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Returns the positive side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the positive side edge intersection shape function values. For non-split edges,
    * the corresponding row is plenty of zeros.
    */
    void ComputeShapeFunctionsOnPositiveEdgeIntersections(
        Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues) override;

    /**
    * Returns the negative side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the negative side edge intersection shape function values. For non-split edges,
    * the corresponding row is plenty of zeros.
    */
    void ComputeShapeFunctionsOnNegativeEdgeIntersections(
        Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues) override;

    /**
    * Returns true if the element is split and false otherwise.
    */
    bool IsSplit() const override;

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix) override;

    void SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix) override;

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    DivideTriangle2D3::Pointer mpTriangleSplitter;

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
    Triangle2D3AusasModifiedShapeFunctions& operator=(Triangle2D3AusasModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    Triangle2D3AusasModifiedShapeFunctions(Triangle2D3AusasModifiedShapeFunctions const& rOther) :
        AusasModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
        mpTriangleSplitter(new DivideTriangle2D3(*rOther.GetInputGeometry(), rOther.GetNodalDistances())) {

        // Perform the element splitting
        mpTriangleSplitter->GenerateDivision();
        mpTriangleSplitter->GenerateIntersectionsSkin();
    };

    ///@}

};// class Triangle2D3AusasModifiedShapeFunctions

}
#endif /* KRATOS_TRIANGLE_2D_3_AUSAS_MODIFIED_SHAPE_FUNCTIONS defined */
