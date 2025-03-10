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

#if !defined(KRATOS_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS)
#define KRATOS_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/modified_shape_functions.h"

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

class KRATOS_API(KRATOS_CORE) Tetrahedra3D4ModifiedShapeFunctions : public ModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of Tetrahedra3D4ModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(Tetrahedra3D4ModifiedShapeFunctions);

    // General type definitions
    typedef ModifiedShapeFunctions                             BaseType;
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
    Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType rpInputGeometry, const Vector& rNodalDistances);

    //
    Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType rpInputGeometry, const Vector& rNodalDistances, const Vector& rStructureNodes);

    /// Destructor
    ~Tetrahedra3D4ModifiedShapeFunctions() override;

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
    const DivideGeometry<Node>::Pointer pGetSplittingUtil() const override;

    /**
    * Returns the shape function values in the positive split element side for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rPositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
        Vector &rPositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
        Vector &rNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod);

    ///@}

    /**
    * Returns the shape function values in the positive split element interface side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfacePositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
        Vector &rInterfacePositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the shape function values in the negative split element interface side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfaceNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
        Vector &rInterfaceNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the shape function values in the negative split element side for a given quadrature on the contact line.
    * @param ContactLineIndices: indices associated with the outer faces that can be considered as contact lines
    * @return rContactLineNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rContactLineNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rContactLineNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
        std::vector<unsigned int>& ContactLineIndices,
        std::vector<Matrix> &rContactLineNegativeSideShapeFunctionsValues,
        std::vector<ShapeFunctionsGradientsType> &rContactLineNegativeSideShapeFunctionsGradientsValues,
        std::vector<Vector> &rContactLineNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) override;

    /**
    * Given a face id, returns the shape function values in the positive split element exterior face side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rPositiveExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
        Vector &rPositiveExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Given a face id, returns the shape function values in the negative split element exterior face side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rNegativeExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
        Vector &rNegativeExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveSideInterfaceAreaNormals: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveSideInterfaceAreaNormals(
        AreaNormalsContainerType &rPositiveSideInterfaceAreaNormals,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeSideInterfaceAreaNormals: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeSideInterfaceAreaNormals(
        AreaNormalsContainerType &rNegativeSideInterfaceAreaNormals,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveExteriorFaceAreaNormals(
        AreaNormalsContainerType &rPositiveExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeExteriorFaceAreaNormals(
        AreaNormalsContainerType &rNegativeExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the positive side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the positive side edge intersection shape function values. For non-split edges,
    * the corresponding row is plenty of zeros.
    */
    void ComputeShapeFunctionsOnPositiveEdgeIntersections(
        Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues);

    /**
    * Returns the negative side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the negative side edge intersection shape function values. For non-split edges,
    * the corresponding row is plenty of zeros.
    */
    void ComputeShapeFunctionsOnNegativeEdgeIntersections(
        Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues);

    /**
    * Returns the negative side tangential vector (in the direction of contact line).
    * @param FaceIndices: vector giving the indices of the contact faces from DivideGeometry::mContactFace
    * @return rNegativeSideContactLineVector: single vector showing the contact line. 
    */
    void ComputeNegativeSideContactLineVector(
        std::vector<unsigned int>& FaceIndices,
        std::vector<Vector> &rNegativeSideContactLineVector) override;
    /* bool ComputeNegativeSideContactLineVector(
        Vector &rNegativeSideContactLineVector) override; */

    /**
    * Returns true if the element is split and false otherwise.
    */
    bool IsSplit();

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

    void SetCondensationMatrix(Matrix& rIntPointCondMatrix) override;

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

    DivideTetrahedra3D4<Node>::Pointer mpTetrahedraSplitter;

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
    Tetrahedra3D4ModifiedShapeFunctions& operator=(Tetrahedra3D4ModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    Tetrahedra3D4ModifiedShapeFunctions(Tetrahedra3D4ModifiedShapeFunctions const& rOther) :
        ModifiedShapeFunctions(rOther.GetInputGeometry(), rOther.GetNodalDistances()),
        mpTetrahedraSplitter(new DivideTetrahedra3D4<Node>(*rOther.GetInputGeometry(), rOther.GetNodalDistances())) {

        // Perform the element splitting
        mpTetrahedraSplitter->GenerateDivision();
        mpTetrahedraSplitter->GenerateIntersectionsSkin();
    };

    ///@}

};// class Tetrahedra3D4ModifiedShapeFunctions

}
#endif /* KRATOS_TETRAHEDRA_3D_4_MODIFIED_SHAPE_FUNCTIONS defined */
