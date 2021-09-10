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

#if !defined(KRATOS_MODIFIED_SHAPE_FUNCTIONS)
#define KRATOS_MODIFIED_SHAPE_FUNCTIONS

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
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

class KRATOS_API(KRATOS_CORE) ModifiedShapeFunctions
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModifiedShapeFunctions
    KRATOS_CLASS_POINTER_DEFINITION(ModifiedShapeFunctions);

    // General type definitions
    typedef Geometry < Node<3> >                                                                    GeometryType;
    typedef GeometryType::Pointer                                                                   GeometryPointerType;
    typedef GeometryType::CoordinatesArrayType                                                      CoordinatesArrayType;
    typedef GeometryData::IntegrationMethod                                                         IntegrationMethodType;
    typedef GeometryData::ShapeFunctionsGradientsType                                               ShapeFunctionsGradientsType;
    typedef std::vector<array_1d<double,3>> AreaNormalsContainerType;

    typedef DivideGeometry::IndexedPointGeometryType                                                IndexedPointGeometryType;
    typedef DivideGeometry::IndexedPointGeometryPointerType                                         IndexedPointGeometryPointerType;

    typedef IntegrationPoint<3>                                                                     IntegrationPointType;
    typedef std::vector<IntegrationPointType>                                                       IntegrationPointsArrayType;
    typedef std::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods>      IntegrationPointsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    virtual ~ModifiedShapeFunctions();

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
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    * Returns the member pointer to the splitting utility.
    */
    virtual const DivideGeometry::Pointer pGetSplittingUtil() const;

    /**
    * Returns a the member pointer to the input geometry.
    */
    const GeometryPointerType GetInputGeometry() const;

    /**
    * Returns a reference to the nodal distances vector member variable.
    */
    const Vector& GetNodalDistances() const;

    /**
     * @brief Compute the positive domain size
     * For the given input geometry, this method calculates the positive side volume (or area in 2D)
     * @return double The positive side domain size
     */
    double ComputePositiveSideDomainSize() const;

    /**
     * @brief Compute the negative domain size
     * For the given input geometry, this method calculates the negative side volume (or area in 2D)
     * @return double The negative side domain size
     */
    double ComputeNegativeSideDomainSize() const;

    void ComputePositiveSideShapeFunctionsAndWeights(
        Matrix &rPositiveSideShapeFunctionsValues,
        Vector &rPositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod);

    void ComputeNegativeSideShapeFunctionsAndWeights(
        Matrix &rNegativeSideShapeFunctionsValues,
        Vector &rNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod);

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
    * Returns the shape function values in the negative split element side for a given quadrature.
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
    * @return rPositiveSideInterfaceAreaNormal: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveSideInterfaceAreaNormals(
        AreaNormalsContainerType& rPositiveSideInterfaceAreaNormal,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeSideInterfaceAreaNormal: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeSideInterfaceAreaNormals(
        AreaNormalsContainerType& rNegativeSideInterfaceAreaNormal,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputePositiveExteriorFaceAreaNormals(
        AreaNormalsContainerType& rPositiveExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeExteriorFaceAreaNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeNegativeExteriorFaceAreaNormals(
        AreaNormalsContainerType& rNegativeExteriorFaceAreaNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the positive side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the positive side edge intersection shape function values.
    */
    void ComputeShapeFunctionsOnPositiveEdgeIntersections(Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues);

    /**
    * Returns the negative side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the negative side edge intersection shape function values.
    */
    void ComputeShapeFunctionsOnNegativeEdgeIntersections(Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues);

    /**
    * Returns true if the element is split and false otherwise.
    */
    bool IsSplit() const;

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

    /**
    * Returns the intersection points condensation matrix.
    * This matrix is used to extrapolate the subdivisions shape funtion values to the
    * original geometry ones. It has size (nnodes+nedges)x(nnodes).
    * @return rIntPointCondMatrix: Reference to the intersection points condensation matrix.
    */
    virtual void SetCondensationMatrix(Matrix& rIntPointCondMatrix)
    {
        KRATOS_ERROR << "Calling base class \'SetCondensationMatrix\'. Call the derived class one instead." << std::endl;
    }

    /**
     * @brief Set the Positive Side Condensation Matrix object
     * This function sets the positive side condensation matrix
     * Note that for the Ausas FE space a different condensation matrix is required for the positive and negative sides
     * @param rPosSideCondMatrix Positive side condensation matrix
     */
    virtual void SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
    {
        KRATOS_ERROR << "Calling base class \'SetPositiveSideCondensationMatrix\'. Call the derived class one instead." << std::endl;
    }

    /**
     * @brief Set the Negative Side Condensation Matrix object
     * This function sets the negative side condensation matrix
     * Note that for the Ausas FE space a different condensation matrix is required for the positive and negative sides
     * @param rNegSideCondMatrix Negative side condensation matrix
     */
    virtual void SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
    {
        KRATOS_ERROR << "Calling base class \'SetNegativeSideCondensationMatrix\'. Call the derived class one instead." << std::endl;
    }

    /**
    * Returns the intersection points condensation matrix.
    * This matrix is used to extrapolate the subdivisions shape funtion values to the
    * original geometry ones. It has size (nnodes+nedges)x(nnodes).
    * @return rIntPointCondMatrix: Reference to the intersection points condensation matrix.
    * @param rEdgeNodeI Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges Integers array containing the original nodes ids and the intersected edges nodes ones.
    */
    void SetCondensationMatrix(
        Matrix& rIntPointCondMatrix,
        const std::vector<int>& rEdgeNodeI,
        const std::vector<int>& rEdgeNodeJ,
        const std::vector<int>& rSplitEdges);

    /**
    * Returns the shape function values in either the positive or negative element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rPmatrix Reference to the condensation matrix.
    * @param rSubdivisionGeom std::vector of subdivisions point based geometries where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeValuesOnOneSide(
        Matrix &rShapeFunctionsValues,
        ShapeFunctionsGradientsType &rShapeFunctionsGradientsValues,
        Vector &rWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
        const Matrix &rPmatrix,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the shape function values in either the positive or negative element interfaces for a given quadrature.
    * @return rInterfaceShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rInterfaceShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rInterfaceWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rInterfacesVector std::vector of intersection point based geometries where the values are to be computed.
    * @param rParentGeometriesVector std::vector of subdivisions point based parent geometries.
    * @param rInterfacesParentIdsVector std::vector containing the parent ids of each interface geometry.
    * @param rPmatrix reference to the interface interpolation matrix
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeFaceValuesOnOneSide(
        Matrix &rInterfaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfaceShapeFunctionsGradientsValues,
        Vector &rInterfaceWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const std::vector<IndexedPointGeometryPointerType> &rParentGeometriesVector,
        const std::vector<unsigned int> &rInterfacesParentIdsVector,
        const Matrix &rPmatrix,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the outwards area normal vector values in either the positive or negative element interfaces for a given quadrature.
    * @return rInterfaceAreaNormalValues: std::vector containing the area normal values for the selected quadrature Gauss pts.
    * @param rInterfacesVector std::vector of intersection point based geometries where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeFaceNormalOnOneSide(
        AreaNormalsContainerType& rInterfaceAreaNormalValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Given a condensation matrix, extracts the edge intersection points shape function values.
    * @param rPmatrix Reference to the condensation matrix.
    * @return rEdgeShapeFunctionValues Reference to the matrix containing the shape function values.
    */
    void ComputeEdgeIntersectionValuesOnOneSide(
        const Matrix &rPmatrix,
        Matrix &rEdgeShapeFunctionValues);

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

    private :
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const GeometryPointerType mpInputGeometry;
    const Vector mNodalDistances;

    ///@}
    ///@name Serialization
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Calculate the domain size in one element side
     * For the given subdivisions vector, this method calculates the sum of their domain sizes (area in 2D and volume in 3D)
     * @param rSubdivisionsVector The side of interest subdivisions vectors
     * @return double The domain size of the side of interest
     */
    double ComputeDomainSizeOnOneSide(const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector) const;

    /**
     * @brief Calculate shape function values and weights in one element side
    * Returns the shape function values and weights in either the positive or negative element subdivision for a given quadrature.
    * @return rShapeFunctionsValues Matrix containing the computed shape function values.
    * @return rWeightsValues Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rSubdivisionsVector std::vector of subdivisions point based geometries where the values are to be computed.
    * @param rPmatrix Reference to the condensation matrix.
    * @param IntegrationMethod Desired integration quadrature.
    */
    void ComputeValuesOnOneSide(
        Matrix &rShapeFunctionsValues,
        Vector &rWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
        const Matrix &rPmatrix,
        const IntegrationMethodType IntegrationMethod);

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
    ModifiedShapeFunctions& operator=(ModifiedShapeFunctions const& rOther);

    /// Copy constructor.
    ModifiedShapeFunctions(ModifiedShapeFunctions const& rOther)
        : mpInputGeometry(rOther.mpInputGeometry) , mNodalDistances(rOther.mNodalDistances) {};

    ///@}

};// class ModifiedShapeFunctions

}
#endif /* KRATOS_MODIFIED_SHAPE_FUNCTIONS defined */
