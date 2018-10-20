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
    const Vector& GetNodalDistances() const ;

    /**
    * Returns the shape function values in the positive split element side for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rPositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
        Vector &rPositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
        Vector &rNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the shape function values in the positive split element interface side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfacePositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
        Vector &rInterfacePositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfaceNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
        Vector &rInterfaceNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Given a face id, returns the shape function values in the positive split element exterior face side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rPositiveExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
        Vector &rPositiveExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Given a face id, returns the shape function values in the negative split element exterior face side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rNegativeExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
        Vector &rNegativeExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveSideInterfaceNormal: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputePositiveSideInterfaceNormals(
        std::vector<Vector> &rPositiveSideInterfaceNormal,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeSideInterfaceNormal: Outwards area normal vector list.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeNegativeSideInterfaceNormals(
        std::vector<Vector> &rNegativeSideInterfaceNormal,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the positive side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveExteriorFaceNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputePositiveExteriorFaceNormals(
        std::vector<Vector> &rPositiveExteriorFaceNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the negative side outwards area normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeExteriorFaceNormal: Outwards area normal vector list.
    * @param FaceId Face local id. in where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeNegativeExteriorFaceNormals(
        std::vector<Vector> &rNegativeExteriorFaceNormal,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the positive side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the positive side edge intersection shape function values.
    */
    virtual void ComputeShapeFunctionsOnPositiveEdgeIntersections(
        Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues) = 0;

    /**
    * Returns the negative side edge intersections shape function values.
    * @return rPositiveEdgeIntersectionsShapeFunctionsValues A matrix, which size is edges x nodes,
    * containing the negative side edge intersection shape function values.
    */
    virtual void ComputeShapeFunctionsOnNegativeEdgeIntersections(
        Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues) = 0;

    /**
    * Returns true if the element is split and false otherwise.
    */
    virtual bool IsSplit() = 0;

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
    virtual void ComputeValuesOnOneSide(
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
    virtual void ComputeFaceValuesOnOneSide(
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
    * @return rInterfaceNormalValues: std::vector containing the area normal values for the selected quadrature Gauss pts.
    * @param rInterfacesVector std::vector of intersection point based geometries where the values are to be computed.
    * @param IntegrationMethod Desired integration quadrature.
    */
    virtual void ComputeFaceNormalOnOneSide(
        std::vector<Vector> &rInterfaceNormalValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Given a condensation matrix, extracts the edge intersection points shape function values.
    * @param rPmatrix Reference to the condensation matrix.
    * @return rEdgeShapeFunctionValues Reference to the matrix containing the shape function values.
    */
    virtual void ComputeEdgeIntersectionValuesOnOneSide(
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
