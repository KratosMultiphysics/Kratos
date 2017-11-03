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
    typedef Geometry < Node<3> >                                                                                      GeometryType;
    typedef typename GeometryType::Pointer                                                                     GeometryPointerType;
    typedef typename GeometryType::CoordinatesArrayType                                                       CoordinatesArrayType;
    typedef GeometryData::IntegrationMethod                                                                  IntegrationMethodType;
    typedef typename GeometryData::ShapeFunctionsGradientsType                                         ShapeFunctionsGradientsType;

    typedef typename DivideGeometry::IndexedPointGeometryType                                             IndexedPointGeometryType;
    typedef typename DivideGeometry::IndexedPointGeometryPointerType                               IndexedPointGeometryPointerType;

    typedef IntegrationPoint<3>                                                                               IntegrationPointType;
    typedef std::vector<IntegrationPointType>                                                           IntegrationPointsArrayType;
    typedef boost::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods>      IntegrationPointsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    ~ModifiedShapeFunctions();

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
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rPositiveSideShapeFunctionsValues,
        std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
        Vector &rPositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rNegativeSideShapeFunctionsValues,
        std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
        Vector &rNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the shape function values in the positive split element interface side for a given quadrature.
    * @return rInterfacePositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rInterfacePositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rInterfacePositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfacePositiveSideShapeFunctionsValues,
        std::vector<Matrix> &rInterfacePositiveSideShapeFunctionsGradientsValues,
        Vector &rInterfacePositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rInterfaceNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rInterfaceNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rInterfaceNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfaceNegativeSideShapeFunctionsValues,
        std::vector<Matrix> &rInterfaceNegativeSideShapeFunctionsGradientsValues,
        Vector &rInterfaceNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod) = 0;
        
    /**
    * Returns the positive side outwards unit normal vector values for the Gauss pts. of given quadrature.
    * @return rPositiveSideInterfaceUnitNormal: Outwards unit normal vector list.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputePositiveSideInterfaceUnitNormals(
        std::vector<Vector> &rPositiveSideInterfaceUnitNormal,
        const IntegrationMethodType IntegrationMethod) = 0;

    /**
    * Returns the positive side outwards unit normal vector values for the Gauss pts. of given quadrature.
    * @return rNegativeSideInterfaceUnitNormal: Outwards unit normal vector list.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeNegativeSideInterfaceUnitNormals(
        std::vector<Vector> &rNegativeSideInterfaceUnitNormal,
        const IntegrationMethodType IntegrationMethod) = 0;

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
    * Returns true if the element is split and false otherwise.
    * @return rIntPointCondMatrix: Reference to the intersection points condensation matrix.
    * @param rEdgeNodeI: Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ: Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges: Integers array containing the original nodes ids and the intersected edges nodes ones.
    * @param splitEdgesNumber: Number of splitted edges.
    */
    void SetIntersectionPointsCondensationMatrix(
        Matrix& rIntPointCondMatrix,
        const std::vector<int>& rEdgeNodeI,
        const std::vector<int>& rEdgeNodeJ,
        const std::vector<int>& rSplitEdges);

    /**
    * Returns the shape function values in either the positive or negative element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rPmatrix: Reference to the condensation matrix.
    * @param rSubdivisionGeom: std::vector of subdivisions point based geometries where the values are to be computed.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeValuesOnOneSide(
        Matrix &rShapeFunctionsValues,
        std::vector<Matrix> &rShapeFunctionsGradientsValues,
        Vector &rWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
        const Matrix &rPmatrix,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the shape function values in either the positive or negative element interfaces for a given quadrature.
    * @return rInterfaceShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rInterfaceShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rInterfaceWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rInterfacesVector: std::vector of intersection point based geometries where the values are to be computed.
    * @param rPmatrix: reference to the interpolation matrix
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeInterfaceValuesOnOneSide(
        Matrix &rInterfaceShapeFunctionsValues,
        std::vector<Matrix> &rInterfaceShapeFunctionsGradientsValues,
        Vector &rInterfaceWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the outwards unit normal vector values in either the positive or negative element interfaces for a given quadrature.
    * @return rInterfaceUnitNormalValues: std::vector of subdivisions point based geometries where the values are to be computed.
    * @param rInterfacesVector: std::vector of intersection point based geometries where the values are to be computed.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeInterfaceNormalOnOneSide(
        std::vector<Vector> &rInterfaceUnitNormalValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const IntegrationMethodType IntegrationMethod);

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
    const Vector& mrNodalDistances;

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
        : mpInputGeometry(rOther.mpInputGeometry) , mrNodalDistances(rOther.mrNodalDistances) {};

    ///@}

};// class ModifiedShapeFunctions

}
#endif /* KRATOS_MODIFIED_SHAPE_FUNCTIONS defined */
