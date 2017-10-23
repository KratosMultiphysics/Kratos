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
    typedef GeometryData::IntegrationMethod                                                                  IntegrationMethodType;
    typedef typename GeometryData::ShapeFunctionsGradientsType                                         ShapeFunctionsGradientsType;
    
    typedef typename DivideGeometry::IndexedPointGeometryType                                             IndexedPointGeometryType;
    typedef typename DivideGeometry::IndexedPointGeometryPointerType                               IndexedPointGeometryPointerType;

    typedef IntegrationPoint<3>                                                                               IntegrationPointType;
    typedef std::vector<IntegrationPointType>                                                           IntegrationPointsArrayType;
    typedef boost::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods>      IntegrationPointsContainerType;

    int mSplitEdgesNumber;  // Number of split edges
    int mDivisionsNumber;   // Number of generated subdivisions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ModifiedShapeFunctions(GeometryPointerType rInputGeometry, Vector& rNodalDistances);

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

    Vector GetNodalDistances() const ;

    GeometryPointerType GetInputGeometry() const;

    /**
    * Returns the shape function values in both the positive or negative split element sides for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void GetShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                     Matrix &rNegativeSideShapeFunctionsValues,
                                                     std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                     std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                     Vector &rPositiveSideWeightsValues,
                                                     Vector &rNegativeSideWeightsValues,
                                                     const IntegrationMethodType IntegrationMethod);

    /**
    * Returns the shape function values in the positive split element side for a given quadrature.
    * @return rPositiveSideShapeFunctionValues: Matrix containing the positive side computed shape function values.
    * @return rPositiveSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the positive side.
    * @return rPositiveSideWeightsValues: Vector containing the Gauss pts. positive side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void GetPositiveSideShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                                 std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                                 Vector &rPositiveSideWeightsValues,
                                                                 const IntegrationMethodType IntegrationMethod);
                                                                 
    /**
    * Returns the shape function values in the negative split element side for a given quadrature.
    * @return rNegativeSideShapeFunctionValues: Matrix containing the negative side computed shape function values.
    * @return rNegativeSideShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values on the negative side.
    * @return rNegativeSideWeightsValues: Vector containing the Gauss pts. negative side weights (already multiplied by the Jacobian).
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void GetNegativeSideShapeFunctionsAndGradientsValues(Matrix &rNegativeSideShapeFunctionsValues,
                                                                 std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                                 Vector &rNegativeSideWeightsValues,
                                                                 const IntegrationMethodType IntegrationMethod);

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
    */
    bool IsSplit();

    /**
    * Returns true if the element is split and false otherwise.
    * @return rIntPointCondMatrix: Reference to the intersection points condensation matrix.
    * @param rEdgeNodeI: Integers array containing the nodes "I" that conform the edges.
    * @param rEdgeNodeJ: Integers array containing the nodes "J" that conform the edges.
    * @param rSplitEdges: Integers array containing the original nodes ids and the intersected edges nodes ones.
    * @param splitEdgesNumber: Number of splitted edges.
    */
    void SetIntersectionPointsCondensationMatrix(Matrix& rIntPointCondMatrix,
                                                 const int rEdgeNodeI[],
                                                 const int rEdgeNodeJ[],
                                                 const int rSplitEdges[],
                                                 const unsigned int splitEdgesNumber);

    /**
    * Returns the shape function values in either the positive or negative element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rSubdivisionGeom: std::vector of subdivisions point based geometries where the values are to be computed.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void ComputeValuesOnOneSide(Matrix &rShapeFunctionsValues,
                                        std::vector<Matrix> &rShapeFunctionsGradientsValues,
                                        Vector &rWeightsValues,
                                        const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
                                        const Matrix &p_matrix,
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

    GeometryPointerType mpInputGeometry;
    Vector& mrNodalDistances;

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
