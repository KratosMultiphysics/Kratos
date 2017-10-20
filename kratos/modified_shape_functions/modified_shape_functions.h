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
    typedef Geometry < Node<3> >                                                  GeometryType;
    typedef GeometryData::IntegrationMethod                              IntegrationMethodType;
    typedef typename GeometryData::ShapeFunctionsGradientsType     ShapeFunctionsGradientsType;
    
    typedef IndexedPoint                                                      IndexedPointType;
    typedef typename IndexedPoint::Pointer                             IndexedPointPointerType;
    typedef Geometry < IndexedPoint >                                 IndexedPointGeometryType;
    typedef Geometry < IndexedPoint >::Pointer                 IndexedPointGeometryPointerType;
    typedef PointerVectorSet<IndexedPointType, IndexedObject>       IndexedPointsContainerType;
    typedef IndexedPointsContainerType::iterator                     IndexedPointsIteratorType;

    typedef IntegrationPoint<3>                                                                          IntegrationPointType;
    typedef std::vector<IntegrationPointType>                                                      IntegrationPointsArrayType;
    typedef boost::array<IntegrationPointsArrayType, GeometryData::NumberOfIntegrationMethods> IntegrationPointsContainerType;

    int mSplitEdgesNumber;  // Number of split edges
    int mDivisionsNumber;   // Number of generated subdivisions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ModifiedShapeFunctions(GeometryType& rInputGeometry, Vector& rNodalDistances);

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

    GeometryType GetInputGeometry() const;

    /**
    * Returns the shape function values in either the positive or negative element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @return rShapeFunctionsGradientsValues: std::vector containing the shape functions gradients values.
    * @return rWeightsValues: Vector containing the Gauss pts. weights (already multiplied by the Jacobian).
    * @param rSubdivisionGeom: std::vector of subdivisions point based geometries where the values are to be computed.
    * @param IntegrationMethod: Desired integration quadrature.
    */
    virtual void GetShapeFunctionsAndGradientsValues(Matrix& rShapeFunctionsValues,
                                                     std::vector<Matrix>& rShapeFunctionsGradientsValues,
                                                     Vector& rWeightsValues,
                                                     const std::vector < IndexedPointGeometryPointerType >& rSubdivisionsVector,
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

    GeometryType& mrInputGeometry;
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
        : mrInputGeometry(rOther.mrInputGeometry) , mrNodalDistances(rOther.mrNodalDistances) {};

    ///@}

};// class ModifiedShapeFunctions

}
#endif /* KRATOS_MODIFIED_SHAPE_FUNCTIONS defined */
