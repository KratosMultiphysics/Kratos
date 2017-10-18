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

#if !defined(KRATOS_GEOMETRY_SPLITTING_UTILS)
#define KRATOS_GEOMETRY_SPLITTING_UTILS

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/indexed_object.h"
#include "containers/pointer_vector_set.h"

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

class KRATOS_API(KRATOS_CORE) IndexedPoint : public Point<3>, public IndexedObject
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IndexedPoint
    KRATOS_CLASS_POINTER_DEFINITION(IndexedPoint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty constructor
    IndexedPoint();

    /// Auxiliar constructor
    IndexedPoint(const unsigned int Id);

    /// Default constructor
    IndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id);

    /// Destructor
    ~IndexedPoint();

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
    ///@name Member variables
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

    friend class Serializer;

    void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point<3>);
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    };

    void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point<3>);
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
    };

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

    ///@}
};

class KRATOS_API(KRATOS_CORE) GeometrySplittingUtils
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometrySplittingUtils
    KRATOS_CLASS_POINTER_DEFINITION(GeometrySplittingUtils);

    // General type definitions
    typedef Geometry < Node<3> >                                                  GeometryType;
    typedef GeometryData::IntegrationMethod                              IntegrationMethodType;
    
    typedef IndexedPoint                                                      IndexedPointType;
    typedef typename IndexedPoint::Pointer                             IndexedPointPointerType;
    typedef Geometry < IndexedPoint >                                 IndexedPointGeometryType;
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
    GeometrySplittingUtils(GeometryType& rInputGeometry, Vector& rNodalDistances);

    /// Destructor
    ~GeometrySplittingUtils();

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
     * Divides the input geometry according to the provided distance data.
     * @return rAuxPoints: Reference to the pointer vector set containing the original nodes plus the intersection ones.
     * @return rPositiveSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the positive subdivisions.
     * @return rNegativeSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the negative subdivisions.
     */
    virtual bool DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                std::vector < IndexedPointGeometryType >& rPositiveSubdivisions,
                                std::vector < IndexedPointGeometryType >& rNegativeSubdivisions);

    /**
    * Returns the shape function values in any element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @param rSubdivisionGeom: Subdivision point based geometry.
    * @param IntegrationMethod: Integration quadrature.
    */
    virtual void GetShapeFunctionValues(Matrix& rShapeFunctionValues,
                                        Vector& rWeightsValues,
                                        const std::vector < IndexedPointGeometryType >& rSubdivisionsVector,
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
    GeometrySplittingUtils& operator=(GeometrySplittingUtils const& rOther);

    /// Copy constructor.
    GeometrySplittingUtils(GeometrySplittingUtils const& rOther)
        : mrInputGeometry(rOther.mrInputGeometry) , mrNodalDistances(rOther.mrNodalDistances) {};

    ///@}

};// class GeometrySplittingUtils

class KRATOS_API(KRATOS_CORE) TriangleSplittingUtils : public GeometrySplittingUtils
{
public:
    ///@name Type Definitions
    ///@{

    typedef Triangle2D3 < IndexedPointType >                        IndexedPointTriangleType;

    /// Pointer definition of TriangleSplittingUtils
    KRATOS_CLASS_POINTER_DEFINITION(TriangleSplittingUtils);

    int mEdgeNodeI[3] = {0, 1, 2};
    int mEdgeNodeJ[3] = {1, 2, 0};
    int mSplitEdges[6] = {0, 1, 2, -1, -1, -1};

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    TriangleSplittingUtils(GeometryType& rInputGeometry, Vector& rNodalDistances);

    /// Destructor
    ~TriangleSplittingUtils();

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
    bool DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                        std::vector < IndexedPointGeometryType >& rPositiveSubdivisions,
                        std::vector < IndexedPointGeometryType >& rNegativeSubdivisions) override;

    /**
    * Returns the shape function values in any element subdivision for a given quadrature.
    * @return rShapeFunctionValues: Matrix containing the computed shape function values.
    * @param rSubdivisionGeom: Subdivision point based geometry.
    * @param IntegrationMethod: Integration quadrature.
    */
    void GetShapeFunctionValues(Matrix& rShapeFunctionValues,
                                Vector& rWeightsValues,
                                const std::vector < IndexedPointGeometryType >& rSubdivisionsVector,
                                const IntegrationMethodType IntegrationMethod) override;

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
    // TriangleSplittingUtils& operator=(TriangleSplittingUtils const& rOther);
    //
    // /// Copy constructor.
    // TriangleSplittingUtils(TriangleSplittingUtils const& rOther)
    //     : GeometrySplittingUtils(rOther.mrInputGeometry, rOther.mrNodalDistances) {};

    ///@}

};// class TriangleSplittingUtils

}
#endif /* KRATOS_GEOMETRY_SPLITTING_UTILS defined */
