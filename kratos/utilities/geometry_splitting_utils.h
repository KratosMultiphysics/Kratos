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
    typedef Node<3>                                                                   NodeType;
    typedef Point<3>                                                                 PointType;
    typedef Geometry < NodeType >                                                 GeometryType;
    typedef Geometry < PointType >                                           PointGeometryType;

    typedef IndexedPoint                                                      IndexedPointType;
    typedef boost::shared_ptr<IndexedPoint>                            IndexedPointPointerType;
    typedef PointerVectorSet<IndexedPointType, IndexedObject>       IndexedPointsContainerType;
    typedef IndexedPointsContainerType::iterator                     IndexedPointsIteratorType;

    int mSplitEdgesNumber;  // Number of split edges
    int mDivisionsNumber;   // Number of generated subdivisions

    int mEdgeNodeI[];
    int mEdgeNodeJ[];
    int mSplitEdges[];

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    GeometrySplittingUtils(const GeometryType& rInputGeometry, const Vector& rNodalDistances);

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

    GeometryType GetInputGeometry() const;

    Vector GetNodalDistances() const;

    /**
     * Divides the input geometry according to the provided distance data.
     * @return rAuxPoints: Reference to the pointer vector set containing the original nodes plus the intersection ones.
     * @return rPositiveSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the positive subdivisions.
     * @return rNegativeSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the negative subdivisions.
     */
    virtual bool DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                std::vector < PointGeometryType >& rPositiveSubdivisions,
                                std::vector < PointGeometryType >& rNegativeSubdivisions);

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

    /**
    * Returns true if the element is split and false otherwise.
    */
    bool IsSplit();

    /**
    * Returns true if the element is split and false otherwise.
    * @return rAuxPoints: Reference to the pointer vector set containing the original nodes plus the intersection ones.
    * @return rPositiveSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the positive subdivisions.
    * @return rNegativeSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the negative subdivisions.
    */
    void SetIntersectionPointsCondensationMatrix(Matrix& rIntPointCondMatrix,
                                                 const unsigned int splitEdgesNumber);

    ///@}
    ///@name Protected Operations
    ///@{

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

    const GeometryType& mrInputGeometry;
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

    /// Pointer definition of TriangleSplittingUtils
    KRATOS_CLASS_POINTER_DEFINITION(TriangleSplittingUtils);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    TriangleSplittingUtils(const GeometryType& rInputGeometry, const array_1d< double, 3 >& rNodalDistances);

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
                        std::vector < PointGeometryType >& rPositiveSubdivisions,
                        std::vector < PointGeometryType >& rNegativeSubdivisions) override;

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
