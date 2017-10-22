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

#if !defined(KRATOS_DIVIDE_GEOMETRY)
#define KRATOS_DIVIDE_GEOMETRY

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

class KRATOS_API(KRATOS_CORE) DivideGeometry
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of DivideGeometry
    KRATOS_CLASS_POINTER_DEFINITION(DivideGeometry);

    // General type definitions
    typedef Geometry < Node<3> >                                                  GeometryType;
    typedef IndexedPoint                                                      IndexedPointType;
    typedef typename IndexedPoint::Pointer                             IndexedPointPointerType;
    typedef Geometry < IndexedPoint >                                 IndexedPointGeometryType;
    typedef Geometry < IndexedPoint >::Pointer                 IndexedPointGeometryPointerType;
    typedef PointerVectorSet<IndexedPointType, IndexedObject>       IndexedPointsContainerType;

    int mSplitEdgesNumber;  // Number of split edges
    int mDivisionsNumber;   // Number of generated subdivisions

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    DivideGeometry(GeometryType& rInputGeometry, Vector& rNodalDistances);

    /// Destructor
    ~DivideGeometry();

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
    virtual bool GenerateDivision(IndexedPointsContainerType& rAuxPoints,
                                  std::vector < IndexedPointGeometryPointerType >& rPositiveSubdivisions,
                                  std::vector < IndexedPointGeometryPointerType >& rNegativeSubdivisions);

    /**
     * Generates a list containing the intersection interface geometries for either the positive or the negative element subdivisions.
     * @return rInterfacesVector: Reference to a std::vector containing the generated interface geometries.
     * @param rSubdivisionsVector: std::vector of subdivisions point based geometries to its intersection geometries.
     */
    virtual void GenerateIntersectionsSkin(std::vector < IndexedPointGeometryPointerType >& rInterfacesVector,
                                           IndexedPointsContainerType& rAuxPoints,
                                           const std::vector < IndexedPointGeometryPointerType >& rSubdivisionsVector);

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
    DivideGeometry& operator=(DivideGeometry const& rOther);

    /// Copy constructor.
    DivideGeometry(DivideGeometry const& rOther)
        : mrInputGeometry(rOther.mrInputGeometry) , mrNodalDistances(rOther.mrNodalDistances) {};

    ///@}

};// class DivideGeometry

}
#endif /* KRATOS_DIVIDE_GEOMETRY defined */
