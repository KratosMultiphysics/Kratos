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

#if !defined(CONTACT_LINE_DIVIDE_GEOMETRY)
#define CONTACT_LINE_DIVIDE_GEOMETRY

// System includes

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/point.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/indexed_object.h"
#include "containers/pointer_vector_set.h"
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

class KRATOS_API(DROPLET_DYNAMICS_APPLICATION) ContactLineIndexedPoint : public Point, public IndexedObject
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of IndexedPoint
    KRATOS_CLASS_POINTER_DEFINITION(ContactLineIndexedPoint);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty constructor
    ContactLineIndexedPoint();

    /// Auxiliar constructor
    explicit ContactLineIndexedPoint(const unsigned int Id);

    /// Default constructor
    ContactLineIndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id);

    /// Destructor
    ~ContactLineIndexedPoint();

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point);
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    };

    void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point);
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

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ContactLineIndexedPoint& rThis) {
    return rIStream;
};

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ContactLineIndexedPoint& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}; 

template<class TPointType>
class KRATOS_API(DROPLET_DYNAMICS_APPLICATION) ContactLineDivideGeometry : public DivideGeometry<TPointType>
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ContactLineDivideGeometry
    KRATOS_CLASS_POINTER_DEFINITION(ContactLineDivideGeometry);

    // General type definitions
    typedef Geometry <TPointType>                                 GeometryType;
    typedef IndexedPoint                                 IndexedPointType;
    typedef IndexedPoint::Pointer                        IndexedPointPointerType;
    typedef Geometry <IndexedPoint>                    IndexedPointGeometryType;//<Node<3>>  
    typedef Geometry <IndexedPoint>::Pointer           IndexedPointGeometryPointerType;
    typedef PointerVectorSet<IndexedPointType, IndexedObject>       IndexedPointsContainerType;

    bool mIsSplit;          // True if the element is split.

    int mSplitEdgesNumber;  // Number of split edges.
    int mDivisionsNumber;   // Number of generated subdivisions.

    IndexedPointsContainerType mAuxPointsContainer;                         // Indexed points container to store the original plus the intersection points.
    // using DivideGeometry<TPointType>::mPositiveSubdivisions;  // Array to store the generated positive subdivisions geometries.
    // using DivideGeometry<TPointType>::mNegativeSubdivisions;  // Array to store the generated negative subdivisions geometries.
    // using DivideGeometry<TPointType>::mPositiveInterfaces;    // Array to store the generated positive interfaces geometries.
    // using DivideGeometry<TPointType>::mNegativeInterfaces;    // Array to store the generated negative interfaces geometries.
    // using DivideGeometry<TPointType>::mPositiveInterfacesParentIds;              // Array to store the parent subgeometries ids of the generated positive interfaces.
    // using DivideGeometry<TPointType>::mNegativeInterfacesParentIds;              // Array to store the parent subgeometries ids of the generated negative interfaces.
    std::vector<IndexedPointGeometryPointerType> mPositiveSubdivisions; // Array to store the generated positive subdivisions geometries.
    std::vector<IndexedPointGeometryPointerType> mNegativeSubdivisions; // Array to store the generated negative subdivisions geometries.
    std::vector<IndexedPointGeometryPointerType> mPositiveInterfaces;   // Array to store the generated positive interfaces geometries.
    std::vector<IndexedPointGeometryPointerType> mNegativeInterfaces;   // Array to store the generated negative interfaces geometries.
    std::vector<unsigned int> mPositiveInterfacesParentIds;             // Array to store the parent subgeometries ids of the generated positive interfaces.
    std::vector<unsigned int> mNegativeInterfacesParentIds;             // Array to store the parent subgeometries ids of the generated negative interfaces.

    std::vector < unsigned int > mContactInterface;     // Zero or One, gives the interface that contacts the solid
    std::vector < unsigned int > mContactEdge;          // Zero, One, or Two, gives the contact edge of the contact interface 
    std::vector < IndexedPointGeometryPointerType > mContactLine; // Object to store the contact line(s) (intersection of the interface with solid).
    //std::vector < unsigned int > mContactLineNodeIds;   // Object to store the contact line(s)' pair of node Ids in order (e.g. 1,2 , 1,3)
    std::vector < unsigned int > mContactFace;          // Object to store the face (local) number corresponding to mContactLine.
    

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ContactLineDivideGeometry(const GeometryType& rInputGeometry, const Vector& rNodalDistances);

    /// Destructor
    virtual ~ContactLineDivideGeometry();

    ///@}
    ///@name Access
    ///@{    

    /**
     * @brief Get the Positive Subdivisions object
     * This method returns the container with the positive side subgeometries
     * @return std::vector<IndexedPointGeometryPointerType> A vector containing pointers to the positive side subgeometries
     */
    std::vector<IndexedPointGeometryPointerType> GetPositiveSubdivisions() const;

    /**
     * @brief Get the Negative Subdivisions object
     * This method returns the container with the negative side subgeometries
     * @return std::vector<IndexedPointGeometryPointerType> A vector containing pointers to the negative side subgeometries
     */
    std::vector<IndexedPointGeometryPointerType> GetNegativeSubdivisions() const;

    /**
     * @brief Get the Positive Interfaces object
     * This method returns the container with the positive side interfaces
     * @return std::vector<IndexedPointGeometryPointerType> A vector containing pointers to the positive side subinterfaces
     */
    std::vector<IndexedPointGeometryPointerType> GetPositiveInterfaces() const;

    /**
     * @brief Get the Negaitive Interfaces object
     * This method returns the container with the negative side interfaces
     * @return std::vector<IndexedPointGeometryPointerType> A vector containing pointers to the negative side subinterfaces
     */
    std::vector<IndexedPointGeometryPointerType> GetNegativeInterfaces() const;

    /**
     * @brief Get the Positive Interfaces Parent Ids object
     * This method returns the container with the positive side interfaces parent ids
     * @return std::vector<IndexedPointGeometryPointerType> A vector containing the ids of the positive side interfaces parents
     */
    std::vector<unsigned int> GetPositiveInterfacesParentIds() const;

    /**
     * @brief Get the Negative Interfaces Parent Ids object
     * This method returns the container with the negative side interfaces parent ids
     * @return std::vector<IndexedPointGeometryPointerType> A vector containing the ids of the negative side interfaces parents
     */
    std::vector<unsigned int> GetNegativeInterfacesParentIds() const;

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
     * Returns the member vector containing the edges node I ids.
     */
    virtual const std::vector<int>& GetEdgeIdsI() const = 0;

    /**
     * Returns the member vector containing the edges node J ids.
     */
    virtual const std::vector<int>& GetEdgeIdsJ() const = 0;

    /**
     * Returns the member vector containing the split edges member vector.
     */
    virtual std::vector<int>& GetSplitEdges() = 0;

    /**
     * Returns the nodal distance values.
     */
    Vector GetNodalDistances() const;

    /**
     * Returns a reference to the input geometry.
     */
    GeometryType GetInputGeometry() const;

    /**
     * Divides the input geometry according to the provided distance data.
     */
    virtual void GenerateDivision() = 0;

    /**
     * Generates a list containing the intersection interface geometries for either the positive or the negative element subdivisions.
     */
    virtual void GenerateIntersectionsSkin() = 0;

    /**
     * Generates a list containing the exterior (boundary) faces geometries for either the positive or the negative element subdivisions.
     * @param rExteriorFacesVector Vector containing the generated exterior subfaces geometries
     * @param rExteriorFacesParentSubdivisionsIdsVector Vector containing the ids of the parent subdivision of each subface
     * @param rSubdivisionsContainer Positive or negative parent geometry subdivisions container
     */
    virtual void GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer)
    {
        KRATOS_ERROR << "Accessing base class \'GenerateExteriorFaces\' method."<< std::endl;
    };

    /**
     * Given a father face id, generates a list containing the exterior (boundary)
     * faces geometries belonging to either the positive or negative side of that that father face.
     * @param rExteriorFacesVector Vector containing the generated exterior subfaces geometries
     * @param rExteriorFacesParentSubdivisionsIdsVector Vector containing the ids of the parent subdivision of each subface
     * @param rSubdivisionsContainer Positive or negative parent geometry subdivisions container
     * @param FatherFaceId Father face in where the positive exterior faces are to be obtained
     */
    virtual void GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer,
        const unsigned int FatherFaceId)
    {
        KRATOS_ERROR << "Accessing base class \'GenerateExteriorFaces\' method."<< std::endl;
    };

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
    ContactLineDivideGeometry& operator=(ContactLineDivideGeometry const& rOther);

    /// Copy constructor.
    ContactLineDivideGeometry(ContactLineDivideGeometry const& rOther)
        : DivideGeometry<TPointType>(rOther.mrInputGeometry, rOther.mrNodalDistances), mrInputGeometry(rOther.mrInputGeometry) , mrNodalDistances(rOther.mrNodalDistances) {};

    ///@}

};// class DivideGeometry

}//namespace Kratos
#endif /* CONTACT_LINE_DIVIDE_GEOMETRY defined */
