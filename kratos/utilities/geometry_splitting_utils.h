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
#include "utilities/openmp_utils.h"
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

class IndexedPoint : public Point<3>, public IndexedObject
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
    IndexedPoint() : Point<3>() , IndexedObject(0) {}

    /// Default constructor
    IndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id) 
        : Point<3>(rCoords) , IndexedObject(Id) {};

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

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
    }

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

class GeometrySplittingUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometrySplittingUtils
    KRATOS_CLASS_POINTER_DEFINITION(GeometrySplittingUtils);
    
    // General type definitions
    typedef Node<3>                                                 NodeType;
    // typedef Point<3>                                            PointType;
    // typedef PointType::CoordinatesArrayType          CoordinatesArrayType;
    typedef Geometry < NodeType >                               GeometryType;

    typedef IndexedPoint                                    IndexedPointType;
    typedef boost::shared_ptr<IndexedPoint>          IndexedPointPointerType;

    typedef PointerVectorSet<IndexedPointType, IndexedObject>    PointsContainerType;
    // typedef Geometry<PointType>                         GeometryPointType;
    // typedef GeometryData::IntegrationMethod             IntegrationMethod;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    GeometrySplittingUtils(const GeometryType& rInputGeometry, const Vector& rNodalDistances) 
        : mrInputGeometry(rInputGeometry), mrNodalDistances(rNodalDistances) {};

    /// Destructor
    ~GeometrySplittingUtils() {};

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
     * Divides the input geometry according to the provided distance data.
     * @return rAuxPoints: Reference to the pointer vector set containing the original nodes plus the intersection ones.
     * @return rPositiveSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the positive subdivisions.
     * @return rNegativeSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the negative subdivisions.
     */
    virtual bool DivideGeometry(PointerVectorSet< IndexedPoint >& rAuxPoints,
                                std::vector < GeometryType >& rPositiveSubdivisions,
                                std::vector < GeometryType >& rNegativeSubdivisions) {    

        KRATOS_ERROR << "Calling the base class geometry splitting DivideGeometry method. Call the specific geometry one.";
    }
    
    ///@}

protected:
    // TODO: DEFINE ALL THE REQUIRED ARRAYS TO BE USED IN DE DERIVED CLASSES
    
    /**
     * Returns true if the element is split and false otherwise
     */
    bool IsSplit() {
        unsigned int n_pos = 0 , n_neg = 0;

        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            if (mrNodalDistances(i) < 0.0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }

        if ((n_pos > 0) && (n_neg > 0)) {
            return true;
        } else {
            return false;
        }
    }

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const GeometryType& mrInputGeometry;
    const Vector& mrNodalDistances;

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

class TriangleSplittingUtils : public GeometrySplittingUtils
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
    TriangleSplittingUtils(const GeometryType& rInputGeometry, const Vector& rNodalDistances) 
        : GeometrySplittingUtils(rInputGeometry, rNodalDistances), mrInputGeometry(rInputGeometry), mrNodalDistances(rNodalDistances) {};

    /// Destructor
    ~TriangleSplittingUtils() {};

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
     * Divides the input geometry according to the provided distance data.
     * @return rAuxPoints: Reference to the pointer vector set containing the original nodes plus the intersection ones.
     * @return rPositiveSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the positive subdivisions.
     * @return rNegativeSubdivisions: Reference to a vector containing the nodal auxiliar ids. that conform the negative subdivisions.
     */
    bool DivideGeometry(PointerVectorSet< IndexedPoint >& rAuxPoints,
                        std::vector < GeometryType >& rPositiveSubdivisions,
                        std::vector < GeometryType >& rNegativeSubdivisions) {    
        // Fill the auxiliar points vector set
        if (this->IsSplit()) {
            // Clear the auxiliar vector points set
            rAuxPoints.clear();

            // Add the original geometry points
            for (unsigned int i = 0; i < mrInputGeometry.size(); ++i) {
                IndexedPointPointerType paux_point(new IndexedPoint(mrInputGeometry[i].Coordinates(), i));
                rAuxPoints.push_back(paux_point);
            }

            // Add the edge intersection points
            unsigned int counter = mrInputGeometry.size();
            for(unsigned int i = 0; i < mrInputGeometry.size(); ++i) {
                for(unsigned int j = i+1; j < mrInputGeometry.size(); ++j) {
                    // Check if the edge is split
                    if(mrNodalDistances(i) * mrNodalDistances(j) < 0) {
                        // Edge nodes coordinates
                        const array_1d<double, 3> i_node_coords = mrInputGeometry[i].Coordinates();
                        const array_1d<double, 3> j_node_coords = mrInputGeometry[j].Coordinates();

                        // Compute the coordinates of the point on the edge
                        const double aux_node_rel_location = std::abs (mrNodalDistances(i)/(mrNodalDistances(j)-mrNodalDistances(i))) ;
                        array_1d<double, 3> aux_point_coords;
                        for (unsigned int comp = 0; comp < 3; ++comp) {
                            aux_point_coords(comp) = i_node_coords(i)*aux_node_rel_location + j_node_coords(j)*(1.0-aux_node_rel_location);
                        }

                        // Add the intersection point the auxiliar points array
                        IndexedPointPointerType paux_point(new IndexedPoint(aux_point_coords, counter));
                        rAuxPoints.push_back(paux_point);
                    }   
                    counter++;
                }   
            }

            // Call the splitting function

            return true;

        } else {

            return false;   
        }
    }
    
    ///@}

protected:
    // TODO: DEFINE ALL THE REQUIRED ARRAYS TO BE USED IN DE DERIVED CLASSES
    
private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const GeometryType& mrInputGeometry;
    const Vector& mrNodalDistances;

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
    TriangleSplittingUtils& operator=(TriangleSplittingUtils const& rOther);

    /// Copy constructor.
    TriangleSplittingUtils(TriangleSplittingUtils const& rOther) 
        : GeometrySplittingUtils(rOther.mrInputGeometry, rOther.mrNodalDistances), 
          mrInputGeometry(rOther.mrInputGeometry), 
          mrNodalDistances(rOther.mrNodalDistances) {};


    ///@}

};// class GeometrySplittingUtils

///@name Explicit Specializations
///@{

}
#endif /* KRATOS_GEOMETRY_SPLITTING_UTILS defined */
