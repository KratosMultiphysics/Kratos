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
#include "utilities/split_triangle.c"
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
    IndexedPoint() 
        : Point<3>() , IndexedObject(0) {};

    /// Default constructor
    IndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id) 
        : Point<3>(rCoords) , IndexedObject(Id) {};

    /// Destructor
    ~IndexedPoint() {};

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
    std::string Info() const override {
        std::stringstream info_string;
        info_string << "Indexed point class created as a combination of the Point<3> and IndexedObject classes.\n";
        info_string << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
        return info_string.str();
    };

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "Indexed point class created as a combination of the Point<3> and IndexedObject classes.\n";
        rOStream << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
    };

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Indexed point object constructed with:\n";
        rOStream << "   Index value: " << this->Id() << "\n";
        const array_1d<double, 3> point_coords = this->Coordinates();
        std::stringstream coordinates_buffer;
        for (unsigned int i = 0; i < 3; ++i) {
            coordinates_buffer << std::to_string(point_coords(i)) << " ";
        }
        rOStream << "   Coordinates: " << coordinates_buffer.str();
    };

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Point<3>);
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Point<3>);
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
    typedef Node<3>                                                                   NodeType;
    typedef Point<3>                                                                 PointType;
    typedef Geometry < NodeType >                                                 GeometryType;
    typedef Geometry < PointType >                                           PointGeometryType;

    typedef IndexedPoint                                                      IndexedPointType;
    typedef boost::shared_ptr<IndexedPoint>                            IndexedPointPointerType;
    typedef PointerVectorSet<IndexedPointType, IndexedObject>       IndexedPointsContainerType;
    typedef IndexedPointsContainerType::iterator                     IndexedPointsIteratorType;
    
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
    virtual std::string Info() const {
        return "Base class for geometries splitting operations.";
    };

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations.";
    };

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations constructed with:\n";
        rOStream << "   Geometry type: " << mrInputGeometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            distances_buffer << std::to_string(mrNodalDistances(i)) << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

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
    virtual bool DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                std::vector < PointGeometryType >& rPositiveSubdivisions,
                                std::vector < PointGeometryType >& rNegativeSubdivisions) {    

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
    ///@name Serialization
    ///@{

    // friend class Serializer;

    // virtual void save(Serializer& rSerializer) const
    // {
    //     KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeometrySplittingUtils);
    // }

    // virtual void load(Serializer& rSerializer)
    // {
    //     KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeometrySplittingUtils);
    // }

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
    std::string Info() const override {
        return "Triangle splitting operations utility.";
    };

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "Triangle splitting operations utility.";
    };

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Triangle splitting operations utility constructed with:\n";
        rOStream << "   Geometry type: " << mrInputGeometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            distances_buffer << std::to_string(mrNodalDistances(i)) << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

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
                        std::vector < PointGeometryType >& rNegativeSubdivisions) override {    
        // Fill the auxiliar points vector set
        if (this->IsSplit()) {

            // Set the triangle geometry parameters
            const unsigned int n_nodes = 3;
            const unsigned int n_edges = 3;

            const int edge_i[] = {0, 0, 1};
            const int edge_j[] = {1, 2, 2};

            int split_edge[] = {0, 1, 2, -1, -1, -1}; // {node_0, node_1, node_2, node_edge_01, node_edge_12, node_edge_20}

            // Clear the auxiliar vector points set
            rAuxPoints.clear();

            // Add the original geometry points
            for (unsigned int i = 0; i < n_nodes; ++i) {
                const array_1d<double, 3> aux_point_coords = mrInputGeometry[i].Coordinates();
                IndexedPointPointerType paux_point(new IndexedPoint(aux_point_coords, i));
                rAuxPoints.push_back(paux_point);
            }
            
            // Decide the splitting pattern
            unsigned int aux_node_id = n_nodes;

            for(unsigned int idedge = 0; idedge < n_edges; ++idedge) {
                const unsigned int edge_node_i = edge_i[idedge];
                const unsigned int edge_node_j = edge_j[idedge];
                
                // Check if the edge is split
                if(mrNodalDistances(edge_node_i) * mrNodalDistances(edge_node_j) < 0) {
                    // Set the new node id. in the split edge array corresponding slot
                    split_edge[idedge + n_nodes] = aux_node_id;

                    // Edge nodes coordinates
                    const array_1d<double, 3> i_node_coords = mrInputGeometry[edge_node_i].Coordinates();
                    const array_1d<double, 3> j_node_coords = mrInputGeometry[edge_node_j].Coordinates();
                    
                    // Compute the coordinates of the point on the edge
                    const double aux_node_rel_location = std::abs (mrNodalDistances(edge_node_i)/(mrNodalDistances(edge_node_j)-mrNodalDistances(edge_node_i))) ;
                    array_1d<double, 3> aux_point_coords;
                    for (unsigned int comp = 0; comp < 3; ++comp) {
                        aux_point_coords(comp) = i_node_coords(comp)*aux_node_rel_location + j_node_coords(comp)*(1.0-aux_node_rel_location);
                    }
                    
                    // Add the intersection point to the auxiliar points array
                    IndexedPointPointerType paux_point(new IndexedPoint(aux_point_coords, aux_node_id));
                    rAuxPoints.push_back(paux_point);
                }

                aux_node_id++;
            }
            
            // Call the splitting mode computation function
            int edge_ids[3];
            TriangleSplitMode(split_edge, edge_ids);

            KRATOS_WATCH(split_edge[0])
            KRATOS_WATCH(split_edge[1])
            KRATOS_WATCH(split_edge[2])
            KRATOS_WATCH(split_edge[3])
            KRATOS_WATCH(split_edge[4])
            KRATOS_WATCH(split_edge[5])
            KRATOS_WATCH(edge_ids[0])
            KRATOS_WATCH(edge_ids[1])
            KRATOS_WATCH(edge_ids[2])

            // Call the splitting function
            int t[12];          // Ids of the generated subdivisions
            int n_divisions;    // Number of generated subdivisions
            int n_split_edges;  // Number of split edges
            int n_int = 0;      // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            Split_Triangle(edge_ids, t, &n_divisions, &n_split_edges, &n_int);

            for (unsigned int i=0; i<12; ++i) {KRATOS_WATCH(t[i]);}
            KRATOS_WATCH(n_divisions)
            KRATOS_WATCH(n_split_edges)
            KRATOS_WATCH(n_int)

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < n_divisions; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2; 
                TriangleGetNewConnectivityGID(idivision, t, split_edge, &i0, &i1, &i2);

                std::cout << "Subdivision " << idivision << " i0: " << i0 << " i1: " << i1 << " i2: " << i2 << std::endl;
                
                // Get a pointer to the point objects corresponding to the indices above
                boost::shared_ptr< Point<3> > p_point_0(new Point<3>((rAuxPoints.find(i0))->Coordinates()));
                boost::shared_ptr< Point<3> > p_point_1(new Point<3>((rAuxPoints.find(i1))->Coordinates()));
                boost::shared_ptr< Point<3> > p_point_2(new Point<3>((rAuxPoints.find(i2))->Coordinates()));
                
                // Generate a triangle (with the same type as the input one but generated with points) with the subdivision pts.
	            Geometry< Point<3> >::PointsArrayType PointsArray;
	            PointsArray.reserve(3);
                PointsArray.push_back(p_point_0);
                PointsArray.push_back(p_point_1);
                PointsArray.push_back(p_point_2);
                PointGeometryType aux_partition(PointsArray);

                // Determine if the subdivision is wether in the negative or the positive side
                bool is_positive;
                if ((i0 == 0) || (i0 == 1) || (i0 == 2)) {
                    is_positive = mrNodalDistances(i0) < 0.0 ? false : true;
                } else if ((i1 == 0) || (i1 == 1) || (i1 == 2)) {
                    is_positive = mrNodalDistances(i1) < 0.0 ? false : true;
                } else {
                    is_positive = mrNodalDistances(i2) < 0.0 ? false : true;
                }

                // Add the generated triangle to its corresponding partition arrays
                if (is_positive) {
                    rPositiveSubdivisions.push_back(aux_partition);
                } else {
                    rNegativeSubdivisions.push_back(aux_partition);
                }   
            }

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
    ///@name Serialization
    ///@{

    // friend class Serializer;

    // void save(Serializer& rSerializer) const 
    // {
    //     KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeometrySplittingUtils);
    // }

    // void load(Serializer& rSerializer) 
    // {
    //     KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeometrySplittingUtils);
    // }

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
