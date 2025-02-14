//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Ruben Zorrilla
//
//  Collaborators:   Franziska Wahl
//

#if !defined(KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "geometries/plane_3d.h"
#include "includes/checks.h"
#include "processes/process.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "utilities/variable_utils.h"
#include "utilities/pointer_communicator.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) CalculateDiscontinuousDistanceToSkinProcessFlags
{
public:
    KRATOS_DEFINE_LOCAL_FLAG(CALCULATE_ELEMENTAL_EDGE_DISTANCES); /// Local flag to switch on/off the elemental edge distances storage
    KRATOS_DEFINE_LOCAL_FLAG(CALCULATE_ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED); /// Local flag to switch on/off the extrapolated elemental edge distances storage
    KRATOS_DEFINE_LOCAL_FLAG(USE_POSITIVE_EPSILON_FOR_ZERO_VALUES); /// Local flag to switch from positive (true) to negative (false) epsilon when replacing zero distance values.
};

/// This only calculates the distance. Calculating the inside outside should be done by a derived class of this.
/** This process takes a volume model part (with tetrahedra mesh) and a skin model part (with triangle mesh) and
     and calculates the distance to the skin for all the elements and nodes of the volume model part.
*/
template<std::size_t TDim = 3>
class KRATOS_API(KRATOS_CORE) CalculateDiscontinuousDistanceToSkinProcess : public Process
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateDiscontinuousDistanceToSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateDiscontinuousDistanceToSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor to be used.
    CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart);

    /// Constructor with option
    CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        const Flags rOptions);

    /// Constructor with parameters
    CalculateDiscontinuousDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        Parameters rParameters);

    /// Destructor.
    ~CalculateDiscontinuousDistanceToSkinProcess() override;

    ///@}
    ///@name Deleted
    ///@{

    /// Default constructor.
    CalculateDiscontinuousDistanceToSkinProcess() = delete;

    /// Copy constructor.
    CalculateDiscontinuousDistanceToSkinProcess(CalculateDiscontinuousDistanceToSkinProcess const& rOther) = delete;

    /// Assignment operator.
    CalculateDiscontinuousDistanceToSkinProcess& operator=(CalculateDiscontinuousDistanceToSkinProcess const& rOther) = delete;

    FindIntersectedGeometricalObjectsProcess mFindIntersectedObjectsProcess;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initializes discontinuous distance computation process
     * This method initializes the TO_SPLIT flag, the DISTANCE and
     * ELEMENTAL_DISTANCES variables as well as the EMBEDDED_VELOCITY
     */
    virtual void Initialize();

    /**
     * @brief Calls the FindIntersectedObjectsProcess to find the intersections
     * This method calls the FindIntersectedObjectsProcess FindIntersections method.
     */
    virtual void FindIntersections();

    /**
     * @brief Get the array containing the intersecting objects
     * This method returns an array containing pointers to the intersecting geometries
     * @return std::vector<PointerVector<GeometricalObject>>&
     */
    virtual std::vector<PointerVector<GeometricalObject>>& GetIntersections();

    /**
     * @brief Computes the elemental distance values
     * Given an intersecting objects vector, this method computes the elemental distance field
     * @param rIntersectedObjects array containing pointers to the intersecting geometries
     */
    virtual void CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects);

    /**
     * @brief Calls the FindIntersectedObjects Clear() method
     * This method calls the FindIntersectedObjects Clear() to empty the intersecting objects geometries array
     */
    void Clear() override;

    /**
     * @brief Executes the CalculateDiscontinuousDistanceToSkinProcess
     * This method automatically does all the calls required to compute the discontinuous distance function.
     */
    void Execute() override;

    /**
     * @brief Calculate embedded variable from skin double specialization
     * This method calls the specialization method for two double variables
     * @param rVariable origin double variable in the skin mesh
     * @param rEmbeddedVariable elemental double variable in the volume mesh to be computed
     */
    void CalculateEmbeddedVariableFromSkin(
        const Variable<double> &rVariable,
        const Variable<double> &rEmbeddedVariable);

    /**
     * @brief Calculate embedded variable from skin array specialization
     * This method calls the specialization method for two double variables
     * @param rVariable origin array variable in the skin mesh
     * @param rEmbeddedVariable elemental array variable in the volume mesh to be computed
     */
    void CalculateEmbeddedVariableFromSkin(
        const Variable<array_1d<double,3>> &rVariable,
        const Variable<array_1d<double,3>> &rEmbeddedVariable);

    /**
     * @brief Obtain the default parameters to construct the class.
     */
    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Access
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
protected:

    const Variable<Vector>* mpElementalDistancesVariable = &ELEMENTAL_DISTANCES;


    ///@name Protected Operations
    ///@{

    /**
     * @brief Set the Intersection Plane object
     * This method returns the plane that defines the element intersection. The 2D
     * case is considered to be a simplification of the 3D one, so a "fake" extra
     * point is created by extruding the first point in the z-direction.
     * @param rIntPtsVector array containing the intersecting points coordinates
     * @return Plane3D the plane defined by the given intersecting points coordinates
     */
    Plane3D SetIntersectionPlane(const std::vector<array_1d<double,3>> &rIntPtsVector);

    /**
     * @brief Calculates the domain characteristic length
     * This method computes the domain characteristic length as the norm of
     * the diagonal vector that joins the maximum and minimum coordinates
     * @return double the calculated characteristic length
     */
    double CalculateCharacteristicLength();

    ///@}
private:
    ///@name Member Variables
    ///@{

    ModelPart& mrSkinPart;
    ModelPart& mrVolumePart;

    Flags mOptions;

    static const std::size_t mNumNodes = TDim + 1;
    static const std::size_t mNumEdges = (TDim == 2) ? 3 : 6;

    const double mZeroToleranceMultiplier = 1e3;
    bool mDetectedZeroDistanceValues = false;
    bool mAreNeighboursComputed = false;
    bool mCalculateElementalEdgeDistances = false;
    bool mCalculateElementalEdgeDistancesExtrapolated = false;
    bool mUsePositiveEpsilonForZeroValues = true;


    const Variable<Vector>* mpElementalEdgeDistancesVariable = &ELEMENTAL_EDGE_DISTANCES;
    const Variable<Vector>* mpElementalEdgeDistancesExtrapolatedVariable = &ELEMENTAL_EDGE_DISTANCES_EXTRAPOLATED;
    const Variable<array_1d<double, 3>>* mpEmbeddedVelocityVariable = &EMBEDDED_VELOCITY;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Computes the discontinuous distance in one element
     * This method computes the discontinuous distance field for a given element
     * @param rElement1 reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     */
    void CalculateElementalDistances(
        Element& rElement1,
        PointerVector<GeometricalObject>& rIntersectedObjects);

    /**
     * @brief Computes the discontinuous edge-based distance in one element
     * This method computes the discontinuous edge-based distance field for a given element
     * @param rElement1 reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     */
    void CalculateElementalAndEdgeDistances(
        Element& rElement1,
        PointerVector<GeometricalObject>& rIntersectedObjects);

    /**
     * @brief Computes the edges intersections in one element
     * Provided a list of elemental intersecting geometries, this
     * method computes the edge intersections for a given element
     * @param rElement1 reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     * @param rEdgesContainer reference to the array containing the edges of the element of interest
     * @param rCutEdgesRatioVector array that stores the relative positions from node zero of the average intersection points
     * @param rCutExtraEdgesRatioVector array that stores the relative positions from node zero of the average intersection points of the extrapolated geometry
     * @param rIntersectionPointsArray array containing the edges intersection points
     * @return unsigned int number of cut edges
     */
    unsigned int ComputeEdgesIntersections(
        Element& rElement1,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        array_1d<double,mNumEdges> &rCutEdgesRatioVector,
        array_1d<double,mNumEdges> &rCutExtraEdgesRatioVector,
        std::vector<array_1d <double,3> > &rIntersectionPointsArray);

    /**
     * @brief Computes the intersection of a single edge
     * This method computes the intersection of a given edge with the candidate
     * intersecting geometry. This operation is performed accordingly to the working
     * space dimension using the intersection utilities implemented in intersection_utilities.h
     * @param rIntObjGeometry candidate intersecting geometry
     * @param rEdgePoint1 edge origin point
     * @param rEdgePoint2 edge end point
     * @param rIntersectionPoint intersection point
     * @return int type of intersection id (see intersection_utilities.h)
     */
    int ComputeEdgeIntersection(
        const Element::GeometryType& rIntObjGeometry,
        const Element::NodeType& rEdgePoint1,
        const Element::NodeType& rEdgePoint2,
        Point& rIntersectionPoint);

    /**
     * @brief Checks if rIntersectionPoint is already present in the
     * intersection point list in rIntersectionPointsVector for the tolerance rTolerance.
     * @param rIntersectionPoint reference to the intersection point
     * @param rIntersectionPointsVector reference to the list of already computed intersected points
     * @param rEdgeTolerance tolerance to compare two points and assess if they are equal
     * @return bool if rIntersectionPoint is present in rIntersectionPointsVector
     */
    bool CheckIfPointIsRepeated(
        const array_1d<double,3>& rIntersectionPoint,
        const std::vector<array_1d<double,3>>&  rIntersectionPointsVector,
        const double& rEdgeTolerance);

    /**
     * @brief Computes the element intersection unit normal
     * This method computes the element intersection unit normal vector using the distance function gradient.
     * @param rGeometry reference to the geometry of the element of interest
     * @param rElementalDistances array containing the ELEMENTAL_DISTANCES values
     * @param rNormal obtained unit normal vector
     */
    void ComputeIntersectionNormal(
        const Element::GeometryType& rGeometry,
        const Vector& rElementalDistances,
        array_1d<double,3> &rNormal);

    /**
     * @brief Computes the nodal distances to the intersection plane
     * This methods creates a plane from the intersection points and then calculates the nodal distances
     * to the intersection plane.
     * In presence of multiple intersections, it performs a least squares approximation of the intersection plane.
     * @param rElement Element to calculate the ELEMENTAL_DISTANCES
     * @param rIntersectedObjects Intersected objects container
     * @param rIntersectionPointsCoordinates The edges intersection points coordinates
     */
    void ComputeIntersectionPlaneElementalDistances(
        Element& rElement,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const std::vector<array_1d<double,3>>& rIntersectionPointsCoordinates);

    /**
     * @brief Computes the intersection plane approximation
     * For complex intersection patterns, this method takes a list containing
     * all the intersecting points and computes the plane that minimizes the
     * distance from all these points in a least squares sense. The approximated
     * plane is defined in terms of an origin point and its normal vector.
     * @param rElement1 reference to the element of interest
     * @param rPointsCoord list containing the coordinates of al the intersecting points
     * @param rPlaneBasePointCoords base point defining the approximated plane
     * @param rPlaneNormal normal vector defining the approximated plane
     */
    void ComputePlaneApproximation(
        const Element& rElement1,
        const std::vector< array_1d<double,3> >& rPointsCoord,
        array_1d<double,3>& rPlaneBasePointCoords,
        array_1d<double,3>& rPlaneNormal);


    /**
     * @brief Computes the elemental distances from the approximation
     * plane defined by the set of points in rPointVector.
     * @param rElement reference to the element of interest
     * @param rElementalDistances reference to the elemental distances container containing the coordinates of al the intersecting points
     * @param rPoitnVector reference to the vector containing the points to define the approximation plane
     */
    void ComputeElementalDistancesFromPlaneApproximation(
        Element& rElement,
        Vector& rElementalDistances,
        const std::vector<array_1d<double,3>>& rPointVector);

    /**
     * @brief Checks and replaces the values of the ELEMENTAL_DISTANCES vector that are
     * zero. The values are replaced by an epsilon (whose sign depends on a flag)
     * that is a fixed factor from the double precision. Can be deactivated by a flag.
     * @param rElementalDistances array containing the ELEMENTAL_DISTANCES values
     */
    void ReplaceZeroDistances(Vector& rElementalDistances);

    /**
     * @brief Checks (and corrects if needed) the intersection normal orientation
     * This method checks the orientation of the previously computed intersection normal.
     * To do that, the normal vector to each one of the intersecting geometries is
     * computed and its directo is compared against the current one. If the negative
     * votes win, the current normal vector orientation is switched.
     * @param rGeometry element of interest geometry
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     * @param rElementalDistances array containing the ELEMENTAL_DISTANCES values
     */
    void CorrectDistanceOrientation(
        const Element::GeometryType& rGeometry,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        Vector& rElementalDistances);

    /**
     * @brief Computes the normal vector to an intersecting object geometry
     * This method computes the normal vector to an intersecting object geometry.
     * @param rGeometry reference to the geometry of the intersecting object
     * @param rIntObjNormal reference to the intersecting object normal vector
     */
    void inline ComputeIntersectionNormalFromGeometry(
        const Element::GeometryType &rGeometry,
        array_1d<double,3> &rIntObjNormal);

    /**
     * @brief Checks if element is incised and then computes the uncut edges intersections of the element
     * with an averaged and extrapolated geometry. Therefore it calls 'ComputeExtrapolatedGeometryIntersections'.
     * Note: for uncut or completely cut elements no ratios of the extrapolated geometry will be calculated.
     * @param rElement reference to the element of interest
     * @param rEdgesContainer reference to the array containing the edges of the element of interest
     * @param rNumCutEdges number of cut edges of the element (by the non-extrapolated geometry)
     * @param rCutEdgesRatioVector array that stores the relative positions from node zero of the average intersection points
     * @param rExtraGeomNormal array as normal vector of the averaged and extrapolated geometry
     * @param rCutExtraEdgesRatioVector array that stores the relative positions from node zero of the additional
     * average intersection points of the extrapolated geometry
     */
    void ComputeExtrapolatedEdgesIntersectionsIfIncised(
        const Element& rElement,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        unsigned int &rNumCutEdges,
        array_1d<double,mNumEdges>& rCutEdgesRatioVector,
        array_1d<double,3> &rExtraGeomNormal,
        array_1d<double,mNumEdges>& rCutExtraEdgesRatioVector);

    /**
     * @brief Computes the uncut edges intersections of one element with an averaged and extrapolated geometry.
     * Therefore it calls 'IntersectionUtilities'.
     * It saves the edge intersections as ratios of the edge's length in rCutExtraEdgesRatioVector.
     * @param rElement reference to the element of interest
     * @param rEdgesContainer reference to the array containing the edges of the element of interest
     * @param rNumCutEdges number of cut edges of the element
     * @param rCutEdgesRatioVector array that stores the relative positions from node zero of the average intersection points
     * @param rExtraGeomNormal normal of the averaged and extrapolated geometry
     * @param rCutExtraEdgesRatioVector array that stores the relative positions from node zero of the additional
     * average intersection points of the extrapolated geometry
     */
    void ComputeExtrapolatedGeometryIntersections(
        const Element& rElement,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        unsigned int& rNumCutEdges,
        array_1d<double,mNumEdges>& rCutEdgesRatioVector,
        array_1d<double,3>& rExtraGeomNormal,
        array_1d<double,mNumEdges>& rCutExtraEdgesRatioVector);

    /**
     * @brief Converts edge ratios and edge ratios of the extrapolated geometry to elemental (node) distances
     * @param rElement reference to the element of interest
     * @param rIntersectedObjects reference to the array containing the element of interest intersecting geometries
     * @param rEdgesContainer reference to the array containing the edges of the element of interest
     * @param rCutEdgesRatioVector array that stores the relative positions from node zero of the average intersection points
     * (ELEMENTAL_EDGE_DISTANCES)
     * @param rCutExtraEdgesRatioVector array that stores the relative positions from node zero of the additional
     * average intersection points of the extrapolated geometry (ELEMENTAL_EXTRA_EDGE_DISTANCES)
     */
    void ComputeElementalDistancesFromEdgeRatios(
        Element& rElement,
        const PointerVector<GeometricalObject>& rIntersectedObjects,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        const array_1d<double,mNumEdges> &rCutEdgesRatioVector,
        const array_1d<double,mNumEdges> &rCutExtraEdgesRatioVector);

    /**
     * @brief Computes the intersection points from the intersection ratios of the edges of the element of interest
     * @param rGeometry reference to geometry of the element of interest
     * @param rEdgesContainer reference to the array containing the edges of the element of interest
     * @param rEdgeRatiosVector array containing the intersection ratios of an element's edges
     * @param rIntersectionPointsVector vector containing the intersection point arrays
     */
    void ConvertRatiosToIntersectionPoints(
        const Element::GeometryType& rGeometry,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        const array_1d<double,mNumEdges> &rEdgeRatiosVector,
        std::vector<array_1d <double,3> > &rIntersectionPointsVector);

    /**
     * @brief Checks whether the edges of an element, which are cut, all share one node
     * @param rEdge reference to the edge of interest
     * @param rIntersectionPoint average intersection point at the edge
     * @return calculated relative positions of the intersection point along the edge from node zero
     */
    double ConvertIntersectionPointToEdgeRatio(
        const Geometry<Node >& rEdge,
        const array_1d<double,3>& rIntersectionPoint);

    /**
     * @brief Checks whether the edges of an element, which are cut, all share one node
     * @param rEdge reference to the edge of interest
     * @param rEdgeRatio relative positions of the intersection point along the edge from node zero
     * @return rIntersectionPoint calculated average intersection point at the edge
     */
    array_1d<double,3> ConvertEdgeRatioToIntersectionPoint(
        const Geometry<Node >& rEdge,
        const double& rEdgeRatio);

    /**
     * @brief Checks whether the edges of an element, which are cut, all share one node
     * @param rElement reference to the element of interest
     * @param rEdgesContainer reference to the array containing the edges of the element of interest
     * @param rCutEdgesRatioVector array that stores the relative positions from node zero of the average intersection points
     * @return boolean true if cut edges share one node
     */
    bool CheckIfCutEdgesShareNode(
        const Element& rElement,
        const Element::GeometryType::GeometriesArrayType& rEdgesContainer,
        const array_1d<double,mNumEdges>& rCutEdgesRatioVector) const;

    /**
     * @brief Computes the value of any embedded variable
     * For a given array variable in the skin mesh, this method calculates the value
     * of such variable in the embedded mesh. This is done in each element of the volume
     * mesh by computing the average value of all the edges intersections. This value
     * is averaged again according to the number of intersected edges.
     * @tparam TVarType variable type
     * @param rVariable origin variable in the skin mesh
     * @param rEmbeddedVariable elemental variable in the volume mesh to be computed
     */
    template<class TVarType>
    void CalculateEmbeddedVariableFromSkinSpecialization(
        const Variable<TVarType> &rVariable,
        const Variable<TVarType> &rEmbeddedVariable)
    {
        const auto &r_int_obj_vect= this->GetIntersections();
        const int n_elems = mrVolumePart.NumberOfElements();

        KRATOS_ERROR_IF((mrSkinPart.NodesBegin())->SolutionStepsDataHas(rVariable) == false)
            << "Skin model part solution step data missing variable: " << rVariable << std::endl;

        // Initialize embedded variable value
        VariableUtils().SetNonHistoricalVariableToZero(rEmbeddedVariable, mrVolumePart.Elements());

        // Compute the embedded variable value for each element
        #pragma omp parallel for schedule(dynamic)
        for (int i_elem = 0; i_elem < n_elems; ++i_elem) {
            // Check if the current element has intersecting entities
            if (r_int_obj_vect[i_elem].size() != 0) {
                // Initialize the element values
                unsigned int n_int_edges = 0;
                auto it_elem = mrVolumePart.ElementsBegin() + i_elem;
                auto &r_geom = it_elem->GetGeometry();
                const auto edges = r_geom.GenerateEdges();

                // Loop the element of interest edges
                for (unsigned int i_edge = 0; i_edge < r_geom.EdgesNumber(); ++i_edge) {
                    // Initialize edge values
                    unsigned int n_int_obj = 0;
                    TVarType i_edge_val = rEmbeddedVariable.Zero();

                    // Check the edge intersection against all the candidates
                    for (auto &r_int_obj : r_int_obj_vect[i_elem]) {
                        Point intersection_point;
                        const int is_intersected = this->ComputeEdgeIntersection(
                            r_int_obj.GetGeometry(),
                            edges[i_edge][0],
                            edges[i_edge][1],
                            intersection_point);

                        // Compute the variable value in the intersection point
                        if (is_intersected == 1) {
                            n_int_obj++;
                            array_1d<double,3> local_coords;
                            r_int_obj.GetGeometry().PointLocalCoordinates(local_coords, intersection_point);
                            Vector int_obj_N;
                            r_int_obj.GetGeometry().ShapeFunctionsValues(int_obj_N, local_coords);
                            for (unsigned int i_node = 0; i_node < r_int_obj.GetGeometry().PointsNumber(); ++i_node) {
                                i_edge_val += r_int_obj.GetGeometry()[i_node].FastGetSolutionStepValue(rVariable) * int_obj_N[i_node];
                            }
                        }
                    }

                    // Check if the edge is intersected
                    if (n_int_obj != 0) {
                        // Update the element intersected edges counter
                        n_int_edges++;
                        // Add the average edge value (there might exist cases in where
                        // more than one geometry intersects the edge of interest).
                        it_elem->GetValue(rEmbeddedVariable) += i_edge_val / n_int_obj;
                    }
                }

                // Average between all the intersected edges
                if (n_int_edges != 0) {
                    it_elem->GetValue(rEmbeddedVariable) /= n_int_edges;
                }
            }
        }
    };

    /**
     * @brief Set the TO_SPLIT Kratos flag
     * This function sets the TO_SPLIT flag in the provided element according to the ELEMENTAL_DISTANCES values
     * Note that the zero distance case is avoided by checking the positiveness and negativeness of the nodal values
     * @param rElement Element to set the TO_SPLIT flag
     * @param ZeroTolerance Tolerance to check the zero distance values
     */
    void SetToSplitFlag(
        Element& rElement,
        const double ZeroTolerance);

    /**
     * @brief Checks the elemental edges distances if zero values of the distance
     * are detected. This ensures that the elementes detected as incised and intersected
     * are consistent with the zero-correction applied by the process.
     */
    void CheckAndCorrectEdgeDistances();

    /**
     * @brief Creates the global pointer communicator that contains all neighbours elements. In MPI, this
     * allows to get information from neighbours elements that are not in the same partition.
     */
    GlobalPointerCommunicator<Element>::Pointer CreatePointerCommunicator();
    ///@}

}; // Class CalculateDiscontinuousDistanceToSkinProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CalculateDiscontinuousDistanceToSkinProcess<>& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CalculateDiscontinuousDistanceToSkinProcess<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISCONTINUOUS_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED  defined
