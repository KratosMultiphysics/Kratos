//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_UTILITIES_H_INCLUDED )
#define  KRATOS_MAPPER_UTILITIES_H_INCLUDED

// System includes

// External includes
#ifdef KRATOS_USING_MPI
#include "mpi.h" // for GetCurrentTime()
#endif

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h" // Cross Product
#include "utilities/openmp_utils.h" // for GetCurrentTime()


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

///@}
///@name Kratos Classes
///@{

/// Auxiliary functions for the MappingApplication
/** This class provides a set of auxiliary functions that are used by several other functions / classes
*/
class MapperUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MapperUtilities);

    ///@}
    ///@name  Enum's
    ///@{

    enum InterfaceObjectConstructionType
    {
        Node_Coords,
        Condition_Center,
        Condition_Gauss_Point,
        // Point or Coordinates or sth => for Contact => create with list of points/coords,
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~MapperUtilities() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static double GetCurrentTime()
    {
        double current_time;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        int mpi_initialized;
        MPI_Initialized(&mpi_initialized);
        if (mpi_initialized)   // parallel execution, i.e. mpi imported in python
        {
            current_time = MPI_Wtime();
        }
        else     // serial execution, i.e. mpi NOT imported in python
        {
            current_time = OpenMPUtils::GetCurrentTime();
        }
#else // serial compilation
        current_time = OpenMPUtils::GetCurrentTime();
#endif

        return current_time;
    }

    static int ComputeNumberOfNodes(ModelPart& rModelPart)
    {
        int num_nodes = rModelPart.GetCommunicator().LocalMesh().NumberOfNodes();
        rModelPart.GetCommunicator().SumAll(num_nodes); // Compute the sum among the partitions
        return num_nodes;
    }

    static int ComputeNumberOfConditions(ModelPart& rModelPart)
    {
        int num_conditions = rModelPart.GetCommunicator().LocalMesh().NumberOfConditions();
        rModelPart.GetCommunicator().SumAll(num_conditions); // Compute the sum among the partitions
        return num_conditions;
    }

    static int ComputeNumberOfElements(ModelPart& rModelPart)
    {
        int num_elements = rModelPart.GetCommunicator().LocalMesh().NumberOfElements();
        rModelPart.GetCommunicator().SumAll(num_elements); // Compute the sum among the partitions
        return num_elements;
    }

    static double ComputeDistance(const array_1d<double, 3>& rCoords1,
                                  const array_1d<double, 3>& rCoords2)
    {
        return sqrt(pow(rCoords1[0] - rCoords2[0] , 2) +
                    pow(rCoords1[1] - rCoords2[1] , 2) +
                    pow(rCoords1[2] - rCoords2[2] , 2));
    }

    template <typename T>
    static double ComputeMaxEdgeLengthLocal(const T& rEntityContainer)
    {
        double max_element_size = 0.0f;
        // Loop through each edge of a geometrical entity ONCE
        for (auto& r_entity : rEntityContainer)
        {
            for (std::size_t i = 0; i < (r_entity.GetGeometry().size() - 1); ++i)
            {
                for (std::size_t j = i + 1; j < r_entity.GetGeometry().size(); ++j)
                {
                    double edge_length = ComputeDistance(r_entity.GetGeometry()[i].Coordinates(),
                                                         r_entity.GetGeometry()[j].Coordinates());
                    max_element_size = std::max(max_element_size, edge_length);
                }
            }
        }
        return max_element_size;
    }

    static double ComputeMaxEdgeLengthLocal(const ModelPart::NodesContainerType& rNodes)
    {
        double max_element_size = 0.0f;
        // TODO modify loop such that it loop only once over the nodes
        for (auto& r_node_1 : rNodes)
        {
            for (auto& r_node_2 : rNodes)
            {
                double edge_length = ComputeDistance(r_node_1.Coordinates(),
                                                     r_node_2.Coordinates());
                max_element_size = std::max(max_element_size, edge_length);
            }
        }
        return max_element_size;
    }

    static double ComputeSearchRadius(ModelPart& rModelPart1, ModelPart& rModelPart2, const int EchoLevel)
    {
        double search_radius = std::max(ComputeSearchRadius(rModelPart1, EchoLevel),
                                        ComputeSearchRadius(rModelPart2, EchoLevel));
        return search_radius;
    }

    static double ComputeSearchRadius(ModelPart& rModelPart, const int EchoLevel)
    {
        double search_safety_factor = 1.2;
        double max_element_size = 0.0;

        int num_conditions_global = ComputeNumberOfConditions(rModelPart);
        int num_elements_global = ComputeNumberOfElements(rModelPart);

        if (num_conditions_global > 0)
        {
            max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Conditions());

        }
        else if (num_elements_global > 0)
        {
            max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Elements());
        }
        else
        {
            if (EchoLevel >= 2 && rModelPart.GetCommunicator().MyPID() == 0)
                std::cout << "MAPPER WARNING, no conditions/elements for search radius "
                          << "computations in ModelPart \"" << rModelPart.Name() << "\" found, "
                          << "using nodes (less efficient, bcs search radius will be larger)" << std::endl;
            max_element_size = ComputeMaxEdgeLengthLocal(rModelPart.GetCommunicator().LocalMesh().Nodes());
        }

        rModelPart.GetCommunicator().MaxAll(max_element_size); // Compute the maximum among the partitions
        return max_element_size * search_safety_factor;
    }

    static double ComputeConservativeFactor(const double NumNodesOrigin,
                                            const double NumNodesDestination)
    {
        // NumNodes* are casted to doubles in order to use the double devision
        // if this function would take ints, then the return value would also be an int!
        KRATOS_ERROR_IF(NumNodesDestination <= 0) << "Division by zero!" << std::endl;

        return NumNodesOrigin / NumNodesDestination;
    }

    // **********************************************************************
    // Functions related to Bounding Boxes **********************************
    // **********************************************************************
    static void ComputeLocalBoundingBox(ModelPart& rModelPart,
                                        double* pLocalBoundingBox)
    {
        // xmax, xmin,  ymax, ymin,  zmax, zmin
        // loop over all nodes (local and ghost(necessary if conditions have only ghost nodes) )
        for (auto &r_node : rModelPart.Nodes())
        {
            pLocalBoundingBox[0] = std::max(r_node.X(), pLocalBoundingBox[0]);
            pLocalBoundingBox[1] = std::min(r_node.X(), pLocalBoundingBox[1]);
            pLocalBoundingBox[2] = std::max(r_node.Y(), pLocalBoundingBox[2]);
            pLocalBoundingBox[3] = std::min(r_node.Y(), pLocalBoundingBox[3]);
            pLocalBoundingBox[4] = std::max(r_node.Z(), pLocalBoundingBox[4]);
            pLocalBoundingBox[5] = std::min(r_node.Z(), pLocalBoundingBox[5]);
        }
    }

    static void ComputeBoundingBoxWithTolerance(double* pLocalBoundingBox,
            const double Tolerance,
            double* pLocalBoundingBoxWithTolerance)
    {
        // xmax, xmin,  ymax, ymin,  zmax, zmin
        pLocalBoundingBoxWithTolerance[0] = pLocalBoundingBox[0] + Tolerance;
        pLocalBoundingBoxWithTolerance[1] = pLocalBoundingBox[1] - Tolerance;
        pLocalBoundingBoxWithTolerance[2] = pLocalBoundingBox[2] + Tolerance;
        pLocalBoundingBoxWithTolerance[3] = pLocalBoundingBox[3] - Tolerance;
        pLocalBoundingBoxWithTolerance[4] = pLocalBoundingBox[4] + Tolerance;
        pLocalBoundingBoxWithTolerance[5] = pLocalBoundingBox[5] - Tolerance;
    }

    static std::vector<double> ComputeModelPartBoundingBox(ModelPart& rModelPart){
        double init_max = std::numeric_limits<double>::max();
        double init_min = std::numeric_limits<double>::lowest();
        
        std::vector<double> model_part_bbox{ init_min, init_max, init_min, init_max, init_min, init_max };
        // xmax, xmin,  ymax, ymin,  zmax, zmin
        MapperUtilities::ComputeLocalBoundingBox(rModelPart, &model_part_bbox[0]);

        // this is not the most efficient way, but it is less implementational 
        // effort and anyway only used for debugging
        rModelPart.GetCommunicator().MaxAll(model_part_bbox[0]);
        rModelPart.GetCommunicator().MinAll(model_part_bbox[1]);
        rModelPart.GetCommunicator().MaxAll(model_part_bbox[2]);
        rModelPart.GetCommunicator().MinAll(model_part_bbox[3]);
        rModelPart.GetCommunicator().MaxAll(model_part_bbox[4]);
        rModelPart.GetCommunicator().MinAll(model_part_bbox[5]);

        return model_part_bbox;
    }

    static bool ComputeBoundingBoxIntersection(std::vector<double>& bbox_origin,
                                               std::vector<double>& bbox_destination)
    {   
        std::vector<double> bbox_origin_tol(6);
        std::vector<double> bbox_destination_tol(6);
        
        const double tolerance = std::numeric_limits<double>::epsilon();

        MapperUtilities::ComputeBoundingBoxWithTolerance(&bbox_origin[0], tolerance, &bbox_origin_tol[0]);
        MapperUtilities::ComputeBoundingBoxWithTolerance(&bbox_destination[0], tolerance, &bbox_destination_tol[0]);

        // xmax, xmin,  ymax, ymin,  zmax, zmin

        // check destination bbox corner points in origin bbox
        // check lower point
        Point point_to_check_1(bbox_destination[1], bbox_destination[3], bbox_destination[5]);
        if (MapperUtilities::PointIsInsideBoundingBox(&bbox_origin_tol[0], point_to_check_1))
            return true;
        // check higher point
        Point point_to_check_2(bbox_destination[0], bbox_destination[2], bbox_destination[4]);
        if (MapperUtilities::PointIsInsideBoundingBox(&bbox_origin_tol[0], point_to_check_2))
            return true;

        // check origin bbox corner points in destination bbox
        // check lower point
        Point point_to_check_3(bbox_origin[1], bbox_origin[3], bbox_origin[5]);
        if (MapperUtilities::PointIsInsideBoundingBox(&bbox_destination_tol[0], point_to_check_3))
            return true;
        // check higher point
        Point point_to_check_4(bbox_origin[0], bbox_origin[2], bbox_origin[4]);
        if (MapperUtilities::PointIsInsideBoundingBox(&bbox_destination_tol[0], point_to_check_4))
            return true;

        return false;
    }

    static bool PointIsInsideBoundingBox(double* BoundingBox,
                                         const Point& rPoint)
    {   // The Bounding Box should have some tolerance already!
        bool is_inside = false;
        if (rPoint.X() < BoundingBox[0] && rPoint.X() > BoundingBox[1])   // check x-direction
        {
            if (rPoint.Y() < BoundingBox[2] && rPoint.Y() > BoundingBox[3])   // check y-direction
            {
                if (rPoint.Z() < BoundingBox[4] && rPoint.Z() > BoundingBox[5])   // check z-direction
                {
                    is_inside = true;
                }
            }
        }
        return is_inside;
    }

    static std::string PrintModelPartBoundingBoxes(std::vector<double>& bbox_origin,
                                          std::vector<double>& bbox_destination)
    {
        std::stringstream buffer;
        buffer << "Bounding Box Origin: "
               << MapperUtilities::BoundingBoxStingStream(&bbox_origin[0]) << ", "
               << "Bounding Box Destination: "
               << MapperUtilities::BoundingBoxStingStream(&bbox_destination[0]) 
               << std::endl;
        return buffer.str();
    }

    static std::string BoundingBoxStingStream(double* pBoundingBox)
    {
        // xmax, xmin,  ymax, ymin,  zmax, zmin
        std::stringstream buffer;
        buffer << "[" << pBoundingBox[1] << " "  // xmin
               << pBoundingBox[3] << " "         // ymin
               << pBoundingBox[5] << "]|["       // zmin
               << pBoundingBox[0] << " "         // xmax
               << pBoundingBox[2] << " "         // ymax
               << pBoundingBox[4] << "]";        // zmax
        return buffer.str();
    }

    // **********************************************************************
    // Functions related to Projection **************************************
    // **********************************************************************
    static bool ProjectPointToLine(const Geometry<Node<3>>* pGeometry,
                                   const array_1d<double, 3>& GlobalCoords,
                                   array_1d<double, 3>& rLocalCoords,
                                   double& rDistance)
    {
        // xi,yi are Nodal Coordinates, n is the destination condition's unit normal
        // and d is the distance along n from the point to its projection in the condition
        // | DestX-0.5(x1+x2) |   | 0.5(x2-x1)  nx |   | Chi |
        // |                  | = |                | . |     |
        // | DestY-0.5(y1+y2) |   | 0.5(y2-y1)  ny |   |  d  |

        Matrix transform_matrix(2, 2, false);
        Matrix inv_transform_matrix(2, 2, false);
        double det;

        array_1d<double, 2> RHS;
        array_1d<double, 2> result_vec;

        RHS[0] = GlobalCoords[0] - 0.5 * (pGeometry->GetPoint(0).X() + pGeometry->GetPoint(1).X());
        RHS[1] = GlobalCoords[1] - 0.5 * (pGeometry->GetPoint(0).Y() + pGeometry->GetPoint(1).Y());

        transform_matrix(0, 0) = 0.5 * (pGeometry->GetPoint(1).X() - pGeometry->GetPoint(0).X());
        transform_matrix(1, 0) = 0.5 * (pGeometry->GetPoint(1).Y() - pGeometry->GetPoint(0).Y());

        array_1d<double, 3> rNormal;
        CalculateLineNormal(pGeometry, rNormal);
        transform_matrix(0, 1) = rNormal[0];
        transform_matrix(1, 1) = rNormal[1];

        MathUtils<double>::InvertMatrix2(transform_matrix, inv_transform_matrix, det);
        noalias(result_vec) = prod(inv_transform_matrix, RHS);

        rLocalCoords[0] = result_vec[0];
        rLocalCoords[1] = 0.0f;
        rLocalCoords[2] = 0.0f;

        rDistance = result_vec[1];

        bool is_inside = false;

        if (fabs(rLocalCoords[0]) <= 1.0f + MapperUtilities::tol_local_coords)
        {
            is_inside = true;
        }
        return is_inside;
    }

    static bool ProjectPointToTriangle(Geometry<Node<3>>* pGeometry,
                                       const array_1d<double, 3>& GlobalCoords,
                                       array_1d<double, 3>& rLocalCoords,
                                       double& rDistance)
    {
        // xi,yi,zi are Nodal Coordinates, n is the destination condition's unit normal
        // and d is the distance along n from the point to its projection in the condition
        // | DestX-x1 |   | x2-x1  x3-x1  nx |   | Chi |
        // | DestY-y1 | = | y2-y1  y3-y1  ny | . | Eta |
        // | DestZ-z1 |   | z2-z1  z3-z1  nz |   |  d  |

        Matrix transform_matrix(3, 3, false);
        Matrix inv_transform_matrix(3, 3, false);
        double det;

        array_1d<double, 3> RHS;
        array_1d<double, 3> result_vec;

        RHS[0] = GlobalCoords[0] - pGeometry->GetPoint(0).X();
        RHS[1] = GlobalCoords[1] - pGeometry->GetPoint(0).Y();
        RHS[2] = GlobalCoords[2] - pGeometry->GetPoint(0).Z();

        transform_matrix(0, 0) = pGeometry->GetPoint(1).X() - pGeometry->GetPoint(0).X();
        transform_matrix(1, 0) = pGeometry->GetPoint(1).Y() - pGeometry->GetPoint(0).Y();
        transform_matrix(2, 0) = pGeometry->GetPoint(1).Z() - pGeometry->GetPoint(0).Z();

        transform_matrix(0, 1) = pGeometry->GetPoint(2).X() - pGeometry->GetPoint(0).X();
        transform_matrix(1, 1) = pGeometry->GetPoint(2).Y() - pGeometry->GetPoint(0).Y();
        transform_matrix(2, 1) = pGeometry->GetPoint(2).Z() - pGeometry->GetPoint(0).Z();

        array_1d<double, 3> rNormal;
        CalculateTriangleNormal(pGeometry, rNormal);
        transform_matrix(0, 2) = rNormal[0];
        transform_matrix(1, 2) = rNormal[1];
        transform_matrix(2, 2) = rNormal[2];

        MathUtils<double>::InvertMatrix3(transform_matrix, inv_transform_matrix, det);
        noalias(result_vec) = prod(inv_transform_matrix, RHS);

        rLocalCoords[0] = result_vec[0];
        rLocalCoords[1] = result_vec[1];
        rLocalCoords[2] = 0.0f;

        rDistance = result_vec[2];

        bool is_inside = false;

        if (1.0f - rLocalCoords[0] - rLocalCoords[1] >= 0.0f - MapperUtilities::tol_local_coords)
        {
            if (rLocalCoords[0] >= 0.0f - MapperUtilities::tol_local_coords)
            {
                if (rLocalCoords[1] >= 0.0f - MapperUtilities::tol_local_coords)
                {
                    is_inside = true;
                }
            }
        }
        return is_inside;
    }

    static bool ProjectPointToQuadrilateral(Geometry<Node<3>>* pGeometry,
                                            const array_1d<double, 3>& GlobalCoords,
                                            array_1d<double, 3>& rLocalCoords,
                                            double& rDistance)
    {
        bool is_inside = pGeometry->IsInside(GlobalCoords, rLocalCoords, tol_local_coords);

        if (is_inside)
        {
            // Calculate Distance
            array_1d<double, 3> projection_global_coords;
            pGeometry->GlobalCoordinates(projection_global_coords, rLocalCoords);
            rDistance = ComputeDistance(GlobalCoords, projection_global_coords);
        }

        return is_inside;
    }

    // static bool ProjectPointToTetrahedra(pGeometry,
    //                                      const array_1d<double, 3>& GlobalCoords,
    //                                      array_1d<double,3>& rLocalCoords,
    //                                      double& rDistance) {
    //     // xi,yi,zi are Nodal Coordinates
    //     // | DestX-x1 |   | x2-x1  x3-x1  x4-x1 |   | Chi1 |
    //     // | DestY-y1 | = | y2-y1  y3-y1  y4-y1 | . | Chi2 |
    //     // | DestZ-z1 |   | z2-z1  z3-z1  z4-z1 |   | Chi3 |

    //     Matrix transform_matrix(3, 3, false);
    //     Matrix inv_transform_matrix(3, 3, false);
    //     double det;

    //     array_1d<double, 3> RHS;

    //     RHS[0] = GlobalCoords[0] - pGeometry->GetPoint(0).X();
    //     RHS[1] = GlobalCoords[1] - pGeometry->GetPoint(0).Y();
    //     RHS[2] = GlobalCoords[2] - pGeometry->GetPoint(0).Z();

    //     transform_matrix(0, 0) = pGeometry->GetPoint(1).X() - pGeometry->GetPoint(0).X();
    //     transform_matrix(1, 0) = pGeometry->GetPoint(1).Y() - pGeometry->GetPoint(0).Y();
    //     transform_matrix(2, 0) = pGeometry->GetPoint(1).Z() - pGeometry->GetPoint(0).Z();

    //     transform_matrix(0, 1) = pGeometry->GetPoint(2).X() - pGeometry->GetPoint(0).X();
    //     transform_matrix(1, 1) = pGeometry->GetPoint(2).Y() - pGeometry->GetPoint(0).Y();
    //     transform_matrix(2, 1) = pGeometry->GetPoint(2).Z() - pGeometry->GetPoint(0).Z();

    //     transform_matrix(0, 2) = pGeometry[3].X() - pGeometry->GetPoint(0).X();
    //     transform_matrix(1, 2) = pGeometry[3].Y() - pGeometry->GetPoint(0).Y();
    //     transform_matrix(2, 2) = pGeometry[3].Z() - pGeometry->GetPoint(0).Z();

    //     MathUtils<double>::InvertMatrix3(transform_matrix,inv_transform_matrix,det);
    //     noalias(rLocalCoords) = prod(inv_transform_matrix, RHS);

    //     bool is_inside = false;

    //     if( rLocalCoords[0] >= 0.0-MapperUtilities::tol_local_coords ) {
    //         if( rLocalCoords[1] >= 0.0-MapperUtilities::tol_local_coords ) {
    //             if( rLocalCoords[2] >= 0.0-MapperUtilities::tol_local_coords ) {
    //                 if( (rLocalCoords[0] + rLocalCoords[1] + rLocalCoords[2]) <= (1.0+MapperUtilities::tol_local_coords)) {
    //                     is_inside = true;

    //                     rDistance = ComputeDistance(GlobalCoords, pGeometry.Center());
    //                     rDistance /= pGeometry.Volume(); // Normalize Distance by Volume
    //                 }
    //             }
    //         }
    //     }

    //     return is_inside;
    // }

    // static bool ProjectPointToCondition(Geometry<Node<3>>* pGeometry,
    //                                     const array_1d<double, 3>& GlobalCoords,
    //                                     array_1d<double,3>& rLocalCoords,
    //                                     double& rDistance) {

    //     Condition::GeometryType& r_condition_geometry = pGeometry;
    //     bool is_inside = r_condition_geometry.IsInside(GlobalCoords, rLocalCoords);

    //     if (is_inside) {
    //         // Calculate Distance
    //         array_1d<double, 3> projection_global_coords;
    //         r_condition_geometry.GlobalCoordinates(projection_global_coords, rLocalCoords);
    //         rDistance = ComputeDistance(GlobalCoords, projection_global_coords);
    //         // KRATOS_WATCH(rDistance)
    //         // std::cout << "Local Coords: [ " << rLocalCoords[0] << " , " << rLocalCoords[1] << " , " << rLocalCoords[2] << " ]" << std::endl;
    //     }

    //     return is_inside;
    // }

    static bool PointLocalCoordinatesInVolume(Geometry<Node<3>>* pGeometry,
            const array_1d<double, 3>& GlobalCoords,
            array_1d<double, 3>& rLocalCoords,
            double& rDistance)
    {
        bool is_inside = pGeometry->IsInside(GlobalCoords, rLocalCoords, tol_local_coords);

        if (is_inside)
        {
            // Calculate Distance
            rDistance = ComputeDistance(GlobalCoords, pGeometry->Center());
            rDistance /= pGeometry->Volume();  // Normalize Distance by Volume

        }

        return is_inside;
    }

    static void CalculateLineNormal(const Geometry<Node<3>>* pGeometry,
                                    array_1d<double, 3>& rNormal)
    {
        // TODO use the normal calculation in geometry.h once it is available
        rNormal.clear();
        array_1d<double, 3> v1, v2;
        v1[0] = pGeometry->GetPoint(1).X() - pGeometry->GetPoint(0).X();
        v1[1] = pGeometry->GetPoint(1).Y() - pGeometry->GetPoint(0).Y();
        v1[2] = pGeometry->GetPoint(1).Z() - pGeometry->GetPoint(0).Z();
        // Assuming plane X-Y in the 2D-case
        v2[0] = 0.0f;
        v2[1] = 0.0f;
        v2[2] = 1.0f;

        // Compute the condition normal
        MathUtils<double>::CrossProduct(rNormal, v1, v2);

        rNormal /= norm_2(rNormal); // normalize the nomal (i.e. length=1)
    }

    static void CalculateTriangleNormal(const Geometry<Node<3>>* pGeometry,
                                        array_1d<double, 3>& rNormal)
    {
        // TODO use the normal calculation in geometry.h once it is available
        rNormal.clear();
        array_1d<double, 3> v1, v2;
        v1[0] = pGeometry->GetPoint(1).X() - pGeometry->GetPoint(0).X();
        v1[1] = pGeometry->GetPoint(1).Y() - pGeometry->GetPoint(0).Y();
        v1[2] = pGeometry->GetPoint(1).Z() - pGeometry->GetPoint(0).Z();

        v2[0] = pGeometry->GetPoint(2).X() - pGeometry->GetPoint(0).X();
        v2[1] = pGeometry->GetPoint(2).Y() - pGeometry->GetPoint(0).Y();
        v2[2] = pGeometry->GetPoint(2).Z() - pGeometry->GetPoint(0).Z();


        // Compute the condition normal
        MathUtils<double>::CrossProduct(rNormal, v1, v2);

        rNormal /= norm_2(rNormal); // normalize the nomal (i.e. length=1)
    }

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MapperUtilities" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


    ///@}
    ///@name Friends
    ///@{


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

    static constexpr double tol_local_coords = 1E-15; // TODO what to use here?
    // static constexpr double tol_local_coords = std::numeric_limits<double>::epsilon(); // gives Problems if nodes are in the same location


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// Default constructor.
    MapperUtilities() { }

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
    MapperUtilities& operator=(MapperUtilities const& rOther);

    //   /// Copy constructor.
    //   MapperUtilities(MapperUtilities const& rOther){}


    ///@}

}; // Class MapperUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MapperUtilities& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MapperUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_H_INCLUDED  defined