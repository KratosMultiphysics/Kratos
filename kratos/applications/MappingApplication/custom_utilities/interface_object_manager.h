//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

// External includes
//@{KRATOS_EXTERNA_INCLUDES}
#include "includes/kratos_flags.h"
#include <boost/python.hpp>

// Project includes

#include "includes/define.h"
// #include "processes/process.h"
#include "includes/kratos_flags.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"

#include "spaces/ublas_space.h"

#include "interface_object.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
	typedef typename SparseSpaceType::MatrixType SparseMatrixType;
	typedef typename SparseSpaceType::VectorType VectorType;

	typedef array_1d<double,3> array_3d;
	typedef Node < 3 > PointType;
	typedef Node < 3 > ::Pointer PointTypePointer;
	typedef std::vector<PointType::Pointer> PointVector;
	typedef std::vector<PointType::Pointer>::iterator PointIterator;
	typedef std::vector<double> DistanceVector;
	typedef std::vector<double>::iterator DistanceIterator;
	typedef ModelPart::ConditionsContainerType ConditionsArrayType;

	typedef Element BaseType;
	typedef BaseType::GeometryType GeometryType;

  typedef IntegrationPoint<3> IntegrationPointType;

  /** A Vector of IntegrationPointType which used to hold
  integration points related to an integration
  method. IntegrationPoints functions used this type to return
  their results.
  */
  typedef std::vector<IntegrationPointType> IntegrationPointsArrayType;


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class InterfaceObjectManager
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManager);
    // InterfaceObjectManager(ModelPart& model_part,int comm_rank, int comm_size) : m_comm_rank(comm_rank), m_comm_size(comm_size) {
    //     m_point_comm_objects.resize(model_part.GetCommunicator().LocalMesh().Nodes().size());
    //
    //     m_candidate_send_points.resize(m_comm_size);
    //     m_send_points.resize(m_comm_size);
    //     m_matching_information.resize(m_comm_size);
    //
    //     int i = 0;
    //     for (auto &local_node : model_part.GetCommunicator().LocalMesh().Nodes()) {
    //         m_point_comm_objects[i] = InterfaceObject::Pointer( new PointCommNode(local_node) );
    //         ++i;
    //     }
    //
    //     ComputeLocalBoundingBox(model_part);
    // }

    ~InterfaceObjectManager() {
        delete [] m_bounding_box;
    }

    static InterfaceObjectManager::Pointer CreateInterfaceNodeManager(ModelPart& model_part, int comm_rank, int comm_size) {
        InterfaceObjectManager::Pointer manager = InterfaceObjectManager::Pointer(new InterfaceObjectManager());

        manager->m_comm_rank = comm_rank;
        manager->m_comm_size = comm_size;

        manager->m_point_comm_objects.resize(model_part.GetCommunicator().LocalMesh().Nodes().size());

        manager->m_candidate_send_points.resize(comm_size);
        manager->m_send_points.resize(comm_size);
        manager->m_matching_information.resize(comm_size);

        manager->m_candidate_receive_points.resize(comm_size);
        manager->m_receive_points.resize(comm_size);

        int i = 0;
        for (auto &local_node : model_part.GetCommunicator().LocalMesh().Nodes()) {
            manager->m_point_comm_objects[i] = InterfaceObject::Pointer( new InterfaceNode(local_node) );
            ++i;
        }

        return manager;
    }

    static InterfaceObjectManager::Pointer CreateInterfaceConditionManager(ModelPart& model_part, int comm_rank, int comm_size,
                                                                           const int construction_type) {
        InterfaceObjectManager::Pointer manager = InterfaceObjectManager::Pointer(new InterfaceObjectManager());

        manager->m_comm_rank = comm_rank;
        manager->m_comm_size = comm_size;

        manager->m_point_comm_objects.resize(3*model_part.GetCommunicator().LocalMesh().Conditions().size()); // TODO segfault expected here if using GPS!

        manager->m_candidate_send_points.resize(comm_size);
        manager->m_send_points.resize(comm_size);
        manager->m_matching_information.resize(comm_size);

        manager->m_candidate_receive_points.resize(comm_size);
        manager->m_receive_points.resize(comm_size);

        // From the fractional step element
        // const GeometryType& rGeom = this->GetGeometry();
        // Vector DetJ;
        // rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_2);
        // NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        // const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);
        //
        // rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2),false);
        //
        // for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
        //     rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();


        int i = 0;
        for (auto& condition : model_part.GetCommunicator().LocalMesh().Conditions()) {
            if(construction_type == 0) { // construct with Gauss Points
                // IntegrationPointsArrayType gauss_points_local_coords = condition.GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_1);
                //
                // for (auto& gp_local_coords : gauss_points_local_coords) {
                //     Node<3>::CoordinatesArrayType sth = condition.GetGeometry().GlobalCoordinates(sth, gp_local_coords);
                //     KRATOS_WATCH(sth.size())
                //     KRATOS_WATCH(sth[0])
                //     KRATOS_WATCH(sth[1])
                //     KRATOS_WATCH(sth[2])
                // }




                // Node<3>::CoordinatesArrayType sth = condition.GetGeometry().GlobalCoordinates(sth, static_cast<Node<3>::CoordinatesArrayType>(condition.GetGeometry().IntegrationPoints(GeometryData::GI_GAUSS_2)));
                // KRATOS_WATCH(arrayIP[0].size())
                // KRATOS_WATCH(arrayIP[0][0])
                // KRATOS_WATCH(arrayIP[0][1])
                // KRATOS_WATCH(arrayIP[0][2])


                const Geometry< Node<3> >& condition_geometry = condition.GetGeometry();

                Matrix shape_functions = condition_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

                const int num_gauss_points = shape_functions.size1();
                const int num_nodes = shape_functions.size2();

                array_1d<double, 3> gauss_point_global_coords;

                for (int g = 0; g < num_gauss_points; ++g) {
                    gauss_point_global_coords[0] = 0;
                    gauss_point_global_coords[1] = 0;
                    gauss_point_global_coords[2] = 0;

                    for (int n = 0; n < num_nodes; ++n) {
                        gauss_point_global_coords[0] += shape_functions(g, n) * condition_geometry[n].X();
                        gauss_point_global_coords[1] += shape_functions(g, n) * condition_geometry[n].Y();
                        gauss_point_global_coords[2] += shape_functions(g, n) * condition_geometry[n].Z();
                    }

                    manager->m_point_comm_objects[i] = InterfaceObject::Pointer( new InterfaceCondition(condition, gauss_point_global_coords) );
                    ++i;
                }
            } else if (construction_type == 1) { // construct with center of condition
                manager->m_point_comm_objects[i] = InterfaceObject::Pointer( new InterfaceCondition(condition, condition.GetGeometry().Center()) );
                ++i;

            // } else if (construction_type == 2) { // construct with nodal coordinates
            //     for (auto& node : condition.GetGeometry().Points()) {
            //         manager->m_point_comm_objects[i] = InterfaceObject::Pointer( new InterfaceCondition(condition, node.Coordinates()) );
            //         ++i;
            //     }
            } else {
                KRATOS_ERROR << "MappingApplication; InterfaceObjectManager; \"CreateInterfaceConditionManager\" unknown type to construct manager with!" << std::endl;
            }
        }
        std::cout << "Initialized Condition Manager" << std::endl;
        return manager;
    }

    void ComputeCandidatePartitions(int* local_comm_list, int* local_memory_size_array, double* global_bounding_boxes) {
        double* bounding_box[6];

        for (auto& point : m_point_comm_objects) {
            for (int partition_index = 0; partition_index < m_comm_size; ++partition_index) { // loop over partitions
                for (int j = 0; j < 6; ++j) { // retrieve bounding box of partition
                    bounding_box[j] = &global_bounding_boxes[(partition_index * 6) + j];
                }
                if (point->IsInBoundingBox(bounding_box)) {
                    m_candidate_send_points[partition_index].push_back(point);
                    local_comm_list[partition_index] = 1;
                    ++local_memory_size_array[partition_index];
                    point->SetIsBeingSent();
                }
            }
            // Robustness check, if point is sent to at least one partition
            if (!point->GetIsBeingSent()) {
                // Send Point to all Partitions
                for (int partition_index = 0; partition_index < m_comm_size; ++partition_index) {
                    std::cout << "MAPPER WARNING, Rank " << m_comm_rank << ", Point [ " <<
                    point->X() << " " << point->Y() << " " << point->Z() <<
                    " ] was not in a bounding box and is sent to all partitions!" << std::endl;

                    m_candidate_send_points[partition_index].push_back(point);
                    local_comm_list[partition_index] = 1;
                    ++local_memory_size_array[partition_index];
                }
            }
        }

    }

    void PrepareMatching(int* local_comm_list, int* local_memory_size_array) {
        for (int partition_index = 0; partition_index < m_comm_size; ++partition_index) {
            for (auto point : m_candidate_send_points[partition_index]) {
                // TODO m_matching_information[partition_index].reserve(m_candidate_send_points[partition_index].size()); ??? Do it
                if (point->HasNeighborInPartition(partition_index)) {
                    m_send_points[partition_index].push_back(point);
                    m_matching_information[partition_index].push_back(1);
                    local_comm_list[partition_index] = 1;
                    ++local_memory_size_array[partition_index];
                } else {
                    m_matching_information[partition_index].push_back(0);
                }
            }
        }
    }

    void FillSendBufferWithMatchInformation(int* send_buffer, int& send_buffer_size, const int comm_partner){
        int i = 0;
        for (auto info : m_matching_information[comm_partner]) {
            send_buffer[i] = info;
            ++i;
        }
        send_buffer_size = i;
    }

    InterfaceObjectConfigure::ContainerType& GetPointList() { // TODO is that correct?
        return m_point_comm_objects;
    }

    InterfaceObjectConfigure::ContainerType& GetPointListSerialSearch() {
        m_candidate_send_points[m_comm_rank] = m_point_comm_objects;
        return m_point_comm_objects;
    }

    void FillBufferLocalSearch(InterfaceObjectConfigure::ContainerType& send_points, int& num_points) {
        int i = 0;
        for (auto point : m_candidate_send_points[m_comm_rank]) {
            send_points[i] = point;
            ++i;
        }
        num_points = i;
    }

    void FillSendBufferRemoteSearch(double* send_buffer, int& send_buffer_size, const int comm_partner) {
        int i = 0;
        for (auto point : m_candidate_send_points[comm_partner]){
            send_buffer[(i*3) + 0] = point->X();
            send_buffer[(i*3) + 1] = point->Y();
            send_buffer[(i*3) + 2] = point->Z();
            ++i;
        }
        send_buffer_size = 3 * i;
    }

    template <typename T> // needed bcs "distances" can be "std::vector<double>" or "double*", local resp. remote search
    void PostProcessReceivedResults(const T& distances, const int comm_partner) {
        int i = 0;
        for (auto point : m_candidate_send_points[comm_partner]) {
            if (distances[i] > -0.5) // failed search has value "-1"
                point->ProcessDistance(distances[i], comm_partner);
            ++i;
        }
    }

    std::vector<InterfaceObject::Pointer> GetLocalMappingVector() {
        return m_send_points[m_comm_rank];
    }

    void CheckResults() {
        for (auto& point : m_point_comm_objects) {
            if (!point->NeighborFound()) {
                std::cout << "\tPoint has not found a neighbor, [ " << point->X() << " | " << point->Y() << " | " << point->Z() << " ]" << std::endl;
            }
        }
    }

// ********** Side that receives the points to find matches, aka Origin ***************************************

    void ProcessReceiveBuffer(InterfaceObjectConfigure::ContainerType& remote_p_point_list,
                              const double* coordinate_list, const int coordinate_list_size,
                              int& num_points){
        num_points = coordinate_list_size / 3;

        for (int i = 0; i < num_points; ++i) { // create Points
            remote_p_point_list[i] = InterfaceObject::Pointer(new InterfaceObject(
                coordinate_list[(i*3)+0], coordinate_list[(i*3)+1], coordinate_list[(i*3)+2]));
        }
    }

    void FillSendBufferWithResults(double* send_buffer, const int send_buffer_size, const std::vector<double>& min_distances){
        for (int i = 0; i < send_buffer_size; ++i)
            send_buffer[i] = min_distances[i];
    }

    void StoreTempSearchResults(const std::vector<InterfaceObject::Pointer> temp_closest_results, const int comm_partner) {
        m_candidate_receive_points[comm_partner] = temp_closest_results;
    }

    void ProcessMatchInformation(int* buffer, const int buffer_size, const int comm_partner){
        for (int i = 0; i < buffer_size; ++i) {
            if (buffer[i] == 1) { // Match
                if (m_candidate_receive_points[comm_partner][i] == nullptr)
                    KRATOS_ERROR << "Point Mismatch" << std::endl;
                m_receive_points[comm_partner].push_back(m_candidate_receive_points[comm_partner][i]);
            }
        }
    }

    std::vector<InterfaceObject::Pointer>& GetRemoteMappingVector() {
        return m_receive_points[m_comm_rank];
    }

    std::vector<InterfaceObject::Pointer>& GetPointListBins() {
        return m_candidate_receive_points[m_comm_rank];
    }

// ************ for Mapping ****************************************************
    void FillSendBufferWithValues(double* send_buffer, int& send_buffer_size, const int comm_partner,
                                  const Variable<double> variable) {
        int i = 0;
        std::vector<InterfaceObject::Pointer> point_list = m_receive_points[comm_partner];

        for (auto point : point_list) {
            send_buffer[i] = point->GetObjectValue(variable);
            ++i;
        }
        send_buffer_size = i;
    }

    void FillSendBufferWithValues(double* send_buffer, int& send_buffer_size, const int comm_partner,
                                  const Variable< array_1d<double,3> > variable) {
        int i = 0;
        std::vector<InterfaceObject::Pointer> point_list = m_receive_points[comm_partner];

        for (auto point : point_list) {
            send_buffer[(i*3) + 0] = point->GetObjectValue(variable)[0];
            send_buffer[(i*3) + 1] = point->GetObjectValue(variable)[1];
            send_buffer[(i*3) + 2] = point->GetObjectValue(variable)[2];
            ++i;
        }
        send_buffer_size = i * 3;
    }

    void ProcessValues(double* buffer, const int buffer_size, const int comm_partner,
                       const Variable<double> variable, const bool add_value) {

        std::vector<InterfaceObject::Pointer> point_list = m_send_points[comm_partner];

        if (static_cast<int>(point_list.size()) != buffer_size) {
            KRATOS_ERROR << "MappingApplication; InterfaceObjectManager; \"ProcessValues, double\": Wrong number of results received!;\
            point_list.size() = " << point_list.size() << ", buffer_size = " << buffer_size << std::endl;
        }

        for (int i = 0; i < buffer_size; ++i) {
            if (add_value) {
                point_list[i]->AddObjectValue(variable, buffer[i]);
            } else {
                point_list[i]->SetObjectValue(variable, buffer[i]);
            }
        }
    }

    void ProcessValues(double* buffer, const int buffer_size, const int comm_partner,
                       const Variable< array_1d<double,3> > variable, const bool add_value) {

        if (buffer_size % 3 != 0) {
            KRATOS_ERROR << "MappingApplication; InterfaceObjectManager; \"ProcessValues, double<3>\": Uneven number of results received!;\
            buffer_size modulo 3 = " << buffer_size % 3 << std::endl;
        }

        const int num_values = buffer_size / 3;

        std::vector<InterfaceObject::Pointer> point_list = m_send_points[comm_partner];

        if (static_cast<int>(point_list.size()) != num_values) {
           KRATOS_ERROR << "MappingApplication; InterfaceObjectManager; \"ProcessValues, double<3>\": Wrong number of results received!;\
           point_list.size() = " << point_list.size() << ", num_values = " << num_values << std::endl;
        }

        array_1d<double,3> value;

        for (int i = 0; i < num_values; ++i) {
            value[0] = buffer[(i*3) + 0];
            value[1] = buffer[(i*3) + 1];
            value[2] = buffer[(i*3) + 2];

            if (add_value) {
                point_list[i]->AddObjectValue(variable, value);
            } else {
                point_list[i]->SetObjectValue(variable, value);
            }
        }
    }

protected:

private:
    InterfaceObjectManager() {} // default constructor
    InterfaceObjectConfigure::ContainerType m_point_comm_objects; // TODO type? => implicit conversions

    double* m_bounding_box = new double[6] {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
                                          // xmax, xmin,  ymax, ymin,  zmax, zmin

    int m_comm_rank = 0;
    int m_comm_size = 0;
    // sending interface
    std::vector<std::vector<InterfaceObject::Pointer>> m_candidate_send_points; // TODO clear after matches were found?
    std::vector<std::vector<InterfaceObject::Pointer>> m_send_points;
    std::vector<std::vector<int>> m_matching_information;

    // receiving interface
    // TODO change to set or map/unordered_map, in case information is not needed from all partitions ?
    std::vector<std::vector<InterfaceObject::Pointer>> m_candidate_receive_points; // TODO clear after matches were found?
    std::vector<std::vector<InterfaceObject::Pointer>> m_receive_points;
};

///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_H_INCLUDED  defined
