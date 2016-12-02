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

#if !defined(KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED

// System includes
#include <string>
#include <iostream>
// #include <algorithm>
// #include <vector>

// External includes
//@{KRATOS_EXTERNA_INCLUDES}
// #include "includes/kratos_flags.h"
// #include <boost/python.hpp>

// Project includes
#include "includes/define.h"
// #include "includes/kratos_flags.h"
// #include "includes/element.h"
// #include "includes/model_part.h"
// #include "geometries/geometry_data.h"

// #include "spaces/ublas_space.h"
// #include "linear_solvers/linear_solver.h"
// #include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
// #include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
// #include "solving_strategies/strategies/residualbased_linear_strategy.h"
// #include "elements/distance_calculation_element_simplex.h"

namespace Kratos {

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
	// typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	// typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
	// typedef typename SparseSpaceType::MatrixType SparseMatrixType;
	// typedef typename SparseSpaceType::VectorType VectorType;
  //
	// typedef array_1d<double,3> array_3d;
	// typedef Node < 3 > PointType;
	// typedef Node < 3 > ::Pointer PointTypePointer;
	// typedef std::vector<PointType::Pointer> PointVector;
	// typedef std::vector<PointType::Pointer>::iterator PointIterator;
	// typedef std::vector<double> DistanceVector;
	// typedef std::vector<double>::iterator DistanceIterator;
	// typedef ModelPart::ConditionsContainerType ConditionsArrayType;
  //
	// typedef Element BaseType;
	// typedef BaseType::GeometryType GeometryType;



///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

class InterfaceObject : public Point<3>
{

public:
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObject);

    // InterfaceObject(TObject& i_object);
    InterfaceObject() : Point<3>(0.0, 0.0, 0.0) {  } // Default Constructor

    InterfaceObject(Node<3>& i_node) : Point<3>(i_node) {  }    // constuct point from node
    InterfaceObject(array_1d<double, 3> coords) : Point<3>(coords) {  }    // constuct point from coordinate-array
    InterfaceObject(double X, double Y, double Z) : Point<3>(X, Y, Z) {  }  // constuct point from coordinates

    bool IsInBoundingBox(double* bounding_box[]){
        // xmax, xmin,  ymax, ymin,  zmax, zmin
        bool is_inside = false;

        if (this->X() < *bounding_box[0] && this->X() > *bounding_box[1]) { // check x-direction
            if (this->Y() < *bounding_box[2] && this->Y() > *bounding_box[3]) { // check y-direction
                if (this->Z() < *bounding_box[4] && this->Z() > *bounding_box[5]) { // check z-direction
                    is_inside = true;
                }
            }
        }
        return is_inside;
    }

    void ProcessDistance(double distance, int rank) {
        m_neighbor_found = true;
        if (distance < m_min_distance_neighbor) {
            m_min_distance_neighbor = distance;
            m_neighbor_rank = rank;
        }
    }

    bool NeighborFound() {
        return m_neighbor_found;
    }

    double MinDistance() {
        return m_min_distance_neighbor;
    }

    bool HasNeighborInPartition(const int partition_index) {
        bool return_value = false;
        if (m_neighbor_found) {
            if (m_neighbor_rank == partition_index)
                return_value = true;
        }
        return return_value;
    }

    void SetIsBeingSent() {
        m_is_being_sent = true;
    }

    bool GetIsBeingSent() {
        return m_is_being_sent;
    }

    virtual int GetObjectId() { // TODO pure virtual ??? => there are some complications
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"GetObjectId\" of the base class called!" << std::endl;
        return -1;
    };

    // These functions have to be duplicated because virtual templates are not possible in C++
    // Scalars
    virtual double GetObjectValue(const Variable<double>& variable) {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"GetObjectValue, double\" of the base class called!" << std::endl;
    }

    virtual void SetObjectValue(const Variable<double>& variable, const double value)  {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"SetObjectValue, double\" of the base class called!" << std::endl;
    }

    virtual void AddObjectValue(const Variable<double>& variable, const double value)  {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"AddObjectValue, double\" of the base class called!" << std::endl;
    }

    // Vectors
    virtual array_1d<double,3> GetObjectValue(const Variable< array_1d<double,3> >& variable) {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"GetObjectValue, double<3>\" of the base class called!" << std::endl;
    }

    virtual void SetObjectValue(const Variable< array_1d<double,3> >& variable, const array_1d<double,3>& value)  {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"SetObjectValue, double<3>\" of the base class called!" << std::endl;
    }

    virtual void AddObjectValue(const Variable< array_1d<double,3> >& variable, const array_1d<double,3>& value)  {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"AddObjectValue, double<3>\" of the base class called!" << std::endl;
    }

    virtual bool CheckResult(array_1d<double, 3> global_coords, double& min_distance, double distance) {
        KRATOS_ERROR << "MappingApplication; InterfaceObject; \"CheckResult\" of the base class called!" << std::endl;
        return false;
    }


    // virtual void SelectBestResult(const ResultContainerType& result_list, const std::vector<double>& distances,
    //                               const std::size_t& num_results, InterfaceObject::Pointer& vec_closest_results, double& closest_distance) const {
    //     KRATOS_ERROR << "MappingApplication; InterfaceObject; \"SelectBestResult\" of the base class called!" << std::endl;
    // };

protected:

private:
    double m_min_distance_neighbor = 1e10;
    bool m_neighbor_found = false;
    bool m_is_being_sent = false;
    int m_neighbor_rank = 0;
};

class InterfaceNode : public InterfaceObject
{
public:
    InterfaceNode(Node<3>& i_node) : InterfaceObject(i_node),  m_p_node(&i_node) {  }

    int GetObjectId() override {
        return m_p_node->Id();
    }

    // Scalars
    double GetObjectValue(const Variable<double>& variable) override {
        return m_p_node->FastGetSolutionStepValue(variable);
    }

    void SetObjectValue(const Variable<double>& variable, const double value) override {
        m_p_node->FastGetSolutionStepValue(variable) = value;
    }

    void AddObjectValue(const Variable<double>& variable, const double value) override {
        m_p_node->FastGetSolutionStepValue(variable) += value;
    }

    // Vectors
    array_1d<double,3> GetObjectValue(const Variable< array_1d<double,3> >& variable) override {
        return m_p_node->FastGetSolutionStepValue(variable);
    }

    void SetObjectValue(const Variable< array_1d<double,3> >& variable, const array_1d<double,3>& value) override {
        m_p_node->FastGetSolutionStepValue(variable) = value;
    }

    void AddObjectValue(const Variable< array_1d<double,3> >& variable, const array_1d<double,3>& value) override {
        m_p_node->FastGetSolutionStepValue(variable) += value;
    }


    bool CheckResult(array_1d<double, 3> global_coords, double& min_distance, double distance) override { // I am an object in the bins
        bool is_closer = false;

        if (distance < min_distance){
            min_distance = distance;
            is_closer = true;
        }

        return is_closer;
    }

protected:

private:
    Node<3>* m_p_node;
};

class InterfaceCondition : public InterfaceObject
{
public:
    InterfaceCondition(Condition& i_condition, array_1d<double, 3> coords) : InterfaceObject(coords), m_p_condition(&i_condition) {
        // InterfaceObject(1.0, 2.0, 3.0); // TODO these coords are the element center coords

        // in case some calculations have to be done in order to construct the base element (e.g. element center),
        // then first call the empty standard constructor, then do the calculations and populate
        // the base class stuff, although populating the derived class (i.e. "this") should also be ok

        // m_num_nodes = i_condition.GetGeometry().PointsNumber();
        m_geometry_family = i_condition.GetGeometry().GetGeometryFamily();
    }

    int GetObjectId() override {
        return m_p_condition->Id();
    }

    bool CheckResult(array_1d<double, 3> global_coords, double& min_distance, double distance) override { // I am an object in the bins
        bool is_inside;
        // TODO do sth with min_distance
        if (m_geometry_family == GeometryData::Kratos_Linear) { // I am a line condition
            is_inside = ProjectPointToLine(global_coords, distance);
        } else if (m_geometry_family == GeometryData::Kratos_Triangle) { // I am a triangular condition
            is_inside = ProjectPointToTriangle(global_coords, distance);
        } else if (m_geometry_family == GeometryData::Kratos_Quadrilateral) { // I am a quadrilateral condition
            is_inside = ProjectPointToQuaddrilateral(global_coords, distance);
        } else {
            KRATOS_ERROR << "MappingApplication; InterfaceCondition; \"ComputeBestResult\" wrong number of points in condition!" << std::endl;
        }

        return is_inside;
    }

protected:

private:
    Condition* m_p_condition;
    int m_geometry_family;

    bool ProjectPointToLine(array_1d<double, 3> global_coords, double& distance) {
        return false;
    }

    bool ProjectPointToTriangle(array_1d<double, 3> global_coords, double& distance) {
        return false;
    }

    bool ProjectPointToQuaddrilateral(array_1d<double, 3> global_coords, double& distance) {
        return false;
    }
};

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_INCLUDED_H_INCLUDED  defined
