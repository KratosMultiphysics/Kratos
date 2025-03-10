//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//
//

#ifndef KRATOS_CONTACT_ANGLE_EVALUATOR_H
#define KRATOS_CONTACT_ANGLE_EVALUATOR_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
//#include "containers/model.h"
//#include "includes/checks.h"
//#include "utilities/openmp_utils.h"
//#include "processes/find_nodal_h_process.h"
#include "utilities/variable_utils.h" //Now necessary!
//#include "processes/compute_nodal_gradient_process.h"
//#include "custom_utilities/element_size_calculator.h"
#include "includes/deprecated_variables.h" //For IS_STRUCTURED
#include "includes/global_pointer_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "droplet_dynamics_application_variables.h"

#define PI 3.14159265358979

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Helper process to record statistics on the integration points of the mesh.
class ContactAngleEvaluator: public Process
{
public:
///@name Type Definitions
///@{

/// Pointer definition of ContactAngleEvaluator
KRATOS_CLASS_POINTER_DEFINITION(ContactAngleEvaluator);

///@}
///@name Life Cycle
///@{

/// Constructor using a ModelPart.
/** @param rModelPart ModelPart the statistics will be calculated on.
 */
ContactAngleEvaluator(ModelPart& rModelPart):
    Process(),
    mrModelPart(rModelPart)
{
    KRATOS_INFO("HERE_CONSTRUCTOR");
}

ContactAngleEvaluator(ModelPart& rModelPart, Parameters &rParameters):
    Process(),
    mrModelPart(rModelPart)
{
    KRATOS_INFO("HERE_CONSTRUCTOR");
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

/// Destructor.
~ContactAngleEvaluator() override
{}

///@}
///@name Operations
///@{

void CheckDefaultsAndProcessSettings(Parameters &rParameters)
    {
        Parameters default_parameters(R"(
    {
        "theta_advancing" : 130,
        "theta_receding" : 130
    }  )");

        rParameters.ValidateAndAssignDefaults(default_parameters);

        theta_advancing = rParameters["theta_advancing"].GetDouble() * PI/180.0;
        theta_receding = rParameters["theta_receding"].GetDouble() * PI/180.0;
    }

void ExecuteInitialize() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void ExecuteBeforeSolutionLoop() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void ExecuteInitializeSolutionStep() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void ExecuteFinalizeSolutionStep() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void Execute() override
{
    KRATOS_TRY;
    KRATOS_INFO("ContactAngleEvaluatorProcess") << "Execute: Start." << std::endl;

    // const double theta_advancing = 130.0*PI/180.0;//180.0*PI/180.0;//149.0*PI/180.0;//129.78*PI/
    // const double theta_receding = 130.0*PI/180.0;//0.0*PI/180.0;//115.0*PI/180.0;//129.78*PI/

    const unsigned int num_nodes = mrModelPart.NumberOfNodes();
    const unsigned int num_elements = mrModelPart.NumberOfElements();

    // Auxiliar containers
    Vector N, distances, solid_normal, gradient;
    Matrix InvJ0, J0, DN_DX;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Initial resize
    const auto& r_first_element_geometry = it_element_begin->GetGeometry();
    const std::size_t number_of_nodes_first_element = r_first_element_geometry.PointsNumber();
    KRATOS_INFO("ContactAngleEvaluatorProcess") << "Number of nodes: " << num_nodes << ", Number of elements: " << num_elements << ", Number of nodes of first element: " << number_of_nodes_first_element << std::endl;
    const std::size_t local_space_dimension_first_element = r_first_element_geometry.LocalSpaceDimension();

    // Resize if needed
    DN_DX.resize(number_of_nodes_first_element, dimension);
    J0.resize(dimension, local_space_dimension_first_element);
    N.resize(number_of_nodes_first_element);
    distances.resize(number_of_nodes_first_element);
    solid_normal.resize(dimension);
    gradient.resize(dimension);

    // First node iterator
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for (unsigned int k = 0; k < num_nodes; ++k) {
        auto it_node = it_node_begin + k;
        it_node->Set(INTERFACE, false);
        it_node->SetValue(CURVATURE, 0.0);
    }

    // Create lookup table of structure nodes at beginning
    std::unordered_map<Kratos::Node::IndexType, bool> is_structure_looKupTable;

    for (Kratos::ModelPart::NodeIterator inode = mrModelPart.NodesBegin();
         inode != mrModelPart.NodesEnd(); ++inode)
    {
        if (inode->GetValue(IS_STRUCTURE))
            is_structure_looKupTable[inode->Id()] = true; // set to true if it is a structure node
        else
            is_structure_looKupTable[inode->Id()] = false; // set to false if it is not a structure node
    }

    // Iterate over the elements
    #pragma omp parallel for firstprivate(N, distances, DN_DX, J0, InvJ0, solid_normal, gradient)
    for(int i_elem=0; i_elem < num_elements; ++i_elem) {
        auto it_elem = it_element_begin + i_elem;

        it_elem->SetValue(CONTACT_ANGLE, 0.0);
        it_elem->SetValue(NORMAL_VECTOR, ZeroVector(3));

        auto& r_geometry = it_elem->GetGeometry();
        unsigned int n_contact_neg = 0;
        unsigned int n_contact_pos = 0;
        //unsigned int n_contact = 0;

        unsigned int n_neg = 0;
        unsigned int n_pos = 0;

        solid_normal = ZeroVector(dimension);

        for (std::size_t i_node = 0; i_node < number_of_nodes_first_element; ++i_node){
            auto& node = r_geometry[i_node];
            const double distance = node.FastGetSolutionStepValue(DISTANCE);
            distances[i_node] = distance;
            if (distance > 0.0){
                n_pos++;
                if (node.GetValue(IS_STRUCTURE) == 1.0){//if (is_structure_looKupTable[node.Id()]) {//
                    n_contact_pos++;
                    solid_normal += node.FastGetSolutionStepValue(NORMAL);
                    //n_contact++;
                }
            } else{
                n_neg++;
                if (node.GetValue(IS_STRUCTURE) == 1.0){//if (is_structure_looKupTable[node.Id()]) {//
                    n_contact_neg++;
                    solid_normal += node.FastGetSolutionStepValue(NORMAL);
                    //n_contact++;
                }
            }
        }

        if (n_pos > 0 && n_neg > 0){
            for (std::size_t i_node = 0; i_node < number_of_nodes_first_element; ++i_node){
                auto& node = r_geometry[i_node];
                node.Set(INTERFACE, true);
            }
        }

        if (n_contact_pos > 0 && n_contact_neg > 0){
            const double normal_norm = Kratos::norm_2(solid_normal);
            solid_normal = (1.0/normal_norm)*solid_normal;

            // // The integration points
            // const auto& r_integration_method = GeometryData::GI_GAUSS_1; //r_geometry.GetDefaultIntegrationMethod();
            // const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            // const std::size_t number_of_integration_points = r_integration_points.size();

            // double detJ0 = 0.0;

            // // The containers of the shape functions and the local gradients
            // // const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
            // const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

            // for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            //     // Getting the shape functions
            //     // noalias(N) = row(rNcontainer, point_number);

            //     // Getting the jacobians and local gradients
            //     GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
            //     MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
            //     const Matrix& rDN_De = rDN_DeContainer[point_number];
            //     GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

            //     gradient += prod(trans(DN_DX), distances);

            //     // const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;
            // }

            for (std::size_t i_node = 0; i_node < number_of_nodes_first_element; ++i_node){
                auto& node = r_geometry[i_node];
                if (node.GetValue(IS_STRUCTURE) == 1.0){//if (is_structure_looKupTable[node.Id()]) {//
                    gradient += node.FastGetSolutionStepValue(DISTANCE_GRADIENT);
                }
            }

            const double gradient_norm = Kratos::norm_2(gradient);
            gradient = (1.0/gradient_norm)*gradient;

            it_elem->GetValue(CONTACT_ANGLE) = PI - std::acos( std::max( std::min( inner_prod(solid_normal, gradient), 1.0 ), -1.0 ) );
            it_elem->GetValue(NORMAL_VECTOR) = gradient;
        }

    }

    #pragma omp parallel for
    for (unsigned int k = 0; k < num_nodes; ++k) {
        auto it_node = it_node_begin + k;

        //it_node->SetValue(CONTACT_ANGLE, 0.0);
        it_node->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
        it_node->SetValue(CONTACT_ANGLE, 0.0);
        it_node->FastGetSolutionStepValue(CONTACT_VELOCITY) = 0.0;
        it_node->Free(DISTANCE);

        double weight = 0.0;
        double avg_contact_angle = 0.0;
        Vector normal_avg = ZeroVector(3);

        auto& neighbour_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);

        for (auto i_element = neighbour_elements.begin(); i_element != neighbour_elements.end(); i_element++){
            const double elemental_contact_angle = i_element->GetValue(CONTACT_ANGLE);

            if (elemental_contact_angle > 1.0e-12){
                avg_contact_angle += elemental_contact_angle;
                normal_avg += i_element->GetValue(NORMAL_VECTOR);
                weight += 1.0;
            }
        }

        if (weight >= 1.0){
            const double contact_angle = avg_contact_angle/weight;
            it_node->FastGetSolutionStepValue(CONTACT_ANGLE) = contact_angle;
            it_node->GetValue(CONTACT_ANGLE) = contact_angle;
            const Vector normal = (1.0/norm_2(normal_avg)) * normal_avg;
            it_node->FastGetSolutionStepValue(NORMAL_VECTOR) = normal;

            //Vector velocity = it_node->FastGetSolutionStepValue(VELOCITY);
            //const double velocity_norm = Kratos::norm_2(velocity);
            //velocity = (1.0/velocity_norm)*velocity;

            const double distance_diff = it_node->FastGetSolutionStepValue(DISTANCE) - it_node->FastGetSolutionStepValue(DISTANCE_AUX);
            const int velocity_direction = (distance_diff < 0.0) - (distance_diff > 0.0);//inner_prod(velocity, normal);
            //it_node->FastGetSolutionStepValue(CONTACT_VELOCITY) = velocity_direction;

            if (it_node->GetValue(IS_STRUCTURE) == 1.0){//if (is_structure_looKupTable[it_node->Id()]) {//
                /* if (velocity_direction > 0.0 && contact_angle < theta_advancing){
                    it_node->Fix(DISTANCE);
                } else if (velocity_direction < 0.0 && contact_angle > theta_receding){
                    it_node->Fix(DISTANCE);
                } else if (velocity_direction == 0.0 && contact_angle < theta_advancing && contact_angle > theta_receding){
                    it_node->Fix(DISTANCE);
                } */

                if ( ( !(velocity_direction <= 0 && contact_angle <= theta_receding) &&
                    !(velocity_direction >= 0 && contact_angle >= theta_advancing) ) /* ||
                    (velocity_direction == 0 && contact_angle < theta_advancing && contact_angle > theta_receding) */ ){ // this last OR condition is unnecessary!
                    it_node->Fix(DISTANCE);
                    it_node->FastGetSolutionStepValue(CONTACT_VELOCITY) = static_cast<double>(velocity_direction);
                }

                /* if (contact_angle < theta_advancing && contact_angle > theta_receding){
                    it_node->Fix(DISTANCE);
                } */
            }

        }

        /* if (norm_2(it_node->FastGetSolutionStepValue(DISTANCE_GRADIENT)) < 0.5){ // Added to prevent a highly irregular level-set
            it_node->Fix(DISTANCE);
        } */
    }

    #pragma omp parallel for
    for (unsigned int i = 0; i < num_nodes; ++i) {

        auto it_node_i = it_node_begin + i;
        auto& node_i_contact_angle = it_node_i->FastGetSolutionStepValue(CONTACT_ANGLE);
        const double node_i_distance = it_node_i->FastGetSolutionStepValue(DISTANCE);
        const auto node_i_coordinates = it_node_i->Coordinates();
        double& node_i_curvature = it_node_i->GetValue(CURVATURE);

        // This part (CURVATURE estimation) should be merged with the next one (CONTACT_ANGLE estimation)
        if (!it_node_i->Is(INTERFACE)){

            double min_dist = 1.0e6;
            unsigned int min_dist_index;

            for (unsigned int j = 0; j < num_nodes; ++j) {
                auto it_node_j = it_node_begin + j;
                if (it_node_j->Is(INTERFACE)){
                    const double nodal_dist = norm_2(node_i_coordinates - it_node_j->Coordinates());
                    if (nodal_dist < min_dist){
                        min_dist = nodal_dist;
                        min_dist_index = j;
                    }
                }
            }
            const double node_j_curvature = (it_node_begin + min_dist_index)->FastGetSolutionStepValue(CURVATURE);
            const double radius_at_nodej = 2.0*node_j_curvature/(node_j_curvature*node_j_curvature + 1.0e-10);
            node_i_curvature = (std::abs(radius_at_nodej)/(std::abs(radius_at_nodej + node_i_distance) + 1.0e-10))*node_j_curvature;

        } else{
            // for interface nodes
            node_i_curvature = it_node_i->FastGetSolutionStepValue(CURVATURE);
        }


        // CONTACT_ANGLE estimation
        if (node_i_contact_angle == 0.0 && node_i_distance > 0.0){
            if (it_node_i->GetValue(IS_STRUCTURE) == 1.0 && it_node_i->Coordinates()[2] == 0.0){    //is_structure_looKupTable[it_node_i->Id()]// This won't work if droplet has no contact with z=0 plane.
                double min_horizontal_dist = 1.0e6;
                double min_dist_contact_angle = 0.0;
                unsigned int min_dist_index;

                for (unsigned int j = 0; j < num_nodes; ++j) {
                    auto it_node_j = it_node_begin + j;

                    const double node_j_contact_angle = it_node_j->GetValue(CONTACT_ANGLE);

                    if (node_j_contact_angle > 1.0e-10){
                        const double nodal_dist = norm_2(node_i_coordinates - it_node_j->Coordinates());
                        if (nodal_dist < min_horizontal_dist){
                            min_horizontal_dist = nodal_dist;
                            min_dist_contact_angle = node_j_contact_angle;
                            min_dist_index = j;
                        }
                    }

                }

                const double node_j_curvature = (it_node_begin + min_dist_index)->FastGetSolutionStepValue(CURVATURE);
                const double radius_at_nodej = 2.0*node_j_curvature/(node_j_curvature*node_j_curvature + 1.0e-10);
                node_i_contact_angle = std::asin( std::max( std::min( radius_at_nodej*std::sin(min_dist_contact_angle - PI/2.0)/(node_i_distance + radius_at_nodej), 1.0 ), -1.0 ) ) + PI/2.0; //min_dist_contact_angle;

            } else if (it_node_i->GetValue(IS_STRUCTURE) == 1.0 || it_node_i->Is(BOUNDARY)){ //is_structure_looKupTable[it_node_i->Id()]// By default contact angle is not set for the NON IS_STRUCTURE nodes

                double min_dist = 1.0e6;
                unsigned int min_dist_index;

                for (unsigned int j = 0; j < num_nodes; ++j) {
                    auto it_node_j = it_node_begin + j;
                    if (it_node_j->Is(INTERFACE)){
                        const double nodal_dist = norm_2(node_i_coordinates - it_node_j->Coordinates());
                        if (nodal_dist < min_dist){
                            min_dist = nodal_dist;
                            min_dist_index = j;
                        }
                    }
                }

                auto min_dist_iter_node = it_node_begin + min_dist_index;
                const double node_j_contact_angle = min_dist_iter_node->GetValue(CONTACT_ANGLE);
                auto min_dist_normal = min_dist_iter_node->FastGetSolutionStepValue(DISTANCE_GRADIENT);
                min_dist_normal /= norm_2(min_dist_normal);
                auto node_i_solid_normal = it_node_i->FastGetSolutionStepValue(NORMAL);
                node_i_solid_normal /= norm_2(node_i_solid_normal);

                if (node_j_contact_angle > 1.0e-10){
                    Vector min_dist_vector = min_dist_iter_node->Coordinates() - node_i_coordinates;
                    min_dist_vector(2) = 0.0;
                    const double min_horizontal_dist = norm_2(min_dist_vector);
                    const double node_j_curvature = min_dist_iter_node->FastGetSolutionStepValue(CURVATURE);
                    const double radius_at_nodej = 2.0*node_j_curvature/(node_j_curvature*node_j_curvature + 1.0e-10);
                    const double gamma = std::acos( std::max( std::min( (min_horizontal_dist + radius_at_nodej*std::sin(node_j_contact_angle))/(
                                        node_i_distance + radius_at_nodej), 1.0 ), -1.0 ) );

                    const double factor = (1.0 - std::sin(gamma)*std::sin(gamma))/(
                                    1.0 - std::sin(node_j_contact_angle)*std::sin(node_j_contact_angle) + 1.0e-10);

                    /* if(factor > 1.0e1){
                        KRATOS_WATCH(gamma)
                        KRATOS_WATCH(node_j_contact_angle)
                        KRATOS_WATCH(node_i_distance)
                        KRATOS_WATCH(min_dist_normal)
                        KRATOS_WATCH(radius_at_nodej)
                    } */

                    Vector node_i_normal = factor*min_dist_normal;
                    node_i_normal(2) = std::sin(gamma);

                    node_i_contact_angle = PI - std::acos( std::max( std::min( inner_prod(node_i_normal, node_i_solid_normal), 1.0 ), -1.0 ) );
                    /* KRATOS_WATCH(factor)
                    KRATOS_WATCH(node_i_contact_angle) */ // A revision is needed to better interpret the cause of nan (out of range argument of acos and asin above)

                } else{
                    node_i_contact_angle = PI - std::acos( std::max( std::min( inner_prod(min_dist_normal, node_i_solid_normal), 1.0 ), -1.0 ) );
                }

            }
        }
    }
    KRATOS_INFO("ContactAngleEvaluatorProcess") << "Execute: End." << std::endl;
    KRATOS_CATCH("");
}

///@}
///@name Input and output
///@{

/// Turn back information as a string.
std::string Info() const override
{
    std::stringstream buffer;
    buffer << "ContactAngleEvaluator";
    return buffer.str();
}

/// Print information about this object.
void PrintInfo(std::ostream &rOStream) const override { rOStream << "ContactAngleEvaluator"; }

/// Print object's data.
void PrintData(std::ostream &rOStream) const override {}

///@}

protected:

///@name Protected Operations
///@{

///@}

private:

///@name Member Variables
///@{

double theta_advancing;
double theta_receding;

ModelPart& mrModelPart;

///@}
///@name Private Operations
///@{

///@name Un accessible methods
///@{

/// Default constructor.
ContactAngleEvaluator() = delete;

/// Assignment operator.
ContactAngleEvaluator &operator=(ContactAngleEvaluator const &rOther) = delete;

/// Copy constructor.
ContactAngleEvaluator(ContactAngleEvaluator const &rOther) = delete;

///@}

}; // Class ContactAngleEvaluator

///@}

///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_CONTACT_ANGLE_EVALUATOR_H  defined