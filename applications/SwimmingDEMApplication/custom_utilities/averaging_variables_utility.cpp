//   $Author: Joaquin Gonzalez-Usua

// System includes
#include <cmath>

// External includes

// Project includes
#include "averaging_variables_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

double AveragingVariablesUtility::AverageVariables(ModelPart& r_spheres_model_part,
                                            const double centroid_current_layer,
                                            const double layer_width,
                                            const double number_of_sublayers,
                                            const double plane_area,
                                            Vector& layer_averaged_velocity,
                                            Vector& layer_averaged_dv_dz,
                                            Matrix& layer_averaged_particle_stress,
                                            Matrix& layer_granular_temperature)
{

    layer_averaged_velocity = ZeroVector(3);
    layer_averaged_dv_dz = ZeroVector(3);
    layer_averaged_particle_stress = ZeroMatrix(3,3);
    layer_granular_temperature = ZeroMatrix(3,3);
    double layer_averaged_packing_density = 0.0;
    Vector first_sublayer_velocity, last_sublayer_velocity;
    for (unsigned int m = 1; m <= number_of_sublayers; ++m){
        double centroid_current_sublayer = centroid_current_layer - layer_width/2.0 + (m - 1.0) * layer_width/10.0;
        double total_contributive_area = 0.0;
        Vector sublayer_averaged_velocity = ZeroVector(3);
        BoundedMatrix<double, 3, 3> sublayer_averaged_particle_stress = ZeroMatrix(3,3);
        double weight = std::cos(1.0/5.0 * Globals::Pi * (m - 1.0) - Globals::Pi);
        #pragma omp for schedule(guided, 1024)
        for (int i = 0; i < (int)r_spheres_model_part.Nodes().size(); ++i){
            NodeIteratorType i_node = r_spheres_model_part.NodesBegin() + i;
            double distance_centroids = std::abs(i_node->Y() - centroid_current_sublayer);
            double particle_radius = i_node->GetSolutionStepValue(RADIUS);
            array_1d<double,3> particle_velocity = i_node->GetSolutionStepValue(VELOCITY);
            double particle_contributive_area;
            if (distance_centroids <= particle_radius){
                particle_contributive_area = Globals::Pi * (std::pow(particle_radius, 2) - std::pow(distance_centroids, 2));
                sublayer_averaged_velocity += particle_contributive_area * particle_velocity;
                sublayer_averaged_particle_stress += particle_contributive_area * i_node->GetSolutionStepValue(DEM_STRESS_TENSOR);
                total_contributive_area += particle_contributive_area;
            }
        }
        if (total_contributive_area < std::numeric_limits<double>::epsilon()){
            total_contributive_area = std::numeric_limits<double>::epsilon();
        }
        sublayer_averaged_velocity *= 1.0/total_contributive_area;
        if (m == 1)
            first_sublayer_velocity = sublayer_averaged_velocity;
        else if (m == number_of_sublayers)
            last_sublayer_velocity = sublayer_averaged_velocity;

        double sublayer_averaged_packing_density = total_contributive_area/plane_area;
        sublayer_averaged_particle_stress *= 1.0/plane_area;
        layer_averaged_velocity -= weight * sublayer_averaged_velocity;
        layer_averaged_particle_stress -= weight * sublayer_averaged_particle_stress;
        layer_averaged_packing_density -= weight * sublayer_averaged_packing_density;
        Matrix layer_velocity_variance = ZeroMatrix(3,3);
        CalculateVelocityVariance(r_spheres_model_part, layer_velocity_variance, sublayer_averaged_velocity, centroid_current_sublayer);
        layer_granular_temperature -= weight * layer_velocity_variance;
    }

    for (unsigned int d = 0; d < 3; ++d)
            layer_averaged_dv_dz[d] += (last_sublayer_velocity[d] - first_sublayer_velocity[d])/layer_width;
    return layer_averaged_packing_density;
}

void AveragingVariablesUtility::CalculateVelocityVariance(ModelPart& r_spheres_model_part,
                                                          Matrix& layer_velocity_variance,
                                                          Vector sublayer_averaged_velocity,
                                                          double& centroid_current_sublayer)
{
    double total_contributive_area = 0.0;
    #pragma omp for schedule(guided, 1024)
    for (int i = 0; i < (int)r_spheres_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_spheres_model_part.NodesBegin() + i;
        double distance_centroids = std::abs(i_node->Y() - centroid_current_sublayer);
        double particle_radius = i_node->GetSolutionStepValue(RADIUS);
        array_1d<double,3> particle_velocity = i_node->GetSolutionStepValue(VELOCITY);
        double particle_contributive_area = 0.0;
        if (distance_centroids <= particle_radius){
            particle_contributive_area = Globals::Pi * (std::pow(particle_radius, 2) - std::pow(distance_centroids, 2));
            array_1d<double,3> distance_to_mean_velocity = particle_velocity - sublayer_averaged_velocity;
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 3; ++j)
                    layer_velocity_variance(i,j) += particle_contributive_area * distance_to_mean_velocity[i] * distance_to_mean_velocity[j];
            total_contributive_area += particle_contributive_area;
        }
    }
    layer_velocity_variance *= 1.0 / total_contributive_area;
}

/// Turn back information as a string.
std::string AveragingVariablesUtility::Info() const {
        return "AveragingVariablesUtility";
}

/// Print information about this object.
void AveragingVariablesUtility::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AveragingVariablesUtility::PrintData(std::ostream& rOStream) const {}


} // namespace Kratos