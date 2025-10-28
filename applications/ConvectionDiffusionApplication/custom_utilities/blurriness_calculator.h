//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala Pascual
//
//

#ifndef KRATOS_BLURRINESS_CALCULATOR
#define KRATOS_BLURRINESS_CALCULATOR

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

#include "includes/node.h"
#include "utilities/math_utils.h"
#include "spatial_containers/spatial_containers.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
    ///@addtogroup ConvectionDiffusionApplication
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

    class BlurrinessCalculator
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of BlurrinessCalculator
        KRATOS_CLASS_POINTER_DEFINITION(BlurrinessCalculator);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Constructor.
        BlurrinessCalculator(
            ModelPart &model_part,
            std::vector<ModelPart *> model_part_vector,
            const std::vector<double> &interfaces_positions,
            const std::vector<double> &interfaces_density_values,
            const double& u0,
            const double& r0) : mMainModelPart(model_part), mModelPartVector(model_part_vector), mIsBlurrinessComputed(false), mU0(u0), mR0(r0)
        {
            // Check that the dimensions are correct
            KRATOS_ERROR_IF(interfaces_positions.size() != interfaces_density_values.size() + 1)
                << "The number of interfaces positions in BlurrinessCalculator do not match the number"
                << "of density values (note that the interfaces positions must include the positions of the walls)" << std::endl;

            // Assign the values for the interfaces positions
            mInterfacesPositions.resize(interfaces_positions.size());
            for (unsigned n = 0; n < interfaces_positions.size() - 1; n++)
            {
                mInterfacesPositions[n] = interfaces_positions[n];
                mInterfacesPositions[n + 1] = interfaces_positions[n + 1];
                KRATOS_ERROR_IF(mInterfacesPositions[n] >= mInterfacesPositions[n + 1])
                    << "The interfaces positions must be sorted in ascending order" << std::endl;
            }

            // Assign the values for the interfaces density values
            mLayersDensityValues.resize(interfaces_positions.size());
            for (unsigned n = 0; n < interfaces_density_values.size(); n++)
            {
                mLayersDensityValues[n] = interfaces_density_values[n];
            }

            // Resize the blurriness matrix
            mNumLayers = mInterfacesPositions.size() - 1; // The interfaces must include the bottom and upper walls!
            mBlurriness.resize(mModelPartVector.size());  // Each model part has their own blurriness
            for (unsigned m = 0; m < mBlurriness.size(); m++)
            {
                mBlurriness[m].resize(mNumLayers - 1); // The blurriness is defined at each interface
            }

            // Check that the main mdel part has the error variable used to compute the blurriness
            // KRATOS_ERROR_IF_NOT(mMainModelPart.Has(STEP_SOLUTION_ERROR))
            // << "The model part " << mMainModelPart.Name() << " does not have the error variable STEP_SOLUTION_ERROR" << std::endl;
        }

        BlurrinessCalculator(
            ModelPart &model_part,
            std::vector<ModelPart *> model_part_vector,
            const std::vector<double> &interfaces_positions,
            const std::vector<double> &interfaces_density_values): BlurrinessCalculator(model_part, model_part_vector, interfaces_positions, interfaces_density_values, 1.0, 1.0){}

        /// Turn back information as a string.
        std::string Info() const
        {
            std::stringstream buffer;
            buffer << "BlurrinessCalculator";
            return buffer.str();
        }

        /// Print information about this object.
        void PrintInfo(std::ostream &rOStream) const { rOStream << "BlurrinessCalculator"; }

        /// Print object's data.
        void PrintData(std::ostream &rOStream) const {}

        /// @brief Compute the blurriness of each interface
        void ComputeBlurriness()
        {
            double scale_factor = mR0 * mR0 * mU0;
            for (unsigned m = 0; m < mModelPartVector.size(); m++)
            {
                std::vector<double> layer_integrals(mNumLayers, 0.); // Defined as \int{ (rho_h - rho_step)^2 dQ } , dQ = u * dS
                std::vector<double> layer_flows(mNumLayers, 0.);     // Velocity flows of each layer

                std::vector<double> areas(mNumLayers, 0.0);
                std::vector<double> element_areas(mNumLayers, 0.0);
                std::vector<unsigned int> n_elements_layers(mNumLayers, 0.0);
                unsigned int n_elements = mModelPartVector[m]->NumberOfConditions();

                ModelPart &r_model_part = *(mModelPartVector[m]);
                const unsigned number_of_conditions = r_model_part.NumberOfConditions();
#pragma omp parallel for schedule(dynamic)
                for (unsigned c = 0; c < number_of_conditions; c++)
                {
                    ModelPart::ConditionsContainerType::iterator it_cond = r_model_part.ConditionsBegin() + c;

                    // Usual chunk of code necessary to integrate on an element
                    Geometry<Node> &r_geometry = it_cond->GetGeometry();
                    unsigned int NumNodes = r_geometry.size();
                    GeometryData::IntegrationMethod integration_method = it_cond->GetIntegrationMethod();
                    const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints(integration_method);
                    unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
                    Vector detJ_vector(r_number_integration_points);
                    r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
                    Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);

                    // // Check at which layer this condition belongs to
                    double condition_center_z = r_geometry.Center()[2];
                    unsigned ith_interface = mNumLayers; // Initialize it with an impossible value
                    for (unsigned i = 0; i < mNumLayers; i++)
                    {
                        if (condition_center_z <= mInterfacesPositions[i + 1])
                        {
                            ith_interface = i;
                            break;
                        }
                    }
                    KRATOS_ERROR_IF(ith_interface == mNumLayers)
                        << "Unable to find layer for condition with Id = " << it_cond->Id() << " with center at z = " << condition_center_z << std::endl;
                    element_areas[ith_interface] += r_geometry.Area();
                    n_elements_layers[ith_interface] += 1;

                    // Perform the layer integral
                    for (unsigned g = 0; g < r_number_integration_points; g++)
                    {
                        // Global position of this gauss point
                        array_1d<double, 3> gauss_point_global = ZeroVector(3);
                        for (unsigned n = 0; n < NumNodes; n++)
                        {
                            Point node_global_pos = r_geometry.GetPoint(n);
                            for (unsigned d = 0; d < 3; d++)
                            {
                                gauss_point_global[d] += node_global_pos[d] * NContainer(g, n);
                            }
                        }

                        // Check at which layer this condition belongs to
                        // double gauss_coord_z = gauss_point_global[2];
                        // unsigned ith_interface = mNumLayers; // Initialize it with an impossible value
                        // for (unsigned i = 0; i < mNumLayers; i++)
                        // {
                        //     if (gauss_coord_z <= mInterfacesPositions[i + 1])
                        //     {
                        //         ith_interface = i;
                        //         break;
                        //     }
                        // }
                        // KRATOS_ERROR_IF(ith_interface == mNumLayers)
                        //     << "Unable to find layer for gauss point with z = " << gauss_coord_z << std::endl;

                        // Find the element containing this gauss point
                        array_1d<double, 3> gauss_point_local;
                        ModelPart::ElementType::Pointer p_elem;
                        FindElementContainingPoint(gauss_point_global, gauss_point_local, p_elem);
                        KRATOS_ERROR_IF(p_elem == nullptr)
                            << "Element with Id = " << p_elem->Id() << " not found in model part " << mMainModelPart.Name() << std::endl;

                        // Interpolate the value of the blurriness at the gauss points
                        array_1d<double, 3> normal_vec = r_geometry.Normal(gauss_point_local);
                        double norm = 0.0;
                        for (unsigned d = 0; d < 3; d++)
                            norm += normal_vec[d] * normal_vec[d];
                        norm = std::sqrt(norm);
                        for (unsigned d = 0; d < 3; d++)
                            normal_vec[d] /= norm;

                        double step_sol_error, normal_vel;
                        InterpolateAtPosition(p_elem, gauss_point_local, normal_vec, step_sol_error, normal_vel);
                        double err_abs = abs(step_sol_error);

                        double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
                        layer_flows[ith_interface] += Weight * normal_vel / scale_factor;
                        layer_integrals[ith_interface] += Weight * step_sol_error * step_sol_error * normal_vel / scale_factor;
                        areas[ith_interface] += Weight;

                        // Divide by its scale factor
                        // layer_flows[ith_interface] /= scale_factor;
                        // layer_integrals[ith_interface] /= scale_factor;
                    }
                }

                // Compute blurriness for this model part
                // double exact_area = M_PI * pow(.5e-3, 2);
                // double error = abs(exact_area - area) / exact_area;

                std::cout << "\nModel part " << m + 1 << std::endl;
                for (unsigned i = 0; i < mNumLayers; i++)
                {
                    layer_integrals[i] = abs(layer_integrals[i]);
                    layer_flows[i] = abs(layer_flows[i]);

                    // Print info
                    std::cout << "Layer " << i << ":" << std::endl;
                    std::cout << "  - ui = " << mLayersDensityValues[i] << std::endl;
                    std::cout << "  - qi = " << layer_flows[i] << std::endl;
                    std::cout << "  - layer integral = " << layer_integrals[i] << std::endl;
                    std::cout << "  - area_i = " << areas[i] << std::endl;
                    std::cout << "  - element_areas = " << element_areas[i] << std::endl;
                    std::cout << "  - n_elements_layers = " << n_elements_layers[i] << " / " << n_elements << std::endl;
                }

                // Compute total velocity flux and normalize the surface integral values
                double total_flow = 0.0;
                for (unsigned i = 0; i < mNumLayers; i++)
                {
                    total_flow += layer_flows[i];
                    layer_integrals[i] = layer_integrals[i] / layer_flows[i];
                }

                // Compute asymptotic solution
                double asymptotic_sol_value = 0.0;
                for (unsigned i = 0; i < mNumLayers; i++)
                {
                    asymptotic_sol_value += mLayersDensityValues[i] * layer_flows[i] / total_flow;
                }
                std::cout << "  qt = " << total_flow << std::endl;
                std::cout << "  u_inf = " << asymptotic_sol_value << std::endl;

                // Set blurriness to 0
                for (unsigned i = 0; i < mNumLayers; i++)
                {
                    mBlurriness[m][i] = 0.0;
                }

                // Compute blurriness
                for (unsigned i = 0; i < mNumLayers - 1; i++)
                {
                    // Norm factor
                    double density_diff_1 = asymptotic_sol_value - mLayersDensityValues[i];
                    double density_diff_2 = asymptotic_sol_value - mLayersDensityValues[i + 1];
                    double norm_factor = density_diff_1 * density_diff_1 + density_diff_2 * density_diff_2;

                    // Layer integrals
                    mBlurriness[m][i] = (layer_integrals[i] + layer_integrals[i + 1]) / norm_factor;
                    mBlurriness[m][i] = std::sqrt(mBlurriness[m][i]);

                    // Print info
                    std::cout << "Interface " << i << ":" << std::endl;
                    std::cout << "  - norm_factor = " << norm_factor << std::endl;
                    std::cout << "  - num / qi = " << layer_integrals[i] + layer_integrals[i + 1] << std::endl;
                    std::cout << "  - blurriness^2 = " << mBlurriness[m][i] * mBlurriness[m][i] << std::endl;
                    std::cout << "  - blurriness = " << mBlurriness[m][i] << std::endl;
                }

            }
            mIsBlurrinessComputed = true;
        }

        /// @brief Return the blurriness at each interface
        /// @return Blurriness
        std::vector<std::vector<double>> GetBlurriness()
        {
            if (!mIsBlurrinessComputed)
            {
                KRATOS_ERROR << "The blurriness is not defined. Method `ComputeBlurriness` must be called before getting the blurriness." << std::endl;
            }
            return mBlurriness;
        }

    private:
        // Member variables
        ModelPart &mMainModelPart;
        std::vector<ModelPart *> mModelPartVector;
        std::vector<double> mInterfacesPositions;
        std::vector<double> mLayersDensityValues;
        std::vector<std::vector<double>> mBlurriness;
        unsigned mNumLayers;
        bool mIsBlurrinessComputed;
        const double mU0; // Velocity scale
        const double mR0; // Length scale

        /// Default constructor.
        BlurrinessCalculator() = delete;

        /// Assignment operator.
        BlurrinessCalculator &operator=(BlurrinessCalculator const &rOther) = delete;

        /// Copy constructor.
        BlurrinessCalculator(BlurrinessCalculator const &rOther) = delete;

        /// @brief Find the element that contains the point
        /// @param point Point inside the element we want to find
        /// @param p_elem Pointer to the element containing the point
        void FindElementContainingPoint(const CoordinateVector &point, array_1d<double, 3> &p_pos_local, ModelPart::ElementType::Pointer &p_elem)
        {
            p_elem = nullptr;
            const unsigned number_of_elements = mMainModelPart.NumberOfElements();
            for (unsigned int e = 0; e < number_of_elements; e++)
            {
                ModelPart::ElementsContainerType::iterator it_elem = mMainModelPart.ElementsBegin() + e;
                Geometry<Node> &r_geometry = it_elem->GetGeometry();
                if (r_geometry.IsInside(point, p_pos_local, 1e-10))
                {
                    p_elem = mMainModelPart.pGetElement(it_elem->Id());
                    break;
                }
            }
        }

        /// @brief Interpolate the step solution error and the normal velocity at global_pos
        /// @param p_elem The element of the main model part
        /// @param global_pos The point in which the variables will be interpolated
        /// @param normal_vec Normal vector of the condition
        /// @param step_function_error Interpolated value of the step solution error, i.e. rho_h - rho_step
        /// @param normal_velocity Interpolated value of the normal velocity, i.e. \vec{u} \cdot \hat{n}, where \hat{n} is the normal vector of the condition
        void InterpolateAtPosition(Element::Pointer p_elem, const array_1d<double, 3> &p_pos_local, const array_1d<double, 3> &normal_vec, double &step_func_error, double &normal_velocity)
        {
            step_func_error = 0.0;
            normal_velocity = 0.0;

            Geometry<Node> &r_geometry = p_elem->GetGeometry();
            unsigned int NumNodes = r_geometry.size();

            for (unsigned n = 0; n < NumNodes; n++)
            {
                double nodal_step_func_error = r_geometry[n].FastGetSolutionStepValue(STEP_SOLUTION_ERROR); // rho_h - rho_step
                double nodal_normal_velocity = 0.0;
                array_1d<double, 3> nodal_vel = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
                for (unsigned d = 0; d < 3; d++)
                {
                    nodal_normal_velocity += (nodal_vel[d]) * normal_vec[d];
                }
                double shape_function_value = r_geometry.ShapeFunctionValue(n, p_pos_local);
                step_func_error += nodal_step_func_error * shape_function_value;
                normal_velocity += nodal_normal_velocity * shape_function_value;
            }
        }
    }; // Class BlurrinessCalculator

}; // namespace Kratos.

#endif // KRATOS_BLURRINESS_CALCULATOR