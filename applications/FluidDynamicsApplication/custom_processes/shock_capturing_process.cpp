//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "processes/calculate_nodal_area_process.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "shock_capturing_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void ShockCapturingProcess::Execute()
    {
        Check();
        ExecuteInitialize();
        CalculatePhysicsBasedShockCapturing();
    }

    void ShockCapturingProcess::ExecuteInitialize()
    {
        // Initialize nodal values
        for (auto& r_node : mrModelPart.Nodes()) {
            r_node.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
            r_node.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);
            r_node.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);
        }

        // Initialize elemental values
        for (auto& r_elem : mrModelPart.Elements()) {
            r_elem.SetValue(SHOCK_SENSOR, 0.0);
            r_elem.SetValue(SHEAR_SENSOR, 0.0);
            r_elem.SetValue(THERMAL_SENSOR, 0.0);
            r_elem.SetValue(ARTIFICIAL_CONDUCTIVITY, 0.0);
            r_elem.SetValue(ARTIFICIAL_BULK_VISCOSITY, 0.0);
            r_elem.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, 0.0);
        }

        // Calculate the NODAL_AREA
        CalculateNodalAreaProcess<false> nodal_area_process(mrModelPart);
        nodal_area_process.Execute();
    }

    void ShockCapturingProcess::ExecuteFinalizeSolutionStep()
    {
        CalculatePhysicsBasedShockCapturing();
    }

    const Parameters ShockCapturingProcess::GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "calculate_nodal_area_at_each_step" : false,
            "shock_sensor" : true,
            "shear_sensor" : true,
            "thermal_sensor" : true
        })");

        return default_parameters;
    }

    int ShockCapturingProcess::Check()
    {
        // Base process check
        int err_code = BaseType::Check();

        // Check that the required variables are in the nodal data
        for (const auto& rNode : mrModelPart.Nodes()) {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DENSITY, rNode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, rNode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TEMPERATURE, rNode)
        }

        // Check that the required material properties are in the elemental properties
        for (const auto& rElement : mrModelPart.Elements()) {
            const auto& r_prop = rElement.GetProperties();
            KRATOS_ERROR_IF_NOT(r_prop.Has(SPECIFIC_HEAT)) << "Element " << rElement.Id() << " properties " << r_prop.Id() <<  " has no SPECIFIC_HEAT." << std::endl;
            KRATOS_ERROR_IF_NOT(r_prop.Has(HEAT_CAPACITY_RATIO)) << "Element " << rElement.Id() << " properties " << r_prop.Id() << " has no HEAT_CAPACITY_RATIO." << std::endl;
        }

        return err_code;
    }

    /* Private functions *****************************************************/

    void ShockCapturingProcess::ValidateAndAssignParameters(Parameters& rParameters)
    {
        // Validate and assign defaults
        rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Assign settings
        mUpdateNodalArea = rParameters["calculate_nodal_area_at_each_step"].GetBool();
        mShockSensor = rParameters["shock_sensor"].GetBool();
        mShearSensor = rParameters["shear_sensor"].GetBool();
        mThermalSensor = rParameters["thermal_sensor"].GetBool();
    }

    ShockCapturingProcess::ElementMetricFunctionType ShockCapturingProcess::SetElementMetricFunction()
    {
        // Get geometry type
        ElementMetricFunctionType elem_metric_function;
        const auto geometry_type = (mrModelPart.ElementsBegin()->GetGeometry()).GetGeometryType();

        // Set the corresponding element metric tensor function
        switch (geometry_type) {
            case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                elem_metric_function = [&](const Geometry<Node<3>> &rGeometry) { return CalculateTriangleMetricTensor(rGeometry); };
                break;
            case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                elem_metric_function = [&](const Geometry<Node<3>> &rGeometry) { return CalculateTetrahedraMetricTensor(rGeometry); };
                break;
            default:
                KRATOS_ERROR << "Asking for a non-implemented geometry.";
        }

        return elem_metric_function;
    }

    /**
     * @brief Physics-based shock capturing
     * This function calculates the artificial magnitudes using a physics-based shock capturing method.
     * References https://arc.aiaa.org/doi/abs/10.2514/6.2018-0062
     */
    void ShockCapturingProcess::CalculatePhysicsBasedShockCapturing()
    {
        // Calculate the model part data
        const auto r_process_info = mrModelPart.GetProcessInfo();
        const int n_nodes = mrModelPart.NumberOfNodes();

        // Initialize the values to zero
        block_for_each(mrModelPart.Nodes(), [](Node<3> &rNode) {
            rNode.GetValue(ARTIFICIAL_CONDUCTIVITY) = 0.0;
            rNode.GetValue(ARTIFICIAL_BULK_VISCOSITY) = 0.0;
            rNode.GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) = 0.0;
        });

        // If required, update the NODAL_AREA
        if (mUpdateNodalArea) {
            CalculateNodalAreaProcess<false> nodal_area_process(mrModelPart);
            nodal_area_process.Execute();
        }

        // Set the functor to calculate the element size
        // Note that this assumes a unique geometry in the computational mesh
        auto elem_metric_function = SetElementMetricFunction();

        // Loop the elements to project the gradients
        // Note that it is assumed that the gradient is constant within the element
        // Hence, only one Gauss point is used
        const double eps = 1.0e-7;

        // Calculate the elemental contributions of the shock capturing
        const auto geometry_type = (mrModelPart.ElementsBegin()->GetGeometry()).GetGeometryType();
        if (geometry_type = GeometryData::KratosGeometryType::Kratos_Triangle2D3) {
            // Set auxiliary TLS container and elemental function
            ShockCapturingTLSType2D3N tls_container_2D3N;
            auto aux_function_2D3N = [&, this] (Element &rElement, ShockCapturingTLSType2D &rShockCapturingTLS) {this->CalculatePhysicsBasedShockCapturingElementContribution(rElement, rShockCapturingTLS);};
            // Perform the elemental loop
            block_for_each(mrModelPart.Elements(), tls_container_2D3N, aux_function_2D3N);
        } else if (geometry_type == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4) {
            // Set auxiliary TLS container and elemental function
            ShockCapturingTLSType3D4N tls_container_3D4N;
            auto aux_function_3D4N = [&, this] (Element &rElement, ShockCapturingTLSType3D &rShockCapturingTLS) {this->CalculatePhysicsBasedShockCapturingElementContribution(rElement, rShockCapturingTLS);};
            // Perform the elemental loop
            block_for_each(mrModelPart.Elements(), tls_container_3D4N, aux_function_3D4N);
        } else {
            KRATOS_ERROR << "Asking for a non-supported geometry. Physics-based shock capturing only supports \'Triangle2D3\' and \'Tetrahedra3D4\' geometries.";
        }

        // Nodal smoothing of the shock capturing magnitudes
        // Note that to avoid calculating the NODAL_AREA we took from the first DOF of the lumped mass vector
        IndexPartition<>(n_nodes).for_each([&](const ModelPart::SizeType iNode) {
            auto it_node = mrModelPart.NodesBegin() + iNode;
            const double nodal_area = it_node->GetValue(NODAL_AREA);
            it_node->GetValue(ARTIFICIAL_CONDUCTIVITY) /= nodal_area;
            it_node->GetValue(ARTIFICIAL_BULK_VISCOSITY) /= nodal_area;
            it_node->GetValue(ARTIFICIAL_DYNAMIC_VISCOSITY) /= nodal_area;
        });
    }

    // TODO: REMOVE AFTER SUNETH'S PR
    template <>
    void ShockCapturingProcess::UpdateValue(double &rOutput, const double &rInput)
    {
        rOutput += rInput;
    }

    // TODO: REMOVE AFTER SUNETH'S PR
    template <>
    void ShockCapturingProcess::UpdateValue(array_1d<double,3> &rOutput, const array_1d<double,3> &rInput)
    {
        noalias(rOutput) += rInput;
    }

    double ShockCapturingProcess::LimitingFunction(
        const double s,
        const double s_0,
        const double s_max,
        const double s_min)
    {
        const double aux_1 = std::max(s - s_0, s_min);
        const double aux_2 = std::min(aux_1 - s_max, s_min);
        return aux_2 + s_max;
    }

    double ShockCapturingProcess::SmoothedLimitingFunction(
        const double s,
        const double s_0,
        const double s_max)
    {
        const double aux_1 = SmoothedMaxFunction(s - s_0);
        const double aux_2 = SmoothedMinFunction(aux_1 - s_max);
        return aux_2 + s_max;
    }

    // Smooth approximation of the max(s,0) function
    double ShockCapturingProcess::SmoothedMaxFunction(const double s)
    {
        const double b = 100;
        const double l_max = (s / Globals::Pi) * std::atan(b * s) + 0.5 * s - (1.0 / Globals::Pi) * std::atan(b) + 0.5;
        return l_max;
    }

    // Smooth approximation of the min(s,0) function
    double ShockCapturingProcess::SmoothedMinFunction(const double s)
    {
        return s - SmoothedMaxFunction(s);
    }

    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner
    std::tuple<double, double, Matrix> ShockCapturingProcess::CalculateTriangleMetricTensor(const Geometry<Node<3>> &rGeometry)
    {
        const array_1d<double, 3> p_1 = rGeometry[0].Coordinates();
        const array_1d<double, 3> p_2 = rGeometry[1].Coordinates();
        const array_1d<double, 3> p_3 = rGeometry[2].Coordinates();

        // Solve the metric problem trans(e)*M*e = 1
        // This means, find the coefficients of the matrix M such that all the edges have unit length
        Vector sol;
        array_1d<double, 3> aux_vect;
        BoundedMatrix<double, 3, 3> aux_mat;
        aux_mat(0, 0) = std::pow(p_1[0] - p_2[0], 2);
        aux_mat(0, 1) = 2.0 * (p_1[0] - p_2[0]) * (p_1[1] - p_2[1]);
        aux_mat(0, 2) = std::pow(p_1[1] - p_2[1], 2);
        aux_mat(1, 0) = std::pow(p_1[0] - p_3[0], 2);
        aux_mat(1, 1) = 2.0 * (p_1[0] - p_3[0]) * (p_1[1] - p_3[1]);
        aux_mat(1, 2) = std::pow(p_1[1] - p_3[1], 2);
        aux_mat(2, 0) = std::pow(p_2[0] - p_3[0], 2);
        aux_mat(2, 1) = 2.0 * (p_2[0] - p_3[0]) * (p_2[1] - p_3[1]);
        aux_mat(2, 2) = std::pow(p_2[1] - p_3[1], 2);
        aux_vect[0] = 1.0;
        aux_vect[1] = 1.0;
        aux_vect[2] = 1.0;
        MathUtils<double>::Solve(aux_mat, sol, aux_vect);

        // Set the metric tensor
        Matrix metric(2, 2);
        metric(0, 0) = sol[0];
        metric(0, 1) = sol[1];
        metric(1, 0) = sol[1];
        metric(1, 1) = sol[2];

        // Calculate the eigenvalues of the metric tensor to obtain the ellipsis of inertia axes lengths
        BoundedMatrix<double, 2, 2> eigenvects, eigenvals;
        MathUtils<double>::GaussSeidelEigenSystem(metric, eigenvects, eigenvals);
        const double h_1 = std::sqrt(1.0 / eigenvals(0, 0));
        const double h_2 = std::sqrt(1.0 / eigenvals(1, 1));

        // Calculate the reference element size as the average of the ellipsis of intertia axes lengths
        const double h_ref = 0.5 * (h_1 + h_2);

        // Make the metric dimensionless
        metric *= std::pow(h_ref, 2);

        // Calculate metric infimum norm
        // TODO: Using h_min should yield a similar behavior
        const double metric_inf = std::min(eigenvals(0, 0), eigenvals(1, 1));

        return std::make_tuple(h_ref, metric_inf, metric);
    }

    // https://es.wikipedia.org/wiki/Circunelipse_de_Steiner --> 3D extension
    std::tuple<double, double, Matrix> ShockCapturingProcess::CalculateTetrahedraMetricTensor(const Geometry<Node<3>> &rGeometry)
    {
        // Solve the metric problem trans(e)*M*e = 1
        // This means, find the coefficients of the matrix M such that all the edges have unit length
        Vector sol;
        array_1d<double, 6> aux_vect;
        BoundedMatrix<double, 6, 6> aux_mat;
        unsigned int row = 0;
        for (unsigned int i = 0; i < 3; ++i)
        {
            const auto &i_coord = rGeometry[i].Coordinates();
            for (unsigned int j = i + 1; j < 4; ++j)
            {
                const auto &j_coord = rGeometry[j].Coordinates();
                aux_mat(row, 0) = std::pow(i_coord[0] - j_coord[0], 2);
                aux_mat(row, 1) = 2.0 * (i_coord[0] - j_coord[0]) * (i_coord[1] - j_coord[1]);
                aux_mat(row, 2) = 2.0 * (i_coord[0] - j_coord[0]) * (i_coord[2] - j_coord[2]);
                aux_mat(row, 3) = std::pow(i_coord[1] - j_coord[1], 2);
                aux_mat(row, 4) = 2.0 * (i_coord[1] - j_coord[1]) * (i_coord[2] - j_coord[2]);
                aux_mat(row, 5) = std::pow(i_coord[2] - j_coord[2], 2);
                aux_vect(row) = 1.0;
                row++;
            }
        }
        MathUtils<double>::Solve(aux_mat, sol, aux_vect);

        // Set the metric tensor
        Matrix metric(3, 3);
        metric(0, 0) = sol[0];
        metric(0, 1) = sol[1];
        metric(0, 2) = sol[2];
        metric(1, 0) = sol[1];
        metric(1, 1) = sol[3];
        metric(1, 2) = sol[4];
        metric(2, 0) = sol[2];
        metric(2, 1) = sol[4];
        metric(2, 2) = sol[5];

        // Calculate the eigenvalues of the metric tensor to obtain the ellipsis of inertia axes lengths
        BoundedMatrix<double, 3, 3> eigenvects, eigenvals;
        MathUtils<double>::GaussSeidelEigenSystem(metric, eigenvects, eigenvals);
        const double h_1 = std::sqrt(1.0 / eigenvals(0, 0));
        const double h_2 = std::sqrt(1.0 / eigenvals(1, 1));
        const double h_3 = std::sqrt(1.0 / eigenvals(2, 2));

        // Calculate the reference element size as the average of the ellipsis of intertia axes lengths
        const double h_ref = (h_1 + h_2 + h_3) / 3.0;

        // Make the metric dimensionless
        metric *= std::pow(h_ref, 2);

        // Calculate metric infimum norm
        // TODO: Using h_min should yield a similar behavior
        const double metric_inf = std::min(eigenvals(0, 0), std::min(eigenvals(1, 1), eigenvals(2, 2)));

        return std::make_tuple(h_ref, metric_inf, metric);
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ShockCapturingProcess& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
