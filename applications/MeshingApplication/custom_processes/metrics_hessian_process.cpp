// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "utilities/math_utils.h"
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/metrics_math_utils.h"
#include "processes/compute_nodal_gradient_process.h"
#include "custom_processes/metrics_hessian_process.h"

namespace Kratos
{
ComputeHessianSolMetricProcess::ComputeHessianSolMetricProcess(
    ModelPart& rThisModelPart,
    Variable<double>& rVariable,
    Parameters ThisParameters
    ):mThisModelPart(rThisModelPart)
{
    // We push the list of double variables
    mrOriginVariableDoubleList.push_back(&rVariable);

    // We check the parameters
    Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    InitializeVariables(ThisParameters);
}

/***********************************************************************************/
/***********************************************************************************/

ComputeHessianSolMetricProcess::ComputeHessianSolMetricProcess(
    ModelPart& rThisModelPart,
    ComponentType& rVariable,
    Parameters ThisParameters
    ):mThisModelPart(rThisModelPart)
{
    // We push the components list
    mrOriginVariableComponentsList.push_back(&rVariable);

    // We check the parameters
    Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    InitializeVariables(ThisParameters);
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeHessianSolMetricProcess::Execute()
{
    // Computing auxiliar Hessian
    CalculateAuxiliarHessian();

    // Some checks
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    if (mrOriginVariableDoubleList.size() > 0) {
        VariableUtils().CheckVariableExists(*mrOriginVariableDoubleList[0], nodes_array);
    } else {
        VariableUtils().CheckVariableExists(*mrOriginVariableComponentsList[0], nodes_array);
    }

    // Checking NODAL_H
    for (const auto& i_node : nodes_array)
        KRATOS_ERROR_IF_NOT(i_node.Has(NODAL_H)) << "NODAL_H must be computed" << std::endl;

    // Getting dimension
    const std::size_t dimension = mThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Computing metric
    if (dimension == 2) { // 2D
        CalculateMetric<2>();
    } else if (dimension == 3) { // 3D
        CalculateMetric<3>();
    } else {
        KRATOS_ERROR << "Dimension can be only 2D or 3D. Dimension: " << dimension << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
array_1d<double, 3 * (TDim - 1)> ComputeHessianSolMetricProcess::ComputeHessianMetricTensor(
    const Vector& rHessian,
    const double AnisotropicRatio,
    const double ElementMinSize, // This way we can impose as minimum as the previous size if we desire
    const double ElementMaxSize, // This way we can impose as maximum as the previous size if we desire
    const double NodalH
    )
{
    /// The type of array considered for the tensor
    typedef typename std::conditional<TDim == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// Matrix type definition
    typedef BoundedMatrix<double, TDim, TDim> MatrixType;

    // We first transform the Hessian into a matrix
    const MatrixType hessian_matrix = MathUtils<double>::VectorToSymmetricTensor<Vector, MatrixType>(rHessian);

    // Calculating Metric parameters (using equation from remark 4.2.2 on Metric-Based Anisotropic Mesh Adaptation)
    double interpolation_error = mInterpError;
    if (mEstimateInterpError) {
        interpolation_error = mMeshConstant * MathUtils<double>::Max(NodalH, NodalH * norm_frobenius(hessian_matrix)); // NOTE: To compute it properly instead of iterating over the nodes you should iterate over the elements and instead of ElementMaxSize you should iterate over the edges, this is equivalent when using nodes and computing NodalH previously
    }

    // Declaring the eigen system
    MatrixType eigen_vector_matrix, eigen_values_matrix;

    MathUtils<double>::EigenSystem<TDim>(hessian_matrix, eigen_vector_matrix, eigen_values_matrix, 1e-18, 20);

    // We check is the interpolation error is near zero. If it is we will correct it
    if (interpolation_error < std::numeric_limits<double>::epsilon()) { // In practice, the Hessian of function u can be 0, e.g. if u is linear, then |Hu| is not definite. In this particular case, the interpolation error is 0 and we want to prescribe a mesh size which is infinite. To solve this issue, this infinite size prescription is truncated by imposing maximal size hmax . This is equivalent to truncate tiny eigenvalues by lambda  = 1/hmax^2 . See [1] pag. 34
        KRATOS_WARNING("ComputeHessianSolMetricProcess") << "WARNING: Your interpolation error is near zero: " << interpolation_error  <<  ". Computing a local L(inf) upper bound of the interpolation error"<< std::endl;
        const double l_square_minus1 = 1.0/std::pow(ElementMaxSize, 2);
        for (IndexType i = 0; i < TDim; ++i) {
            eigen_values_matrix(i, i) = l_square_minus1;
        }
    } else { // Equation 4.4 from Metric-Based Anisotropic Mesh Adaptation
        const double c_epsilon = mMeshConstant/interpolation_error;
        const double min_ratio = 1.0/std::pow(ElementMinSize, 2);
        const double max_ratio = 1.0/std::pow(ElementMaxSize, 2);

        // Recalculate the Metric eigen values
        for (IndexType i = 0; i < TDim; ++i) {
            eigen_values_matrix(i, i) = MathUtils<double>::Min(MathUtils<double>::Max(c_epsilon * std::abs(eigen_values_matrix(i, i)), max_ratio), min_ratio);
        }
    }

    // Considering anisotropic
    if (AnisotropicRatio < 1.0) {
        double eigen_max = eigen_values_matrix(0, 0);
        double eigen_min = eigen_values_matrix(0, 0);
        for (IndexType i = 1; i < TDim; ++i) {
            eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
            eigen_min = MathUtils<double>::Min(eigen_min, eigen_values_matrix(i, i));
        }

        const double eigen_radius = std::abs(eigen_max - eigen_min) * (1.0 - AnisotropicRatio);
        const double relative_eigen_radius = std::abs(eigen_max - eigen_radius);

        for (IndexType i = 0; i < TDim; ++i)
            eigen_values_matrix(i, i) = MathUtils<double>::Max(MathUtils<double>::Min(eigen_values_matrix(i, i), eigen_max), relative_eigen_radius);
    } else { // NOTE: For isotropic we should consider the maximum of the eigenvalues
        double eigen_max = eigen_values_matrix(0, 0);
        for (IndexType i = 1; i < TDim; ++i)
            eigen_max = MathUtils<double>::Max(eigen_max, eigen_values_matrix(i, i));
        for (IndexType i = 0; i < TDim; ++i)
            eigen_values_matrix(i, i) = eigen_max;
        eigen_vector_matrix = IdentityMatrix(TDim, TDim);
    }

    // We compute the product
    const MatrixType& metric_matrix =  prod(trans(eigen_vector_matrix), prod<MatrixType>(eigen_values_matrix, eigen_vector_matrix));

    // Finally we transform to a vector
    const TensorArrayType& metric = MathUtils<double>::StressTensorToVector<MatrixType, TensorArrayType>(metric_matrix);

    return metric;
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeHessianSolMetricProcess::CalculateAuxiliarHessian()
{
    // Iterate in the elements
    ElementsArrayType& elements_array = mThisModelPart.Elements();
    const int num_elements = static_cast<int>(elements_array.size());
    const auto& it_element_begin = elements_array.begin();

    // Geometry information
    const std::size_t dimension = mThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Declaring auxiliar vector
    const Vector aux_zero_hessian = ZeroVector(3 * (dimension - 1));
    const array_1d<double, 3> aux_zero_vector = ZeroVector(3);

    // Iterate in the nodes
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    // Initialize auxiliar variables
    const auto& it_nodes_begin = nodes_array.begin();
    #pragma omp parallel for
    for(int i_node = 0; i_node < num_nodes; ++i_node) {
        auto it_node = it_nodes_begin + i_node;
        it_node->SetValue(NODAL_AREA, 0.0);
        it_node->SetValue(AUXILIAR_HESSIAN, aux_zero_hessian);
        it_node->SetValue(AUXILIAR_GRADIENT, aux_zero_vector);
    }

    // Compute auxiliar gradient
    if (mrOriginVariableDoubleList.size() > 0) {
        auto gradient_process = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mThisModelPart, *mrOriginVariableDoubleList[0], AUXILIAR_GRADIENT, NODAL_AREA);
        gradient_process.Execute();
    } else {
        auto gradient_process = ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>(mThisModelPart, *mrOriginVariableComponentsList[0], AUXILIAR_GRADIENT, NODAL_AREA);
        gradient_process.Execute();
    }

    // Auxiliar containers
    Matrix DN_DX, J0;
    Vector N;

    #pragma omp parallel for firstprivate(DN_DX,  N, J0)
    for(int i_elem = 0; i_elem < num_elements; ++i_elem) {
        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_geometry.PointsNumber();

        // Resize if needed
        if (DN_DX.size1() != number_of_nodes || DN_DX.size2() != dimension)
            DN_DX.resize(number_of_nodes, dimension);
        if (N.size() != number_of_nodes)
            N.resize(number_of_nodes);
        if (J0.size1() != dimension || J0.size2() != local_space_dimension)
            J0.resize(dimension, local_space_dimension);

        // The integration points
        const auto& integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const std::size_t number_of_integration_points = integration_points.size();

        // The containers of the shape functions and the local gradients
        const auto& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);
        const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(integration_method);

        // 2D case
        if (dimension == 2) {
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                noalias(N) = row(rNcontainer, point_number);

                // Getting the jacobians and local gradients
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
                double detJ0;
                Matrix InvJ0;
                MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
                const Matrix& rDN_De = rDN_DeContainer[point_number];
                GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;

                Matrix values(number_of_nodes, 2);
                for(IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    const array_1d<double, 3>& aux_grad = r_geometry[i_node].GetValue(AUXILIAR_GRADIENT);
                    for (IndexType i_dim = 0; i_dim < 2; ++i_dim)
                        values(i_node, i_dim) = aux_grad[i_dim];
                }

                const BoundedMatrix<double,2, 2>& hessian = prod(trans(DN_DX), values);
                const array_1d<double, 3>& hessian_cond = MathUtils<double>::StressTensorToVector<BoundedMatrix<double, 2, 2>, array_1d<double, 3>>(hessian);

                for(IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    auto& aux_hessian = r_geometry[i_node].GetValue(AUXILIAR_HESSIAN);
                    for(IndexType k = 0; k < 3; ++k) {
                        #pragma omp atomic
                        aux_hessian[k] += N[i_node] * gauss_point_volume * hessian_cond[k];
                    }
                }
            }
        } else { // 3D case
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                noalias(N) = row(rNcontainer, point_number);

                // Getting the jacobians and local gradients
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
                double detJ0;
                Matrix InvJ0;
                MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
                const Matrix& rDN_De = rDN_DeContainer[point_number];
                GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;

                Matrix values(number_of_nodes, 3);
                for(IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    const array_1d<double, 3>& aux_grad = r_geometry[i_node].GetValue(AUXILIAR_GRADIENT);
                    for (IndexType i_dim = 0; i_dim < 3; ++i_dim)
                        values(i_node, i_dim) = aux_grad[i_dim];
                }

                const BoundedMatrix<double, 3, 3> hessian = prod(trans(DN_DX), values);
                const array_1d<double, 6>& hessian_cond = MathUtils<double>::StressTensorToVector<BoundedMatrix<double, 3, 3>, array_1d<double, 6>>(hessian);

                for(IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                    auto& aux_hessian = r_geometry[i_node].GetValue(AUXILIAR_HESSIAN);
                    for(IndexType k = 0; k < 6; ++k) {
                        #pragma omp atomic
                        aux_hessian[k] += N[i_node] * gauss_point_volume * hessian_cond[k];
                    }
                }
            }
        }
    }

    #pragma omp parallel for
    for(int i_node = 0; i_node < num_nodes; ++i_node) {
        auto it_node = nodes_array.begin() + i_node;
        const double nodal_area = it_node->GetValue(NODAL_AREA);
        if (nodal_area > std::numeric_limits<double>::epsilon()) {
            it_node->GetValue(AUXILIAR_HESSIAN) /= nodal_area;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

double ComputeHessianSolMetricProcess::CalculateAnisotropicRatio(
    const double Distance,
    const double AnisotropicRatio,
    const double BoundLayer,
    const Interpolation rInterpolation
    )
{
    const double tolerance = 1.0e-12;
    double ratio = 1.0; // NOTE: Isotropic mesh
    if (AnisotropicRatio < 1.0) {
        if (std::abs(Distance) <= BoundLayer) {
            if (rInterpolation == Interpolation::CONSTANT)
                ratio = AnisotropicRatio;
            else if (rInterpolation == Interpolation::LINEAR)
                ratio = AnisotropicRatio + (std::abs(Distance)/BoundLayer) * (1.0 - AnisotropicRatio);
            else if (rInterpolation == Interpolation::EXPONENTIAL) {
                ratio = - std::log(std::abs(Distance)/BoundLayer) * AnisotropicRatio + tolerance;
                if (ratio > 1.0) ratio = 1.0;
            }
        }
    }

    return ratio;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void ComputeHessianSolMetricProcess::CalculateMetric()
{
    /// The type of array considered for the tensor
    typedef typename std::conditional<TDim == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    // Iterate in the nodes
    NodesArrayType& nodes_array = mThisModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    // Tensor variable definition
    const Variable<TensorArrayType>& tensor_variable = KratosComponents<Variable<TensorArrayType>>::Get("METRIC_TENSOR_"+std::to_string(TDim)+"D");

    // Setting metric in case not defined
    const auto it_node_begin = nodes_array.begin();
    if (!it_node_begin->Has(tensor_variable)) {
        // Declaring auxiliar vector
        const TensorArrayType aux_zero_vector = ZeroVector(3 * (TDim - 1));
        VariableUtils().SetNonHistoricalVariable(tensor_variable, aux_zero_vector, nodes_array);
    }

    // Ratio reference variable
    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mRatioReferenceVariable)) << "Variable " << mRatioReferenceVariable << " is not a double variable" << std::endl;
    const auto& reference_var = KratosComponents<Variable<double>>::Get(mRatioReferenceVariable);

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = it_node_begin + i;

        const Vector& hessian = it_node->GetValue(AUXILIAR_HESSIAN);

        const double nodal_h = it_node->GetValue(NODAL_H);

        double element_min_size = mMinSize;
        if ((element_min_size > nodal_h) && mEnforceCurrent) element_min_size = nodal_h;
        double element_max_size = mMaxSize;
        if ((element_max_size > nodal_h) && mEnforceCurrent) element_max_size = nodal_h;

        // Isotropic by default
        double ratio = 1.0;

        if (it_node->SolutionStepsDataHas(reference_var)) {
            const double ratio_reference = it_node->FastGetSolutionStepValue(reference_var);
            ratio = CalculateAnisotropicRatio(ratio_reference, mAnisotropicRatio, mBoundLayer, mInterpolation);
        }

        // For postprocess pourposes
        it_node->SetValue(ANISOTROPIC_RATIO, ratio);

        // We compute the metric
        KRATOS_DEBUG_ERROR_IF_NOT(it_node->Has(tensor_variable)) << "METRIC_TENSOR_" + std::to_string(TDim) + "D  not defined for node " << it_node->Id() << std::endl;
        TensorArrayType& metric = it_node->GetValue(tensor_variable);

        const double norm_metric = norm_2(metric);
        if (norm_metric > 0.0) {// NOTE: This means we combine differents metrics, at the same time means that the metric should be reseted each time
            const TensorArrayType& old_metric = it_node->GetValue(tensor_variable);
            const TensorArrayType& new_metric = ComputeHessianMetricTensor<TDim>(hessian, ratio, element_min_size, element_max_size, nodal_h);

            metric = MetricsMathUtils<TDim>::IntersectMetrics(old_metric, new_metric);
        } else {
            metric = ComputeHessianMetricTensor<TDim>(hessian, ratio, element_min_size, element_max_size, nodal_h);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters ComputeHessianSolMetricProcess::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size"                        : 0.1,
        "maximal_size"                        : 10.0,
        "enforce_current"                     : true,
        "hessian_strategy_parameters":
        {
            "metric_variable"                  : ["DISTANCE"],
            "estimate_interpolation_error"         : false,
            "interpolation_error"                  : 1.0e-6,
            "mesh_dependent_constant"              : 0.28125
        },
        "anisotropy_remeshing"                : true,
        "anisotropy_parameters":
        {
            "reference_variable_name"              : "DISTANCE",
            "hmin_over_hmax_anisotropic_ratio"     : 1.0,
            "boundary_layer_max_distance"          : 1.0,
            "interpolation"                        : "Linear"
        }
    })" );

    // Identify the dimension first
    const SizeType dimension = mThisModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // The mesh dependent constant depends on dimension
    if (dimension == 2) {
        default_parameters["hessian_strategy_parameters"]["mesh_dependent_constant"].SetDouble(2.0/9.0);
    } else if (dimension == 3) {
        default_parameters["hessian_strategy_parameters"]["mesh_dependent_constant"].SetDouble(9.0/32.0);
    } else {
        KRATOS_ERROR << "Dimension can be only 2D or 3D. Dimension: " << dimension << std::endl;
    }

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeHessianSolMetricProcess::InitializeVariables(Parameters ThisParameters)
{
    // Get default variables
    Parameters default_parameters = GetDefaultParameters();

    // Set variables
    mMinSize = ThisParameters["minimal_size"].GetDouble();
    mMaxSize = ThisParameters["maximal_size"].GetDouble();
    mEnforceCurrent = ThisParameters["enforce_current"].GetBool();

    // In case we have isotropic remeshing (default values)
    if (ThisParameters["anisotropy_remeshing"].GetBool() == false) {
        mEstimateInterpError = default_parameters["hessian_strategy_parameters"]["estimate_interpolation_error"].GetBool();
        mInterpError = default_parameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble();
        mMeshConstant = default_parameters["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble();
        mRatioReferenceVariable = default_parameters["anisotropy_parameters"]["reference_variable_name"].GetString();
        mAnisotropicRatio = default_parameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = default_parameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(default_parameters["anisotropy_parameters"]["interpolation"].GetString());
    } else {
        mEstimateInterpError = ThisParameters["hessian_strategy_parameters"]["estimate_interpolation_error"].GetBool();
        mInterpError = ThisParameters["hessian_strategy_parameters"]["interpolation_error"].GetDouble();
        mMeshConstant = ThisParameters["hessian_strategy_parameters"]["mesh_dependent_constant"].GetDouble();
        mRatioReferenceVariable = ThisParameters["anisotropy_parameters"]["reference_variable_name"].GetString();
        mAnisotropicRatio = ThisParameters["anisotropy_parameters"]["hmin_over_hmax_anisotropic_ratio"].GetDouble();
        mBoundLayer = ThisParameters["anisotropy_parameters"]["boundary_layer_max_distance"].GetDouble();
        mInterpolation = ConvertInter(ThisParameters["anisotropy_parameters"]["interpolation"].GetString());
    }
}

};// namespace Kratos.
