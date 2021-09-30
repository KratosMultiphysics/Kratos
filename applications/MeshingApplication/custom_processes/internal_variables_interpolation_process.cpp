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

// System includes

// External includes

// Project includes
#include "utilities/variable_utils.h"
#include "custom_processes/internal_variables_interpolation_process.h"

namespace Kratos
{
InternalVariablesInterpolationProcess::InternalVariablesInterpolationProcess(
    ModelPart& rOriginMainModelPart,
    ModelPart& rDestinationMainModelPart,
    Parameters ThisParameters
    ):mrOriginMainModelPart(rOriginMainModelPart),
    mrDestinationMainModelPart(rDestinationMainModelPart),
    mDimension(rDestinationMainModelPart.GetProcessInfo()[DOMAIN_SIZE])
{
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mAllocationSize = ThisParameters["allocation_size"].GetInt();
    mBucketSize = ThisParameters["bucket_size"].GetInt();
    mSearchFactor = ThisParameters["search_factor"].GetDouble();
    mThisInterpolationType = ConvertInter(ThisParameters["interpolation_type"].GetString());

    if (ThisParameters["internal_variable_interpolation_list"].IsArray() == true) {
        auto variable_array_list = ThisParameters["internal_variable_interpolation_list"];

        for (IndexType i_var = 0; i_var < variable_array_list.size(); ++i_var) {
            const std::string& r_variable_name = variable_array_list[i_var].GetString();
            mInternalVariableList.push_back(r_variable_name);
        }
    } else {
        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: No variables to interpolate, look that internal_variable_interpolation_list is correctly defined in your parameters" << std::endl;
        mInternalVariableList.clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::Execute()
{
    if (mThisInterpolationType == InterpolationTypes::CLOSEST_POINT_TRANSFER && ComputeTotalNumberOfVariables() > 0) {
        InterpolateGaussPointsClosestPointTransfer();
    } else if (mThisInterpolationType == InterpolationTypes::LEAST_SQUARE_TRANSFER && ComputeTotalNumberOfVariables() > 0) {
        InterpolateGaussPointsLeastSquareTransfer();
    } else if (mThisInterpolationType == InterpolationTypes::SHAPE_FUNCTION_TRANSFER && ComputeTotalNumberOfVariables() > 0) {
//         InterpolateGaussPointsShapeFunctionTransfer();
        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: SHAPE FUNCTION TRANSFER THIS DOESN'T WORK, AND REQUIRES EXTRA STORE. PLEASE COOSE ANY OTHER ALTERNATIVE" << std::endl;
    } else
        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: INTERPOLATION TYPE NOT AVALAIBLE OR EMPTY LIST" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

PointVector InternalVariablesInterpolationProcess::CreateGaussPointList(ModelPart& ThisModelPart)
{
    PointVector this_point_vector;

    GeometryData::IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_1;

    // Iterate in the elements
    ElementsArrayType& r_elements_array = ThisModelPart.Elements();
    const int num_elements = static_cast<int>(r_elements_array.size());

    const ProcessInfo& r_current_process_info = ThisModelPart.GetProcessInfo();
    const auto it_elem_begin = r_elements_array.begin();

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        PointVector points_buffer;

        #pragma omp for firstprivate(this_integration_method)
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
            if (!old_entity) { // We don't interpolate from preserved meshes
                // Getting the geometry
                GeometryType& r_this_geometry = it_elem->GetGeometry();

                // Getting the integration points
                this_integration_method = it_elem->GetIntegrationMethod();
                const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Computing the Jacobian
                Vector vector_det_j(integration_points_number);
                r_this_geometry.DeterminantOfJacobian(vector_det_j,this_integration_method);

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                it_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,r_current_process_info);

                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point ) {
                    const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();

                    // We compute the corresponding weight
                    const double weight = vector_det_j[i_gauss_point] * integration_points[i_gauss_point].Weight();

                    // We compute the global coordinates
                    array_1d<double, 3> global_coordinates;
                    global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                    // We create the respective GP
                    ConstitutiveLaw::Pointer p_origin_cl = constitutive_law_vector[i_gauss_point];
                    PointTypePointer p_point = PointTypePointer(new PointType(global_coordinates, p_origin_cl, weight));

                    // We save the values not accesible from the CL (it consummes memory, so preferably use variable from the CL)
                    for (auto& variable_name : mInternalVariableList) {
                        if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                            const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                            if (!p_origin_cl->Has(r_variable)) {
                                SaveValuesOnGaussPoint(r_variable, p_point, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                            const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                            if (!p_origin_cl->Has(r_variable)) {
                                SaveValuesOnGaussPoint(r_variable, p_point, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                            const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                            if (!p_origin_cl->Has(r_variable)) {
                                SaveValuesOnGaussPoint(r_variable, p_point, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                            const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                            if (!p_origin_cl->Has(r_variable)) {
                                SaveValuesOnGaussPoint(r_variable, p_point, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else {
                            KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                        }
                    }

                    // Finally we push over the the buffer vector
                    points_buffer.push_back(p_point);
                }
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(points_buffer.begin(),points_buffer.end(),back_inserter(this_point_vector));
        }
    }

    return this_point_vector;
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsClosestPointTransfer()
{
    // We Initialize the process info
    const ProcessInfo& r_current_process_info = mrDestinationMainModelPart.GetProcessInfo();

    // We update the list of points
    mPointListOrigin.clear();
    mPointListOrigin = CreateGaussPointList(mrOriginMainModelPart);

    //#pragma omp parallel firstprivate(mPointListOrigin)
    //{
        // We initialize the intergration method
        GeometryData::IntegrationMethod this_integration_method;

        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTree tree_points(mPointListOrigin.begin(), mPointListOrigin.end(), mBucketSize);

        // Iterate over the destination elements
        ElementsArrayType& r_elements_array = mrDestinationMainModelPart.Elements();
        const auto it_elem_begin = r_elements_array.begin();
        const int num_elements = static_cast<int>(r_elements_array.size());

//         #pragma omp for
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
            if (!old_entity) { // We don't interpolate from preserved meshes
                // Getting the geometry
                GeometryType& r_this_geometry = it_elem->GetGeometry();

                // Getting the integration points
                this_integration_method = it_elem->GetIntegrationMethod();
                const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                it_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,r_current_process_info);

                for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point ) {
                    // We compute the global coordinates
                    const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();
                    array_1d<double, 3> global_coordinates;
                    global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                    PointTypePointer p_gp_origin = tree_points.SearchNearestPoint(global_coordinates);

                    ConstitutiveLaw::Pointer p_origin_cl = p_gp_origin->GetConstitutiveLaw();
                    ConstitutiveLaw::Pointer p_destination_cl = constitutive_law_vector[i_gauss_point];

                    // Get and set variable
                    for (auto& variable_name : mInternalVariableList) {
                        if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                            const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                            if (p_destination_cl->Has(r_variable)) {
                                GetAndSetDirectVariableOnConstitutiveLaw(r_variable, p_origin_cl, p_destination_cl, r_current_process_info);
                            } else {
                                GetAndSetDirectVariableOnElements(r_variable, p_gp_origin, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                            const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                            if (p_destination_cl->Has(r_variable)) {
                                GetAndSetDirectVariableOnConstitutiveLaw(r_variable, p_origin_cl, p_destination_cl, r_current_process_info);
                            } else {
                                GetAndSetDirectVariableOnElements(r_variable, p_gp_origin, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                            const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                            if (p_destination_cl->Has(r_variable)) {
                                GetAndSetDirectVariableOnConstitutiveLaw(r_variable, p_origin_cl, p_destination_cl, r_current_process_info);
                            } else {
                                GetAndSetDirectVariableOnElements(r_variable, p_gp_origin, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                            const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                            if (p_destination_cl->Has(r_variable)) {
                                GetAndSetDirectVariableOnConstitutiveLaw(r_variable, p_origin_cl, p_destination_cl, r_current_process_info);
                            } else {
                                GetAndSetDirectVariableOnElements(r_variable, p_gp_origin, it_elem, i_gauss_point, r_current_process_info);
                            }
                        } else {
                            KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                        }
                    }
                }
            }
        }
    //}
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsLeastSquareTransfer()
{
    // We Initialize the process info
    const ProcessInfo& r_current_process_info = mrDestinationMainModelPart.GetProcessInfo();

    // We update the list of points
    mPointListOrigin.clear();
    mPointListOrigin = CreateGaussPointList(mrOriginMainModelPart);

    // Check the NODAL_H
    NodesArrayType& r_nodes_array = mrDestinationMainModelPart.Nodes();
    for (auto& i_node : r_nodes_array)
        KRATOS_ERROR_IF_NOT(i_node.Has(NODAL_H)) << "NODAL_H must be computed" << std::endl;

    //#pragma omp parallel firstprivate(mPointListOrigin)
    //{
        // We initialize the intergration method
        GeometryData::IntegrationMethod this_integration_method;

        // Initialize values
        PointVector points_found(mAllocationSize);
        std::vector<double> point_distances(mAllocationSize);
        std::size_t number_points_found = 0;

        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTree tree_points(mPointListOrigin.begin(), mPointListOrigin.end(), mBucketSize);

        // Iterate over the destination elements
        ElementsArrayType& r_elements_array = mrDestinationMainModelPart.Elements();
        const auto it_elem_begin = r_elements_array.begin();
        const int num_elements = static_cast<int>(r_elements_array.size());

//         #pragma omp for firstprivate(this_integration_method)
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            const bool old_entity = it_elem->IsDefined(OLD_ENTITY) ? it_elem->Is(OLD_ENTITY) : false;
            if (!old_entity) { // We don't interpolate from preserved meshes
                // Getting the geometry
                GeometryType& r_this_geometry = it_elem->GetGeometry();

                // Getting the integration points
                this_integration_method = it_elem->GetIntegrationMethod();
                const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
                const std::size_t integration_points_number = integration_points.size();

                // Getting the CL
                std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
                it_elem->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,r_current_process_info);

                // Computing the radius
                const double radius = mSearchFactor *  (mDimension == 2 ? std::sqrt(r_this_geometry.Area()) : std::cbrt(r_this_geometry.Volume()));

                // We get the NODAL_H vector
                Vector nodal_h_vector(r_this_geometry.size());
                for (std::size_t i_node = 0; i_node < r_this_geometry.size(); ++i_node)
                    nodal_h_vector[i_node] = r_this_geometry[i_node].GetValue(NODAL_H);

                for (std::size_t i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point ) {
                    // We compute the global coordinates
                    const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();
                    array_1d<double, 3> global_coordinates;
                    global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                    // We compute the pondered characteristic length
                    Vector N( r_this_geometry.size() );
                    r_this_geometry.ShapeFunctionsValues( N, local_coordinates );
                    const double characteristic_length = inner_prod(N, nodal_h_vector);

                    number_points_found = tree_points.SearchInRadius(global_coordinates, radius, points_found.begin(), point_distances.begin(), mAllocationSize);

                    // The destination CL
                    ConstitutiveLaw::Pointer p_destination_cl = constitutive_law_vector[i_gauss_point];

                    if (number_points_found > 0) {
                        for (auto& variable_name : mInternalVariableList) {
                            if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                                const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                                if (p_destination_cl->Has(r_variable)) {
                                    GetAndSetWeightedVariableOnConstitutiveLaw(r_variable, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, r_current_process_info);
                                } else {
                                    GetAndSetWeightedVariableOnElements(r_variable, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, r_current_process_info);
                                }
                            } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                                const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                                if (p_destination_cl->Has(r_variable)) {
                                    GetAndSetWeightedVariableOnConstitutiveLaw(r_variable, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, r_current_process_info);
                                } else {
                                    GetAndSetWeightedVariableOnElements(r_variable, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, r_current_process_info);
                                }
                            } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                                const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                                if (p_destination_cl->Has(r_variable)) {
                                    GetAndSetWeightedVariableOnConstitutiveLaw(r_variable, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, r_current_process_info);
                                } else {
                                    GetAndSetWeightedVariableOnElements(r_variable, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, r_current_process_info);
                                }
                            } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                                const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                                if (p_destination_cl->Has(r_variable)) {
                                    GetAndSetWeightedVariableOnConstitutiveLaw(r_variable, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, r_current_process_info);
                                } else {
                                    GetAndSetWeightedVariableOnElements(r_variable, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, r_current_process_info);
                                }
                            } else {
                                KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                            }
                        }
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: It wasn't impossible to find any Gauss Point from where interpolate the internal variables" << std::endl;
                    }
                }
            }
        }
    //}
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsShapeFunctionTransfer()
{
    // Initialize some values
    GeometryData::IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_1;

    // Iterate in the nodes to initialize the values
    NodesArrayType& r_nodes_array = mrOriginMainModelPart.Nodes();

    /* Nodes */
    for (auto& r_variable_name : mInternalVariableList) {
        if (KratosComponents<DoubleVarType>::Has(r_variable_name)) {
            const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(r_variable_name);
            VariableUtils().SetNonHistoricalVariable(r_variable, r_variable.Zero(), r_nodes_array);
        } else if (KratosComponents<ArrayVarType>::Has(r_variable_name)) {
            const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(r_variable_name);
            VariableUtils().SetNonHistoricalVariable(r_variable, r_variable.Zero(), r_nodes_array);
        } else if (KratosComponents<VectorVarType>::Has(r_variable_name)) {
            const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(r_variable_name);
            VariableUtils().SetNonHistoricalVariable(r_variable, r_variable.Zero(), r_nodes_array);
        } else if (KratosComponents<MatrixVarType>::Has(r_variable_name)) {
            const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(r_variable_name);
            VariableUtils().SetNonHistoricalVariable(r_variable, r_variable.Zero(), r_nodes_array);
        } else {
            KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << r_variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
        }
    }

    // Iterate in the elements to ponderate the values
    ElementsArrayType& r_elements_array = mrOriginMainModelPart.Elements();

    const ProcessInfo& r_origin_process_info = mrOriginMainModelPart.GetProcessInfo();

    /* Elements */
    block_for_each(r_elements_array, this_integration_method,
        [&r_origin_process_info, this](Element& rElement, GeometryData::IntegrationMethod& this_integration_method) {

        const bool old_entity = rElement.IsDefined(OLD_ENTITY) ? rElement.Is(OLD_ENTITY) : false;
        if (!old_entity) { // We don't interpolate from preserved meshes
            // Getting the geometry
            GeometryType& r_this_geometry = rElement.GetGeometry();

            // Getting the integration points
            this_integration_method = rElement.GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const std::size_t integration_points_number = integration_points.size();

            // Computing the Jacobian
            Vector vector_det_j(integration_points_number);
            r_this_geometry.DeterminantOfJacobian(vector_det_j,this_integration_method);

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,r_origin_process_info);

            // We initialize the total weigth
            double total_weight = 0.0;

            for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point ) {
                const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();

                // We compute the corresponding weight
                const double weight = vector_det_j[i_gauss_point] * integration_points[i_gauss_point].Weight();
                total_weight += weight;

                // We compute the pondered characteristic length
                Vector N( r_this_geometry.size() );
                r_this_geometry.ShapeFunctionsValues( N, local_coordinates );

                // We compute the global coordinates
                array_1d<double, 3> global_coordinates;
                global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                // The origin CL
                ConstitutiveLaw::Pointer p_origin_cl = constitutive_law_vector[i_gauss_point];

                // We interpolate and add the variable
                for (auto& variable_name : mInternalVariableList) {
                    if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                        const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                        if (p_origin_cl->Has(r_variable)) {
                            InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, r_variable, N, p_origin_cl, weight);
                        } else {
                            InterpolateAddVariableOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, weight, r_origin_process_info);
                        }
                    } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                        const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                        if (p_origin_cl->Has(r_variable)) {
                            InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, r_variable, N, p_origin_cl, weight);
                        } else {
                            InterpolateAddVariableOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, weight, r_origin_process_info);
                        }
                    } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                        const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                        if (p_origin_cl->Has(r_variable)) {
                            InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, r_variable, N, p_origin_cl, weight);
                        } else {
                            InterpolateAddVariableOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, weight, r_origin_process_info);
                        }
                    } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                        const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                        if (p_origin_cl->Has(r_variable)) {
                            InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, r_variable, N, p_origin_cl, weight);
                        } else {
                            InterpolateAddVariableOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, weight, r_origin_process_info);
                        }
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                    }
                }
            }

            // We divide by the total weight
            for (auto& variable_name : mInternalVariableList) {
                if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                    const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                    PonderateVariable(r_this_geometry, r_variable, total_weight);
                } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                    const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                    PonderateVariable(r_this_geometry, r_variable, total_weight);
                } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                    const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                    PonderateVariable(r_this_geometry, r_variable, total_weight);
                } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                    const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                    PonderateVariable(r_this_geometry, r_variable, total_weight);
                } else {
                    KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                }
            }
        }
    });

    // We interpolate to the new nodes
    if (mDimension == 2) {
        InterpolateToNodes<2>();
    } else {
        InterpolateToNodes<3>();
    }

    // Finally we interpolate to the new GP
    ElementsArrayType& r_elements_array_destination = mrDestinationMainModelPart.Elements();

    // The destination process info
    const ProcessInfo& r_destination_process_info = mrOriginMainModelPart.GetProcessInfo();

    /* Elements */
    block_for_each(r_elements_array_destination, this_integration_method,
        [&r_destination_process_info, this](Element& rElement, GeometryData::IntegrationMethod& this_integration_method) {

        const bool old_entity = rElement.IsDefined(OLD_ENTITY) ? rElement.Is(OLD_ENTITY) : false;
        if (!old_entity) { // We don't interpolate from preserved meshes
            // Getting the geometry
            GeometryType& r_this_geometry = rElement.GetGeometry();

            // Getting the integration points
            this_integration_method = rElement.GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const std::size_t integration_points_number = integration_points.size();

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector,r_destination_process_info);

            for (IndexType i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point ) {
                const array_1d<double, 3>& local_coordinates = integration_points[i_gauss_point].Coordinates();

                // We compute the pondered characteristic length
                Vector N( r_this_geometry.size() );
                r_this_geometry.ShapeFunctionsValues( N, local_coordinates );

                // We compute the global coordinates
                array_1d<double, 3> global_coordinates;
                global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                Vector values(r_this_geometry.size() );

                // The destination CL
                ConstitutiveLaw::Pointer p_destination_cl = constitutive_law_vector[i_gauss_point];

                for (auto& variable_name : mInternalVariableList) {
                    if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                        const DoubleVarType& r_variable = KratosComponents<DoubleVarType>::Get(variable_name);
                        if (p_destination_cl->Has(r_variable)) {
                            SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, r_variable, N, p_destination_cl, r_destination_process_info);
                        } else {
                            SetInterpolatedValueOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, r_destination_process_info);
                        }
                    } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                        const ArrayVarType& r_variable = KratosComponents<ArrayVarType>::Get(variable_name);
                        if (p_destination_cl->Has(r_variable)) {
                            SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, r_variable, N, p_destination_cl, r_destination_process_info);
                        } else {
                            SetInterpolatedValueOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, r_destination_process_info);
                        }
                    } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                        const VectorVarType& r_variable = KratosComponents<VectorVarType>::Get(variable_name);
                        if (p_destination_cl->Has(r_variable)) {
                            SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, r_variable, N, p_destination_cl, r_destination_process_info);
                        } else {
                            SetInterpolatedValueOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, r_destination_process_info);
                        }
                    } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                        const MatrixVarType& r_variable = KratosComponents<MatrixVarType>::Get(variable_name);
                        if (p_destination_cl->Has(r_variable)) {
                            SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, r_variable, N, p_destination_cl, r_destination_process_info);
                        } else {
                            SetInterpolatedValueOnElement(r_this_geometry, r_variable, N, rElement, i_gauss_point, r_destination_process_info);
                        }
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                    }
                }
            }
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t InternalVariablesInterpolationProcess::ComputeTotalNumberOfVariables()
{
    return mInternalVariableList.size();
}

/***********************************************************************************/
/***********************************************************************************/

InternalVariablesInterpolationProcess::InterpolationTypes InternalVariablesInterpolationProcess::ConvertInter(const std::string& Str)
{
    if(Str == "CPT" || Str == "CLOSEST_POINT_TRANSFER")
        return InterpolationTypes::CLOSEST_POINT_TRANSFER;
    else if(Str == "LST" || Str == "LEAST_SQUARE_TRANSFER")
        return InterpolationTypes::LEAST_SQUARE_TRANSFER;
    else if(Str == "SFT" || Str == "SHAPE_FUNCTION_TRANSFER")
        return InterpolationTypes::SHAPE_FUNCTION_TRANSFER;
    else
        return InterpolationTypes::LEAST_SQUARE_TRANSFER;
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters InternalVariablesInterpolationProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "allocation_size"                      : 1000,
        "bucket_size"                          : 4,
        "search_factor"                        : 2,
        "interpolation_type"                   : "LST",
        "internal_variable_interpolation_list" :[]
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<double>(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    )
{
    double origin_value;
    origin_value = pConstitutiveLaw->GetValue(rThisVar, origin_value);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAdd(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_value * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<array_1d<double, 3>>(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    )
{
    array_1d<double, 3> origin_value;
    origin_value = pConstitutiveLaw->GetValue(rThisVar, origin_value);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAddVector(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_value * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<Vector>(
    GeometryType& rThisGeometry,
    const Variable<Vector>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    )
{
    Vector origin_value;
    origin_value = pConstitutiveLaw->GetValue(rThisVar, origin_value);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAddVector(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_value * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnConstitutiveLaw<Matrix>(
    GeometryType& rThisGeometry,
    const Variable<Matrix>& rThisVar,
    const Vector& N,
    ConstitutiveLaw::Pointer& pConstitutiveLaw,
    const double Weight
    )
{
    Matrix origin_value;
    origin_value = pConstitutiveLaw->GetValue(rThisVar, origin_value);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAddMatrix(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_value * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    std::vector<double> origin_values;
    rElement.CalculateOnIntegrationPoints(rThisVar, origin_values, rCurrentProcessInfo);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAdd(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_values[GaussPointId] * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    std::vector<array_1d<double, 3>> origin_values;
    rElement.CalculateOnIntegrationPoints(rThisVar, origin_values, rCurrentProcessInfo);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAddVector(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_values[GaussPointId] * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<Vector>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    std::vector<Vector> origin_values;
    rElement.CalculateOnIntegrationPoints(rThisVar, origin_values, rCurrentProcessInfo);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAddVector(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_values[GaussPointId] * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::InterpolateAddVariableOnElement(
    GeometryType& rThisGeometry,
    const Variable<Matrix>& rThisVar,
    const Vector& N,
    Element& rElement,
    const IndexType GaussPointId,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    std::vector<Matrix> origin_values;
    rElement.CalculateOnIntegrationPoints(rThisVar, origin_values, rCurrentProcessInfo);

    // We sum all the contributions
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicAddMatrix(rThisGeometry[i_node].GetValue(rThisVar), N[i_node] * origin_values[GaussPointId] * Weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<double>& rThisVar,
    const double TotalWeight
    )
{
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicDiv(rThisGeometry[i_node].GetValue(rThisVar), TotalWeight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<array_1d<double, 3>>& rThisVar,
    const double TotalWeight
    )
{
    array_1d<double, 3> aux;
    for (std::size_t i = 0; i < 3; ++i) {
        aux[i] = TotalWeight;
    }
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicDivVector(rThisGeometry[i_node].GetValue(rThisVar), aux);
    }
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<Vector>& rThisVar,
    const double TotalWeight
    )
{
    Vector aux(rThisGeometry[0].GetValue(rThisVar).size());
    for (std::size_t i = 0; i < aux.size(); ++i) {
        aux[i] = TotalWeight;
    }
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicDivVector(rThisGeometry[i_node].GetValue(rThisVar), aux);
    }
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void InternalVariablesInterpolationProcess::PonderateVariable(
    GeometryType& rThisGeometry,
    const Variable<Matrix>& rThisVar,
    const double TotalWeight
    )
{
    const std::size_t size_1 = rThisGeometry[0].GetValue(rThisVar).size1();
    const std::size_t size_2 = rThisGeometry[0].GetValue(rThisVar).size2();
    Matrix aux(size_1, size_2);
    for (std::size_t i = 0; i < size_1; ++i) {
        for (std::size_t j = 0; j < size_2; ++j) {
            aux(i, j) = TotalWeight;
        }
    }
    for (std::size_t i_node = 0; i_node < rThisGeometry.size(); ++i_node) {
        AtomicDivMatrix(rThisGeometry[i_node].GetValue(rThisVar), aux);
    }
}

}  // namespace Kratos.
