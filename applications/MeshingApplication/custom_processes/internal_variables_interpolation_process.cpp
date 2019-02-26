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
    Parameters default_parameters = Parameters(R"(
    {
        "allocation_size"                      : 1000,
        "bucket_size"                          : 4,
        "search_factor"                        : 2,
        "interpolation_type"                   : "LST",
        "internal_variable_interpolation_list" :[]
    })" );

    ThisParameters.ValidateAndAssignDefaults(default_parameters);

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
    ElementsArrayType& elements_array = ThisModelPart.Elements();
    const int num_elements = static_cast<int>(elements_array.size());

    const ProcessInfo& current_process_info = ThisModelPart.GetProcessInfo();
    const auto it_elem_begin = elements_array.begin();

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        PointVector points_buffer;

        #pragma omp for firstprivate(this_integration_method)
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;

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
            it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

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
                        const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                        if (!p_origin_cl->Has(this_var)) {
                            SaveValuesOnGaussPoint(this_var, p_point, it_elem, i_gauss_point, current_process_info);
                        }
                    } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                        const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                        if (!p_origin_cl->Has(this_var)) {
                            SaveValuesOnGaussPoint(this_var, p_point, it_elem, i_gauss_point, current_process_info);
                        }
                    } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                        const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                        if (!p_origin_cl->Has(this_var)) {
                            SaveValuesOnGaussPoint(this_var, p_point, it_elem, i_gauss_point, current_process_info);
                        }
                    } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                        const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                        if (!p_origin_cl->Has(this_var)) {
                            SaveValuesOnGaussPoint(this_var, p_point, it_elem, i_gauss_point, current_process_info);
                        }
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                    }
                }

                // Finally we push over the the buffer vector
                points_buffer.push_back(p_point);
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
    const ProcessInfo& current_process_info = mrDestinationMainModelPart.GetProcessInfo();

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
        ElementsArrayType& elements_array = mrDestinationMainModelPart.Elements();
        auto num_elements = elements_array.end() - elements_array.begin();

        //#pragma omp for
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = elements_array.begin() + i;

            // Getting the geometry
            GeometryType& r_this_geometry = it_elem->GetGeometry();

            // Getting the integration points
            this_integration_method = it_elem->GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const std::size_t integration_points_number = integration_points.size();

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

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
                        const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                        if (p_destination_cl->Has(this_var)) {
                            GetAndSetDirectVariableOnConstitutiveLaw(this_var, p_origin_cl, p_destination_cl, current_process_info);
                        } else {
                            GetAndSetDirectVariableOnElements(this_var, p_gp_origin, it_elem, i_gauss_point, current_process_info);
                        }
                    } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                        const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                        if (p_destination_cl->Has(this_var)) {
                            GetAndSetDirectVariableOnConstitutiveLaw(this_var, p_origin_cl, p_destination_cl, current_process_info);
                        } else {
                            GetAndSetDirectVariableOnElements(this_var, p_gp_origin, it_elem, i_gauss_point, current_process_info);
                        }
                    } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                        const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                        if (p_destination_cl->Has(this_var)) {
                            GetAndSetDirectVariableOnConstitutiveLaw(this_var, p_origin_cl, p_destination_cl, current_process_info);
                        } else {
                            GetAndSetDirectVariableOnElements(this_var, p_gp_origin, it_elem, i_gauss_point, current_process_info);
                        }
                    } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                        const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                        if (p_destination_cl->Has(this_var)) {
                            GetAndSetDirectVariableOnConstitutiveLaw(this_var, p_origin_cl, p_destination_cl, current_process_info);
                        } else {
                            GetAndSetDirectVariableOnElements(this_var, p_gp_origin, it_elem, i_gauss_point, current_process_info);
                        }
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
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
    const ProcessInfo& current_process_info = mrDestinationMainModelPart.GetProcessInfo();

    // We update the list of points
    mPointListOrigin.clear();
    mPointListOrigin = CreateGaussPointList(mrOriginMainModelPart);

    // Check the NODAL_H
    NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
    for (auto& i_node : nodes_array)
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
        ElementsArrayType& elements_array = mrDestinationMainModelPart.Elements();
        const int num_elements = static_cast<int>(elements_array.size());

//         #pragma omp for firstprivate(this_integration_method)
        for(int i = 0; i < num_elements; ++i)
        {
            auto it_elem = elements_array.begin() + i;

            // Getting the geometry
            GeometryType& r_this_geometry = it_elem->GetGeometry();

            // Getting the integration points
            this_integration_method = it_elem->GetIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const std::size_t integration_points_number = integration_points.size();

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

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
                            const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                            if (p_destination_cl->Has(this_var)) {
                                GetAndSetWeightedVariableOnConstitutiveLaw(this_var, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, current_process_info);
                            } else {
                                GetAndSetWeightedVariableOnElements(this_var, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, current_process_info);
                            }
                        } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                            const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                            if (p_destination_cl->Has(this_var)) {
                                GetAndSetWeightedVariableOnConstitutiveLaw(this_var, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, current_process_info);
                            } else {
                                GetAndSetWeightedVariableOnElements(this_var, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, current_process_info);
                            }
                        } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                            const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                            if (p_destination_cl->Has(this_var)) {
                                GetAndSetWeightedVariableOnConstitutiveLaw(this_var, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, current_process_info);
                            } else {
                                GetAndSetWeightedVariableOnElements(this_var, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, current_process_info);
                            }
                        } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                            const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                            if (p_destination_cl->Has(this_var)) {
                                GetAndSetWeightedVariableOnConstitutiveLaw(this_var, number_points_found, points_found, point_distances, characteristic_length, p_destination_cl, current_process_info);
                            } else {
                                GetAndSetWeightedVariableOnElements(this_var, number_points_found, points_found, point_distances, characteristic_length, it_elem, i_gauss_point, current_process_info);
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
    //}
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsShapeFunctionTransfer()
{
    // Initialize some values
    GeometryData::IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_1;

    // Iterate in the nodes to initialize the values
    NodesArrayType& nodes_array = mrOriginMainModelPart.Nodes();

    /* Nodes */
    for (auto& variable_name : mInternalVariableList) {
        if (KratosComponents<DoubleVarType>::Has(variable_name)) {
            const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
            VariableUtils().SetNonHistoricalVariable(this_var, this_var.Zero(), nodes_array);
        } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
            const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
            VariableUtils().SetNonHistoricalVariable(this_var, this_var.Zero(), nodes_array);
        } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
            const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
            VariableUtils().SetNonHistoricalVariable(this_var, this_var.Zero(), nodes_array);
        } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
            const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
            VariableUtils().SetNonHistoricalVariable(this_var, this_var.Zero(), nodes_array);
        } else {
            KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
        }
    }

    // Iterate in the elements to ponderate the values
    ElementsArrayType& elements_array = mrOriginMainModelPart.Elements();
    int num_elements = static_cast<int>(elements_array.size());

    const ProcessInfo& origin_process_info = mrOriginMainModelPart.GetProcessInfo();

    /* Elements */
    #pragma omp parallel for firstprivate(this_integration_method)
    for(int i = 0; i < num_elements; ++i) {
        auto it_elem = elements_array.begin() + i;

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
        it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,origin_process_info);

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
                    const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                    if (p_origin_cl->Has(this_var)) {
                        InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, this_var, N, p_origin_cl, weight);
                    } else {
                        InterpolateAddVariableOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, weight, origin_process_info);
                    }
                } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                    const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                    if (p_origin_cl->Has(this_var)) {
                        InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, this_var, N, p_origin_cl, weight);
                    } else {
                        InterpolateAddVariableOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, weight, origin_process_info);
                    }
                } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                    const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                    if (p_origin_cl->Has(this_var)) {
                        InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, this_var, N, p_origin_cl, weight);
                    } else {
                        InterpolateAddVariableOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, weight, origin_process_info);
                    }
                } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                    const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                    if (p_origin_cl->Has(this_var)) {
                        InterpolateAddVariableOnConstitutiveLaw(r_this_geometry, this_var, N, p_origin_cl, weight);
                    } else {
                        InterpolateAddVariableOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, weight, origin_process_info);
                    }
                } else {
                    KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                }
            }
        }

        // We divide by the total weight
        for (auto& variable_name : mInternalVariableList) {
            if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                PonderateVariable(r_this_geometry, this_var, total_weight);
            } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                PonderateVariable(r_this_geometry, this_var, total_weight);
            } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                PonderateVariable(r_this_geometry, this_var, total_weight);
            } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                PonderateVariable(r_this_geometry, this_var, total_weight);
            } else {
                KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
            }
        }
    }

    // We interpolate to the new nodes
    if (mDimension == 2) {
        // We create the locator
        BinBasedFastPointLocator<2> point_locator(mrOriginMainModelPart);
        point_locator.UpdateSearchDatabase();

        // Iterate over nodes
        NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        const auto it_node_begin = nodes_array.begin();

        /* Nodes */
        #pragma omp parallel for firstprivate(point_locator)
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;

            Vector N;
            Element::Pointer p_element;

            const bool found = point_locator.FindPointOnMeshSimplified(it_node->Coordinates(), N, p_element, mAllocationSize);

            if (found == false) {
                KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING: GP not found (interpolation not posible)" << "\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << std::endl;
            } else {
                for (auto& variable_name : mInternalVariableList) {
                    if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                        const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                        const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                        const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                        const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                    }
                }
            }
        }
    } else {
        // We create the locator
        BinBasedFastPointLocator<3> point_locator(mrOriginMainModelPart);
        point_locator.UpdateSearchDatabase();

        // Iterate over nodes
        NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());
        const auto it_node_begin = nodes_array.begin();

        /* Nodes */
        #pragma omp parallel for firstprivate(point_locator)
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = it_node_begin + i;

            Vector N;
            Element::Pointer p_element;

            const bool found = point_locator.FindPointOnMeshSimplified(it_node->Coordinates(), N, p_element, mAllocationSize);

            if (found == false) {
                KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING: Node "<< it_node->Id() << " not found (interpolation not posible)" <<  "\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << "\t Z:"<< it_node->Z() << std::endl;
            } else {
                for (auto& variable_name : mInternalVariableList) {
                    if (KratosComponents<DoubleVarType>::Has(variable_name)) {
                        const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                        const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                        const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                        const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                        InterpolateToNode(this_var, N, (*it_node.base()), p_element);
                    } else {
                        KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                    }
                }
            }
        }
    }

    // Finally we interpolate to the new GP
    ElementsArrayType& elements_array_destination = mrDestinationMainModelPart.Elements();
    num_elements = static_cast<int>(elements_array_destination.size());

    const ProcessInfo& destination_process_info = mrOriginMainModelPart.GetProcessInfo();

    /* Elements */
    #pragma omp parallel for firstprivate(this_integration_method)
    for(int i = 0; i < num_elements; ++i) {
        auto it_elem = elements_array_destination.begin() + i;

        // Getting the geometry
        GeometryType& r_this_geometry = it_elem->GetGeometry();

        // Getting the integration points
        this_integration_method = it_elem->GetIntegrationMethod();
        const GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
        const std::size_t integration_points_number = integration_points.size();

        // Getting the CL
        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
        it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_law_vector,destination_process_info);

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
                    const DoubleVarType& this_var = KratosComponents<DoubleVarType>::Get(variable_name);
                    if (p_destination_cl->Has(this_var)) {
                        SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, this_var, N, p_destination_cl, destination_process_info);
                    } else {
                        SetInterpolatedValueOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, destination_process_info);
                    }
                } else if (KratosComponents<ArrayVarType>::Has(variable_name)) {
                    const ArrayVarType& this_var = KratosComponents<ArrayVarType>::Get(variable_name);
                    if (p_destination_cl->Has(this_var)) {
                        SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, this_var, N, p_destination_cl, destination_process_info);
                    } else {
                        SetInterpolatedValueOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, destination_process_info);
                    }
                } else if (KratosComponents<VectorVarType>::Has(variable_name)) {
                    const VectorVarType& this_var = KratosComponents<VectorVarType>::Get(variable_name);
                    if (p_destination_cl->Has(this_var)) {
                        SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, this_var, N, p_destination_cl, destination_process_info);
                    } else {
                        SetInterpolatedValueOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, destination_process_info);
                    }
                } else if (KratosComponents<MatrixVarType>::Has(variable_name)) {
                    const MatrixVarType& this_var = KratosComponents<MatrixVarType>::Get(variable_name);
                    if (p_destination_cl->Has(this_var)) {
                        SetInterpolatedValueOnConstitutiveLaw(r_this_geometry, this_var, N, p_destination_cl, destination_process_info);
                    } else {
                        SetInterpolatedValueOnElement(r_this_geometry, this_var, N, it_elem, i_gauss_point, destination_process_info);
                    }
                } else {
                    KRATOS_WARNING("InternalVariablesInterpolationProcess") << "WARNING:: " << variable_name << " is not registered as any type of compatible variable: DOUBLE or ARRAY_1D or VECTOR or Matrix" << std::endl;
                }
            }
        }
    }
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

}  // namespace Kratos.
