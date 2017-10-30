// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes

// Project includes
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
    Parameters DefaultParameters = Parameters(R"(
        {
            "allocation_size"                      : 1000,
            "bucket_size"                          : 4,
            "search_factor"                        : 2,
            "interpolation_type"                   : "LST",
            "internal_variable_interpolation_list" :[]
        })" );

    ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

    mAllocationSize = ThisParameters["allocation_size"].GetInt();
    mBucketSize = ThisParameters["bucket_size"].GetInt();
    mSearchFactor = ThisParameters["search_factor"].GetDouble();
    mThisInterpolationType = ConvertInter(ThisParameters["interpolation_type"].GetString());

    if (ThisParameters["internal_variable_interpolation_list"].IsArray() == true)
    {
        auto variable_array_list = ThisParameters["internal_variable_interpolation_list"];

        for (unsigned int i_var = 0; i_var < variable_array_list.size(); ++i_var)
        {
            mInternalVariableList.push_back(KratosComponents<Variable<double>>::Get(variable_array_list[i_var].GetString()));
        }
    }
    else
    {
        std::cout << "WARNING:: No variables to interpolate, look that internal_variable_interpolation_list is correctly defined in your parameters" << std::endl;
        mInternalVariableList.clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::Execute()
{
    /** NOTE: There are mainly two ways to interpolate the internal variables (there are three, but just two are behave correctly)
    * CPT: Closest point transfer. It transfer the values from the closest GP
    * LST: Least-square projection transfer. It transfers from the closest GP from the old mesh
    * SFT: It transfer GP values to the nodes in the old mesh and then interpolate to the new mesh using the sahpe functions all the time (NOTE: THIS DOESN"T WORK, AND REQUIRES EXTRA STORE)
    */

    if (mThisInterpolationType == CPT && mInternalVariableList.size() > 0)
    {
        InterpolateGaussPointsCPT();
    }
    else if (mThisInterpolationType == LST && mInternalVariableList.size() > 0)
    {
        InterpolateGaussPointsLST();
    }
    else if (mThisInterpolationType == SFT && mInternalVariableList.size() > 0)
    {
        InterpolateGaussPointsSFT();
    }
}

/***********************************************************************************/
/***********************************************************************************/

PointVector InternalVariablesInterpolationProcess::CreateGaussPointList(ModelPart& ThisModelPart)
{
    PointVector this_point_vector;

    GeometryData::IntegrationMethod this_integration_method;

    // Iterate in the elements
    ElementsArrayType& elements_array = ThisModelPart.Elements();
    int num_elements = ThisModelPart.NumberOfElements();

    const ProcessInfo& current_process_info = ThisModelPart.GetProcessInfo();

    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<PointVector> points_buffer(num_threads);

    #pragma omp parallel
    {
        const int id = OpenMPUtils::ThisThread();

        #pragma omp for
        for(int i = 0; i < num_elements; ++i)
        {
            auto it_elem = elements_array.begin() + i;

            // Getting the geometry
            Element::GeometryType& r_this_geometry = it_elem->GetGeometry();

            // Getting the integration points
            this_integration_method = it_elem->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const unsigned int integration_points_number = integration_points.size();

            // Computing the Jacobian
            Vector vector_det_j(integration_points_number);
            r_this_geometry.DeterminantOfJacobian(vector_det_j,this_integration_method);

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

            for (unsigned int i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point )
            {
                const array_1d<double, 3> local_coordinates = integration_points[i_gauss_point].Coordinates();

                // We compute the corresponding weight
                const double weight = vector_det_j[i_gauss_point] * integration_points[i_gauss_point].Weight();

                // We compute the global coordinates
                array_1d<double, 3> global_coordinates;
                global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                // We create the respective GP
                PointTypePointer p_point = PointTypePointer(new PointType(global_coordinates, constitutive_law_vector[i_gauss_point], weight));
                (points_buffer[id]).push_back(p_point);
            }
        }

        // Combine buffers together
        #pragma omp single
        {
            for( auto& point_buffer : points_buffer)
            {
                std::move(point_buffer.begin(),point_buffer.end(),back_inserter(this_point_vector));
            }
        }
    }

    return this_point_vector;
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsCPT()
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
        for(int i = 0; i < num_elements; ++i)
        {
            auto it_elem = elements_array.begin() + i;

            // Getting the geometry
            Element::GeometryType& r_this_geometry = it_elem->GetGeometry();

            // Getting the integration points
            this_integration_method = it_elem->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const unsigned int integration_points_number = integration_points.size();

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

            for (unsigned int i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point )
            {
                // We compute the global coordinates
                const array_1d<double, 3> local_coordinates = integration_points[i_gauss_point].Coordinates();
                array_1d<double, 3> global_coordinates;
                global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                PointTypePointer p_gp_origin = tree_points.SearchNearestPoint(global_coordinates);

                for (auto this_var : mInternalVariableList)
                {
                    double origin_value;
                    origin_value = (p_gp_origin->GetConstitutiveLaw())->GetValue(this_var, origin_value);

                    (constitutive_law_vector[i_gauss_point])->SetValue(this_var, origin_value, current_process_info);
                }
            }
        }
    //}
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsLST()
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

        // Initialize values
        PointVector points_found(mAllocationSize);
        std::vector<double> point_distnaces(mAllocationSize);
        unsigned int number_points_found = 0;

        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTree tree_points(mPointListOrigin.begin(), mPointListOrigin.end(), mBucketSize);

        // Iterate over the destination elements
        ElementsArrayType& elements_array = mrDestinationMainModelPart.Elements();
        auto num_elements = elements_array.end() - elements_array.begin();

        //#pragma omp for
        for(int i = 0; i < num_elements; ++i)
        {
            auto it_elem = elements_array.begin() + i;

            // Getting the geometry
            Element::GeometryType& r_this_geometry = it_elem->GetGeometry();

            // Getting the integration points
            this_integration_method = it_elem->GetIntegrationMethod();
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
            const unsigned int integration_points_number = integration_points.size();

            // Getting the CL
            std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
            it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,current_process_info);

            // Computing the radius
            const double radius = mSearchFactor *  (mDimension == 2 ? std::sqrt(r_this_geometry.Area()) : std::cbrt(r_this_geometry.Volume()));

            // We get the NODAL_H vector
            Vector nodal_h_vector(r_this_geometry.size());
            for (unsigned int i_node = 0; i_node < r_this_geometry.size(); ++i_node)
            {
                if ( r_this_geometry[i_node].SolutionStepsDataHas( NODAL_H ) == false )
                {
                    KRATOS_ERROR << "NODAL_H is not defined in the node ID: " << r_this_geometry[i_node].Id() << std::endl;
                }

                nodal_h_vector[i_node] = r_this_geometry[i_node].FastGetSolutionStepValue(NODAL_H);
            }

            for (unsigned int i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point )
            {
                // We compute the global coordinates
                const array_1d<double, 3> local_coordinates = integration_points[i_gauss_point].Coordinates();
                array_1d<double, 3> global_coordinates;
                global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

                // We compute the pondered characteristic length
                Vector N( r_this_geometry.size() );
                r_this_geometry.ShapeFunctionsValues( N, local_coordinates );
                const double characteristic_length = inner_prod(N, nodal_h_vector);

                number_points_found = tree_points.SearchInRadius(global_coordinates, radius, points_found.begin(), point_distnaces.begin(), mAllocationSize);

                if (number_points_found > 0)
                {
                    for (auto this_var : mInternalVariableList)
                    {
                        double weighting_function_numerator   = 0.0;
                        double weighting_function_denominator = 0.0;
                        double origin_value;

                        for (unsigned int i_point_found = 0; i_point_found < number_points_found; ++i_point_found)
                        {
                            PointTypePointer p_gp_origin = points_found[i_point_found];

                            const double distance = point_distnaces[i_point_found];

                            origin_value = (p_gp_origin->GetConstitutiveLaw())->GetValue(this_var, origin_value);

                            const double ponderated_weight = p_gp_origin->GetWeight() * std::exp( -4.0 * distance * distance /(characteristic_length * characteristic_length));

                            weighting_function_numerator   += ponderated_weight * origin_value;
                            weighting_function_denominator += ponderated_weight;
                        }

                        const double destination_value = weighting_function_numerator/weighting_function_denominator;

                        (constitutive_law_vector[i_gauss_point])->SetValue(this_var, destination_value, current_process_info);
                    }
                }
                else
                {
                    std::cout << "WARNING:: It wasn't impossible to find any Gauss Point from where interpolate the internal variables" << std::endl;
                }
            }
        }
    //}
}

/***********************************************************************************/
/***********************************************************************************/

void InternalVariablesInterpolationProcess::InterpolateGaussPointsSFT()
{
    // Initialize some values
    GeometryData::IntegrationMethod this_integration_method;

    // Iterate in the nodes to initialize the values
    NodesArrayType& nodes_array = mrOriginMainModelPart.Nodes();
    auto num_nodes = nodes_array.end() - nodes_array.begin();

    /* Nodes */
    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)
    {
        auto it_node = nodes_array.begin() + i;

        for (auto this_var : mInternalVariableList)
        {
            it_node->SetValue(this_var, 0.0);
        }
    }

    // Iterate in the elements to ponderate the values
    ElementsArrayType& elements_array = mrOriginMainModelPart.Elements();
    auto num_elements = elements_array.end() - elements_array.begin();

    const ProcessInfo& origin_process_info = mrOriginMainModelPart.GetProcessInfo();

    /* Elements */
    #pragma omp parallel for
    for(int i = 0; i < num_elements; ++i)
    {
        auto it_elem = elements_array.begin() + i;

        // Getting the geometry
        Element::GeometryType& r_this_geometry = it_elem->GetGeometry();

        // Getting the integration points
        this_integration_method = it_elem->GetIntegrationMethod();
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
        const unsigned int integration_points_number = integration_points.size();

        // Computing the Jacobian
        Vector vector_det_j(integration_points_number);
        r_this_geometry.DeterminantOfJacobian(vector_det_j,this_integration_method);

        // Getting the CL
        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
        it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,origin_process_info);

        // We initialize the total weigth
        double total_weight = 0.0;

        for (unsigned int i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point )
        {
            const array_1d<double, 3> local_coordinates = integration_points[i_gauss_point].Coordinates();

            // We compute the corresponding weight
            const double weight = vector_det_j[i_gauss_point] * integration_points[i_gauss_point].Weight();
            total_weight += weight;

            // We compute the pondered characteristic length
            Vector N( r_this_geometry.size() );
            r_this_geometry.ShapeFunctionsValues( N, local_coordinates );

            // We compute the global coordinates
            array_1d<double, 3> global_coordinates;
            global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

            for (auto this_var : mInternalVariableList)
            {
                double origin_value;
                origin_value = constitutive_law_vector[i_gauss_point]->GetValue(this_var, origin_value);

                // We sum all the contributions
                for (unsigned int i_node = 0; i_node < r_this_geometry.size(); ++i_node)
                {
                    #pragma omp atomic
                    r_this_geometry[i_node].GetValue(this_var) += N[i_node] * origin_value * weight;
                }
            }
        }

        // We divide by the total weight
        for (auto this_var : mInternalVariableList)
        {
            for (unsigned int i_node = 0; i_node < r_this_geometry.size(); ++i_node)
            {
                #pragma omp critical
                r_this_geometry[i_node].GetValue(this_var) /= total_weight;
            }
        }
    }

    // We interpolate to the new nodes
    if (mDimension == 2)
    {
        // We create the locator
        BinBasedFastPointLocator<2> point_locator = BinBasedFastPointLocator<2>(mrOriginMainModelPart);
        point_locator.UpdateSearchDatabase();

        // Iterate in the nodes
        NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
        auto num_nodes = nodes_array.end() - nodes_array.begin();

        /* Nodes */
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = nodes_array.begin() + i;

            Vector N;
            Element::Pointer p_element;

            const bool found = point_locator.FindPointOnMeshSimplified(it_node->Coordinates(), N, p_element, mAllocationSize);

            if (found == false)
            {
                std::cout << "WARNING: GP not found (interpolation not posible)" << std::endl;
                std::cout << "\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << std::endl;
            }
            else
            {
                for (auto this_var : mInternalVariableList)
                {
                    Vector values(p_element->GetGeometry().size());

                    for (unsigned int i_node = 0; i_node < p_element->GetGeometry().size(); ++i_node)
                    {
                        values[i_node] = p_element->GetGeometry()[i_node].GetValue(this_var);
                    }

                    it_node->GetValue(this_var) = inner_prod(values, N);
                }
            }
        }
    }
    else
    {
        // We create the locator
        BinBasedFastPointLocator<3> point_locator = BinBasedFastPointLocator<3>(mrOriginMainModelPart);
        point_locator.UpdateSearchDatabase();

        // Iterate in the nodes
        NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
        auto num_nodes = nodes_array.end() - nodes_array.begin();

        /* Nodes */
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i)
        {
            auto it_node = nodes_array.begin() + i;

            Vector N;
            Element::Pointer p_element;

            const bool found = point_locator.FindPointOnMeshSimplified(it_node->Coordinates(), N, p_element, mAllocationSize);

            if (found == false)
            {
                std::cout << "WARNING: Node "<< it_node->Id() << " not found (interpolation not posible)" << std::endl;
                std::cout << "\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << "\t Z:"<< it_node->Z() << std::endl;
            }
            else
            {
                for (auto this_var : mInternalVariableList)
                {
                    Vector values(p_element->GetGeometry().size());

                    for (unsigned int i_node = 0; i_node < p_element->GetGeometry().size(); ++i_node)
                    {
                        values[i_node] = p_element->GetGeometry()[i_node].GetValue(this_var);
                    }

                    it_node->GetValue(this_var) = inner_prod(values, N);
                }
            }
        }
    }

    // Finally we interpolate to the new GP
    ElementsArrayType& elements_arrayDestination = mrDestinationMainModelPart.Elements();
    num_elements = elements_arrayDestination.end() - elements_arrayDestination.begin();

    const ProcessInfo& destination_process_info = mrOriginMainModelPart.GetProcessInfo();

    /* Elements */
    #pragma omp parallel for
    for(int i = 0; i < num_elements; ++i)
    {
        auto it_elem = elements_arrayDestination.begin() + i;

        // Getting the geometry
        Element::GeometryType& r_this_geometry = it_elem->GetGeometry();

        // Getting the integration points
        this_integration_method = it_elem->GetIntegrationMethod();
        const Element::GeometryType::IntegrationPointsArrayType& integration_points = r_this_geometry.IntegrationPoints(this_integration_method);
        const unsigned int integration_points_number = integration_points.size();

        // Getting the CL
        std::vector<ConstitutiveLaw::Pointer> constitutive_law_vector(integration_points_number);
        it_elem->GetValueOnIntegrationPoints(CONSTITUTIVE_LAW,constitutive_law_vector,destination_process_info);

        for (unsigned int i_gauss_point = 0; i_gauss_point < integration_points_number; ++i_gauss_point )
        {
            array_1d<double, 3> local_coordinates = integration_points[i_gauss_point].Coordinates();

            // We compute the pondered characteristic length
            Vector N( r_this_geometry.size() );
            r_this_geometry.ShapeFunctionsValues( N, local_coordinates );

            // We compute the global coordinates
            array_1d<double, 3> global_coordinates;
            global_coordinates = r_this_geometry.GlobalCoordinates( global_coordinates, local_coordinates );

            Vector values(r_this_geometry.size() );

            for (auto this_var : mInternalVariableList)
            {
                for (unsigned int i_node = 0; i_node < r_this_geometry.size(); ++i_node)
                {
                    values[i_node] = r_this_geometry[i_node].GetValue(this_var);
                }

                const double destination_value = inner_prod(values, N);

                constitutive_law_vector[i_gauss_point]->SetValue(this_var, destination_value, destination_process_info);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

InterpolationTypes InternalVariablesInterpolationProcess::ConvertInter(const std::string& str)
{
    if(str == "CPT")
    {
        return CPT;
    }
    else if(str == "LST")
    {
        return LST;
    }
    else if(str == "SFT")
    {
        return SFT;
    }
    else
    {
        return LST;
    }
}

}  // namespace Kratos.
