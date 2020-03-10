//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/simple_mortar_mapper_process.h"

/* Custom utilities */
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::AVERAGE_NORMAL(Kratos::Flags::Create(0));
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::DISCONTINOUS_INTERFACE(Kratos::Flags::Create(1));
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ORIGIN_IS_HISTORICAL(Kratos::Flags::Create(2));
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::DESTINATION_IS_HISTORICAL(Kratos::Flags::Create(3));
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ORIGIN_SKIN_IS_CONDITION_BASED(Kratos::Flags::Create(4));
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::DESTINATION_SKIN_IS_CONDITION_BASED(Kratos::Flags::Create(5));

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::SimpleMortarMapperProcess(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    Parameters ThisParameters,
    LinearSolverType::Pointer pThisLinearSolver
    ):  mOriginModelPart(rOriginModelPart),
        mDestinationModelPart(rDestinationModelPart),
        mpOriginVariable(&(KratosComponents<TVarType>::Get(ThisParameters["origin_variable"].GetString()))),
        mpDestinationVariable((ThisParameters["destination_variable"].GetString() == "") ? mpOriginVariable : &(KratosComponents<TVarType>::Get(ThisParameters["destination_variable"].GetString()))),
        mThisParameters(ThisParameters),
        mpThisLinearSolver(pThisLinearSolver)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Setting flags
    mOptions.Set(AVERAGE_NORMAL, mThisParameters["using_average_nodal_normal"].GetBool());
    mOptions.Set(DISCONTINOUS_INTERFACE, mThisParameters["discontinuous_interface"].GetBool());
    mOptions.Set(ORIGIN_IS_HISTORICAL, mThisParameters["origin_variable_historical"].GetBool());
    mOptions.Set(DESTINATION_IS_HISTORICAL, mThisParameters["destination_variable_historical"].GetBool());
    mOptions.Set(ORIGIN_SKIN_IS_CONDITION_BASED, mThisParameters["origin_are_conditions"].GetBool());
    mOptions.Set(DESTINATION_SKIN_IS_CONDITION_BASED, mThisParameters["destination_are_conditions"].GetBool());

    // We set some values
    mMappingCoefficient = mThisParameters["mapping_coefficient"].GetDouble();
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::SimpleMortarMapperProcess(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    TVarType& rThisVariable,
    Parameters ThisParameters,
    LinearSolverType::Pointer pThisLinearSolver
    ):  mOriginModelPart(rOriginModelPart),
        mDestinationModelPart(rDestinationModelPart),
        mpOriginVariable(&rThisVariable),
        mpDestinationVariable(&rThisVariable),
        mThisParameters(ThisParameters),
        mpThisLinearSolver(pThisLinearSolver)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Setting flags
    mOptions.Set(AVERAGE_NORMAL, mThisParameters["using_average_nodal_normal"].GetBool());
    mOptions.Set(DISCONTINOUS_INTERFACE, mThisParameters["discontinuous_interface"].GetBool());
    mOptions.Set(ORIGIN_IS_HISTORICAL, mThisParameters["origin_variable_historical"].GetBool());
    mOptions.Set(DESTINATION_IS_HISTORICAL, mThisParameters["destination_variable_historical"].GetBool());
    mOptions.Set(ORIGIN_SKIN_IS_CONDITION_BASED, mThisParameters["origin_are_conditions"].GetBool());
    mOptions.Set(DESTINATION_SKIN_IS_CONDITION_BASED, mThisParameters["destination_are_conditions"].GetBool());

    // We set some values
    mMappingCoefficient = mThisParameters["mapping_coefficient"].GetDouble();
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::SimpleMortarMapperProcess(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    TVarType& rOriginVariable,
    TVarType& rDestinationVariable,
    Parameters ThisParameters,
    LinearSolverType::Pointer pThisLinearSolver
    ): mOriginModelPart(rOriginModelPart),
       mDestinationModelPart(rDestinationModelPart),
       mpOriginVariable(&rOriginVariable),
       mpDestinationVariable(&rDestinationVariable),
       mThisParameters(ThisParameters),
       mpThisLinearSolver(pThisLinearSolver)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Setting flags
    mOptions.Set(AVERAGE_NORMAL, mThisParameters["using_average_nodal_normal"].GetBool());
    mOptions.Set(DISCONTINOUS_INTERFACE, mThisParameters["discontinuous_interface"].GetBool());
    mOptions.Set(ORIGIN_IS_HISTORICAL, mThisParameters["origin_variable_historical"].GetBool());
    mOptions.Set(DESTINATION_IS_HISTORICAL, mThisParameters["destination_variable_historical"].GetBool());
    mOptions.Set(ORIGIN_SKIN_IS_CONDITION_BASED, mThisParameters["origin_are_conditions"].GetBool());
    mOptions.Set(DESTINATION_SKIN_IS_CONDITION_BASED, mThisParameters["destination_are_conditions"].GetBool());

    // We set some values
    mMappingCoefficient = mThisParameters["mapping_coefficient"].GetDouble();
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>:: Execute()
{
    KRATOS_TRY;

    ExecuteInitializeSolutionStep();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>:: ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;

    // We reset the database if needed
    const bool update_interface = mThisParameters["update_interface"].GetBool();
    if (update_interface) {
        this->UpdateInterface();
    }

    if (mpThisLinearSolver == nullptr) {
        ExecuteExplicitMapping();
    } else {
        ExecuteImplicitMapping();
    }

    // We apply the coeffcient if different of one
    if (mMappingCoefficient != 1.0) {
        auto& r_nodes_array = mDestinationModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();

        if (mOptions.Is(DESTINATION_IS_HISTORICAL)) {
            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (r_nodes_array.size()); ++k) {
                auto it_node = it_node_begin + k;
                it_node->FastGetSolutionStepValue(*mpDestinationVariable) *= mMappingCoefficient;
            }
        } else {
            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (r_nodes_array.size()); ++k) {
                auto it_node = it_node_begin + k;
                it_node->GetValue(*mpDestinationVariable) *= mMappingCoefficient;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CheckAndPerformSearch()
{
    // If we froce to update the database
    const bool update_interface = mThisParameters["update_interface"].GetBool();

    // First we check if search already exists
    bool search_exists = true;

    // Iterate in the conditions and elements
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        auto& r_destination_conditions_array = mDestinationModelPart.Conditions();
        const auto it_cond_begin = r_destination_conditions_array.begin();
        for(IndexType i = 0; i < r_destination_conditions_array.size(); ++i) {
            auto it_cond = it_cond_begin + i;

            if (!(it_cond->Has( INDEX_SET ) || it_cond->Has( INDEX_MAP ))) {
                search_exists = false;
                break;
            }
        }

        if (!search_exists) {
            // We create the variable INDEX_SET
            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(r_destination_conditions_array.size()); ++i) {
                auto it_cond = it_cond_begin + i;
                if (!it_cond->Has(INDEX_SET)) {
                    it_cond->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
                }
            }
        }
    } else {
        auto& r_destination_elements_array = mDestinationModelPart.Elements();
        const auto it_elem_begin = r_destination_elements_array.begin();
        for(IndexType i = 0; i < r_destination_elements_array.size(); ++i) {
            auto it_elem = it_elem_begin + i;

            if (!(it_elem->Has( INDEX_SET ) || it_elem->Has( INDEX_MAP ))) {
                search_exists = false;
                break;
            }
        }

        if (!search_exists) {
            // We create the variable INDEX_SET
            #pragma omp parallel for
            for(int i = 0; i < static_cast<int>(r_destination_elements_array.size()); ++i) {
                auto it_elem = it_elem_begin + i;
                if (!it_elem->Has(INDEX_SET))
                    it_elem->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
            }
        }
    }

    // Now we perform the corresponding search
    if (!search_exists || update_interface) {
        // A list that contents the all the points (from nodes) from the modelpart
        PointVector point_list_destination;
        point_list_destination.clear();

        if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
            // Iterate in the conditions
            auto& r_origin_conditions_array = mOriginModelPart.Conditions();
            const auto it_cond_begin = r_origin_conditions_array.begin();

            #pragma omp parallel
            {
                PointVector points_buffer;

                #pragma omp for schedule(guided, 512) nowait
                for(int i = 0; i < static_cast<int>(r_origin_conditions_array.size()); ++i) {
                    auto it_cond = it_cond_begin + i;

                    const PointTypePointer p_point = Kratos::make_unique<PointMapperType>(*it_cond.base());
                    (points_buffer).push_back(p_point);
                }

                // Combine buffers together
                #pragma omp critical
                {
                    std::move(points_buffer.begin(),points_buffer.end(),back_inserter(point_list_destination));
                }
            }
        } else {
            // Iterate in the elements
            auto& r_origin_elements_array = mOriginModelPart.Elements();
            const auto it_elem_begin = r_origin_elements_array.begin();

            #pragma omp parallel
            {
                PointVector points_buffer;

                #pragma omp for schedule(guided, 512) nowait
                for(int i = 0; i < static_cast<int>(r_origin_elements_array.size()); ++i) {
                    auto it_elem = it_elem_begin + i;

                    const PointTypePointer p_point = Kratos::make_unique<PointMapperType>(*it_elem.base());
                    (points_buffer).push_back(p_point);
                }

                // Combine buffers together
                #pragma omp critical
                {
                    std::move(points_buffer.begin(),points_buffer.end(),back_inserter(point_list_destination));
                }
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(point_list_destination.size()); ++i) {
            point_list_destination[i]->UpdatePoint();
        }

        // Some auxiliar values
        const SizeType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt(); // Allocation size for the vectors and max number of potential results
        const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble(); // The search factor to be considered
        SizeType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt(); // Bucket size for kd-tree

        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTreeType tree_points(point_list_destination.begin(), point_list_destination.end(), bucket_size);

        if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
            auto& r_destination_conditions_array = mDestinationModelPart.Conditions();
            const auto it_cond_begin = r_destination_conditions_array.begin();

            // Iterate over conditions
            for(IndexType i = 0; i < r_destination_conditions_array.size(); ++i) {
                auto it_cond = it_cond_begin + i;

                FillDatabase<Condition>((*it_cond.base()), tree_points, allocation_size, search_factor);
            }
        } else {
            auto& r_destination_elements_array = mDestinationModelPart.Elements();
            const auto it_elem_begin = r_destination_elements_array.begin();

            // Iterate over elements
            for(IndexType i = 0; i < r_destination_elements_array.size(); ++i) {
                auto it_elem = it_elem_begin + i;

                FillDatabase<Element>((*it_elem.base()), tree_points, allocation_size, search_factor);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ResetNodalArea()
{
    // Iterating over nodes
    auto& r_nodes_array = mDestinationModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_nodes_array.size()); ++i) {
        auto it_node = it_node_begin + i;
        it_node->SetValue(NODAL_AREA, 0.0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
double SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetReferenceArea()
{
    double ref_area = 0.0;

    // We look for the max area in the origin model part
    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
        auto& r_conditions_array_origin = mOriginModelPart.Conditions();

        double current_area;
        for(int i = 0; i < static_cast<int>(r_conditions_array_origin.size()); ++i) {
            auto it_cond = r_conditions_array_origin.begin() + i;
            current_area = it_cond->GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    } else {
        auto& r_elements_array_origin = mOriginModelPart.Elements();

        double current_area;
        for(int i = 0; i < static_cast<int>(r_elements_array_origin.size()); ++i) {
            auto it_elem = r_elements_array_origin.begin() + i;
            current_area = it_elem->GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    }

    // We look for the max area in the destination model part
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        auto& r_conditions_array_destination = mDestinationModelPart.Conditions();

        double current_area;
        for(int i = 0; i < static_cast<int>(r_conditions_array_destination.size()); ++i) {
            auto it_cond = r_conditions_array_destination.begin() + i;
            current_area = it_cond->GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    } else {
        auto& r_elements_array_destination = mDestinationModelPart.Elements();

        double current_area;
        for(int i = 0; i < static_cast<int>(r_elements_array_destination.size()); ++i) {
            auto it_elem = r_elements_array_destination.begin() + i;
            current_area = it_elem->GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    }

    return ref_area;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::AssemblyMortarOperators(
    const std::vector<array_1d<PointType,TDim>>& rGeometricalObjectsPointSlave,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const array_1d<double, 3>& rMasterNormal,
    MortarKinematicVariablesType& rThisKinematicVariables,
    MortarOperatorType& rThisMortarOperators,
    const IntegrationMethod& rThisIntegrationMethod,
    const BoundedMatrixType Ae
    )
{
    for (IndexType i_geom = 0; i_geom < rGeometricalObjectsPointSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            rSlaveGeometry.GlobalCoordinates(global_point, rGeometricalObjectsPointSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompositionType decomp_geom( PointerVector<PointType>{points_array} );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, rSlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (!bad_shape) {
            const GeometryType::IntegrationPointsArrayType& r_integration_points_slave = decomp_geom.IntegrationPoints( rThisIntegrationMethod );

            // Integrating the mortar operators
            PointType local_point_parent, gp_global, projected_gp_global;
            array_1d<double,3> gp_normal;
            GeometryType::CoordinatesArrayType projected_gp_local;
            for ( IndexType point_number = 0; point_number < r_integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp{r_integration_points_slave[point_number].Coordinates()};

                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                rSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

                /// SLAVE CONDITION ///
                rSlaveGeometry.ShapeFunctionsValues( rThisKinematicVariables.NSlave, local_point_parent.Coordinates() );
                rThisKinematicVariables.PhiLagrangeMultipliers = prod(Ae, rThisKinematicVariables.NSlave);

                rThisKinematicVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                /// MASTER CONDITION ///
                if (mOptions.Is(AVERAGE_NORMAL)) {
                    noalias(gp_normal) = MortarUtilities::GaussPointUnitNormal(rThisKinematicVariables.NSlave, rSlaveGeometry);
                } else {
                    noalias(gp_normal) = rSlaveGeometry.UnitNormal(local_point_decomp);
                }

                GeometryType::CoordinatesArrayType slave_gp_global;
                rSlaveGeometry.GlobalCoordinates( slave_gp_global, local_point_parent );
                GeometricalProjectionUtilities::FastProjectDirection( rMasterGeometry, gp_global, projected_gp_global, rMasterNormal, -gp_normal ); // The opposite direction

                rMasterGeometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                // SHAPE FUNCTIONS
                rMasterGeometry.ShapeFunctionsValues( rThisKinematicVariables.NMaster, projected_gp_local );

                const double integration_weight = r_integration_points_slave[point_number].Weight();

                rThisMortarOperators.CalculateMortarOperators(rThisKinematicVariables, integration_weight);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
inline BoundedMatrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CalculateAe(
    GeometryType& rSlaveGeometry,
    MortarKinematicVariablesType& rThisKinematicVariables,
    std::vector<array_1d<PointType,TDim>>& rGeometricalObjectsPointSlave,
    const IntegrationMethod& rThisIntegrationMethod
    )
{
    // We initilize the Ae components
    DualLagrangeMultiplierOperatorsType Ae_data;
    Ae_data.Initialize();

    // Initialize general variables for the current master element
    rThisKinematicVariables.Initialize();

    PointType global_point;
    for (IndexType i_geom = 0; i_geom < rGeometricalObjectsPointSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            rSlaveGeometry.GlobalCoordinates(global_point, rGeometricalObjectsPointSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompositionType decomp_geom( PointerVector<PointType>{points_array} );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, rSlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (!bad_shape) {
            const GeometryType::IntegrationPointsArrayType& r_integration_points_slave = decomp_geom.IntegrationPoints( rThisIntegrationMethod );

            // Integrating the mortar operators
            PointType local_point_parent, gp_global;
            for ( IndexType point_number = 0; point_number < r_integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp{r_integration_points_slave[point_number].Coordinates()};
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                rSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

                rSlaveGeometry.ShapeFunctionsValues( rThisKinematicVariables.NSlave, local_point_parent.Coordinates() );
                rThisKinematicVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                // Integrate
                const double integration_weight = r_integration_points_slave[point_number].Weight();

                Ae_data.CalculateAeComponents(rThisKinematicVariables, integration_weight);
            }
        }
    }

    BoundedMatrixType Ae;
    Ae_data.CalculateAe(Ae);

    return Ae;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
inline BoundedMatrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::InvertDiagonalMatrix(const BoundedMatrixType& rInputMatrix)
{
    BoundedMatrix<double, TNumNodes, TNumNodes> inv_matrix;

    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            if (i == j) {
                inv_matrix(i, i) = 1.0/rInputMatrix(i, i);
            } else {
                inv_matrix(i, i) = 0.0;
            }
        }
    }

    return inv_matrix;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
inline void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::InvertDiagonalMatrix(
    const BoundedMatrixType& rInputMatrix,
    BoundedMatrixType& rInvertedMatrix
    )
{
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            if (i == j) {
                rInvertedMatrix(i, i) = 1.0/rInputMatrix(i, i);
            } else {
                rInvertedMatrix(i, j) = 0.0;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::LumpMatrix(BoundedMatrixType& rInputMatrix)
{
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            if (i != j) {
                rInputMatrix(i, i) += rInputMatrix(i, j);
                rInputMatrix(i, j) = 0.0;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetSystemSize(std::size_t& rSizeSystem)
{
    rSizeSystem = mDestinationModelPart.Nodes().size();
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CreateSlaveConectivityDatabase(
    std::size_t& rSizeSystem,
    IntMap& rConectivityDatabase,
    IntMap& rInverseConectivityDatabase
)
{
    // Initialize the value
    rSizeSystem = 0;

    // Iterating over nodes
    auto& r_nodes_array = mDestinationModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // We create the database
    for(IndexType i = 0; i < r_nodes_array.size(); ++i) {
        auto it_node = it_node_begin + i;
        rConectivityDatabase[rSizeSystem] = it_node->Id();
        rInverseConectivityDatabase[it_node->Id()] = rSizeSystem;
        rSizeSystem += 1;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
IntegrationMethod SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetIntegrationMethod()
{
    const int integration_order = mThisParameters["integration_order"].GetInt();
    switch ( integration_order )
    {
        case 1: return GeometryData::GI_GAUSS_1;
        case 2: return GeometryData::GI_GAUSS_2;
        case 3: return GeometryData::GI_GAUSS_3;
        case 4: return GeometryData::GI_GAUSS_4;
        case 5: return GeometryData::GI_GAUSS_5;
        default: return GeometryData::GI_GAUSS_2;
    }

    return GeometryData::GI_GAUSS_2;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
bool SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CheckWholeVector(std::vector<bool>& rVectorToCheck)
{
    bool result = true;

    for(IndexType i = 0; i < rVectorToCheck.size(); ++i) {
        if (rVectorToCheck[i] == false) return false;
    }

    return result;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ComputeResidualMatrix(
    Matrix& rResidualMatrix,
    const GeometryType& rSlaveGeometry,
    const GeometryType& rMasterGeometry,
    const MortarOperatorType& rThisMortarOperators
    )
{
    const SizeType size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
    Matrix var_origin_matrix(TNumNodesMaster, size_to_compute);
    if (mOptions.Is(ORIGIN_IS_HISTORICAL)) {
        MortarUtilities::MatrixValue<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(rMasterGeometry, *mpOriginVariable, var_origin_matrix);
    } else {
        MortarUtilities::MatrixValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(rMasterGeometry, *mpOriginVariable, var_origin_matrix);
    }
    Matrix var_destination_matrix(TNumNodes, size_to_compute);
    if (mOptions.Is(DESTINATION_IS_HISTORICAL)) {
        MortarUtilities::MatrixValue<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(rSlaveGeometry, *mpDestinationVariable, var_destination_matrix);
    } else {
        MortarUtilities::MatrixValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(rSlaveGeometry, *mpDestinationVariable, var_destination_matrix);
    }

    const SizeType size_1 = var_destination_matrix.size1();
    const SizeType size_2 = var_destination_matrix.size2();
    if (rResidualMatrix.size1() != size_1  || rResidualMatrix.size2() !=  size_2) {
        rResidualMatrix.resize(size_1, size_2, false);
    }

    noalias(rResidualMatrix) = prod(rThisMortarOperators.MOperator, var_origin_matrix) - prod(rThisMortarOperators.DOperator, var_destination_matrix);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::AssembleRHSAndLHS(
    MatrixType& rA,
    std::vector<VectorType>& rb,
    const SizeType VariableSize,
    const Matrix& rResidualMatrix,
    const GeometryType& rSlaveGeometry,
    IntMap& rInverseConectivityDatabase,
    const MortarOperatorType& rThisMortarOperators
    )
{
    double* values_vector = rA.value_data().begin();
    SizeType* index1_vector = rA.index1_data().begin();
    SizeType* index2_vector = rA.index2_data().begin();

    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const SizeType node_i_id = rSlaveGeometry[i_node].Id();
        const SizeType pos_i_id = static_cast<SizeType>(rInverseConectivityDatabase[node_i_id]);

        // First we assemble the RHS
        for (IndexType i_var = 0; i_var < VariableSize; ++i_var) {
            double& aux_b = rb[i_var][pos_i_id];
            #pragma omp atomic
            aux_b += rResidualMatrix(i_node, i_var);
        }

        SizeType left_limit = index1_vector[pos_i_id];
        SizeType last_pos = left_limit;
        while(pos_i_id != index2_vector[last_pos]) last_pos++;
        double& a_value = values_vector[last_pos];

        #pragma omp atomic
        a_value += rThisMortarOperators.DOperator(i_node, i_node);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::AssembleRHS(
    std::vector<VectorType>& rb,
    const SizeType VariableSize,
    const Matrix& rResidualMatrix,
    const GeometryType& rSlaveGeometry,
    IntMap& rInverseConectivityDatabase
    )
{
    // First we assemble the RHS
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const int node_i_id = rSlaveGeometry[i_node].Id();
        const int pos_i_id = rInverseConectivityDatabase[node_i_id];
        for (IndexType i_var = 0; i_var < VariableSize; ++i_var) {
            double& aux_b = rb[i_var][pos_i_id];
            #pragma omp atomic
            aux_b += rResidualMatrix(i_node, i_var);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ExecuteExplicitMapping()
{
    KRATOS_TRY;

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mOriginModelPart, mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED));
    MortarUtilities::ComputeNodesMeanNormalModelPart(mDestinationModelPart, mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED));

    // Defining tolerance
    const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
    const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
    const double distance_threshold = mThisParameters["distance_threshold"].GetDouble();
    const double zero_tolerance_factor = mThisParameters["zero_tolerance_factor"].GetDouble();
    const bool remove_isolated_conditions = mThisParameters["remove_isolated_conditions"].GetBool();
    const SizeType max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
    IndexType iteration = 0;

    // We set to zero the variables
    if (mOptions.Is(DESTINATION_IS_HISTORICAL)) {
        MortarUtilities::ResetValue<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(mDestinationModelPart, *mpDestinationVariable);
    } else {
        MortarUtilities::ResetValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(mDestinationModelPart, *mpDestinationVariable);
    }

    // Declaring auxiliar values
    IntMap inverse_conectivity_database;
    MatrixType A;
    std::vector<VectorType> b;

    // Creating the assemble database
    SizeType system_size;
    GetSystemSize(system_size);

    // Defining the operators
    const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
    std::vector<bool> is_converged(variable_size, false);
    std::vector<double> residual_norm(variable_size, 0.0);
    std::vector<double> norm_b0(variable_size, 0.0);
    std::vector<double> norm_bi(variable_size, 0.0);
    double increment_residual_norm = 0.0;

    // Create and initialize condition variables:
    MortarKinematicVariablesType this_kinematic_variables;

    // Create the mortar operators
    MortarOperatorType this_mortar_operators;

    // We call the exact integration utility
    ExactMortarIntegrationUtilityType integration_utility = ExactMortarIntegrationUtilityType(TDim, distance_threshold, 0, zero_tolerance_factor);

    // We reset the nodal area
    ResetNodalArea();

    // We get the reference area
    const double ref_area = GetReferenceArea();

    // Check if the pairs has been created
    CheckAndPerformSearch();

    // We clear unused conditions before compute
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // Iterate over conditions
        auto& r_conditions_array = mDestinationModelPart.Conditions();
        const int num_conditions = static_cast<int>(r_conditions_array.size());
        const auto it_cond_begin = r_conditions_array.begin();

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            if (it_cond->Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                ClearIndexes<IndexSet>(p_indexes_pairs, (*it_cond.base()), integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                ClearIndexes<IndexMap>(p_indexes_pairs, (*it_cond.base()), integration_utility);
            }
        }
    } else {
        // Iterate over elements
        auto& r_elements_array = mDestinationModelPart.Elements();
        const int num_elements = static_cast<int>(r_elements_array.size());
        const auto it_elem_begin = r_elements_array.begin();

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            if (it_elem->Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_SET ); // These are the master elements
                ClearIndexes<IndexSet>(p_indexes_pairs, (*it_elem.base()), integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_MAP ); // These are the master elements
                ClearIndexes<IndexMap>(p_indexes_pairs, (*it_elem.base()), integration_utility);
            }
        }
    }

    // In case of discontinous interface we create an inverse mapping
    if (mOptions.Is(DISCONTINOUS_INTERFACE)) {
        CreateInverseDatabase();
    }

    // Compute until convergence
    while (!CheckWholeVector(is_converged) && iteration < max_number_iterations) {
        // We reset the auxiliar variable
        MortarUtilities::ResetAuxiliarValue<TVarType>(mOriginModelPart);
        MortarUtilities::ResetAuxiliarValue<TVarType>(mDestinationModelPart);

        if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
            // Iterate over conditions
            auto& r_conditions_array = mDestinationModelPart.Conditions();
            const int num_conditions = static_cast<int>(r_conditions_array.size());
            const auto it_cond_begin = r_conditions_array.begin();

            // We map the values from one side to the other
            #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = it_cond_begin + i;

                if (it_cond->Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                    PerformMortarOperations<IndexSet>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_cond.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                    PerformMortarOperations<IndexMap>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_cond.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                }
            }
        } else {
            // Iterate over elements
            auto& r_elements_array = mDestinationModelPart.Elements();
            const int num_elements = static_cast<int>(r_elements_array.size());
            const auto it_elem_begin = r_elements_array.begin();

            // We map the values from one side to the other
            #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = it_elem_begin + i;

                if (it_elem->Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_SET ); // These are the master elements
                    PerformMortarOperations<IndexSet>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_elem.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_MAP ); // These are the master elements
                    PerformMortarOperations<IndexMap>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_elem.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                }
            }
        }

        // We remove the not used conditions
        if (remove_isolated_conditions) {
            ModelPart& r_root_model_part = mOriginModelPart.GetRootModelPart();
            if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                r_root_model_part.RemoveConditionsFromAllLevels(TO_ERASE);
            } else {
                r_root_model_part.RemoveElementsFromAllLevels(TO_ERASE);
            }
        }

        auto& r_nodes_array = mDestinationModelPart.Nodes();
        const auto it_node_begin = r_nodes_array.begin();
        const int num_nodes = static_cast<int>(r_nodes_array.size());

        // We compute the residual norm
        for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
            residual_norm[i_size] = 0.0;
        }
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = it_node_begin + i;
            NodeType::Pointer pnode = *(it_node.base());
            if (mOptions.Is(DESTINATION_IS_HISTORICAL)) {
                MortarUtilities::AddAreaWeightedNodalValue<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(pnode, *mpDestinationVariable, ref_area);
            } else {
                MortarUtilities::AddAreaWeightedNodalValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(pnode, *mpDestinationVariable, ref_area);
            }
            for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
                const double value = MortarUtilities::GetAuxiliarValue<TVarType>(pnode, i_size);
                #pragma omp atomic
                residual_norm[i_size] += std::pow(value, 2);
            }
        }

        // Finally we compute the residual
        for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
            residual_norm[i_size] = std::sqrt(residual_norm[i_size])/system_size;
            if (iteration == 0) norm_b0[i_size] = residual_norm[i_size];
            if (residual_norm[i_size] < absolute_convergence_tolerance) is_converged[i_size] = true;
            if (residual_norm[i_size]/norm_b0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;
            if (iteration > 0) {
                increment_residual_norm = std::abs((residual_norm[i_size] - norm_bi[i_size])/norm_bi[i_size]);
                if (increment_residual_norm < relative_convergence_tolerance) is_converged[i_size] = true;
            }
            norm_bi[i_size] = residual_norm[i_size];
            KRATOS_INFO_IF("Mortar mapper ", mEchoLevel > 0) << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm[i_size] << "\tRELATIVE: " << residual_norm[i_size]/norm_b0[i_size] << "\tINCREMENT: " << increment_residual_norm << std::endl;
        }

        iteration += 1;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ExecuteImplicitMapping()
{
    KRATOS_TRY;

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mOriginModelPart, mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED));
    MortarUtilities::ComputeNodesMeanNormalModelPart(mDestinationModelPart, mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED));

    // Defining tolerance
    const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
    const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
    const double distance_threshold = mThisParameters["distance_threshold"].GetDouble();
    const double zero_tolerance_factor = mThisParameters["zero_tolerance_factor"].GetDouble();
    const bool remove_isolated_conditions = mThisParameters["remove_isolated_conditions"].GetBool();
    const SizeType max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
    IndexType iteration = 0;

    // We set to zero the variables
    if (mOptions.Is(DESTINATION_IS_HISTORICAL)) {
        MortarUtilities::ResetValue<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(mDestinationModelPart,  *mpDestinationVariable);
    } else {
        MortarUtilities::ResetValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(mDestinationModelPart, *mpDestinationVariable);
    }

    // Creating the assemble database
    SizeType system_size;
    IntMap conectivity_database, inverse_conectivity_database;
    CreateSlaveConectivityDatabase(system_size, conectivity_database, inverse_conectivity_database);

    // Defining the operators
    const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
    std::vector<bool> is_converged(variable_size, false);
    MatrixType A(system_size, system_size);
    for (IndexType i = 0; i < system_size; ++i) {
        A.push_back(i, i, 0.0);
    }
    const VectorType zero_vector = ZeroVector(system_size);
    std::vector<VectorType> b(variable_size, zero_vector);
    std::vector<double> norm_b0(variable_size, 0.0);
    std::vector<double> norm_Dx0(variable_size, 0.0);
    VectorType Dx(system_size);

    // Create and initialize condition variables:
    MortarKinematicVariablesType this_kinematic_variables;

    // Create the mortar operators
    MortarOperatorType this_mortar_operators;

    // We call the exact integration utility
    ExactMortarIntegrationUtilityType integration_utility = ExactMortarIntegrationUtilityType(TDim, distance_threshold, 0, zero_tolerance_factor);

    // Check if the pairs has been created
    CheckAndPerformSearch();

    // We clear unused conditions before compute
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // Iterate over conditions
        auto& r_conditions_array = mDestinationModelPart.Conditions();
        const int num_conditions = static_cast<int>(r_conditions_array.size());
        const auto it_cond_begin = r_conditions_array.begin();

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = it_cond_begin + i;

            if (it_cond->Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                ClearIndexes<IndexSet>(p_indexes_pairs, (*it_cond.base()), integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                ClearIndexes<IndexMap>(p_indexes_pairs, (*it_cond.base()), integration_utility);
            }
        }
    } else {
        // Iterate over elements
        auto& r_elements_array = mDestinationModelPart.Elements();
        const int num_elements = static_cast<int>(r_elements_array.size());
        const auto it_elem_begin = r_elements_array.begin();

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;

            if (it_elem->Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_SET ); // These are the master elements
                ClearIndexes<IndexSet>(p_indexes_pairs, (*it_elem.base()), integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_MAP ); // These are the master elements
                ClearIndexes<IndexMap>(p_indexes_pairs, (*it_elem.base()), integration_utility);
            }
        }
    }

    // In case of discontinous interface we create an inverse mapping
    if (mOptions.Is(DISCONTINOUS_INTERFACE)) {
        CreateInverseDatabase();
    }

    // Compute until convergence
    while (!CheckWholeVector(is_converged) && iteration < max_number_iterations) {
        // We reset the RHS
        if (iteration > 0)
            for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
                b[i_size] = zero_vector;
            }

        if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
            // Iterate over conditions
            auto& r_conditions_array = mDestinationModelPart.Conditions();
            const int num_conditions = static_cast<int>(r_conditions_array.size());
            const auto it_cond_begin = r_conditions_array.begin();

            // We map the values from one side to the other
            #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
            for(int i = 0; i < num_conditions; ++i) {
                auto it_cond = it_cond_begin + i;
                if (it_cond->Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                    PerformMortarOperations<IndexSet, true>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_cond.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                    PerformMortarOperations<IndexMap, true>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_cond.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                }
            }
        } else {
            // Iterate over elements
            auto& r_elements_array = mDestinationModelPart.Elements();
            const int num_elements = static_cast<int>(r_elements_array.size());
            const auto it_elem_begin = r_elements_array.begin();

            // We map the values from one side to the other
            #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
            for(int i = 0; i < num_elements; ++i) {
                auto it_elem = it_elem_begin + i;
                if (it_elem->Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_SET ); // These are the master elements
                    PerformMortarOperations<IndexSet, true>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_elem.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_MAP ); // These are the master elements
                    PerformMortarOperations<IndexMap, true>(A, b, inverse_conectivity_database, p_indexes_pairs, (*it_elem.base()), integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
                }
            }
        }

        // We remove the not used conditions
        if (remove_isolated_conditions) {
            ModelPart& r_root_model_part = mOriginModelPart.GetRootModelPart();
            if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                r_root_model_part.RemoveConditionsFromAllLevels(TO_ERASE);
            } else {
                r_root_model_part.RemoveElementsFromAllLevels(TO_ERASE);
            }
        }

        // Finally we solve the system
        for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
            mpThisLinearSolver->Solve(A, Dx, b[i_size]);
            if (mOptions.Is(DESTINATION_IS_HISTORICAL)) {
                MortarUtilities::UpdateDatabase<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(mDestinationModelPart, *mpDestinationVariable, Dx, i_size, conectivity_database);
            } else {
                MortarUtilities::UpdateDatabase<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(mDestinationModelPart, *mpDestinationVariable, Dx, i_size, conectivity_database);
            }
            const double residual_norm = norm_2(b[i_size])/system_size;
            if (iteration == 0) norm_b0[i_size] = residual_norm;
            const double increment_norm = norm_2(Dx)/system_size;
            if (iteration == 0) norm_Dx0[i_size] = increment_norm;
            if (residual_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
            if (residual_norm/norm_b0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;
            if (increment_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
            if (increment_norm/norm_Dx0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;

            KRATOS_INFO_IF("Mortar mapper ", mEchoLevel > 0) << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm << "\tRELATIVE: " << residual_norm/norm_b0[i_size] << "\tINCREMENT::\tABS" << increment_norm << "\tRELATIVE: " << increment_norm/norm_Dx0[i_size] << std::endl;
        }

        iteration += 1;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CreateInverseDatabase()
{
    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
        /* Conditions */
        // First we clear database
        auto& r_master_conditions_array = mOriginModelPart.Conditions();
        const int num_origin_conditions = static_cast<int>(r_master_conditions_array.size());
        const auto it_cond_origin_begin = r_master_conditions_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < num_origin_conditions; ++i) {
            auto it_cond = it_cond_origin_begin + i;
            if (!it_cond->Has(INDEX_SET)) {
                it_cond->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
            } else {
                it_cond->GetValue(INDEX_SET)->clear();
            }
        }
    } else {
        /* Elements */
        // First we clear database
        auto& r_master_elements_array = mOriginModelPart.Elements();
        const int num_origin_elements = static_cast<int>(r_master_elements_array.size());
        const auto it_elem_origin_begin = r_master_elements_array.begin();

        #pragma omp parallel for
        for(int i = 0; i < num_origin_elements; ++i) {
            auto it_elem = it_elem_origin_begin + i;
            if (!it_elem->Has(INDEX_SET)) {
                it_elem->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
            } else {
                it_elem->GetValue(INDEX_SET)->clear();
            }
        }
    }

    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // Conditions array
        auto& r_conditions_array = mDestinationModelPart.Conditions();
        const int num_conditions = static_cast<int>(r_conditions_array.size());
        const auto it_cond_begin = r_conditions_array.begin();

        // Create an inverted database
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = it_cond_begin + i;
            if (it_cond->Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto p_cond_master = mOriginModelPart.pGetCondition(master_id); // MASTER
                        (p_cond_master->GetValue(INDEX_SET))->AddId(it_cond->Id());
                    } else {
                        auto p_elem_master = mOriginModelPart.pGetElement(master_id); // MASTER
                        (p_elem_master->GetValue(INDEX_SET))->AddId(it_cond->Id());
                    }
                }
            } else {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto p_cond_master = mOriginModelPart.pGetCondition(master_id); // MASTER
                        (p_cond_master->GetValue(INDEX_SET))->AddId(it_cond->Id());
                    } else {
                        auto p_elem_master = mOriginModelPart.pGetElement(master_id); // MASTER
                        (p_elem_master->GetValue(INDEX_SET))->AddId(it_cond->Id());
                    }
                }
            }
        }
    } else {
        // Elements array
        auto& r_elements_array = mDestinationModelPart.Elements();
        const int num_elements = static_cast<int>(r_elements_array.size());
        const auto it_elem_begin = r_elements_array.begin();

        // Create an inverted database
        for(int i = 0; i < num_elements; ++i) {
            auto it_elem = it_elem_begin + i;
            if (it_elem->Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_SET ); // These are the master elements
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto p_cond_master = mOriginModelPart.pGetCondition(master_id); // MASTER
                        (p_cond_master->GetValue(INDEX_SET))->AddId(it_elem->Id());
                    } else {
                        auto p_elem_master = mOriginModelPart.pGetElement(master_id); // MASTER
                        (p_elem_master->GetValue(INDEX_SET))->AddId(it_elem->Id());
                    }

                }
            } else {
                IndexMap::Pointer p_indexes_pairs = it_elem->GetValue( INDEX_MAP ); // These are the master elements
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto p_cond_master = mOriginModelPart.pGetCondition(master_id); // MASTER
                        (p_cond_master->GetValue(INDEX_SET))->AddId(it_elem->Id());
                    } else {
                        auto p_elem_master = mOriginModelPart.pGetElement(master_id); // MASTER
                        (p_elem_master->GetValue(INDEX_SET))->AddId(it_elem->Id());
                    }
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::UpdateInterface()
{
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // Iterate in the conditions
        auto& r_destination_conditions_array = mDestinationModelPart.Conditions();
        const auto it_cond_begin = r_destination_conditions_array.begin();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_destination_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;
            // Reset the index set
            if (it_cond->Has(INDEX_SET)) {
                (it_cond->GetValue(INDEX_SET))->clear();
            }
            // Reset the index map
            if (it_cond->Has(INDEX_MAP)) {
                (it_cond->GetValue(INDEX_MAP))->clear();
            }
        }
    } else {
        // Iterate in the elements
        auto& r_destination_elements_array = mDestinationModelPart.Elements();
        const auto it_elem_begin = r_destination_elements_array.begin();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_destination_elements_array.size()); ++i) {
            auto it_elem = it_elem_begin + i;
            // Reset the index set
            if (it_elem->Has(INDEX_SET)) {
                (it_elem->GetValue(INDEX_SET))->clear();
            }
            // Reset the index map
            if (it_elem->Has(INDEX_MAP)) {
                (it_elem->GetValue(INDEX_MAP))->clear();
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
Parameters SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "echo_level"                       : 0,
        "using_average_nodal_normal"       : true,
        "discontinuous_interface"          : false,
        "discontinous_interface_factor"    : 1.0e-4,
        "absolute_convergence_tolerance"   : 1.0e-9,
        "relative_convergence_tolerance"   : 1.0e-4,
        "max_number_iterations"            : 10,
        "integration_order"                : 2,
        "distance_threshold"               : 1.0e24,
        "zero_tolerance_factor"            : 1.0e0,
        "remove_isolated_conditions"       : false,
        "mapping_coefficient"              : 1.0e0,
        "origin_variable"                  : "TEMPERATURE",
        "destination_variable"             : "",
        "origin_variable_historical"       : true,
        "origin_are_conditions"            : true,
        "destination_variable_historical"  : true,
        "destination_are_conditions"       : true,
        "update_interface"                 : true,
        "search_parameters"                : {
            "allocation_size"                  : 1000,
            "bucket_size"                      : 4,
            "search_factor"                    : 3.5
        }
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class SimpleMortarMapperProcess<2, 2, Variable<double>>;
template class SimpleMortarMapperProcess<3, 3, Variable<double>>;
template class SimpleMortarMapperProcess<3, 4, Variable<double>>;
template class SimpleMortarMapperProcess<3, 3, Variable<double>, 4>;
template class SimpleMortarMapperProcess<3, 4, Variable<double>, 3>;

template class SimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>, 4>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>, 3>;

}  // namespace Kratos.
