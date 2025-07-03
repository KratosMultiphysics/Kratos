//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "processes/simple_mortar_mapper_process.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/mortar_utilities.h"
#include "utilities/normal_calculation_utils.h"
#include "utilities/math_utils.h"
#include "utilities/atomic_utilities.h"

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
template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Kratos::Flags SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CONSIDER_TESELLATION(Kratos::Flags::Create(6));

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
    const Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Setting flags
    mOptions.Set(AVERAGE_NORMAL, mThisParameters["using_average_nodal_normal"].GetBool());
    mOptions.Set(DISCONTINOUS_INTERFACE, mThisParameters["discontinuous_interface"].GetBool());
    mOptions.Set(ORIGIN_IS_HISTORICAL, mThisParameters["origin_variable_historical"].GetBool());
    mOptions.Set(DESTINATION_IS_HISTORICAL, mThisParameters["destination_variable_historical"].GetBool());
    mOptions.Set(ORIGIN_SKIN_IS_CONDITION_BASED, mThisParameters["origin_are_conditions"].GetBool());
    mOptions.Set(DESTINATION_SKIN_IS_CONDITION_BASED, mThisParameters["destination_are_conditions"].GetBool());
    mOptions.Set(CONSIDER_TESELLATION, mThisParameters["consider_tessellation"].GetBool());

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
    const Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Setting flags
    mOptions.Set(AVERAGE_NORMAL, mThisParameters["using_average_nodal_normal"].GetBool());
    mOptions.Set(DISCONTINOUS_INTERFACE, mThisParameters["discontinuous_interface"].GetBool());
    mOptions.Set(ORIGIN_IS_HISTORICAL, mThisParameters["origin_variable_historical"].GetBool());
    mOptions.Set(DESTINATION_IS_HISTORICAL, mThisParameters["destination_variable_historical"].GetBool());
    mOptions.Set(ORIGIN_SKIN_IS_CONDITION_BASED, mThisParameters["origin_are_conditions"].GetBool());
    mOptions.Set(DESTINATION_SKIN_IS_CONDITION_BASED, mThisParameters["destination_are_conditions"].GetBool());
    mOptions.Set(CONSIDER_TESELLATION, mThisParameters["consider_tessellation"].GetBool());

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
    const Parameters default_parameters = GetDefaultParameters();
    mThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Setting flags
    mOptions.Set(AVERAGE_NORMAL, mThisParameters["using_average_nodal_normal"].GetBool());
    mOptions.Set(DISCONTINOUS_INTERFACE, mThisParameters["discontinuous_interface"].GetBool());
    mOptions.Set(ORIGIN_IS_HISTORICAL, mThisParameters["origin_variable_historical"].GetBool());
    mOptions.Set(DESTINATION_IS_HISTORICAL, mThisParameters["destination_variable_historical"].GetBool());
    mOptions.Set(ORIGIN_SKIN_IS_CONDITION_BASED, mThisParameters["origin_are_conditions"].GetBool());
    mOptions.Set(DESTINATION_SKIN_IS_CONDITION_BASED, mThisParameters["destination_are_conditions"].GetBool());
    mOptions.Set(CONSIDER_TESELLATION, mThisParameters["consider_tessellation"].GetBool());

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

        if (mOptions.Is(DESTINATION_IS_HISTORICAL)){
            block_for_each(r_nodes_array, [this](Node& rNode){
                rNode.FastGetSolutionStepValue(*mpDestinationVariable) *= mMappingCoefficient;
            });
        } else {
            block_for_each(r_nodes_array, [this](Node& rNode){
                rNode.GetValue(*mpDestinationVariable) *= mMappingCoefficient;
            });
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::CheckAndPerformSearch()
{
    // If we force to update the database
    const bool update_interface = mThisParameters["update_interface"].GetBool();

    // First we check if search already exists
    bool search_exists = true;

    // Iterate in the conditions and elements
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        auto& r_destination_conditions_array = mDestinationModelPart.Conditions();
        for(auto& r_cond : r_destination_conditions_array) {
            if (!(r_cond.Has( INDEX_SET ) || r_cond.Has( INDEX_MAP ))) {
                search_exists = false;
                break;
            }
        }

        if (!search_exists) {
            // We create the variable INDEX_SET
            block_for_each(r_destination_conditions_array, [&](Condition& rCond){
                if(!rCond.Has(INDEX_SET)){
                    rCond.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
                }
            });
        }
    } else {
        auto& r_destination_elements_array = mDestinationModelPart.Elements();
        for(auto& r_elem : r_destination_elements_array) {
            if (!(r_elem.Has( INDEX_SET ) || r_elem.Has( INDEX_MAP ))) {
                search_exists = false;
                break;
            }
        }

        if (!search_exists) {
            // We create the variable INDEX_SET
            block_for_each(r_destination_elements_array, [&](Element& rElem){
                if(!rElem.Has(INDEX_SET))
                    rElem.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
            });
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

        IndexPartition<std::size_t>(point_list_destination.size()).for_each([point_list_destination](std::size_t Index){
            point_list_destination[Index]->UpdatePoint();
        });

        // Some auxiliar values
        const SizeType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt(); // Allocation size for the vectors and max number of potential results
        const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble(); // The search factor to be considered
        SizeType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt(); // Bucket size for kd-tree

        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTreeType tree_points(point_list_destination.begin(), point_list_destination.end(), bucket_size);

        if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
            // Iterate over conditions
            for(auto& r_cond : mDestinationModelPart.Conditions()) {
                FillDatabase<Condition>(r_cond, tree_points, allocation_size, search_factor);
            }
        } else {
            // Iterate over elements
            for(auto& r_elem : mDestinationModelPart.Elements()) {
                FillDatabase<Element>(r_elem, tree_points, allocation_size, search_factor);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::ResetNodalArea()
{
    // We set to zero
    VariableUtils().SetNonHistoricalVariableToZero(NODAL_AREA, mDestinationModelPart.Nodes());
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
double SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetReferenceArea()
{
    double ref_area = 0.0;
    double current_area;

    // We look for the max area in the origin model part
    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
        for(auto& r_cond : mOriginModelPart.Conditions()) {
            current_area = r_cond.GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    } else {
        for(auto& r_elem : mOriginModelPart.Elements()) {
            current_area = r_elem.GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    }

    // We look for the max area in the destination model part
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        for(auto& r_cond : mDestinationModelPart.Conditions()) {
            current_area = r_cond.GetGeometry().Area();
            if (current_area > ref_area) ref_area = current_area;
        }
    } else {
        for(auto& r_elem : mDestinationModelPart.Elements()) {
            current_area = r_elem.GetGeometry().Area();
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
    // Auxiliar value
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
    // We initialize the Ae components
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

    // We create the database
    IndexType aux_index;
    for(auto& r_node : mDestinationModelPart.Nodes()) {
        aux_index = r_node.Id();
        rConectivityDatabase[rSizeSystem] = aux_index;
        rInverseConectivityDatabase[aux_index] = rSizeSystem;
        rSizeSystem += 1;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
GeometryData::IntegrationMethod SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetIntegrationMethod()
{
    const int integration_order = mThisParameters["integration_order"].GetInt();
    switch ( integration_order ) {
        case 1: return GeometryData::IntegrationMethod::GI_GAUSS_1;
        case 2: return GeometryData::IntegrationMethod::GI_GAUSS_2;
        case 3: return GeometryData::IntegrationMethod::GI_GAUSS_3;
        case 4: return GeometryData::IntegrationMethod::GI_GAUSS_4;
        case 5: return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default: return GeometryData::IntegrationMethod::GI_GAUSS_2;
    }

    return GeometryData::IntegrationMethod::GI_GAUSS_2;
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
            AtomicAdd(aux_b, rResidualMatrix(i_node, i_var));
        }

        SizeType left_limit = index1_vector[pos_i_id];
        SizeType last_pos = left_limit;
        while(pos_i_id != index2_vector[last_pos]) last_pos++;
        double& a_value = values_vector[last_pos];

        AtomicAdd(a_value, rThisMortarOperators.DOperator(i_node, i_node));
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
        const SizeType node_i_id = rSlaveGeometry[i_node].Id();
        const SizeType pos_i_id = rInverseConectivityDatabase[node_i_id];
        for (IndexType i_var = 0; i_var < VariableSize; ++i_var) {
            double& aux_b = rb[i_var][pos_i_id];
            AtomicAdd(aux_b, rResidualMatrix(i_node, i_var));
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
    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(mOriginModelPart, true);
    } else {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ElementsContainerType>(mOriginModelPart, true);
    }
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(mDestinationModelPart, true);
    } else {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ElementsContainerType>(mDestinationModelPart, true);
    }

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

    // We reset the nodal area
    ResetNodalArea();

    // We get the reference area
    const double ref_area = GetReferenceArea();

    TLS tls;
    tls.integration_utility = ExactMortarIntegrationUtilityType(TDim, distance_threshold, 0, zero_tolerance_factor, mOptions.Is(CONSIDER_TESELLATION));

    // Check if the pairs has been created
    CheckAndPerformSearch();

    // We clear unused conditions before compute
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // We map the values from one side to the other
        block_for_each(mDestinationModelPart.Conditions(), tls, [&](Condition& rCond, TLS& rTLS) {
            if (rCond.Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = rCond.GetValue( INDEX_SET ); // These are the master conditions
                ClearIndexes<IndexSet>(p_indexes_pairs, rCond, rTLS.integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = rCond.GetValue( INDEX_MAP ); // These are the master conditions
                ClearIndexes<IndexMap>(p_indexes_pairs, rCond, rTLS.integration_utility);
            }
        });
    } else {
        // We map the values from one side to the other
        block_for_each(mDestinationModelPart.Elements(), tls, [&](Element& rElem, TLS& rTLS) {
            if (rElem.Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = rElem.GetValue( INDEX_SET ); // These are the master elements
                ClearIndexes<IndexSet>(p_indexes_pairs, rElem, rTLS.integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = rElem.GetValue( INDEX_MAP ); // These are the master elements
                ClearIndexes<IndexMap>(p_indexes_pairs, rElem, rTLS.integration_utility);
            }
        });
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
            // We map the values from one side to the other
            block_for_each(mDestinationModelPart.Conditions(), tls, [&](Condition& rCond, TLS& rTLS) {
                if (rCond.Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = rCond.GetValue( INDEX_SET ); // These are the master conditions
                    PerformMortarOperations<IndexSet>(A, b, inverse_conectivity_database, p_indexes_pairs, rCond, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = rCond.GetValue( INDEX_MAP ); // These are the master conditions
                    PerformMortarOperations<IndexMap>(A, b, inverse_conectivity_database, p_indexes_pairs, rCond, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                }
            });
        } else {
            // We map the values from one side to the other
            block_for_each(mDestinationModelPart.Elements(), tls, [&](Element& rElem, TLS& rTLS) {
                if (rElem.Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = rElem.GetValue( INDEX_SET ); // These are the master elements
                    PerformMortarOperations<IndexSet>(A, b, inverse_conectivity_database, p_indexes_pairs, rElem, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = rElem.GetValue( INDEX_MAP ); // These are the master elements
                    PerformMortarOperations<IndexMap>(A, b, inverse_conectivity_database, p_indexes_pairs, rElem, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                }
            });
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

        // We compute the residual norm
        for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
            residual_norm[i_size] = 0.0;
        }

        block_for_each(mDestinationModelPart.Nodes(), [&](Node& rNode) {
            if(mOptions.Is(DESTINATION_IS_HISTORICAL)) {
                MortarUtilities::AddAreaWeightedNodalValue<TVarType, MortarUtilitiesSettings::SaveAsHistoricalVariable>(rNode, *mpDestinationVariable, ref_area);
            } else {
                MortarUtilities::AddAreaWeightedNodalValue<TVarType, MortarUtilitiesSettings::SaveAsNonHistoricalVariable>(rNode, *mpDestinationVariable, ref_area);
            }
            for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
                const double value = MortarUtilities::GetAuxiliarValue<TVarType>(rNode, i_size);
                AtomicAdd(residual_norm[i_size], std::pow(value, 2));
            }
        });

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
    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(mOriginModelPart, true);
    } else {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ElementsContainerType>(mOriginModelPart, true);
    }
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ConditionsContainerType>(mDestinationModelPart, true);
    } else {
        NormalCalculationUtils().CalculateUnitNormals<ModelPart::ElementsContainerType>(mDestinationModelPart, true);
    }

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

    TLS tls;
    tls.integration_utility = ExactMortarIntegrationUtilityType(TDim, distance_threshold, 0, zero_tolerance_factor, mOptions.Is(CONSIDER_TESELLATION));

    // Check if the pairs has been created
    CheckAndPerformSearch();

    // We clear unused conditions before compute
    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // We map the values from one side to the other
        block_for_each(mDestinationModelPart.Conditions(), tls, [&](Condition& rCond, TLS& rTLS) {
            if (rCond.Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = rCond.GetValue( INDEX_SET ); // These are the master conditions
                ClearIndexes<IndexSet>(p_indexes_pairs, rCond, rTLS.integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = rCond.GetValue( INDEX_MAP ); // These are the master conditions
                ClearIndexes<IndexMap>(p_indexes_pairs, rCond, rTLS.integration_utility);
            }
        });
    } else {
        // We map the values from one side to the other
        block_for_each(mDestinationModelPart.Elements(), tls, [&](Element& rElem, TLS& rTLS) {
            if (rElem.Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = rElem.GetValue( INDEX_SET ); // These are the master elements
                ClearIndexes<IndexSet>(p_indexes_pairs, rElem, rTLS.integration_utility);
            } else {
                IndexMap::Pointer p_indexes_pairs = rElem.GetValue( INDEX_MAP ); // These are the master elements
                ClearIndexes<IndexMap>(p_indexes_pairs, rElem, rTLS.integration_utility);
            }
        });
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
            // We map the values from one side to the other
            block_for_each(mDestinationModelPart.Conditions(), tls, [&](Condition& rCond, TLS& rTLS) {
                if (rCond.Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = rCond.GetValue( INDEX_SET ); // These are the master conditions
                    PerformMortarOperations<IndexSet, true>(A, b, inverse_conectivity_database, p_indexes_pairs, rCond, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = rCond.GetValue( INDEX_MAP ); // These are the master conditions
                    PerformMortarOperations<IndexMap, true>(A, b, inverse_conectivity_database, p_indexes_pairs, rCond, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                }
            });
        } else {
            // We map the values from one side to the other
            block_for_each(mDestinationModelPart.Elements(), tls, [&](Element& rElem, TLS& rTLS) {
                if (rElem.Has( INDEX_SET )) {
                    IndexSet::Pointer p_indexes_pairs = rElem.GetValue( INDEX_SET ); // These are the master elements
                    PerformMortarOperations<IndexSet, true>(A, b, inverse_conectivity_database, p_indexes_pairs, rElem, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                } else {
                    IndexMap::Pointer p_indexes_pairs = rElem.GetValue( INDEX_MAP ); // These are the master elements
                    PerformMortarOperations<IndexMap, true>(A, b, inverse_conectivity_database, p_indexes_pairs, rElem, rTLS.integration_utility, rTLS.this_kinematic_variables, rTLS.this_mortar_operators, iteration);
                }
            });
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
        block_for_each(mOriginModelPart.Conditions(), [&](Condition& rCond){
            if (!rCond.Has(INDEX_SET)) {
                rCond.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
            } else {
                rCond.GetValue(INDEX_SET)->clear();
            }
        });
    } else {
        /* Elements */
        block_for_each(mOriginModelPart.Elements(), [&](Element& rElem){
            if (!rElem.Has(INDEX_SET)) {
                rElem.SetValue(INDEX_SET, Kratos::make_shared<IndexSet>());
            } else {
                rElem.GetValue(INDEX_SET)->clear();
            }
        });

    }

    if (mOptions.Is(DESTINATION_SKIN_IS_CONDITION_BASED)) {
        // Conditions array: Create an inverted database
        for(auto& r_cond : mDestinationModelPart.Conditions()) {
            if (r_cond.Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = r_cond.GetValue( INDEX_SET ); // These are the master conditions
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto& r_cond_master = mOriginModelPart.GetCondition(master_id); // MASTER
                        (r_cond_master.GetValue(INDEX_SET))->AddId(r_cond.Id());
                    } else {
                        auto& r_elem_master = mOriginModelPart.GetElement(master_id); // MASTER
                        (r_elem_master.GetValue(INDEX_SET))->AddId(r_cond.Id());
                    }
                }
            } else {
                IndexMap::Pointer p_indexes_pairs = r_cond.GetValue( INDEX_MAP ); // These are the master conditions
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto& r_cond_master = mOriginModelPart.GetCondition(master_id); // MASTER
                        (r_cond_master.GetValue(INDEX_SET))->AddId(r_cond.Id());
                    } else {
                        auto& r_elem_master = mOriginModelPart.GetElement(master_id); // MASTER
                        (r_elem_master.GetValue(INDEX_SET))->AddId(r_cond.Id());
                    }
                }
            }
        }
    } else {
        // Elements array. Create an inverted database
        for(auto& r_elem : mDestinationModelPart.Elements()) {
            if (r_elem.Has( INDEX_SET )) {
                IndexSet::Pointer p_indexes_pairs = r_elem.GetValue( INDEX_SET ); // These are the master elements
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto& r_cond_master = mOriginModelPart.GetCondition(master_id); // MASTER
                        (r_cond_master.GetValue(INDEX_SET))->AddId(r_elem.Id());
                    } else {
                        auto& r_elem_master = mOriginModelPart.GetElement(master_id); // MASTER
                        (r_elem_master.GetValue(INDEX_SET))->AddId(r_elem.Id());
                    }

                }
            } else {
                IndexMap::Pointer p_indexes_pairs = r_elem.GetValue( INDEX_MAP ); // These are the master elements
                for (auto it_pair = p_indexes_pairs->begin(); it_pair != p_indexes_pairs->end(); ++it_pair ) {
                    const IndexType master_id = p_indexes_pairs->GetId(it_pair);
                    if (mOptions.Is(ORIGIN_SKIN_IS_CONDITION_BASED)) {
                        auto& r_cond_master = mOriginModelPart.GetCondition(master_id); // MASTER
                        (r_cond_master.GetValue(INDEX_SET))->AddId(r_elem.Id());
                    } else {
                        auto& r_elem_master = mOriginModelPart.GetElement(master_id); // MASTER
                        (r_elem_master.GetValue(INDEX_SET))->AddId(r_elem.Id());
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
        block_for_each(mDestinationModelPart.Conditions(), [&](Condition& rCond){
            // Reset the index set
            if (rCond.Has(INDEX_SET)) {
                (rCond.GetValue(INDEX_SET))->clear();
            }
            // Reset the index map
            if (rCond.Has(INDEX_MAP)) {
                (rCond.GetValue(INDEX_MAP))->clear();
            }
        });
    } else {
        block_for_each(mDestinationModelPart.Elements(), [&](Element& rElem){
            // Reset the index set
            if (rElem.Has(INDEX_SET)) {
                (rElem.GetValue(INDEX_SET))->clear();
            }
            // Reset the index map
            if (rElem.Has(INDEX_MAP)) {
                (rElem.GetValue(INDEX_MAP))->clear();
            }
        });
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, class TVarType, const SizeType TNumNodesMaster>
const Parameters SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, TNumNodesMaster>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "echo_level"                       : 0,
        "consider_tessellation"            : false,
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