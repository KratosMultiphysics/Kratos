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
template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::SimpleMortarMapperProcess(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    TVarType& OriginVariable,
    TVarType& DestinationVariable,
    Parameters ThisParameters,
    LinearSolverType::Pointer pThisLinearSolver
    ): mOriginModelPart(rOriginModelPart),
       mDestinationModelPart(rDestinationModelPart),
       mOriginVariable(OriginVariable),
       mDestinationVariable(DestinationVariable),
       mThisParameters(ThisParameters),
       mpThisLinearSolver(pThisLinearSolver)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    // We set some values
    mOriginHistorical = mThisParameters["origin_variable_historical"].GetBool();
    mDestinationHistorical = mThisParameters["destination_variable_historical"].GetBool();
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::SimpleMortarMapperProcess(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    TVarType& ThisVariable,
    Parameters ThisParameters,
    LinearSolverType::Pointer pThisLinearSolver
):  mOriginModelPart(rOriginModelPart),
    mDestinationModelPart(rDestinationModelPart),
    mOriginVariable(ThisVariable),
    mDestinationVariable(ThisVariable),
    mThisParameters(ThisParameters),
    mpThisLinearSolver(pThisLinearSolver)
{
    // The default parameters
    Parameters default_parameters = GetDefaultParameters();
    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    // We set some values
    mOriginHistorical = mThisParameters["origin_variable_historical"].GetBool();
    mDestinationHistorical = mThisParameters["destination_variable_historical"].GetBool();
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>:: Execute()
{
    KRATOS_TRY;

    if (mpThisLinearSolver == nullptr)
        ExecuteExplicitMapping();
    else
        ExecuteImplicitMapping();

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::CheckAndPerformSearch()
{
    // First we check if search already exists
    bool search_exists = true;
    // Iterate in the conditions
    ConditionsArrayType& destination_conditions_array = mDestinationModelPart.Conditions();
    for(IndexType i = 0; i < destination_conditions_array.size(); ++i) {
        auto it_cond = destination_conditions_array.begin() + i;

        if (!(it_cond->Has( INDEX_SET ) || it_cond->Has( INDEX_MAP ))) {
            search_exists = false; 
            break;
        }
    }
        
    // Now we perform the corresponding search
    if (search_exists == false) {
        
        // We create the variable INDEX_SET
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(destination_conditions_array.size()); ++i) {
            auto it_cond = destination_conditions_array.begin() + i;
            if (it_cond->Has(INDEX_SET) == false) 
                it_cond->SetValue(INDEX_SET, Kratos::make_shared<IndexSet>()); 
        }
        
        // A list that contents the all the points (from nodes) from the modelpart 
        PointVector point_list_destination; 
        
        point_list_destination.clear();
        
        // Iterate in the conditions
        ConditionsArrayType& origin_conditions_array = mOriginModelPart.Conditions();
        
        // Creating a buffer for parallel vector fill
        const int num_threads = OpenMPUtils::GetNumThreads();
        std::vector<PointVector> points_buffer(num_threads);

        #pragma omp parallel
        {
            const int thread_id = OpenMPUtils::ThisThread();

            #pragma omp for
            for(int i = 0; i < static_cast<int>(origin_conditions_array.size()); ++i) {
                auto it_cond = origin_conditions_array.begin() + i;
                
                const PointTypePointer& p_point = PointTypePointer(new PointMapperType((*it_cond.base())));
                (points_buffer[thread_id]).push_back(p_point);
            }
            
            // Combine buffers together
            #pragma omp single
            {
                for( auto& point_buffer : points_buffer)
                    std::move(point_buffer.begin(),point_buffer.end(),back_inserter(point_list_destination));
            }
        }
        
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(point_list_destination.size()); ++i) 
            point_list_destination[i]->UpdatePoint();
        
        // Some auxiliar values
        const SizeType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt(); // Allocation size for the vectors and max number of potential results
        const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble(); // The search factor to be considered 
        SizeType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt(); // Bucket size for kd-tree
        
        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTreeType tree_points(point_list_destination.begin(), point_list_destination.end(), bucket_size);
        
        for(IndexType i = 0; i < destination_conditions_array.size(); ++i) {
            auto it_cond = destination_conditions_array.begin() + i;
            
            // Initialize values
            PointVector points_found(allocation_size);
            
            GeometryType& geometry = it_cond->GetGeometry();
            const Point& center = geometry.Center();
            
            double radius = 0.0; 
            for(IndexType i_node = 0; i_node < it_cond->GetGeometry().PointsNumber(); ++i_node)  {
                const array_1d<double, 3> aux_vector = center.Coordinates() - it_cond->GetGeometry()[i_node].Coordinates();
                const double aux_value = inner_prod(aux_vector, aux_vector);
                if(aux_value > radius) radius = aux_value; 
            } 
            
            const double search_radius = search_factor * std::sqrt(radius);

            SizeType number_points_found = tree_points.SearchInRadius(center, search_radius, points_found.begin(), allocation_size);
            
            if (number_points_found > 0) {  
                IndexSet::Pointer indexes_set = it_cond->GetValue(INDEX_SET);
                
                for (IndexType i_point = 0; i_point < number_points_found; ++i_point ) {
                    Condition::Pointer p_cond_master = points_found[i_point]->GetCondition();
                    indexes_set->AddId(p_cond_master->Id());
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::ResetNodalArea()
{
    NodesArrayType& nodes_array = mDestinationModelPart.Nodes();

    // We set to zero
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        it_node->SetValue(NODAL_AREA, 0.0);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
double SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::GetReferenceArea()
{
    double ref_area = 0.0;
    ConditionsArrayType& conditions_array_origin = mOriginModelPart.Conditions();

    // We look for the max area in the origin model part
    for(int i = 0; i < static_cast<int>(conditions_array_origin.size()); ++i) {
        auto it_cond = conditions_array_origin.begin() + i;
        const double current_area = it_cond->GetGeometry().Area();
        if (current_area > ref_area) ref_area = current_area;
    }
    
    ConditionsArrayType& conditions_array_destination = mDestinationModelPart.Conditions();
    
    // We look for the max area in the destination model part
    for(int i = 0; i < static_cast<int>(conditions_array_destination.size()); ++i) {
        auto it_cond = conditions_array_destination.begin() + i;
        const double current_area = it_cond->GetGeometry().Area();
        if (current_area > ref_area) ref_area = current_area;
    }
    
    return ref_area;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::AssemblyMortarOperators(
    const std::vector<array_1d<PointType,TDim>>& ConditionsPointSlave,
    GeometryType& SlaveGeometry,
    GeometryType& MasterGeometry,
    const array_1d<double, 3>& MasterNormal,
    MortarKinematicVariables<TNumNodes>& ThisKinematicVariables,
    MortarOperator<TNumNodes>& ThisMortarOperators,
    const IntegrationMethod& ThisIntegrationMethod,
    const BoundedMatrixType Ae
)
{
    for (IndexType i_geom = 0; i_geom < ConditionsPointSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompType decomp_geom( PointerVector<PointType>{points_array} );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (bad_shape == false) {
            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp{integration_points_slave[point_number].Coordinates()};
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                SlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

                /// SLAVE CONDITION ///
                SlaveGeometry.ShapeFunctionsValues( ThisKinematicVariables.NSlave, local_point_parent.Coordinates() );
                ThisKinematicVariables.PhiLagrangeMultipliers = prod(Ae, ThisKinematicVariables.NSlave);

                ThisKinematicVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                /// MASTER CONDITION ///
                PointType projected_gp_global;
                const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(ThisKinematicVariables.NSlave, SlaveGeometry);

                GeometryType::CoordinatesArrayType slave_gp_global;
                SlaveGeometry.GlobalCoordinates( slave_gp_global, local_point_parent );
                GeometricalProjectionUtilities::FastProjectDirection( MasterGeometry, gp_global, projected_gp_global, MasterNormal, -gp_normal ); // The opposite direction

                GeometryType::CoordinatesArrayType projected_gp_local;

                MasterGeometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;

                // SHAPE FUNCTIONS
                MasterGeometry.ShapeFunctionsValues( ThisKinematicVariables.NMaster, projected_gp_local );

                const double integration_weight = integration_points_slave[point_number].Weight();

                ThisMortarOperators.CalculateMortarOperators(ThisKinematicVariables, integration_weight);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
inline BoundedMatrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::CalculateAe(
    GeometryType& SlaveGeometry,
    MortarKinematicVariables<TNumNodes>& ThisKinematicVariables,
    std::vector<array_1d<PointType,TDim>>& ConditionsPointsSlave,
    const IntegrationMethod& ThisIntegrationMethod
)
{
    // We initilize the Ae components
    DualLagrangeMultiplierOperators<TNumNodes> Ae_data;
    Ae_data.Initialize();

    // Initialize general variables for the current master element
    ThisKinematicVariables.Initialize();

    for (IndexType i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompType decomp_geom( PointerVector<PointType>{points_array} );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (bad_shape == false) {
            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp{integration_points_slave[point_number].Coordinates()};
                PointType local_point_parent;
                PointType gp_global;
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                SlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

                SlaveGeometry.ShapeFunctionsValues( ThisKinematicVariables.NSlave, local_point_parent.Coordinates() );
                ThisKinematicVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                // Integrate
                const double integration_weight = integration_points_slave[point_number].Weight();

                Ae_data.CalculateAeComponents(ThisKinematicVariables, integration_weight);
            }
        }
    }

    BoundedMatrixType Ae;
    Ae_data.CalculateAe(Ae);

    return Ae;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
inline BoundedMatrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::InvertDiagonalMatrix(const BoundedMatrixType& InputMatrix)
{
    BoundedMatrixType inv_matrix = ZeroMatrix(TNumNodes);

    for (IndexType i = 0; i < TNumNodes; ++i)
        inv_matrix(i, i) = 1.0/InputMatrix(i, i);

    return inv_matrix;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
inline void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::InvertDiagonalMatrix(
    const BoundedMatrixType& InputMatrix,
    BoundedMatrixType& InvertedMatrix
)
{
    InvertedMatrix = ZeroMatrix(TNumNodes);

    for (IndexType i = 0; i < TNumNodes; ++i)
        InvertedMatrix(i, i) = 1.0/InputMatrix(i, i);
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::LumpMatrix(BoundedMatrixType& InputMatrix)
{
    for (IndexType i = 0; i < TNumNodes; ++i) {
        for (IndexType j = 0; j < TNumNodes; ++j) {
            if (i != j) {
                InputMatrix(i, i) += InputMatrix(i, j);
                InputMatrix(i, j) = 0.0;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::GetSystemSize(std::size_t& SizeSystem)
{
    SizeSystem = mDestinationModelPart.Nodes().size();
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::CreateSlaveConectivityDatabase(
    std::size_t& SizeSystem,
    IntMap& ConectivityDatabase,
    IntMap& InverseConectivityDatabase
)
{
    // Initialize the value
    SizeSystem = 0;

    NodesArrayType& nodes_array = mDestinationModelPart.Nodes();

    // We create the database
    for(IndexType i = 0; i < nodes_array.size(); ++i) {
        auto it_node = nodes_array.begin() + i;
        ConectivityDatabase[SizeSystem] = it_node->Id();
        InverseConectivityDatabase[it_node->Id()] = SizeSystem;
        SizeSystem += 1;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
IntegrationMethod SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::GetIntegrationMethod()
{
    const int& integration_order = mThisParameters["integration_order"].GetInt();
    switch ( integration_order )
    {
    case 1:
        return GeometryData::GI_GAUSS_1;
    case 2:
        return GeometryData::GI_GAUSS_2;
    case 3:
        return GeometryData::GI_GAUSS_3;
    case 4:
        return GeometryData::GI_GAUSS_4;
    case 5:
        return GeometryData::GI_GAUSS_5;
    default:
        return GeometryData::GI_GAUSS_2;
    }

    return GeometryData::GI_GAUSS_2;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
bool SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::CheckWholeVector(std::vector<bool> VectorToCheck)
{
    bool result = true;

    for(IndexType i = 0; i < VectorToCheck.size(); ++i)
        if (VectorToCheck[i] == false) return false;

    return result;
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::ComputeResidualMatrix(
    Matrix& ResidualMatrix,
    GeometryType& SlaveGeometry,
    GeometryType& MasterGeometry,
    const MortarOperator<TNumNodes>& ThisMortarOperators
)
{
    const SizeType size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
    Matrix var_origin_matrix(TNumNodes, size_to_compute);
    if (mOriginHistorical)
        MortarUtilities::MatrixValue<TVarType, Historical>(MasterGeometry, mOriginVariable, var_origin_matrix);
    else
        MortarUtilities::MatrixValue<TVarType, NonHistorical>(MasterGeometry, mOriginVariable, var_origin_matrix);
    Matrix var_destination_matrix(TNumNodes, size_to_compute);
    if (mDestinationHistorical)
        MortarUtilities::MatrixValue<TVarType, Historical>(SlaveGeometry, mDestinationVariable, var_destination_matrix);
    else
        MortarUtilities::MatrixValue<TVarType, NonHistorical>(SlaveGeometry, mDestinationVariable, var_destination_matrix);

    const SizeType size_1 = var_origin_matrix.size1();
    const SizeType size_2 = var_origin_matrix.size2();
    if (ResidualMatrix.size1() != size_1  || ResidualMatrix.size2() !=  size_2)
        ResidualMatrix.resize(size_1, size_2, false);

    noalias(ResidualMatrix) = prod(ThisMortarOperators.MOperator, var_origin_matrix) - prod(ThisMortarOperators.DOperator, var_destination_matrix);
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::AssembleRHSAndLHS(
    MatrixType& A,
    std::vector<VectorType>& b,
    const SizeType& VariableSize,
    const Matrix& ResidualMatrix,
    GeometryType& SlaveGeometry,
    IntMap& InverseConectivityDatabase,
    const MortarOperator<TNumNodes>& ThisMortarOperators
)
{
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const std::size_t& node_i_id = SlaveGeometry[i_node].Id();
        const std::size_t& pos_i_id = static_cast<std::size_t>(InverseConectivityDatabase[node_i_id]);

        // First we assemble the RHS
        for (IndexType i_var = 0; i_var < VariableSize; ++i_var) {
            double& aux_b = b[i_var][pos_i_id];
            #pragma omp atomic
            aux_b += ResidualMatrix(i_node, i_var);
        }

        double* values_vector = A.value_data().begin();
        std::size_t* index1_vector = A.index1_data().begin();
        std::size_t* index2_vector = A.index2_data().begin();
        std::size_t left_limit = index1_vector[pos_i_id];
        std::size_t last_pos = left_limit;
        while(pos_i_id != index2_vector[last_pos]) last_pos++;
        double& a_value = values_vector[last_pos];

        #pragma omp atomic
        a_value += ThisMortarOperators.DOperator(i_node, i_node);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::AssembleRHS(
    std::vector<VectorType>& b,
    const SizeType& VariableSize,
    const Matrix& ResidualMatrix,
    GeometryType& SlaveGeometry,
    IntMap& InverseConectivityDatabase
)
{
    // First we assemble the RHS
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const int& node_i_id = SlaveGeometry[i_node].Id();
        const int& pos_i_id = InverseConectivityDatabase[node_i_id];
        for (IndexType i_var = 0; i_var < VariableSize; ++i_var) {
            double& aux_b = b[i_var][pos_i_id];
            #pragma omp atomic
            aux_b += ResidualMatrix(i_node, i_var);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::ExecuteExplicitMapping()
{
    KRATOS_TRY;

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mOriginModelPart);
    MortarUtilities::ComputeNodesMeanNormalModelPart(mDestinationModelPart);

    // Defining tolerance
    const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
    const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
    const double distance_threshold = mThisParameters["distance_threshold"].GetDouble();
    const SizeType max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
    IndexType iteration = 0;

    // We set to zero the variables
    if (mDestinationHistorical)
        MortarUtilities::ResetValue<TVarType, Historical>(mDestinationModelPart, mDestinationVariable);
    else
        MortarUtilities::ResetValue<TVarType, NonHistorical>(mDestinationModelPart, mDestinationVariable);

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
    MortarKinematicVariables<TNumNodes> this_kinematic_variables;

    // Create the mortar operators
    MortarOperator<TNumNodes> this_mortar_operators;

    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim, distance_threshold);
    
    // We reset the nodal area
    ResetNodalArea();
    
    // We get the reference area
    const double ref_area = GetReferenceArea();

    // Check if the pairs has been created
    CheckAndPerformSearch();
    
    while (CheckWholeVector(is_converged) == false && iteration < max_number_iterations) {
        // We reset the auxiliar variable
        MortarUtilities::ResetAuxiliarValue<TVarType>(mOriginModelPart);
        MortarUtilities::ResetAuxiliarValue<TVarType>(mDestinationModelPart);

        ConditionsArrayType& conditions_array = mDestinationModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = conditions_array.begin() + i;

            if (it_cond->Has( INDEX_MAP )) {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                PerformMortarOperations<IndexMap>(A, b, inverse_conectivity_database, p_indexes_pairs, it_cond, integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
            } else {
                IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                PerformMortarOperations<IndexSet>(A, b, inverse_conectivity_database, p_indexes_pairs, it_cond, integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
            }
        }

        // We remove the not used conditions
        ModelPart& root_model_part = mOriginModelPart.GetRootModelPart();
        root_model_part.RemoveConditionsFromAllLevels(TO_ERASE);

        NodesArrayType& nodes_array = mDestinationModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());

        // We compute the residual norm
        for (IndexType i_size = 0; i_size < variable_size; ++i_size) residual_norm[i_size] = 0.0;
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;
            NodeType::Pointer pnode = *(it_node.base());
            if (mDestinationHistorical)
                MortarUtilities::AddAreaWeightedNodalValue<TVarType, Historical>(pnode, mDestinationVariable, ref_area);
            else
                MortarUtilities::AddAreaWeightedNodalValue<TVarType, NonHistorical>(pnode, mDestinationVariable, ref_area);
            for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
                const double& value = MortarUtilities::GetAuxiliarValue<TVarType>(pnode, i_size);
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
            if (mEchoLevel > 0)
                KRATOS_INFO("Mortar mapper ")  << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm[i_size] << "\tRELATIVE: " << residual_norm[i_size]/norm_b0[i_size] << "\tINCREMENT: " << increment_residual_norm << std::endl;
        }

        iteration += 1;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::ExecuteImplicitMapping()
{
    KRATOS_TRY;

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mOriginModelPart);
    MortarUtilities::ComputeNodesMeanNormalModelPart(mDestinationModelPart);

    // Defining tolerance
    const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
    const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
    const double distance_threshold = mThisParameters["distance_threshold"].GetDouble();
    const SizeType max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
    IndexType iteration = 0;

    // We set to zero the variables
    if (mDestinationHistorical)
        MortarUtilities::ResetValue<TVarType, Historical>(mDestinationModelPart,  mDestinationVariable);
    else
        MortarUtilities::ResetValue<TVarType, NonHistorical>(mDestinationModelPart, mDestinationVariable);

    // Creating the assemble database
    SizeType system_size;
    IntMap conectivity_database, inverse_conectivity_database;
    CreateSlaveConectivityDatabase(system_size, conectivity_database, inverse_conectivity_database);

    // Defining the operators
    const SizeType variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
    std::vector<bool> is_converged(variable_size, false);
    MatrixType A(system_size, system_size);
    for (IndexType i = 0; i < system_size; ++i)
        A.push_back(i, i, 0.0);
    const VectorType zero_vector = ZeroVector(system_size);
    std::vector<VectorType> b(variable_size, zero_vector);
    std::vector<double> norm_b0(variable_size, 0.0);
    std::vector<double> norm_Dx0(variable_size, 0.0);
    VectorType Dx(system_size);

    // Create and initialize condition variables:
    MortarKinematicVariables<TNumNodes> this_kinematic_variables;

    // Create the mortar operators
    MortarOperator<TNumNodes> this_mortar_operators;

    // We call the exact integration utility
    ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim, distance_threshold);

    // Check if the pairs has been created
    CheckAndPerformSearch();
    
    while (CheckWholeVector(is_converged) == false && iteration < max_number_iterations) {
        // We reset the RHS
        if (iteration > 0)
            for (IndexType i_size = 0; i_size < variable_size; ++i_size)
                b[i_size] = zero_vector;

        ConditionsArrayType& conditions_array = mDestinationModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = conditions_array.begin() + i;
            if (it_cond->Has( INDEX_MAP )) {
                IndexMap::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_MAP ); // These are the master conditions
                PerformMortarOperations<IndexMap, true>(A, b, inverse_conectivity_database, p_indexes_pairs, it_cond, integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
            } else {
                IndexSet::Pointer p_indexes_pairs = it_cond->GetValue( INDEX_SET ); // These are the master conditions
                PerformMortarOperations<IndexSet, true>(A, b, inverse_conectivity_database, p_indexes_pairs, it_cond, integration_utility, this_kinematic_variables, this_mortar_operators, iteration);
            }
        }

        // We remove the not used conditions
        ModelPart& root_model_part = mOriginModelPart.GetRootModelPart();
        root_model_part.RemoveConditionsFromAllLevels(TO_ERASE);

        // Finally we solve the system
        for (IndexType i_size = 0; i_size < variable_size; ++i_size) {
            mpThisLinearSolver->Solve(A, Dx, b[i_size]);
            if (mDestinationHistorical)
                MortarUtilities::UpdateDatabase<TVarType, Historical>(mDestinationModelPart, mDestinationVariable, Dx, i_size, conectivity_database);
            else
                MortarUtilities::UpdateDatabase<TVarType, NonHistorical>(mDestinationModelPart, mDestinationVariable, Dx, i_size, conectivity_database);
            const double residual_norm = norm_2(b[i_size])/system_size;
            if (iteration == 0) norm_b0[i_size] = residual_norm;
            const double increment_norm = norm_2(Dx)/system_size;
            if (iteration == 0) norm_Dx0[i_size] = increment_norm;
            if (residual_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
            if (residual_norm/norm_b0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;
            if (increment_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
            if (increment_norm/norm_Dx0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;

            if (mEchoLevel > 0)
                KRATOS_INFO("Mortar mapper ") << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm << "\tRELATIVE: " << residual_norm/norm_b0[i_size] << "\tINCREMENT::\tABS" << increment_norm << "\tRELATIVE: " << increment_norm/norm_Dx0[i_size] << std::endl;
        }

        iteration += 1;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< std::size_t TDim, std::size_t TNumNodes, class TVarType>
Parameters SimpleMortarMapperProcess<TDim, TNumNodes, TVarType>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "echo_level"                       : 0,
        "absolute_convergence_tolerance"   : 1.0e-9,
        "relative_convergence_tolerance"   : 1.0e-4,
        "max_number_iterations"            : 10,
        "integration_order"                : 2,
        "distance_threshold"               : 1.0e24,
        "origin_variable_historical"       : true,
        "destination_variable_historical"  : true,
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

template class SimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>>;

}  // namespace Kratos.
