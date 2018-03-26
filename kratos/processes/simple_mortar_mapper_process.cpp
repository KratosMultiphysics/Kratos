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
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::SimpleMortarMapperProcess(
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
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::SimpleMortarMapperProcess(
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
    mEchoLevel = mThisParameters["echo_level"].GetInt();
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>:: Execute()
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::CheckAndPerformSearch()
{
    // First we check if search already exists
    bool search_exists = true;
    // Iterate in the conditions
    ConditionsArrayType& destination_conditions_array = mDestinationModelPart.Conditions();
    for(std::size_t i = 0; i < destination_conditions_array.size(); ++i) {
        auto it_cond = destination_conditions_array.begin() + i;

        if (!it_cond->Has( INDEX_SET )) {
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
        const unsigned int allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt(); // Allocation size for the vectors and max number of potential results 
        const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble(); // The search factor to be considered 
        unsigned int bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt(); // Bucket size for kd-tree
        
        // Create a tree
        // It will use a copy of mNodeList (a std::vector which contains pointers)
        // Copying the list is required because the tree will reorder it for efficiency
        KDTreeType tree_points(point_list_destination.begin(), point_list_destination.end(), bucket_size);
        
        for(std::size_t i = 0; i < destination_conditions_array.size(); ++i) {
            auto it_cond = destination_conditions_array.begin() + i;
            
            // Initialize values
            PointVector points_found(allocation_size);
            
            GeometryType& geometry = it_cond->GetGeometry();
            const Point& center = geometry.Center();
            
            double radius = 0.0; 
            for(std::size_t i_node = 0; i_node < it_cond->GetGeometry().PointsNumber(); ++i_node)  { 
                const array_1d<double, 3>& aux_vector = center.Coordinates() - it_cond->GetGeometry()[i_node].Coordinates();
                const double aux_value = inner_prod(aux_vector, aux_vector);
                if(aux_value > radius) radius = aux_value; 
            } 
            
            const double search_radius = search_factor * std::sqrt(radius);

            std::size_t number_points_found = tree_points.SearchInRadius(center, search_radius, points_found.begin(), allocation_size);
            
            if (number_points_found > 0) {  
                IndexSet::Pointer indexes_set = it_cond->GetValue(INDEX_SET);
                
                for (std::size_t i_point = 0; i_point < number_points_found; ++i_point ) {
                    Condition::Pointer p_cond_master = points_found[i_point]->GetCondition();
                    indexes_set->AddId(p_cond_master->Id());
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::ResetNodalArea()
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
double SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::GetReferenceArea()
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::AssemblyMortarOperators(
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
    for (unsigned int i_geom = 0; i_geom < ConditionsPointSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (unsigned int i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompType decomp_geom( points_array );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (bad_shape == false) {
            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
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
                MortarUtilities::FastProjectDirection( MasterGeometry, gp_global, projected_gp_global, MasterNormal, -gp_normal ); // The opposite direction

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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
inline bounded_matrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::CalculateAe(
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

    for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom) {
        std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (unsigned int i_node = 0; i_node < TDim; ++i_node) {
            PointType global_point;
            SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompType decomp_geom( points_array );

        const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);

        if (bad_shape == false) {
            const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );

            // Integrating the mortar operators
            for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
inline bounded_matrix<double, TNumNodes, TNumNodes> SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::InvertDiagonalMatrix(const BoundedMatrixType& InputMatrix)
{
    BoundedMatrixType inv_matrix = ZeroMatrix(TNumNodes);

    for (std::size_t i = 0; i < TNumNodes; ++i)
        inv_matrix(i, i) = 1.0/InputMatrix(i, i);

    return inv_matrix;
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
inline void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::InvertDiagonalMatrix(
    const BoundedMatrixType& InputMatrix,
    BoundedMatrixType& InvertedMatrix
)
{
    InvertedMatrix = ZeroMatrix(TNumNodes);

    for (std::size_t i = 0; i < TNumNodes; ++i)
        InvertedMatrix(i, i) = 1.0/InputMatrix(i, i);
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::LumpMatrix(BoundedMatrixType& InputMatrix)
{
    for (std::size_t i = 0; i < TNumNodes; ++i) {
        for (std::size_t j = 0; j < TNumNodes; ++j) {
            if (i != j) {
                InputMatrix(i, i) += InputMatrix(i, j);
                InputMatrix(i, j) = 0.0;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::GetSystemSize(std::size_t& SizeSystem)
{
    SizeSystem = mDestinationModelPart.Nodes().size();
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::CreateSlaveConectivityDatabase(
    std::size_t& SizeSystem,
    IntMap& ConectivityDatabase,
    IntMap& InverseConectivityDatabase
)
{
    // Initialize the value
    SizeSystem = 0;

    NodesArrayType& nodes_array = mDestinationModelPart.Nodes();

    // We create the database
    for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
        auto it_node = nodes_array.begin() + i;
        ConectivityDatabase[SizeSystem] = it_node->Id();
        InverseConectivityDatabase[it_node->Id()] = SizeSystem;
        SizeSystem += 1;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
IntegrationMethod SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::GetIntegrationMethod()
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
bool SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::CheckWholeVector(std::vector<bool> VectorToCheck)
{
    bool result = true;

    for(std::size_t i = 0; i < VectorToCheck.size(); ++i)
        if (VectorToCheck[i] == false) return false;

    return result;
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::ComputeResidualMatrix(
    Matrix& ResidualMatrix,
    GeometryType& SlaveGeometry,
    GeometryType& MasterGeometry,
    const MortarOperator<TNumNodes>& ThisMortarOperators
)
{
    const unsigned int size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
    Matrix var_origin_matrix(TNumNodes, size_to_compute);
    MortarUtilities::MatrixValue<TVarType, THistOrigin>(MasterGeometry, mOriginVariable, var_origin_matrix);
    Matrix var_destination_matrix(TNumNodes, size_to_compute);
    MortarUtilities::MatrixValue<TVarType, THistDestination>(SlaveGeometry, mDestinationVariable, var_destination_matrix);

    const std::size_t size_1 = var_origin_matrix.size1();
    const std::size_t size_2 = var_origin_matrix.size2();
    if (ResidualMatrix.size1() != size_1  || ResidualMatrix.size2() !=  size_2)
        ResidualMatrix.resize(size_1, size_2, false);

    noalias(ResidualMatrix) = prod(ThisMortarOperators.MOperator, var_origin_matrix) - prod(ThisMortarOperators.DOperator, var_destination_matrix);
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::AssembleRHSAndLHS(
    MatrixType& A,
    std::vector<VectorType>& b,
    const unsigned int& VariableSize,
    const Matrix& ResidualMatrix,
    GeometryType& SlaveGeometry,
    IntMap& InverseConectivityDatabase,
    const MortarOperator<TNumNodes>& ThisMortarOperators
)
{
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
        const std::size_t& node_i_id = SlaveGeometry[i_node].Id();
        const std::size_t& pos_i_id = static_cast<std::size_t>(InverseConectivityDatabase[node_i_id]);

        // First we assemble the RHS
        for (unsigned int i_var = 0; i_var < VariableSize; ++i_var) {
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::AssembleRHS(
    std::vector<VectorType>& b,
    const unsigned int& VariableSize,
    const Matrix& ResidualMatrix,
    GeometryType& SlaveGeometry,
    IntMap& InverseConectivityDatabase
)
{
    // First we assemble the RHS
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node) {
        const int& node_i_id = SlaveGeometry[i_node].Id();
        const int& pos_i_id = InverseConectivityDatabase[node_i_id];
        for (unsigned int i_var = 0; i_var < VariableSize; ++i_var) {
            double& aux_b = b[i_var][pos_i_id];
            #pragma omp atomic
            aux_b += ResidualMatrix(i_node, i_var);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::ExecuteExplicitMapping()
{
    KRATOS_TRY;

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mOriginModelPart);
    MortarUtilities::ComputeNodesMeanNormalModelPart(mDestinationModelPart);

    // Defining tolerance
    const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
    const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
    const double distance_threshold = mThisParameters["distance_threshold"].GetDouble();
    const unsigned int max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
    unsigned int iteration = 0;

    // We set to zero the variables
    MortarUtilities::ResetValue<TVarType, THistDestination>(mDestinationModelPart, mDestinationVariable);

    // Getting the auxiliar variable
    TVarType aux_variable = MortarUtilities::GetAuxiliarVariable<TVarType>();

    // Creating the assemble database
    std::size_t system_size;
    GetSystemSize(system_size);

    // Defining the operators
    const unsigned int variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
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

            const array_1d<double, 3>& slave_normal = it_cond->GetValue(NORMAL);
            GeometryType& slave_geometry = it_cond->GetGeometry();

            IndexSet::Pointer& indexes_set = it_cond->GetValue( INDEX_SET ); // These are the master conditions
            std::vector<std::size_t> indexes_to_remove;

            for (auto it_pair = indexes_set->begin(); it_pair != indexes_set->end(); ++it_pair ) {
                Condition::Pointer p_cond_master = mOriginModelPart.pGetCondition(*it_pair); // MASTER
                const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL);
                GeometryType& master_geometry = p_cond_master->GetGeometry();

                const IntegrationMethod& this_integration_method = GetIntegrationMethod();

                // Reading integration points
                std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping
                const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);

                if (is_inside == true) {
                    // Initialize general variables for the current master element
                    this_kinematic_variables.Initialize();

                    // Initialize the mortar operators
                    this_mortar_operators.Initialize();

                    const BoundedMatrixType& Ae = CalculateAe(slave_geometry, this_kinematic_variables, conditions_points_slave, this_integration_method);

                    AssemblyMortarOperators( conditions_points_slave, slave_geometry, master_geometry,master_normal,this_kinematic_variables, this_mortar_operators, this_integration_method, Ae);

                    /* We compute the residual */
                    const unsigned int size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
                    Matrix residual_matrix(TNumNodes, size_to_compute);
                    ComputeResidualMatrix(residual_matrix, slave_geometry, master_geometry, this_mortar_operators);

                    MortarUtilities::AddValue<TVarType, NonHistorical>(slave_geometry, aux_variable, residual_matrix);

                    // We check if DOperator is diagonal
                    if (mEchoLevel > 1) {
                        BoundedMatrixType aux_copy_D = this_mortar_operators.DOperator;
                        LumpMatrix(aux_copy_D);
                        const BoundedMatrixType aux_diff = aux_copy_D - this_mortar_operators.DOperator;
                        const double norm_diff = norm_frobenius(aux_diff);
                        if (norm_diff > 1.0e-4)
                            KRATOS_WARNING("D OPERATOR") << " THE MORTAR OPERATOR D IS NOT DIAGONAL" << std::endl;
                        if (mEchoLevel == 3) {
                            KRATOS_WATCH(norm_diff);
                            KRATOS_WATCH(this_mortar_operators.DOperator);
                        }
                    }

                    if (iteration == 0) // Just assembled the first iteration
                        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
                            slave_geometry[i_node].GetValue(NODAL_AREA) += this_mortar_operators.DOperator(i_node, i_node);
                } else
                    indexes_to_remove.push_back(*it_pair);
            }

            // Clear indexes
            for (unsigned int i_to_remove = 0; i_to_remove < indexes_to_remove.size(); ++i_to_remove)
                indexes_set->RemoveId(indexes_to_remove[i_to_remove]);
        }

        NodesArrayType& nodes_array = mDestinationModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size());

        // We compute the residual norm
        for (unsigned int i_size = 0; i_size < variable_size; ++i_size) residual_norm[i_size] = 0.0;
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) {
            auto it_node = nodes_array.begin() + i;
            NodeType::Pointer pnode = *(it_node.base());
            MortarUtilities::AddAreaWeightedNodalValue<TVarType, THistDestination>(pnode, mDestinationVariable, ref_area);
            for (unsigned int i_size = 0; i_size < variable_size; ++i_size) {
                const double& value = MortarUtilities::GetAuxiliarValue<TVarType>(pnode, i_size);
                #pragma omp atomic
                residual_norm[i_size] += std::pow(value, 2);
            }
        }

        // Finally we compute the residual
        for (unsigned int i_size = 0; i_size < variable_size; ++i_size) {
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
                KRATOS_INFO("Mortar mapper ")  << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm[i_size] << "\tRELATIVE: " << residual_norm[i_size]/norm_b0[i_size] << "\tINCREMENT: " << increment_residual_norm
                          << std::endl;
        }

        iteration += 1;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
void SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::ExecuteImplicitMapping()
{
    KRATOS_TRY;

    // Calculate the mean of the normal in all the nodes
    MortarUtilities::ComputeNodesMeanNormalModelPart(mOriginModelPart);
    MortarUtilities::ComputeNodesMeanNormalModelPart(mDestinationModelPart);

    // Defining tolerance
    const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
    const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
    const double distance_threshold = mThisParameters["distance_threshold"].GetDouble();
    const unsigned int max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
    unsigned int iteration = 0;

    // We set to zero the variables
    MortarUtilities::ResetValue<TVarType, THistDestination>(mDestinationModelPart, mDestinationVariable);

    // Creating the assemble database
    std::size_t system_size;
    IntMap conectivity_database, inverse_conectivity_database;
    CreateSlaveConectivityDatabase(system_size, conectivity_database, inverse_conectivity_database);

    // Defining the operators
    const unsigned int variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
    std::vector<bool> is_converged(variable_size, false);
    MatrixType A(system_size, system_size);
    for (std::size_t i = 0; i < system_size; ++i)
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
            for (unsigned int i_size = 0; i_size < variable_size; ++i_size)
                b[i_size] = zero_vector;

        ConditionsArrayType& conditions_array = mDestinationModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());

        // We map the values from one side to the other
        #pragma omp parallel for firstprivate(this_kinematic_variables, this_mortar_operators, integration_utility)
        for(int i = 0; i < num_conditions; ++i) {
            auto it_cond = conditions_array.begin() + i;
            const array_1d<double, 3>& slave_normal = it_cond->GetValue(NORMAL);
            GeometryType& slave_geometry = it_cond->GetGeometry();

            IndexSet::Pointer& indexes_set = it_cond->GetValue( INDEX_SET ); // These are the master conditions
            std::vector<std::size_t> indexes_to_remove;

            for (auto it_pair = indexes_set->begin(); it_pair != indexes_set->end(); ++it_pair ) {
                Condition::Pointer p_cond_master = mOriginModelPart.pGetCondition(*it_pair); // MASTER
                const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL);
                GeometryType& master_geometry = p_cond_master->GetGeometry();

                const IntegrationMethod& this_integration_method = GetIntegrationMethod();

                // Reading integration points
                std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping
                const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);

                if (is_inside == true) {
                    // Initialize general variables for the current master element
                    this_kinematic_variables.Initialize();

                    // Initialize the mortar operators
                    this_mortar_operators.Initialize();

                    const BoundedMatrixType& Ae = CalculateAe(slave_geometry, this_kinematic_variables, conditions_points_slave, this_integration_method);

                    AssemblyMortarOperators( conditions_points_slave, slave_geometry, master_geometry,master_normal,this_kinematic_variables, this_mortar_operators, this_integration_method, Ae);

                    /* We compute the residual */
                    const unsigned int size_to_compute = MortarUtilities::SizeToCompute<TDim, TVarType>();
                    Matrix residual_matrix(TNumNodes, size_to_compute);
                    ComputeResidualMatrix(residual_matrix, slave_geometry, master_geometry, this_mortar_operators);

                    // We check if DOperator is diagonal
                    if (mEchoLevel > 1) {
                        BoundedMatrixType aux_copy_D = this_mortar_operators.DOperator;
                        LumpMatrix(aux_copy_D);
                        const BoundedMatrixType aux_diff = aux_copy_D - this_mortar_operators.DOperator;
                        const double norm_diff = norm_frobenius(aux_diff);
                        if (norm_diff > 1.0e-4)
                            KRATOS_WARNING("D OPERATOR") << " THE MORTAR OPERATOR D IS NOT DIAGONAL" << std::endl;
                        if (mEchoLevel == 3) {
                            KRATOS_WATCH(norm_diff);
                            KRATOS_WATCH(this_mortar_operators.DOperator);
                        }
                    }

                    /* We compute the residual and assemble */
                    if (iteration == 0)
                        AssembleRHSAndLHS(A, b, variable_size, residual_matrix, slave_geometry, inverse_conectivity_database, this_mortar_operators);
                    else
                        AssembleRHS(b, variable_size, residual_matrix, slave_geometry, inverse_conectivity_database);
                }
                else
                    indexes_to_remove.push_back(*it_pair);
            }

            // Clear indexes
            for (unsigned int i_to_remove = 0; i_to_remove < indexes_to_remove.size(); ++i_to_remove)
                indexes_set->RemoveId(indexes_to_remove[i_to_remove]);
        }

        // Finally we solve the system
        for (unsigned int i_size = 0; i_size < variable_size; ++i_size)
        {
            mpThisLinearSolver->Solve(A, Dx, b[i_size]);
            MortarUtilities::UpdateDatabase<TVarType, THistDestination>(mDestinationModelPart, mDestinationVariable, Dx, i_size, conectivity_database);
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

template< int TDim, int TNumNodes, class TVarType, HistoricalValues THistOrigin, HistoricalValues THistDestination>
Parameters SimpleMortarMapperProcess<TDim, TNumNodes, TVarType, THistOrigin, THistDestination>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "echo_level"                       : 0,
        "absolute_convergence_tolerance"   : 1.0e-9,
        "relative_convergence_tolerance"   : 1.0e-4,
        "max_number_iterations"            : 10,
        "integration_order"                : 2,
        "distance_threshold"               : 1.0e24,
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

template class SimpleMortarMapperProcess<2, 2, Variable<double>, Historical>;
template class SimpleMortarMapperProcess<3, 3, Variable<double>, Historical>;
template class SimpleMortarMapperProcess<3, 4, Variable<double>, Historical>;

template class SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical>;
template class SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical>;
template class SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical>;

template class SimpleMortarMapperProcess<2, 2, ComponentType, Historical>;
template class SimpleMortarMapperProcess<3, 3, ComponentType, Historical>;
template class SimpleMortarMapperProcess<3, 4, ComponentType, Historical>;

template class SimpleMortarMapperProcess<2, 2, ComponentType, NonHistorical>;
template class SimpleMortarMapperProcess<3, 3, ComponentType, NonHistorical>;
template class SimpleMortarMapperProcess<3, 4, ComponentType, NonHistorical>;

template class SimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>, Historical>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>, Historical>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>, Historical>;

template class SimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>, NonHistorical>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>, NonHistorical>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>, NonHistorical>;

template class SimpleMortarMapperProcess<2, 2, Variable<double>, Historical, NonHistorical>;
template class SimpleMortarMapperProcess<3, 3, Variable<double>, Historical, NonHistorical>;
template class SimpleMortarMapperProcess<3, 4, Variable<double>, Historical, NonHistorical>;

template class SimpleMortarMapperProcess<2, 2, Variable<double>, NonHistorical, Historical>;
template class SimpleMortarMapperProcess<3, 3, Variable<double>, NonHistorical, Historical>;
template class SimpleMortarMapperProcess<3, 4, Variable<double>, NonHistorical, Historical>;

template class SimpleMortarMapperProcess<2, 2, ComponentType, Historical, NonHistorical>;
template class SimpleMortarMapperProcess<3, 3, ComponentType, Historical, NonHistorical>;
template class SimpleMortarMapperProcess<3, 4, ComponentType, Historical, NonHistorical>;

template class SimpleMortarMapperProcess<2, 2, ComponentType, NonHistorical, Historical>;
template class SimpleMortarMapperProcess<3, 3, ComponentType, NonHistorical, Historical>;
template class SimpleMortarMapperProcess<3, 4, ComponentType, NonHistorical, Historical>;

template class SimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>, Historical, NonHistorical>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>, Historical, NonHistorical>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>, Historical, NonHistorical>;

template class SimpleMortarMapperProcess<2, 2, Variable<array_1d<double, 3>>, NonHistorical, Historical>;
template class SimpleMortarMapperProcess<3, 3, Variable<array_1d<double, 3>>, NonHistorical, Historical>;
template class SimpleMortarMapperProcess<3, 4, Variable<array_1d<double, 3>>, NonHistorical, Historical>;

}  // namespace Kratos.
