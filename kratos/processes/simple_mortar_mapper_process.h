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

#if !defined(KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS)
#define KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/point.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

/* Custom includes */
#include "includes/mortar_classes.h"

/* Custom utilities */
#include "utilities/exact_mortar_segmentation_utility.h"
#include "utilities/mortar_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
///@}
///@name  Enum's
///@{

#if !defined(HISTORICAL_VALUES)
#define HISTORICAL_VALUES
    enum HistoricalValues {Historical = 0, NonHistorical = 1};
#endif
    
///@}
///@name  Functions
///@{
    
template< int TDim, int TNumNodes, class TVarType, HistoricalValues THist> 
class SimpleMortarMapperProcess
        : public Process
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of SimpleMortarMapperProcess
    KRATOS_CLASS_POINTER_DEFINITION(SimpleMortarMapperProcess);
    
    typedef Point<3>                                     PointType;
    typedef Node<3>                                       NodeType;
    typedef Geometry<NodeType>                        GeometryType;
    typedef Geometry<PointType>                  GeometryPointType;
    typedef ModelPart::NodesContainerType           NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
    // Type definition for integration methods
    typedef GeometryData::IntegrationMethod      IntegrationMethod;
    
    // Auxiliar geometries
    typedef Line2D2<PointType>                            LineType;
    typedef Triangle3D3<PointType>                    TriangleType;
    
    // Component type
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;

    // Linear solver
    typedef UblasSpace<double, CompressedMatrix, Vector>    SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector>               LocalSpaceType;
    typedef typename SparseSpaceType::MatrixType                 MatrixType;
    typedef typename SparseSpaceType::VectorType                 VectorType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& ThisVariable, 
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(ThisVariable),
           mDestinationVariable(ThisVariable),
           mpThisLinearSolver(pThisLinearSolver)
    {
    }
    
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& OriginVariable,
        TVarType& DestinationVariable,
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(OriginVariable),
           mDestinationVariable(DestinationVariable),
           mpThisLinearSolver(pThisLinearSolver)
    {
    }

    /// Destructor.
    ~SimpleMortarMapperProcess() override
    = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }
    
    ///@}
    ///@name Operations
    ///@{
    
    void Execute() override
    {
        KRATOS_TRY;

        if (mpThisLinearSolver == nullptr)
        {
            ExecuteExplicitMapping();
        }
        else
        {
            ExecuteImplicitMapping();
        }
        
        KRATOS_CATCH("");
    }

    /**
     * This method sets both variables (origin and destination) with the same variable
     */
    void SetVariable(TVarType ThisVariable)
    {
        mOriginVariable = ThisVariable;
        mDestinationVariable = ThisVariable;
    }
    
    /**
     * This method sets both variables (origin and destination) in a separated way
     */
    void SetVariables(
        TVarType OriginVariable,
        TVarType DestinationVariable
        )
    {
        mOriginVariable = OriginVariable;
        mDestinationVariable = DestinationVariable;
    }
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SimpleMortarMapperProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimpleMortarMapperProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{
    
    ModelPart& mrThisModelPart;                   // The model part to compute
    TVarType mOriginVariable;                     // The origin variable to map
    TVarType mDestinationVariable;                // The destiny variable to map
    LinearSolverType::Pointer mpThisLinearSolver; // The linear solver used to compute the solution

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * This method computes the nodal area
     */
    void ComputeNodalArea()
    {
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        // We set to zero
        #pragma omp parallel for
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(NODAL_AREA, 0.0);
        }
        
        // Sum all the nodes areas
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < num_conditions; ++i) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            const unsigned int number_nodes = it_cond->GetGeometry().PointsNumber();
            const double& rArea = it_cond->GetGeometry().Area()/number_nodes;
            
            for (unsigned int i = 0; i < number_nodes; ++i)
            {
                #pragma omp atomic
                it_cond->GetGeometry()[i].GetValue(NODAL_AREA) += rArea;
            }
        }
    }
    
    /**
     * This method computes the Ae matrix
     * @param SlaveGeometry: The slave geometry
     * @param ThisKinematicVariables: The kinematic variables
     * @param ConditionsPointsSlave: The list of decomposed triangles
     * @param ThisIntegrationMethod: The integration method considered
     * @return Ae: The matrix of dual LM
     */
    bounded_matrix<double, TNumNodes, TNumNodes> CalculateAe(
        GeometryType& SlaveGeometry,
        MortarKinematicVariables<TNumNodes>& ThisKinematicVariables,
        std::vector<array_1d<PointType,TDim>>& ConditionsPointsSlave,
        IntegrationMethod ThisIntegrationMethod
        )
    {
        // We initilize the Ae components
        DualLagrangeMultiplierOperators<TNumNodes> Ae_data;
        Ae_data.Initialize();

        // Initialize general variables for the current master element
        ThisKinematicVariables.Initialize();
        
        for (unsigned int i_geom = 0; i_geom < ConditionsPointsSlave.size(); ++i_geom)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
            {
                PointType global_point;
                SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointsSlave[i_geom][i_node]);
                points_array[i_node] = boost::make_shared<PointType>(global_point);
            }
            
            typename std::conditional<TDim == 2, LineType, TriangleType >::type decomp_geom( points_array );
            
            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, SlaveGeometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
            
            if (bad_shape == false)
            { 
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( ThisIntegrationMethod );
                
                // Integrating the mortar operators
                for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                {
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
        
        return Ae_data.CalculateAe();
    }
        
    /**
     * This method inverts a diagonal matrix
     * @param InputMatrix: The matrix to InvertMatrix
     * @return The matrix inverted
     */
    bounded_matrix<double, TNumNodes, TNumNodes> FastInverse(bounded_matrix<double, TNumNodes, TNumNodes> InputMatrix)
    {
        bounded_matrix<double, TNumNodes, TNumNodes> inv_matrix = ZeroMatrix(TNumNodes);
        
        for (std::size_t i = 0; i < TNumNodes; ++i)
        {
            inv_matrix(i, i) = 1.0/InputMatrix(i, i);
        }
        
        return inv_matrix;
    }
        
    /**
     * This method creates a master database needed to assemble the system
     * @param SizeSystem: The size of the system
     * @param ConectivityDatabase: The database that will be used to assemble the system
     * @param InverseConectivityDatabase: The inverse database that will be used to assemble the system
     */
        
    void CreateMasterConectivityDatabase(
        std::size_t& SizeSystem,
        std::unordered_map<int, int>& ConectivityDatabase,
        std::unordered_map<int, int>& InverseConectivityDatabase
        )
    {
        // Initialize the value
        SizeSystem = 0;
        
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        // We create the database
//         #pragma omp parallel for // TODO: Think if it is possible to use OMP
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(MASTER) == true)
            {
                ConectivityDatabase[SizeSystem] = it_node->Id();
                InverseConectivityDatabase[it_node->Id()] = SizeSystem;
                SizeSystem += 1;
            }
        }
    }
    
    /**
     * This method checks if all components of a vector are true
     * @param VectorToCheck: The vector to check
     * @return result: True if all componets are true
     */
    bool CheckWholeVector(std::vector<bool> VectorToCheck)
    {
        bool result = true;

        for(std::size_t i = 0; i < VectorToCheck.size(); ++i)
        {
            if (VectorToCheck[i] == false) return false;
        }
        
        return result;
    }
    
    /**
     * This method executes the explicit mapping (when no linear solver is avalaible)
     */
    void ExecuteExplicitMapping()
    {
        KRATOS_TRY;

        // Defining tolerance
        const double tolerance = 1.0e-15;
//         const double tolerance = std::numeric_limits<double>::epsilon();
        
        // We compute the nodal area
        ComputeNodalArea();
        
        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> this_kinematic_variables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> this_mortar_condition_matrices;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        // Initialize to zero // TODO: Add OMP
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); ++i)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                boost::shared_ptr<ConditionMap>& all_conditions_maps = it_cond->GetValue( MAPPING_PAIRS ); // These are the master conditions
                
                for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
                {
                    Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                    GeometryType& master_geometry = p_cond_master->GetGeometry();
                    
                    MortarUtilities::ResetValue<TVarType, THist>(master_geometry, mDestinationVariable);
                }
            }
        }
        
        // We map the values from one side to the other // TODO: Add OMP
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); ++i)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                const array_1d<double, 3>& slave_normal = it_cond->GetValue(NORMAL);
                GeometryType& slave_geometry = it_cond->GetGeometry();
                
                boost::shared_ptr<ConditionMap>& all_conditions_maps = it_cond->GetValue( MAPPING_PAIRS ); // These are the master conditions
                
                for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
                {
                    // Initialize area
                    double area = 0.0;
                    
                    Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                    const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                    GeometryType& master_geometry = p_cond_master->GetGeometry();
                    
                    IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_3;
                    
                    // Reading integration points
                    std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping  
                    const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);
                    
                    if (is_inside == true)
                    {
                        // Initialize general variables for the current master element
                        this_kinematic_variables.Initialize();
                
                        // Initialize the mortar operators
                        this_mortar_condition_matrices.Initialize();
                        
                        const bounded_matrix<double, TNumNodes, TNumNodes> Ae = CalculateAe(slave_geometry, this_kinematic_variables, conditions_points_slave, this_integration_method);
                        
//                         const bounded_matrix<double, TNumNodes, TNumNodes> aux_matrix = Ae - IdentityMatrix(TNumNodes);
//                         const double norm_aux = norm_frobenius(aux_matrix);
//                         
//                         const bool is_dual = (norm_aux > 0) ? true : false;
                        
                        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
                        {
                            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
                            {
                                PointType global_point;
                                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                points_array[i_node] = boost::make_shared<PointType>(global_point);
                            }
                            
                            typename std::conditional<TDim == 2, LineType, TriangleType >::type decomp_geom( points_array );
                            
                            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                            
                            if (bad_shape == false)
                            {
                                area += decomp_geom.Area();
                                
                                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                                
                                // Integrating the mortar operators
                                for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                                {
                                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                    PointType local_point_parent;
                                    PointType gp_global;
                                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                    slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
           
                                    /// SLAVE CONDITION ///
                                    slave_geometry.ShapeFunctionsValues( this_kinematic_variables.NSlave, local_point_parent.Coordinates() );
//                                     this_kinematic_variables.PhiLagrangeMultipliers = this_kinematic_variables.NSlave;
                                    this_kinematic_variables.PhiLagrangeMultipliers = prod(Ae, this_kinematic_variables.NSlave);
//                                     this_kinematic_variables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
                                    this_kinematic_variables.DetjSlave = slave_geometry.DeterminantOfJacobian( local_point_parent );
                                    
                                    /// MASTER CONDITION ///
                                    PointType projected_gp_global;
                                    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(this_kinematic_variables.NSlave, slave_geometry);
                                    
                                    GeometryType::CoordinatesArrayType slave_gp_global;
                                    slave_geometry.GlobalCoordinates( slave_gp_global, local_point_parent );
                                    MortarUtilities::FastProjectDirection( master_geometry, gp_global, projected_gp_global, master_normal, -gp_normal ); // The opposite direction
                                    
                                    GeometryType::CoordinatesArrayType projected_gp_local;
                                    
                                    master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
                                    
                                    // SHAPE FUNCTIONS 
                                    master_geometry.ShapeFunctionsValues( this_kinematic_variables.NMaster, projected_gp_local );    
                                    
                                    const double integration_weight = integration_points_slave[point_number].Weight();
                                    
                                    this_mortar_condition_matrices.CalculateMortarOperators(this_kinematic_variables, integration_weight);   
                                }
                            }
                        }
                        
                        // We compute the P operator
                        double aux_det = MathUtils<double>::DetMat<TNumNodes>(this_mortar_condition_matrices.DOperator);
                        const bounded_matrix<double, TNumNodes, TNumNodes> inv_d_operator = (std::abs(aux_det) > tolerance) ? MathUtils<double>::InvertMatrix<TNumNodes>(this_mortar_condition_matrices.DOperator, aux_det, tolerance) : ZeroMatrix(TNumNodes);
//                         const bounded_matrix<double, TNumNodes, TNumNodes> inv_d_operator = (is_dual == true) ? FastInverse(this_mortar_condition_matrices.DOperator) : (std::abs(aux_det) > tolerance) ? MathUtils<double>::InvertMatrix<TNumNodes>(this_mortar_condition_matrices.DOperator, aux_det, tolerance) : ZeroMatrix(TNumNodes);
                        const bounded_matrix<double, TNumNodes, TNumNodes> p_operator =prod(inv_d_operator, this_mortar_condition_matrices.MOperator); 
                        
                        Matrix var_origin_matrix;
                        MortarUtilities::MatrixValue<TVarType, THist>(slave_geometry, mOriginVariable, var_origin_matrix);
        
                        const Matrix var_destiny_matrix = (area/master_geometry.Area()) * prod(p_operator, var_origin_matrix);
                        MortarUtilities::AddAreaWeightedValue<TVarType, THist>(master_geometry, mDestinationVariable, var_destiny_matrix);
                    }
                }
            }
        }
        
        KRATOS_CATCH("");
    }
    
    /**
     * This method executes the mapping when a linear solver is avalaible and a system of equations can be solved
     */
    void ExecuteImplicitMapping()
    {
        KRATOS_TRY;

        // Defining tolerance
//         const double relative_convergence_tolerance = 1.0e-4;
        const double absolute_convergence_tolerance = 1.0e-9;
        const unsigned int max_number_iterations = 10;
        unsigned int iteration = 0;
        
        // Creating the assemble database
        std::size_t size_system;
        std::unordered_map<int, int> conectivity_database, inverse_conectivity_database;
        CreateMasterConectivityDatabase(size_system, conectivity_database, inverse_conectivity_database);
        
        // Defining the operators
        const unsigned int variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
        std::vector<bool> is_converged(variable_size, false);
        MatrixType A = ZeroMatrix(size_system, size_system);
        const VectorType zero_vector = ZeroVector(size_system);
        std::vector<VectorType> b(variable_size, zero_vector);
        VectorType Dx(size_system);
        
        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> this_kinematic_variables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> this_mortar_condition_matrices;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        while (CheckWholeVector(is_converged) == false && iteration < max_number_iterations)
        {
            // We map the values from one side to the other // TODO: Add OMP
            for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); ++i)
            {
                auto it_cond = mrThisModelPart.ConditionsBegin() + i;
                
                if (it_cond->Is(SLAVE) == true)
                {
                    const array_1d<double, 3>& slave_normal = it_cond->GetValue(NORMAL);
                    GeometryType& slave_geometry = it_cond->GetGeometry();
                    
                    boost::shared_ptr<ConditionMap>& all_conditions_maps = it_cond->GetValue( MAPPING_PAIRS ); // These are the master conditions
                    
                    for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
                    {
                        Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                        const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                        GeometryType& master_geometry = p_cond_master->GetGeometry();
                        
                        IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_3;
                        
                        // Reading integration points
                        std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping  
                        const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);
                        
                        if (is_inside == true)
                        {
                            // Initialize general variables for the current master element
                            this_kinematic_variables.Initialize();
                    
                            // Initialize the mortar operators
                            this_mortar_condition_matrices.Initialize();
                            
                            for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom)
                            {
                                std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                                for (unsigned int i_node = 0; i_node < TDim; ++i_node)
                                {
                                    PointType global_point;
                                    slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                    points_array[i_node] = boost::make_shared<PointType>(global_point);
                                }
                                
                                typename std::conditional<TDim == 2, LineType, TriangleType >::type decomp_geom( points_array );
                                
                                const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                                
                                if (bad_shape == false)
                                {
                                    const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                                    
                                    // Integrating the mortar operators
                                    for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); ++point_number )
                                    {
                                        const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                        PointType local_point_parent;
                                        PointType gp_global;
                                        decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                        slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
            
                                        /// SLAVE CONDITION ///
                                        slave_geometry.ShapeFunctionsValues( this_kinematic_variables.NSlave, local_point_parent.Coordinates() );
                                        this_kinematic_variables.PhiLagrangeMultipliers = this_kinematic_variables.NSlave;

                                        this_kinematic_variables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );
                                        
                                        /// MASTER CONDITION ///
                                        PointType projected_gp_global;
                                        const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(this_kinematic_variables.NSlave, slave_geometry);
                                        
                                        GeometryType::CoordinatesArrayType slave_gp_global;
                                        slave_geometry.GlobalCoordinates( slave_gp_global, local_point_parent );
                                        MortarUtilities::FastProjectDirection( master_geometry, gp_global, projected_gp_global, master_normal, -gp_normal ); // The opposite direction
                                        
                                        GeometryType::CoordinatesArrayType projected_gp_local;
                                        
                                        master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
                                        
                                        // SHAPE FUNCTIONS 
                                        master_geometry.ShapeFunctionsValues( this_kinematic_variables.NMaster, projected_gp_local );    
                                        
                                        const double integration_weight = integration_points_slave[point_number].Weight();
                                        
                                        this_mortar_condition_matrices.CalculateMortarOperators(this_kinematic_variables, integration_weight);   
                                    }
                                }
                            }
                            
                            Matrix var_origin_matrix;
                            MortarUtilities::MatrixValue<TVarType, THist>(slave_geometry, mOriginVariable, var_origin_matrix);
                            Matrix var_destination_matrix;
                            MortarUtilities::MatrixValue<TVarType, THist>(master_geometry, mDestinationVariable, var_destination_matrix);
                            
                            Matrix residual_matrix = prod(this_mortar_condition_matrices.DOperator, var_origin_matrix) - prod(this_mortar_condition_matrices.MOperator, var_destination_matrix);
                            
                            /* We compute the residual and assemble */
                            // First we assemble the RHS
                            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
                            {
                                const int& node_i_id = master_geometry[i_node].Id();
                                const int& pos_i_id = inverse_conectivity_database[node_i_id];
                                for (unsigned int i_var = 0; i_var < variable_size; ++i_var)
                                {
                                    b[i_var][pos_i_id] += residual_matrix(i_node, i_var);
                                }
                                // We assemble the LHS
                                for (unsigned int j_node = 0; j_node < TNumNodes; ++j_node)
                                {
                                    const int& node_j_id = master_geometry[j_node].Id();
                                    const int& pos_j_id = inverse_conectivity_database[node_j_id];
                                    A(pos_i_id, pos_j_id) += this_mortar_condition_matrices.MOperator(i_node, j_node);
                                }
                            }
                        
                        }
                    }
                }
                
                // Finally we solve the system # TODO: Use iterations to minimize the residual
                for (unsigned int i_size = 0; i_size < variable_size; ++i_size)
                {
                    mpThisLinearSolver->Solve(A, Dx, b[i_size]);
                    MortarUtilities::UpdateDatabase<TVarType, THist>(mrThisModelPart, mDestinationVariable, Dx, i_size, conectivity_database);
                    const double residual_norm = norm_2(b[i_size])/size_system;
//                     const double increment_norm = norm_2(Dx)/size_system;
                    if (residual_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
//                     if (increment_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
                }
            }
            
            iteration += 1;
        }
        
        KRATOS_CATCH(""); 
    }
        
    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    
    /// Assignment operator.
    SimpleMortarMapperProcess& operator=(SimpleMortarMapperProcess const& rOther) = delete;

    /// Copy constructor.
    //SimpleMortarMapperProcess(SimpleMortarMapperProcess const& rOther);
    
    ///@}
};// class SimpleMortarMapperProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SIMPLE_MORTAR_MAPPER_PROCESS defined */
