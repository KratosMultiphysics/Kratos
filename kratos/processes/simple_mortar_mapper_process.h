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
#include "includes/kratos_parameters.h"
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
    
    typedef Point                                        PointType;
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

    typedef bounded_matrix<double, TNumNodes, TNumNodes>  BoundedMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& ThisVariable, 
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(ThisVariable),
           mDestinationVariable(ThisVariable),
           mThisParameters(ThisParameters),
           mpThisLinearSolver(pThisLinearSolver)
    {
        Parameters DefaultParameters = Parameters(R"(
        {
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2
        })" );
        
        mThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mEchoLevel = mThisParameters["echo_level"].GetInt();
    }
    
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& OriginVariable,
        TVarType& DestinationVariable,
        Parameters ThisParameters = Parameters(R"({})" ),
        LinearSolverType::Pointer pThisLinearSolver = nullptr
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(OriginVariable),
           mDestinationVariable(DestinationVariable),
           mThisParameters(ThisParameters),
           mpThisLinearSolver(pThisLinearSolver)
    {
        Parameters DefaultParameters = Parameters(R"(
        {
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2
        })" );
        
        mThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mEchoLevel = mThisParameters["echo_level"].GetInt();
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
    
    unsigned int mEchoLevel;                      // The verbosity level
    Parameters mThisParameters;                   // The configuration parameters
    
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
     * This method assemble locally the mortar operators
     * @param ConditionsPointSlave: The list of points that form the triangle decomposition
     * @param SlaveGeometry: The slave geometry
     * @param MasterGeometry: The master geometry
     * @param MasterNormal: The normal vector of the master geometry
     * @param ThisKinematicVariables: The kinematic variables of the geometries, needed to integrate the mortar operators
     * @param ThisMortarOperators: The mortar operators
     * @param ThisIntegrationMethod: The integration method used, determines the integration order
     * @param Ae: The dual lagrange multiplier operator
     */
    void AssemblyMortarOperators(
        const std::vector<array_1d<PointType,TDim>>& ConditionsPointSlave,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const array_1d<double, 3>& MasterNormal,
        MortarKinematicVariables<TNumNodes>& ThisKinematicVariables,
        MortarOperator<TNumNodes>& ThisMortarOperators,
        const IntegrationMethod& ThisIntegrationMethod,
        const BoundedMatrixType Ae = IdentityMatrix(TNumNodes)
        )
    {
        for (unsigned int i_geom = 0; i_geom < ConditionsPointSlave.size(); ++i_geom)
        {
            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
            for (unsigned int i_node = 0; i_node < TDim; ++i_node)
            {
                PointType global_point;
                SlaveGeometry.GlobalCoordinates(global_point, ConditionsPointSlave[i_geom][i_node]);
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
    
    /**
     * This method computes the Ae matrix
     * @param SlaveGeometry: The slave geometry
     * @param ThisKinematicVariables: The kinematic variables
     * @param ConditionsPointsSlave: The list of decomposed triangles
     * @param ThisIntegrationMethod: The integration method considered
     * @return Ae: The matrix of dual LM
     */
    static inline BoundedMatrixType CalculateAe(
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
     * @param InputMatrix: The matrix to invert
     * @return The matrix inverted
     */
    static inline BoundedMatrixType FastInverse(const BoundedMatrixType& InputMatrix)
    {
        BoundedMatrixType inv_matrix = ZeroMatrix(TNumNodes);
        
        for (std::size_t i = 0; i < TNumNodes; ++i)
        {
            inv_matrix(i, i) = 1.0/InputMatrix(i, i);
        }
        
        return inv_matrix;
    }
    
    /**
     * This method inverts a diagonal matrix
     * @param InputMatrix: The matrix to invert
     * @param InvertedMatrix: The matrix inverted
     */
    static inline void FastInverse(
        const BoundedMatrixType& InputMatrix,
        BoundedMatrixType& InvertedMatrix
        )
    {
        InvertedMatrix = ZeroMatrix(TNumNodes);
        
        for (std::size_t i = 0; i < TNumNodes; ++i)
        {
            InvertedMatrix(i, i) = 1.0/InputMatrix(i, i);
        }
    }
    
    /**
     * This method lumps a matrix
     * @param InputMatrix: The matrix to lump
     */
    void LumpMatrix(BoundedMatrixType& InputMatrix)
    {        
        for (std::size_t i = 0; i < TNumNodes; ++i)
        {
            for (std::size_t j = 0; j < TNumNodes; ++j)
            {
                if (i != j) 
                {
                    InputMatrix(i, i) += InputMatrix(i, j);
                    InputMatrix(i, j) = 0.0;
                }
            }
        }
    }
    
    /**
     * This method creates a slave database needed to assemble the system
     * @param SizeSystem: The size of the system
     * @param ConectivityDatabase: The database that will be used to assemble the system
     * @param InverseConectivityDatabase: The inverse database that will be used to assemble the system
     */
        
    void CreateSlaveConectivityDatabase(
        std::size_t& SizeSystem,
        std::unordered_map<int, int>& InverseConectivityDatabase
        )
    {
        // Initialize the value
        SizeSystem = 0;
        
        NodesArrayType& nodes_array = mrThisModelPart.Nodes();
        const int num_nodes = static_cast<int>(nodes_array.size()); 
        
        // We create the database
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(SLAVE) == true)
            {
                InverseConectivityDatabase[it_node->Id()] = SizeSystem;
                SizeSystem += 1;
            }
        }
    }
    
    /**
     * This method creates a slave database needed to assemble the system
     * @param SizeSystem: The size of the system
     * @param ConectivityDatabase: The database that will be used to assemble the system
     * @param InverseConectivityDatabase: The inverse database that will be used to assemble the system
     */
        
    void CreateSlaveConectivityDatabase(
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
        for(int i = 0; i < num_nodes; ++i) 
        {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(SLAVE) == true)
            {
                ConectivityDatabase[SizeSystem] = it_node->Id();
                InverseConectivityDatabase[it_node->Id()] = SizeSystem;
                SizeSystem += 1;
            }
        }
    }
    
    /**
     * This method returns the corresponding integration order considered
     * @return The integration order considered
     */
    IntegrationMethod GetIntegrationMethod()
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
     * This method computes the residual matrix of the mapping
     * @param ResidualMatrix: The matrix containing the residual of the mappping
     * @param SlaveGeometry: The slave geometry
     * @param MasterGeometry: The master geometry
     * @param ThisMortarOperators: The mortar operators
     */
    void ComputeResidualMatrix(       
        Matrix& ResidualMatrix,
        GeometryType& SlaveGeometry,
        GeometryType& MasterGeometry,
        const MortarOperator<TNumNodes>& ThisMortarOperators
        )
    {
        Matrix var_origin_matrix;
        MortarUtilities::MatrixValue<TVarType, THist>(MasterGeometry, mOriginVariable, var_origin_matrix);
        Matrix var_destination_matrix;
        MortarUtilities::MatrixValue<TVarType, THist>(SlaveGeometry, mDestinationVariable, var_destination_matrix);
     
        const std::size_t size_1 = var_origin_matrix.size1();
        const std::size_t size_2 = var_origin_matrix.size2();
        if (ResidualMatrix.size1() != size_1  || ResidualMatrix.size2() !=  size_2)
        {
            ResidualMatrix.resize(size_1, size_2, false);
        }
        
        noalias(ResidualMatrix) = prod(ThisMortarOperators.MOperator, var_origin_matrix) - prod(ThisMortarOperators.DOperator, var_destination_matrix);
    }
    
    /**
     * This method assembles the LHS and the RHS
     * @param A: The LHS of the system
     * @param b: The RHS of the system
     * @param VariableSize: The size of the variable
     * @param ResidualMatrix: The matrix containing the residual of the mappping
     * @param SlaveGeometry: The slave geometry
     * @param InverseConectivityDatabase: The inverse database that will be used to assemble the system
     * @param ThisMortarOperators: The mortar operators
     */
    void AssembleRHSAndLHS(
        MatrixType& A,
        std::vector<VectorType>& b,
        const unsigned int& VariableSize,
        const Matrix& ResidualMatrix,
        GeometryType& SlaveGeometry,
        std::unordered_map<int, int>& InverseConectivityDatabase,
        const MortarOperator<TNumNodes>& ThisMortarOperators
        )
    {
        // First we assemble the RHS
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const int& node_i_id = SlaveGeometry[i_node].Id();
            const int& pos_i_id = InverseConectivityDatabase[node_i_id];
            for (unsigned int i_var = 0; i_var < VariableSize; ++i_var)
            {
                b[i_var][pos_i_id] += ResidualMatrix(i_node, i_var);
            }

            // We assemble the LHS
            for (unsigned int j_node = 0; j_node < TNumNodes; ++j_node)
            {
                const int& node_j_id = SlaveGeometry[j_node].Id();
                const int& pos_j_id = InverseConectivityDatabase[node_j_id];
                A(pos_i_id, pos_j_id) += ThisMortarOperators.DOperator(i_node, j_node);
            }
        }
    }
    
    /**
     * This method assembles the RHS
     * @param b: The RHS of the system
     * @param VariableSize: The size of the variable
     * @param ResidualMatrix: The matrix containing the residual of the mappping
     * @param SlaveGeometry: The slave geometry
     * @param InverseConectivityDatabase: The inverse database that will be used to assemble the system
     */
    void AssembleRHS(
        std::vector<VectorType>& b,
        const unsigned int& VariableSize,
        const Matrix& ResidualMatrix,
        GeometryType& SlaveGeometry,
        std::unordered_map<int, int>& InverseConectivityDatabase
        )
    {
        // First we assemble the RHS
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const int& node_i_id = SlaveGeometry[i_node].Id();
            const int& pos_i_id = InverseConectivityDatabase[node_i_id];
            for (unsigned int i_var = 0; i_var < VariableSize; ++i_var)
            {
                b[i_var][pos_i_id] += ResidualMatrix(i_node, i_var);
            }
        }
    }
    
    /**
     * This method executes the explicit mapping (when no linear solver is avalaible)
     */
    void ExecuteExplicitMapping() // TODO: Correct, many mistakes
    {
        KRATOS_TRY;

        // Defining tolerance
        const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
        const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
        const unsigned int max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
        unsigned int iteration = 0;
        
        // We compute the nodal area
        ComputeNodalArea();
        
        // We set to zero the variables
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); ++i)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                GeometryType& slave_geometry = it_cond->GetGeometry();
                MortarUtilities::ResetValue<TVarType, THist>(slave_geometry, mDestinationVariable);
            }
        }
        
        // Creating the assemble database
        std::size_t size_system;
        std::unordered_map<int, int> inverse_conectivity_database;
        CreateSlaveConectivityDatabase(size_system, inverse_conectivity_database);
        
        // Defining the operators
        const unsigned int variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
        std::vector<bool> is_converged(variable_size, false);
        const VectorType zero_vector = ZeroVector(size_system);
        std::vector<VectorType> b(variable_size, zero_vector);
        std::vector<double> norm_b0(variable_size, 0.0);
        std::vector<double> norm_bi(variable_size, 0.0);
        double increment_residual_norm = 0.0;
        
        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> this_kinematic_variables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> this_mortar_operators;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        while (CheckWholeVector(is_converged) == false && iteration < max_number_iterations)
        {
            // We reset the RHS
            if (iteration > 0)
            {
                for (unsigned int i_size = 0; i_size < variable_size; ++i_size) 
                {
                    b[i_size] = zero_vector;
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
                        Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                        const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                        GeometryType& master_geometry = p_cond_master->GetGeometry();
                        
                        const IntegrationMethod& this_integration_method = GetIntegrationMethod();
                        
                        // Reading integration points
                        std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping  
                        const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);
                        
                        if (is_inside == true)
                        {
                            // Initialize general variables for the current master element
                            this_kinematic_variables.Initialize();
                    
                            // Initialize the mortar operators
                            this_mortar_operators.Initialize();
                            
                            const BoundedMatrixType& Ae = CalculateAe(slave_geometry, this_kinematic_variables, conditions_points_slave, this_integration_method);
                            
                            AssemblyMortarOperators( conditions_points_slave, slave_geometry, master_geometry,master_normal ,this_kinematic_variables, this_mortar_operators, this_integration_method, Ae);
                            
                            Matrix residual_matrix;
                            ComputeResidualMatrix(residual_matrix, slave_geometry, master_geometry, this_mortar_operators);
                            
                            /* We compute the residual and assemble */
                            AssembleRHS(b, variable_size, residual_matrix, slave_geometry, inverse_conectivity_database);
                            
                            BoundedMatrixType inverse_D_operator;
                            FastInverse(this_mortar_operators.DOperator, inverse_D_operator);
                            
                            if (mEchoLevel > 1)
                            {
                                BoundedMatrixType aux_copy_D = this_mortar_operators.DOperator;
                                LumpMatrix(aux_copy_D);
                                const BoundedMatrixType aux_diff = aux_copy_D - this_mortar_operators.DOperator;
                                const double norm_diff = norm_frobenius(aux_diff);
                                if (norm_diff > 1.0e-4) 
                                {
                                    std::cout << "WARNING: THE MORTAR OPERATOR D IS NOT DIAGONAL" << std::endl;
                                }
                                if (mEchoLevel == 3) 
                                {
                                    KRATOS_WATCH(norm_diff);
                                    KRATOS_WATCH(this_mortar_operators.DOperator);
                                }
                            }
//                             double aux_det;
//                             BoundedMatrixType inverse_D_operator = MathUtils<double>::InvertMatrix<TNumNodes>(this_mortar_operators.DOperator, aux_det);
                            const Matrix solution_matrix = prod(inverse_D_operator, residual_matrix);
                            MortarUtilities::AddAreaWeightedValue<TVarType, THist>(slave_geometry, mDestinationVariable, solution_matrix);
                        }
                    }
                }
            }
            
            // Finally we compute the residual
            for (unsigned int i_size = 0; i_size < variable_size; ++i_size)
            {
                if (mEchoLevel == 4)
                {
                    KRATOS_WATCH(b[i_size])
                }
                const double residual_norm = norm_2(b[i_size])/size_system;
                if (iteration == 0) norm_b0[i_size] = residual_norm;
                if (residual_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
                if (residual_norm/norm_b0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;
                if (iteration > 0) 
                {   
                    increment_residual_norm = std::abs((residual_norm - norm_bi[i_size])/norm_bi[i_size]);
                    if (increment_residual_norm < relative_convergence_tolerance) is_converged[i_size] = true;
                }
                norm_bi[i_size] = residual_norm;
                if (mEchoLevel > 0)
                {
                    std::cout << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm << "\tRELATIVE: " << residual_norm/norm_b0[i_size] << "\tINCREMENT: " << increment_residual_norm << std::endl;
                }
            }
            
            iteration += 1;
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
        const double relative_convergence_tolerance = mThisParameters["relative_convergence_tolerance"].GetDouble();
        const double absolute_convergence_tolerance = mThisParameters["absolute_convergence_tolerance"].GetDouble();
        const unsigned int max_number_iterations = mThisParameters["max_number_iterations"].GetInt();
        unsigned int iteration = 0;
        
        // We set to zero the variables
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); ++i)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                GeometryType& slave_geometry = it_cond->GetGeometry();
                MortarUtilities::ResetValue<TVarType, THist>(slave_geometry, mDestinationVariable);
            }
        }
        
        // Creating the assemble database
        std::size_t size_system;
        std::unordered_map<int, int> conectivity_database, inverse_conectivity_database;
        CreateSlaveConectivityDatabase(size_system, conectivity_database, inverse_conectivity_database);
        
        // Defining the operators
        const unsigned int variable_size = MortarUtilities::SizeToCompute<TDim, TVarType>();
        std::vector<bool> is_converged(variable_size, false);
        MatrixType A = ZeroMatrix(size_system, size_system);
        const VectorType zero_vector = ZeroVector(size_system);
        std::vector<VectorType> b(variable_size, zero_vector);
        std::vector<double> norm_b0(variable_size, 0.0);
        std::vector<double> norm_Dx0(variable_size, 0.0);
        VectorType Dx(size_system);
        
        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> this_kinematic_variables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> this_mortar_operators;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        while (CheckWholeVector(is_converged) == false && iteration < max_number_iterations)
        {
            // We reset the RHS
            if (iteration > 0)
            {
                for (unsigned int i_size = 0; i_size < variable_size; ++i_size) 
                {
                    b[i_size] = zero_vector;
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
                        Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                        const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                        GeometryType& master_geometry = p_cond_master->GetGeometry();
                        
                        const IntegrationMethod& this_integration_method = GetIntegrationMethod();
                        
                        // Reading integration points
                        std::vector<array_1d<PointType,TDim>> conditions_points_slave; // These are the segmentation points, with this points it is possible to create the lines or triangles used on the mapping  
                        const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);
                        
                        if (is_inside == true)
                        {
                            // Initialize general variables for the current master element
                            this_kinematic_variables.Initialize();
                    
                            // Initialize the mortar operators
                            this_mortar_operators.Initialize();
                            
//                             const BoundedMatrixType Ae = CalculateAe(slave_geometry, this_kinematic_variables, conditions_points_slave, this_integration_method);
                            
                            AssemblyMortarOperators( conditions_points_slave, slave_geometry, master_geometry,master_normal ,this_kinematic_variables, this_mortar_operators, this_integration_method);
                            
                            Matrix residual_matrix;
                            ComputeResidualMatrix(residual_matrix, slave_geometry, master_geometry, this_mortar_operators);
                            
                            /* We compute the residual and assemble */
                        
                            if (iteration == 0)
                            {
                                AssembleRHSAndLHS(A, b, variable_size, residual_matrix, slave_geometry, inverse_conectivity_database, this_mortar_operators);
                            }
                            else
                            {
                                AssembleRHS(b, variable_size, residual_matrix, slave_geometry, inverse_conectivity_database);
                            }
                        }
                    }
                }
            }
            
            // Finally we solve the system
            for (unsigned int i_size = 0; i_size < variable_size; ++i_size)
            {
                mpThisLinearSolver->Solve(A, Dx, b[i_size]);
                MortarUtilities::UpdateDatabase<TVarType, THist>(mrThisModelPart, mDestinationVariable, Dx, i_size, conectivity_database);
                const double residual_norm = norm_2(b[i_size])/size_system;
                if (iteration == 0) norm_b0[i_size] = residual_norm;
                const double increment_norm = norm_2(Dx)/size_system;
                if (iteration == 0) norm_Dx0[i_size] = increment_norm;
                if (residual_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
                if (residual_norm/norm_b0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;
                if (increment_norm < absolute_convergence_tolerance) is_converged[i_size] = true;
                if (increment_norm/norm_Dx0[i_size] < relative_convergence_tolerance) is_converged[i_size] = true;
                
                if (mEchoLevel > 0)
                {
                    std::cout << "Iteration: " << iteration + 1 << "\tRESISUAL::\tABS: " << residual_norm << "\tRELATIVE: " << residual_norm/norm_b0[i_size] << "\tINCREMENT::\tABS" << increment_norm << "\tRELATIVE: " << increment_norm/norm_Dx0[i_size] << std::endl;
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
