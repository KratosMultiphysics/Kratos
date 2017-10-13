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

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/point.h"

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& ThisVariable
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(ThisVariable),
           mDestinyVariable(ThisVariable)
    {
    }
    
    SimpleMortarMapperProcess( 
        ModelPart& rThisModelPart,
        TVarType& OriginVariable,
        TVarType& DestinyVariable
        ): mrThisModelPart(rThisModelPart),
           mOriginVariable(OriginVariable),
           mDestinyVariable(DestinyVariable)
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

        // Defining tolerance
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        // We compute the nodal area
        ComputeNodalArea();
        
        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> this_kinematic_variables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> this_mortar_condition_matrices;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(TDim);
        
        // Initialize to zero // TODO: Add OMP
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); i++)
        {
            auto it_cond = mrThisModelPart.ConditionsBegin() + i;
            
            if (it_cond->Is(SLAVE) == true)
            {
                boost::shared_ptr<ConditionMap>& all_conditions_maps = it_cond->GetValue( MAPPING_PAIRS ); // These are the master conditions
                
                for (auto it_pair = all_conditions_maps->begin(); it_pair != all_conditions_maps->end(); ++it_pair )
                {
                    Condition::Pointer p_cond_master = (it_pair->first); // MASTER
                    GeometryType& master_geometry = p_cond_master->GetGeometry();
                    
                    MortarUtilities::ResetValue<TVarType, THist>(master_geometry, mDestinyVariable);
                }
            }
        }
        
        // We map the values from one side to the other // TODO: Add OMP
        for(int i=0; i<static_cast<int>(mrThisModelPart.Conditions().size()); i++)
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
                        
                    const double total_master_area = master_geometry.Area();
                    double master_area = 0.0;
                    
                    IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;
                    
                    // Reading integration points
                    std::vector<array_1d<PointType,TDim>> conditions_points_slave;
                    const bool is_inside = integration_utility.GetExactIntegration(slave_geometry, slave_normal, master_geometry, master_normal, conditions_points_slave);
                    
                    if (is_inside == true)
                    {
                        // Initialize general variables for the current master element
                        this_kinematic_variables.Initialize();
                
                        // Initialize the mortar operators
                        this_mortar_condition_matrices.Initialize();
                        
                        for (unsigned int i_geom = 0; i_geom < conditions_points_slave.size(); i_geom++)
                        {
                            std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
                            for (unsigned int i_node = 0; i_node < TDim; i_node++)
                            {
                                PointType global_point;
                                slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                                points_array[i_node] = boost::make_shared<PointType>(global_point);
                            }
                            
                            typename std::conditional<TDim == 2, LineType, TriangleType >::type decomp_geom( points_array );
                            
                            const bool bad_shape = (TDim == 2) ? MortarUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : MortarUtilities::HeronCheck(decomp_geom);
                            
                            if (bad_shape == false)
                            {
                                master_area += decomp_geom.Area();
                                
                                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                                
                                // Integrating the mortar operators
                                for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                                {
                                    const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                                    PointType local_point_parent;
                                    PointType gp_global;
                                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                                    slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
           
                                    /// SLAVE CONDITION ///
                                    slave_geometry.ShapeFunctionsValues( this_kinematic_variables.NSlave, local_point_parent.Coordinates() );
                                    this_kinematic_variables.PhiLagrangeMultipliers = this_kinematic_variables.NSlave;
                                    this_kinematic_variables.DetjSlave = slave_geometry.DeterminantOfJacobian( local_point_parent );
                                    
                                    /// MASTER CONDITION ///
                                    PointType projected_gp_global;
                                    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointNormal(this_kinematic_variables.NSlave, slave_geometry);
                                    
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
                        
                        // We compute the norm
                        const double norm_D = norm_frobenius(this_mortar_condition_matrices.DOperator);
                        
                        // Now we normalize the matrix
                        const bounded_matrix<double, TNumNodes, TNumNodes> normalized_D = this_mortar_condition_matrices.DOperator/norm_D;
                        
                        double aux_det = MathUtils<double>::DetMat<TNumNodes>(normalized_D);
                        const bounded_matrix<double, TNumNodes, TNumNodes> inv_d_operator = (std::abs(aux_det) < tolerance) ? ZeroMatrix(TNumNodes) : MathUtils<double>::InvertMatrix<TNumNodes>(normalized_D, aux_det, tolerance);
                        const bounded_matrix<double, TNumNodes, TNumNodes> p_operator = 1.0/norm_D * prod(inv_d_operator, this_mortar_condition_matrices.MOperator); 
                        
                        Matrix var_origin_matrix;
                        MortarUtilities::MatrixValue<TVarType, THist>(slave_geometry, mDestinyVariable, var_origin_matrix);
        
                        const Matrix var_destiny_matrix = (master_area/total_master_area) * prod(p_operator, var_origin_matrix);
                        MortarUtilities::AddValue<TVarType, THist>(master_geometry, mDestinyVariable, var_destiny_matrix);
                    }
                }
            }
        }
        
        KRATOS_CATCH("");
    }

    /**
     * This method sets both variables (origin and destination) with the same variable
     */
    void SetVariable(TVarType ThisVariable)
    {
        mOriginVariable = ThisVariable;
        mDestinyVariable = ThisVariable;
    }
    
    /**
     * This method sets both variables (origin and destination) in a separated way
     */
    void SetVariables(
        TVarType OriginVariable,
        TVarType DestinyVariable
        )
    {
        mOriginVariable = OriginVariable;
        mDestinyVariable = DestinyVariable;
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
    
    ModelPart& mrThisModelPart; // The model part to compute
    TVarType mOriginVariable;   // The origin variable to map
    TVarType mDestinyVariable;  // The destiny variable to map

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
        for(int i = 0; i < num_nodes; i++) 
        {
            auto it_node = nodes_array.begin() + i;
            it_node->SetValue(NODAL_AREA, 0.0);
        }
        
        // Sum all the nodes areas
        ConditionsArrayType& conditions_array = mrThisModelPart.Conditions();
        const int num_conditions = static_cast<int>(conditions_array.size());
        
        #pragma omp parallel for
        for(int i = 0; i < num_conditions; i++) 
        {
            auto it_cond = conditions_array.begin() + i;
            
            const unsigned int number_nodes = it_cond->GetGeometry().PointsNumber();
            const double& rArea = it_cond->GetGeometry().Area()/number_nodes;
            
            for (unsigned int i = 0; i < number_nodes; i++)
            {
                #pragma omp atomic
                it_cond->GetGeometry()[i].GetValue(NODAL_AREA) += rArea;
            }
        }
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
