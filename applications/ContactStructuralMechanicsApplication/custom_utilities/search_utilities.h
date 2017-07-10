// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_SEARCH_UTILITIES)
#define KRATOS_SEARCH_UTILITIES

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "geometries/point.h"

/* Custom includes */
#include "custom_includes/mortar_operator.h"
// #include "custom_includes/dual_LM_operators.h" // NOTE: This is expensive, so we wil, use the standard shape function (this is to compute the gap, not to assemble the system, so it is acceptable)
#include "custom_includes/mortar_kinematic_variables.h"
#include "custom_includes/point_item.h"

/* Custom utilities */
#include "custom_utilities/exact_mortar_segmentation_utility.h"
#include "custom_utilities/contact_utilities.h"

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
    
class SearchUtilities
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef Point<3>                                     PointType;
    typedef Node<3>                                       NodeType;
    typedef Geometry<NodeType>                        GeometryType;
    typedef Geometry<PointType>                  GeometryPointType;
    typedef ModelPart::NodesContainerType           NodesArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;
    
    // Type definition for integration methods
    typedef GeometryData::IntegrationMethod      IntegrationMethod;
    
    // Type definitions for the point item
    typedef PointItem                                PointItemType;
    typedef PointItemType::Pointer            PointItemTypePointer;
    typedef std::vector<PointItemTypePointer>      PointItemVector;
    
    // Auxiliar geometries
    typedef Line2D2<PointType>                            LineType;
    typedef Triangle3D3<PointType>                    TriangleType;
    
    ///@}
    ///@name Life Cycle
    ///@{

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
    ///@name Operations
    ///@{
    
    /**
     * This function fills the ConditionSet for the Mortar condition // LEGACY WAY
     * @param ConditionPointers: The map storing the potential contact conditions
     * @param pCond1: The condition pointer of the slave 
     * @param pCond2: The  condition pointer of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @param DualCheck: The threshold distance to check the potential contact
     * @param StrictCheck: If the node must be inside or not
     */
    
    template< const bool TFill>
    static inline void ContactContainerFiller(
        boost::shared_ptr<ConditionSet>& ConditionPointers,
        Condition::Pointer & pCond1,       // SLAVE
        const Condition::Pointer & pCond2, // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckLength,
        const bool DualCheck = false, 
        const bool StrictCheck = true
        )
    {
        // Define the basic information
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        // Initialize geometries
        GeometryType& geom_1 = pCond1->GetGeometry(); // SLAVE
        GeometryType& geom_2 = pCond2->GetGeometry(); // MASTER
        
        const bool condition_is_active = CheckGeometryNodes( geom_1, geom_2, ContactNormal1, ContactNormal2, ActiveCheckLength, DualCheck, StrictCheck, tolerance);
        
        // If condition is active we add
        if (condition_is_active == true && TFill == true)
        {
            ConditionPointers->AddNewCondition(pCond2);
        }
        if (condition_is_active == false && TFill == false)
        {
            ConditionPointers->RemoveCondition(pCond2);
        }
    }
    
    /**
     * This function checks if the geometry can be potentially in contact using exact integration
     * @param rVariables: the kinematic variables
     * @param rThisMortarConditionMatrices: The mortar operators
     * @param IntegrationUtility: The exact integration utility
     * @param pCondSlave: The pointer of the slave
     * @param pCondMaster: The pointer of the master
     * @param SlaveNormal: The normals of the slave
     * @param MasterNormal: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @return condition_is_active: True if at least one node is active, false otherwise
     */
    template< const unsigned int TDim, const unsigned int TNumNodes >
    static inline bool CheckExactIntegration(
        MortarKinematicVariables<TNumNodes> rVariables,
        MortarOperator<TNumNodes> rThisMortarConditionMatrices,   
        ExactMortarIntegrationUtility<TDim, TNumNodes> IntegrationUtility,
        Condition::Pointer& pCondSlave,          // SLAVE
        Condition::Pointer& pCondMaster,         // MASTER
        const array_1d<double, 3>& SlaveNormal,  // SLAVE
        const array_1d<double, 3>& MasterNormal, // SLAVE
        const double ActiveCheckLength
        )
    {
        // Define the basic information
        bool condition_is_active = false;
        IntegrationMethod this_integration_method = GeometryData::GI_GAUSS_2;

        GeometryType& slave_geometry = pCondSlave->GetGeometry();
        GeometryType& master_geometry = pCondMaster->GetGeometry();

        // Reading integration points
        std::vector<array_1d<PointType,TDim>> conditions_points_slave;
        const bool is_inside = IntegrationUtility.GetExactIntegration(slave_geometry, SlaveNormal, master_geometry, MasterNormal, conditions_points_slave);
        
        if (is_inside == true)
        {
            // Initialize general variables for the current master element
            rVariables.Initialize();
    
            // Initialize the mortar operators
            rThisMortarConditionMatrices.Initialize();
            
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
                
                const bool bad_shape = (TDim == 2) ? ContactUtilities::LengthCheck(decomp_geom, slave_geometry.Length() * 1.0e-6) : ContactUtilities::HeronCheck(decomp_geom);
                
                if (bad_shape == false)
                {
                    const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );
                    
                    // Integrating the mortar operators
                    for ( unsigned int point_number = 0; point_number < integration_points_slave.size(); point_number++ )
                    {
                        const PointType local_point_decomp = integration_points_slave[point_number].Coordinates();
                        PointType local_point_parent;
                        PointType gp_global;
                        decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                        slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);
                        
                        // Calculate the kinematic variables
                        const PointType& LocalPoint = integration_points_slave[point_number].Coordinates();

                        /// SLAVE CONDITION ///
                        slave_geometry.ShapeFunctionsValues( rVariables.NSlave, LocalPoint.Coordinates() );
                        rVariables.PhiLagrangeMultipliers = rVariables.NSlave;
                        rVariables.DetjSlave = slave_geometry.DeterminantOfJacobian( LocalPoint );
                        
                        /// MASTER CONDITION ///
                        PointType projected_gp_global;
                        const array_1d<double,3> gp_normal = ContactUtilities::GaussPointNormal(rVariables.NSlave, slave_geometry);
                        
                        GeometryType::CoordinatesArrayType slave_gp_global;
                        slave_geometry.GlobalCoordinates( slave_gp_global, LocalPoint );
                        ContactUtilities::FastProjectDirection( master_geometry, slave_gp_global, projected_gp_global, MasterNormal, -gp_normal ); // The opposite direction
                        
                        GeometryType::CoordinatesArrayType projected_gp_local;
                        
                        master_geometry.PointLocalCoordinates(projected_gp_local, projected_gp_global.Coordinates( ) ) ;
                        
                        // SHAPE FUNCTIONS 
                        master_geometry.ShapeFunctionsValues( rVariables.NMaster, projected_gp_local );    
                        
                        const double integration_weight = integration_points_slave[point_number].Weight(); // NOTE: Error in axisymmetric
                        
                        rThisMortarConditionMatrices.CalculateMortarOperators(rVariables, integration_weight);   
                    }

                    // Setting the gap

                    // Current coordinates 
                    const bounded_matrix<double, TNumNodes, TDim> x1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(slave_geometry);
                    const bounded_matrix<double, TNumNodes, TDim> x2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(master_geometry);
            
                    const bounded_matrix<double, TNumNodes, TDim> Dx1Mx2 = prod(rThisMortarConditionMatrices.DOperator, x1) - prod(rThisMortarConditionMatrices.MOperator, x2); 
                    
                    array_1d<double, TDim> aux_slave_normal;
                    for (unsigned int i_dim = 0; i_dim < TDim; i_dim++)
                    {
                        aux_slave_normal[i_dim] = SlaveNormal[i_dim];
                    }
                    
                    for (unsigned int i_node = 0; i_node < TNumNodes; i_node++)
                    {
                        if (slave_geometry[i_node].Is(ACTIVE) == false)
                        {
                            const array_1d<double, TDim> aux_array = row(Dx1Mx2, i_node);
                            
                            const double nodal_gap = inner_prod(aux_array, - aux_slave_normal); 
                            
                            if (nodal_gap < ActiveCheckLength)
                            {
                                slave_geometry[i_node].Set(ACTIVE, true);
                                condition_is_active = true;
                            }
                        }
                        else
                        {
                            condition_is_active = true;
                        }
                    }
                }
            }
        } 
        
        return condition_is_active;
    }
    
    /**
     * This function checks the ConditionSet for the Mortar condition 
     * @param ConditionPointers: The map storing the potential contact conditions
     * @param pCondSlave: The condition pointer of the slave 
     * @param SlaveNormal: The normals of the slave
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     */
    
    template< const unsigned int TDim, const unsigned int TNumNodes >
    static inline void ExactContactContainerChecker(
        boost::shared_ptr<ConditionSet>& ConditionPointers,
        Condition::Pointer& pCondSlave,         // SLAVE
        const array_1d<double, 3>& SlaveNormal, // SLAVE
        const double ActiveCheckLength
        )
    {
        // Define the basic information
        bool condition_is_active = false;
        
        // Create and initialize condition variables:
        MortarKinematicVariables<TNumNodes> rVariables;
    
        // Create the mortar operators
        MortarOperator<TNumNodes> rThisMortarConditionMatrices;
        
        // We call the exact integration utility
        ExactMortarIntegrationUtility<TDim, TNumNodes> integration_utility = ExactMortarIntegrationUtility<TDim, TNumNodes>(2);
        
        for (auto it_pair = ConditionPointers->begin(); it_pair != ConditionPointers->end(); ++it_pair )
        {  
            Condition::Pointer p_cond_master = *(it_pair); // MASTER
            const array_1d<double, 3>& master_normal = p_cond_master->GetValue(NORMAL); 
                    
            condition_is_active = CheckExactIntegration<TDim, TNumNodes>(rVariables, rThisMortarConditionMatrices, integration_utility, pCondSlave, p_cond_master, SlaveNormal, master_normal, ActiveCheckLength);
            
            // If condition is active we add
            if (condition_is_active == false)
            {
                ConditionPointers->RemoveCondition(p_cond_master);
            }
        }
    }
    
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    /**
     * This function checks the nodes of the geometry if they can be potentially in contact
     * @param Geom1: The geometry of the slave 
     * @param Geom1: The geometry of the master 
     * @param ContactNormal1: The normals of the slave
     * @param ContactNormal2: The normals of the master
     * @param ActiveCheckLength: The threshold distance to check the potential contact
     * @param DualCheck: The threshold distance to check the potential contact
     * @param StrictCheck: If the node must be inside or not
     * @param tolerance: The tolerance considered in the calculations
     * @return at_least_one_node_potential_contact: True if at least one node is active, false otherwise
     */
    static inline bool CheckGeometryNodes(
        GeometryType& Geom1, // SLAVE
        GeometryType& Geom2,  // MASTER
        const array_1d<double, 3> & ContactNormal1, // SLAVE
        const array_1d<double, 3> & ContactNormal2, // MASTER
        const double ActiveCheckLength,
        const bool DualCheck = false, 
        const bool StrictCheck = true,
        const double tolerance = std::numeric_limits<double>::epsilon()
        )
    {
        bool at_least_one_node_potential_contact = false;
        for (unsigned int i_node = 0; i_node < Geom1.size(); i_node++)
        {
            if (Geom1[i_node].Is(ACTIVE) == false)
            {
                Point<3> projected_point;
                double aux_distance = 0.0;
                const array_1d<double, 3> normal = Geom1[i_node].GetValue(NORMAL);
                if (norm_2(normal) < tolerance)
                {
                    aux_distance = ContactUtilities::FastProjectDirection(Geom2, Geom1[i_node], projected_point, ContactNormal2, ContactNormal1);
                }
                else
                {
                    aux_distance = ContactUtilities::FastProjectDirection(Geom2, Geom1[i_node], projected_point, ContactNormal2, normal);
                }  
                
                array_1d<double, 3> result;
                if (aux_distance <= ActiveCheckLength && (StrictCheck == true ? Geom2.IsInside(projected_point, result, tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                {
                    at_least_one_node_potential_contact = true;
                    
                    Geom1[i_node].Set(ACTIVE, true);
                }
                
                if (DualCheck == true)
                {
                    aux_distance = ContactUtilities::FastProjectDirection(Geom2, Geom1[i_node], projected_point, ContactNormal2, -ContactNormal2);
                    if (aux_distance <= ActiveCheckLength && (StrictCheck == true ? Geom2.IsInside(projected_point, result, tolerance) : true)) // NOTE: This can be problematic (It depends the way IsInside() and the LocalPointCoordinates() are implemented)
                    {
                        at_least_one_node_potential_contact = true;
                        
                        Geom1[i_node].Set(ACTIVE, true);
                    }
                }
             }
             else
             {
                 at_least_one_node_potential_contact = true;
             }
         }
         
         return at_least_one_node_potential_contact;
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
    ///@}
};// class SearchUtilities

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
#endif /* KRATOS_SEARCH_UTILITIES defined */
