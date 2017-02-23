// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED )
#define  KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "contact_structural_mechanics_application_variables.h"
/* Utilities */
#include "custom_utilities/contact_utilities.h"
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{
    
///@}
///@name Type Definitions
///@{
    
    typedef Point<3>                                             PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;
    
    ///Type definition for integration methods
    typedef GeometryType::IntegrationPointsArrayType IntegrationPointsType;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

/** \brief ExactMortarIntegrationUtility 
 * This utility calculates the exact integration necessary for the Mortar Conditions
 */
template< unsigned int TDim, unsigned int TNumNodes>
class ExactMortarIntegrationUtility
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ExactMortarIntegrationUtility
    KRATOS_CLASS_POINTER_DEFINITION(ExactMortarIntegrationUtility);
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Constructor
    /**
     * @param SlaveGeometry: The geometry of the slave condition
     * @param IntegrationOrder: The integration order to consider
     */
    ExactMortarIntegrationUtility(
        const GeometryType& SlaveGeometry,
        unsigned int IntegrationOrder = 0
    )
    :mSlaveGeometry(SlaveGeometry),
    mIntegrationOrder(IntegrationOrder)
    {
        // TODO: Add more if necessary
    }
    
    /// Destructor.
    virtual ~ExactMortarIntegrationUtility() {}
    
    ///@}
    ///@name Operators
    ///@{
    
    ///@}
    ///@name Operations
    ///@{    
    
    /**
     * This utility computes the exact integration of the mortar condition
     * @param MasterGeometry: The geometry of the master condition
     * @param IntegrationPointsSlave: The integrations points that belong to the slave
     * @return True if there is a coomon area (the geometries intersect), false otherwise
     */
    
    bool GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        );
    
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

    GeometryType mSlaveGeometry;      // The geometry of the slave condition
    unsigned int mIntegrationOrder;   // The integration order to consider
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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
}; // Class ExactMortarIntegrationUtility

///@name Explicit Specializations
///@{

    template<>  
    bool ExactMortarIntegrationUtility<2,2>::GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {             
        // Using exact integration
        const double tol = 1.0e-4; 
        IntegrationMethod AuxIntegrationMethod;
        
        if (mIntegrationOrder == 1)
        {
            AuxIntegrationMethod = GeometryData::GI_GAUSS_1;
        }
        else if (mIntegrationOrder == 2)
        {
            AuxIntegrationMethod = GeometryData::GI_GAUSS_2;
        }
        else if (mIntegrationOrder == 3)
        {
            AuxIntegrationMethod = GeometryData::GI_GAUSS_3;
        }
        else if (mIntegrationOrder == 4)
        {
            AuxIntegrationMethod = GeometryData::GI_GAUSS_4;
        }
        else if (mIntegrationOrder == 5)
        {
            AuxIntegrationMethod = GeometryData::GI_GAUSS_5;
        }
        else
        {
            AuxIntegrationMethod = mSlaveGeometry.GetDefaultIntegrationMethod();
        }
        
        double total_weight = 0.0;
        array_1d<double,2> coor_aux = ZeroVector(2);
        
        // Declaring auxiliar values
        PointType projected_gp_global;
        GeometryType::CoordinatesArrayType projected_gp_local;
        double aux_dist = 0.0;
        
        // Declare the boolean of full integral
        bool full_int = true;
        
        // First look if the edges of the slave are inside of the master, if not check if the opposite is true, if not then the element is not in contact
        for (unsigned int i_slave = 0; i_slave < 2; i_slave++)
        {
            const array_1d<double, 3> normal = mSlaveGeometry[i_slave].GetValue(NORMAL);
            
            ContactUtilities::ProjectDirection(MasterGeometry, mSlaveGeometry[i_slave].Coordinates(), projected_gp_global, aux_dist, -normal ); // The opposite direction
            
            const bool inside = MasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local );
            
            if (inside == false)
            {
                full_int = false;
            }
            else
            {
                if (i_slave == 0)
                {
                    coor_aux[0] = - 1.0;
                }
                else if (i_slave == 1)
                {
                    coor_aux[1] =   1.0;
                }
            }
        }
        
        if (full_int == true)
        {
            total_weight = 2.0;
        }
        else
        {
            std::vector<double> aux_xi;
            for (unsigned int i_master = 0; i_master < 2; i_master++)
            {
                // NOTE: Check this!!!!
                const bool inside = ContactUtilities::ProjectIterative(mSlaveGeometry, MasterGeometry[i_master].Coordinates(), projected_gp_global);
                
                if (inside == true)
                {
                    aux_xi.push_back(projected_gp_local[0]);
                }
            }
            
            if (aux_xi.size() == 1)
            {
                if (coor_aux[0] == - 1.0)
                {
                    coor_aux[1] = aux_xi[0];
                }
                else if (coor_aux[1] == 1.0)
                {
                    coor_aux[0] = aux_xi[0];
                }
                else
                {
                    KRATOS_WATCH("WARNING: THIS IS NOT SUPPOSED TO HAPPEN!!!!");
                }
            }
            else  if (aux_xi.size() == 2)
            {
                if (aux_xi[0] < aux_xi[1])
                {
                    coor_aux[0] = aux_xi[0];
                    coor_aux[1] = aux_xi[1];
                }
                else
                {
                    coor_aux[1] = aux_xi[0];
                    coor_aux[0] = aux_xi[1];
                }
            }
            
            total_weight = coor_aux[1] - coor_aux[0];
        }
        
        if(total_weight < 0.0)
        {
            KRATOS_ERROR << "WAAAAAAAAAAAAARNING!!!!!!!!, wrong order of the coordinates: "<< coor_aux << std::endl;
        }
        
        if (total_weight > tol)
//             if (total_weight > 0.0)
        {
            // With the proportion of the weigth you recalculate the integration weight. Change the coordinates of the integration to accomodate
            const GeometryType::IntegrationPointsArrayType& integration_points = mSlaveGeometry.IntegrationPoints(AuxIntegrationMethod);
            IntegrationPointsSlave.resize(integration_points.size());
            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
            {
                const double weight = integration_points[PointNumber].Weight() * total_weight/2.0;
                const double xi = 0.5 * (1.0 - integration_points[PointNumber].Coordinate(1)) * coor_aux[0] 
                                + 0.5 * (1.0 + integration_points[PointNumber].Coordinate(1)) * coor_aux[1];
                
                IntegrationPointsSlave[PointNumber] = IntegrationPoint<2>( xi, weight );
            }
        }
        else
        {
//                 IntegrationPointsSlave.resize(0); // An empty std::vector
            IntegrationPointsSlave.clear(); // An empty std::vector
        }
        
//             if (IntegrationPointsSlave.size() > 0)
//             {
//                 std::cout <<  GetGeometry()[0].X() << " " << GetGeometry()[0].Y() << " " << GetGeometry()[1].X() << " " << GetGeometry()[1].Y() << std::endl;
//                 std::cout <<  MasterGeometry[0].X() << " " << MasterGeometry[0].Y() << " " << MasterGeometry[1].X() << " " << MasterGeometry[1].Y() << std::endl;
//                 KRATOS_WATCH(coor_aux);

//                 std::cout << "IntegrationPoints : " << IntegrationPointsSlave.size( ) << std::endl;
//                 for ( unsigned int i_vec = 0; i_vec < IntegrationPointsSlave.size( ); ++i_vec )
//                 {
//                     KRATOS_WATCH( IntegrationPointsSlave[i_vec] );
//                 }
//             }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    bool ExactMortarIntegrationUtility<3,3>::GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        // TODO: Finish this
        // Compute the local Coordinates of the master condition
        PointType projected_gp_global;
        
        GeometryType::CoordinatesArrayType slave_gp_global;
        double aux_dist = 0.0;
        
        for (unsigned int i = 0; i < 3; i++)
        {
            const array_1d<double,3> normal = mSlaveGeometry.GetValue(NORMAL);
            
            mSlaveGeometry.GlobalCoordinates( slave_gp_global, local_point );
            ContactUtilities::ProjectDirection( MasterGeometry, slave_gp_global, projected_gp_global, aux_dist, -normal ); // The opposite direction
            
            GeometryType::CoordinatesArrayType projected_gp_local;
            
            const bool inside = MasterGeometry.IsInside( projected_gp_global.Coordinates( ), projected_gp_local ) ;
        }
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    template<>  
    bool ExactMortarIntegrationUtility<3,4>::GetExactIntegration(         
        GeometryType& MasterGeometry,
        IntegrationPointsType& IntegrationPointsSlave
        )
    {
        // TODO: Finish this!!!
    }
}

#endif  /* KRATOS_EXACT_MORTAR_INTEGRATION_UTILITY_H_INCLUDED defined */
