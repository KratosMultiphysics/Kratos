 
// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_MORTAR_KINEMATIC_VARIABLES )
#define  KRATOS_MORTAR_KINEMATIC_VARIABLES

// System includes

// External includes

// Project includes
// #include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"

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
    
    // Type definition for integration methods
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
    
/** \brief MortarKinematicVariables
 * This is the definition of the kinematic variables
 */

template< const unsigned int TNumNodes>
class MortarKinematicVariables
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of MortarKinematicVariables
    KRATOS_CLASS_POINTER_DEFINITION( MortarKinematicVariables );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarKinematicVariables(){}
    
    virtual ~MortarKinematicVariables(){}
    
    // Shape functions for contact pair
    Vector NMaster;
    Vector NSlave;
    Vector PhiLagrangeMultipliers;

    // Determinant of slave cell's jacobian
    double DetJSlave;
        
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This method initialized the operators
     */
    void Initialize()
    {
        // Shape functions
        NMaster                = ZeroVector(TNumNodes);
        NSlave                 = ZeroVector(TNumNodes);
        PhiLagrangeMultipliers = ZeroVector(TNumNodes);
        
        // Jacobian of slave
        DetJSlave = 0.0;
    }
    
    /**
     * This method prints the current operators
     */
    void print( )
    {
        KRATOS_WATCH( NSlave );
        KRATOS_WATCH( NMaster );
        KRATOS_WATCH( PhiLagrangeMultipliers );
        KRATOS_WATCH( DetJSlave );
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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class MortarKinematicVariables

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MORTAR_KINEMATIC_VARIABLES  defined 
