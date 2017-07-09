 
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

#if !defined(KRATOS_MORTAR_OPERATOR )
#define  KRATOS_MORTAR_OPERATOR

// System includes

// External includes

// Project includes
// #include "contact_structural_mechanics_application.h"
#include "contact_structural_mechanics_application_variables.h"
#include "custom_includes/mortar_kinematic_variables.h"

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
    
/** \brief MortarOperator
 * This is the definition of the mortar operator
 */

template< const unsigned int TNumNodes>
class MortarOperator
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarKinematicVariables<TNumNodes> KinematicVariables;
    
    /// Counted pointer of MortarOperator
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperator );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarOperator(){}
    
    virtual ~MortarOperator(){}
    
    // Mortar condition matrices - DOperator and MOperator
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> DOperator;
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> MOperator;

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
        // We initialize the D and M operators
        DOperator = ZeroMatrix(TNumNodes, TNumNodes);
        MOperator = ZeroMatrix(TNumNodes, TNumNodes);
    }
    
    /**
     * It calculates the mortar operators
     * @param rKinematicVariables: Corresponds with the kinematic variables
     * @param rIntegrationWeight: The corresponding integration weight
     */
    void CalculateMortarOperators(
        KinematicVariables& rKinematicVariables,
        const double& rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double DetjSlave = rKinematicVariables.DetjSlave; 
        const Vector Phi = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector N1  = rKinematicVariables.NSlave;
        const Vector N2  = rKinematicVariables.NMaster;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
            {
                const double phi = Phi[i_slave];
                
                DOperator(i_slave, j_slave) += DetjSlave * rIntegrationWeight * phi * N1[j_slave];
                MOperator(i_slave, j_slave) += DetjSlave * rIntegrationWeight * phi * N2[j_slave];
            }
        }
    }
    
    /**
     * It calculates the POperator (Inverse(D x M))
     */
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> ComputePOperator()
    {
        // We calculate the inverse of D operator
        double auxdet;
        const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> InvDOperator = MathUtils<double>::InvertMatrix<TNumNodes>(DOperator, auxdet);
        
        // We calculate the P operator
        const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> POperator = prod(InvDOperator, MOperator);
        
        return POperator;
    }
    
    /**
     * This method prints the current operators
     */
    void print() 
    {
        KRATOS_WATCH(DOperator);
        KRATOS_WATCH(MOperator);
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

}; // Class MortarOperatorWithDerivatives

template< const unsigned int TDim, const unsigned int TNumNodes, bool TFrictional>
class MortarOperatorWithDerivatives : MortarOperator<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarOperator<TNumNodes>                         BaseClassType;
    
    typedef MortarKinematicVariables<TNumNodes>          KinematicVariables;
    
    typedef DerivativeData<TDim, TNumNodes, TFrictional> DerivativeDataType;
    
    /// Counted pointer of MortarOperatorWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperatorWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarOperatorWithDerivatives(){}
    
    virtual ~MortarOperatorWithDerivatives(){}
    
    static const unsigned int Size2 = 2 * (TNumNodes * TDim);
       
    // D and M directional derivatives
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2> DeltaDOperator;
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, Size2> DeltaMOperator;

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
        BaseClassType::Initialize();
        
        // We initialize the D and M derivatives operators 
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaDOperator[i] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaDOperator[i + TNumNodes * TDim] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaMOperator[i] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaMOperator[i + TNumNodes * TDim] = ZeroMatrix(TNumNodes, TNumNodes);
        }
    }
    
    /**
     * It calculates the mortar operators
     * @param rKinematicVariables: Corresponds with the kinematic variables
     * @param rIntegrationWeight: The corresponding integration weight
     */
    void CalculateMortarOperators(
        KinematicVariables& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double& rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double DetjSlave = rKinematicVariables.DetjSlave; 
        const Vector Phi = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector N1  = rKinematicVariables.NSlave;
        const Vector N2  = rKinematicVariables.NMaster;
        
        // Derivatives
        constexpr unsigned int Size1 =     (TNumNodes * TDim);
        constexpr unsigned int Size2 = 2 * (TNumNodes * TDim);

        const array_1d<double, Size1> DeltaJSlave  = rDerivativeData.DeltaJSlave;
        const array_1d<array_1d<double, TNumNodes >, Size1> DeltaPhi = rDerivativeData.DeltaPhi;
        const array_1d<array_1d<double, TNumNodes >, Size2> DeltaN1  = rDerivativeData.DeltaN1;
        const array_1d<array_1d<double, TNumNodes >, Size2> DeltaN2  = rDerivativeData.DeltaN2;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            const double phi = Phi[i_slave];
            
            for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
            {
                const double n1  = N1[j_slave];
                const double n2  = N2[j_slave];
                
                BaseClassType::DOperator(i_slave, j_slave) += DetjSlave * rIntegrationWeight * phi * n1;
                BaseClassType::MOperator(i_slave, j_slave) += DetjSlave * rIntegrationWeight * phi * n2;
                
                for (unsigned int i = 0; i < TDim * TNumNodes; i++)
                {
                    DeltaDOperator[i](i_slave, j_slave) += DeltaJSlave[i] * rIntegrationWeight * phi* n1        
                                                    + DetjSlave * rIntegrationWeight * DeltaPhi[i][i_slave] * n1
                                                    + DetjSlave * rIntegrationWeight * phi* DeltaN1[i][j_slave];
                                                                                
                    DeltaMOperator[i](i_slave, j_slave) += DeltaJSlave[i] * rIntegrationWeight * phi* n2        
                                                    + DetjSlave * rIntegrationWeight * DeltaPhi[i][i_slave] * n2
                                                    + DetjSlave * rIntegrationWeight * phi* DeltaN2[i][j_slave];
                }
                for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; i++)
                {
                    DeltaDOperator[i](i_slave, j_slave) += DetjSlave * rIntegrationWeight * phi * DeltaN1[i][j_slave];
                                                                                
                    DeltaMOperator[i](i_slave, j_slave) += DetjSlave * rIntegrationWeight * phi * DeltaN2[i][j_slave];
                }
            }
        }
    }
    
    /**
     * This method prints the current operators
     */
    void print() 
    {
        BaseClassType::print();
        
//             for (unsigned int i = 0; i < TNumNodes * TDim; i++)
//             {
//                 KRATOS_WATCH(DeltaDOperator[i]);
//                 KRATOS_WATCH(DeltaMOperator[i]);
//             }
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

}; // Class MortarOperatorWithDerivatives

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MORTAR_OPERATOR  defined 
