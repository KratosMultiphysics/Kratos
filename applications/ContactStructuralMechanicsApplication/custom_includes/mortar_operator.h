 
// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MORTAR_OPERATOR )
#define  KRATOS_MORTAR_OPERATOR

// System includes

// External includes

// Project includes
// #include "contact_classural_mechanics_application.h"
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

    MortarOperator()= default;
    
    virtual ~MortarOperator()= default;
    
    // Mortar condition matrices - DOperator and MOperator
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> DOperator, MOperator;

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
        const double det_j_slave = rKinematicVariables.DetjSlave; 
        const Vector phi_vector = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector n1_vector  = rKinematicVariables.NSlave;
        const Vector n2_vector  = rKinematicVariables.NMaster;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
            {
                const double phi = phi_vector[i_slave];
                
                DOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n1_vector[j_slave];
                MOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n2_vector[j_slave];
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
class MortarOperatorWithDerivatives : public MortarOperator<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarOperator<TNumNodes>                         BaseClassType;
    
    typedef MortarKinematicVariables<TNumNodes>          KinematicVariables;
    
    typedef typename std::conditional<TFrictional == true, DerivativeDataFrictional<TDim, TNumNodes>, DerivativeData<TDim, TNumNodes> >::type DerivativeDataType;
    
    /// Counted pointer of MortarOperatorWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperatorWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarOperatorWithDerivatives()= default;
    
    ~MortarOperatorWithDerivatives() override= default;
    
    static const unsigned int size_2 = 2 * (TNumNodes * TDim);
    
    // D and M directional derivatives
    array_1d<bounded_matrix<double, TNumNodes, TNumNodes>, size_2> DeltaDOperator, DeltaMOperator;

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
    void CalculateDeltaMortarOperators(
        KinematicVariables& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double& rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const double det_j_slave = rKinematicVariables.DetjSlave; 
        const Vector vector_phi = rKinematicVariables.PhiLagrangeMultipliers;
        const Vector vector_n1  = rKinematicVariables.NSlave;
        const Vector vector_n2  = rKinematicVariables.NMaster;
        
        // Derivatives
        constexpr unsigned int size_1 =     (TNumNodes * TDim);
        constexpr unsigned int size_2 = 2 * (TNumNodes * TDim);

        const array_1d<double, size_1> delta_j_slave  = rDerivativeData.DeltaDetjSlave;
        const array_1d<array_1d<double, TNumNodes >, size_1> delta_phi = rDerivativeData.DeltaPhi;
        const array_1d<array_1d<double, TNumNodes >, size_2> delta_n1  = rDerivativeData.DeltaN1;
        const array_1d<array_1d<double, TNumNodes >, size_2> delta_n2  = rDerivativeData.DeltaN2;
        
        for (unsigned int i_slave = 0; i_slave < TNumNodes; i_slave++)
        {
            const double phi = vector_phi[i_slave];
            
            for (unsigned int j_slave = 0; j_slave < TNumNodes; j_slave++)
            {
                const double n1 = vector_n1[j_slave];
                const double n2 = vector_n2[j_slave];
                
                BaseClassType::DOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n1;
                BaseClassType::MOperator(i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * n2;
                
                for (unsigned int i = 0; i < TDim * TNumNodes; i++)
                {
                    DeltaDOperator[i](i_slave, j_slave) += delta_j_slave[i] * rIntegrationWeight * phi* n1        
                                                    + det_j_slave * rIntegrationWeight * delta_phi[i][i_slave] * n1
                                                    + det_j_slave * rIntegrationWeight * phi* delta_n1[i][j_slave];
                                                                                
                    DeltaMOperator[i](i_slave, j_slave) += delta_j_slave[i] * rIntegrationWeight * phi* n2        
                                                    + det_j_slave * rIntegrationWeight * delta_phi[i][i_slave] * n2
                                                    + det_j_slave * rIntegrationWeight * phi* delta_n2[i][j_slave];
                }
                for (unsigned int i = TDim * TNumNodes; i < 2 * TDim * TNumNodes; i++)
                {
                    DeltaDOperator[i](i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * delta_n1[i][j_slave];
                                                                                
                    DeltaMOperator[i](i_slave, j_slave) += det_j_slave * rIntegrationWeight * phi * delta_n2[i][j_slave];
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
