 
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

#if !defined(KRATOS_DUAL_LM_OPERATORS )
#define  KRATOS_DUAL_LM_OPERATORS

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
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
    
/** \brief DualLagrangeMultiplierOperators
 * This is the definition dual lagrange multiplier operators
 */

template< const unsigned int TNumNodes>
class DualLagrangeMultiplierOperators
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef MortarKinematicVariables<TNumNodes> KinematicVariables;
    
    /// Counted pointer of DualLagrangeMultiplierOperators
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperators );
         
    ///@}
    ///@name Life Cycle
    ///@{

    DualLagrangeMultiplierOperators(){}
    
    virtual ~DualLagrangeMultiplierOperators(){}
    
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> Me;
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> De;
        
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
        // We initialize the De and Me operators
        Me = ZeroMatrix(TNumNodes, TNumNodes);
        De = ZeroMatrix(TNumNodes, TNumNodes);
    }
    
    /**
     * Calculates the Ae components necessary to compute the Phi_LagrangeMultipliers shape functions
     * @param rKinematicVariables: The kinematic variables
     * @param rIntegrationWeight: The integration weight considered
     */
    void CalculateAeComponents(
        KinematicVariables& rKinematicVariables,
        const double& rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const Vector N1 = rKinematicVariables.NSlave;
        const double Detj = rKinematicVariables.DetjSlave; 
        
        De += rIntegrationWeight * (ComputeDe(N1, Detj));
        Me += rIntegrationWeight * Detj * outer_prod(N1, N1);
    }
    
    /**
     * Calculates the matrix Ae
     */
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> CalculateAe()
    {        
        // We compute the norm
        const double norm_me = norm_frobenius(Me);
        
        // Now we normalize the matrix
        const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> normalized_me = Me/norm_me;
        
        // We compute the normalized inverse
        double aux_det;
        const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> normalized_inv_me = MathUtils<double>::InvertMatrix<TNumNodes>(normalized_me, aux_det, std::numeric_limits<double>::epsilon()); 
        
        // Now we compute the inverse
        const boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> inv_me = normalized_inv_me/norm_me;

        return prod(De, inv_me);
    }   
    
    /**
     * Calculates the matrix De
     * @param N1: The shape function 
     * @param detJ: The jacobian of the geometry 
     */
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> ComputeDe(        
        const Vector N1, 
        const double detJ 
        )
    {
        boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> De;
    
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            for (unsigned int j = 0; j < TNumNodes; j++)
            {
                if (i == j)
                {
                    De(i,i) = detJ * N1[i];
                }
                else
                {
                    De(i,j) = 0.0;
                }
            }
        }
        
        return De;
    }
    
    /**
     * This method prints the current operators
     */
    void print( )
    {
        KRATOS_WATCH( Me );
        KRATOS_WATCH( De );
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

}; // Class DualLagrangeMultiplierOperators

/** \brief DualLagrangeMultiplierOperatorsWithDertivatives
 * This is the definition dual lagrange multiplier operators with derivatives
 */

template< const unsigned int TDim, const unsigned int TNumNodes, bool TFrictional>
class DualLagrangeMultiplierOperatorsWithDertivatives : DualLagrangeMultiplierOperators<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef DualLagrangeMultiplierOperators<TNumNodes> BaseClassType;  
    
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes> KinematicVariables;
    
    typedef DerivativeData<TDim, TNumNodes, TFrictional> DerivativeDataType;
    
    /// Counted pointer of DualLagrangeMultiplierOperatorsWithDertivatives
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperatorsWithDertivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    DualLagrangeMultiplierOperatorsWithDertivatives(){}
    
    virtual ~DualLagrangeMultiplierOperatorsWithDertivatives(){}
    
    // Auxiliar sizes
    static const unsigned int Size1 = (TNumNodes * TDim);
    
    // Derivatives matrices
    array_1d<boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>, Size1> DeltaMe;
    array_1d<boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>, Size1> DeltaDe;
        
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
        
        // Derivatives matrices
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaMe[i] = ZeroMatrix(TNumNodes, TNumNodes);
            DeltaDe[i] = ZeroMatrix(TNumNodes, TNumNodes);
        }
    }
    
    /**
     * Calculates the Ae components and its derivatives necessary to compute the Phi_LagrangeMultipliers shape functions
     * @param rKinematicVariables: The kinematic variables
     * @param rIntegrationWeight: The integration weight considered
     */
    void CalculateDeltaAeComponents(
        KinematicVariables& rKinematicVariables,
        DerivativeDataType& rDerivativeData,
        const double& rIntegrationWeight
        )
    {
        /* DEFINITIONS */
        const Vector N1 = rKinematicVariables.NSlave;
        
        BaseClassType::CalculateAeComponents(rKinematicVariables, rIntegrationWeight);
        
        for (unsigned int i = 0; i < TDim * TNumNodes; i++)
        {
            const double DeltaDetJ = rDerivativeData.DeltaJSlave[i];
            
            DeltaDe[i] += rIntegrationWeight * this->ComputeDe( N1, DeltaDetJ );
            DeltaMe[i] += rIntegrationWeight * DeltaDetJ * outer_prod(N1, N1);
        }
    }
    
    /**
     * This method prints the current operators
     */
    void print( )
    {
        BaseClassType::print();
        
//         // Derivatives matrices
//         for (unsigned int i = 0; i < TNumNodes * TDim; i++)
//         {
//             KRATOS_WATCH( DeltaMe[i] );
//             KRATOS_WATCH( DeltaDe[i] );
//         }
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

}; // Class DualLagrangeMultiplierOperatorsWithDertivatives

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_DUAL_LM_OPERATORS  defined 
