 
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
    
    boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes> Me, De;
        
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
        const Vector n1 = rKinematicVariables.NSlave;
        const double det_j = rKinematicVariables.DetjSlave; 
        
        De += rIntegrationWeight * (ComputeDe(n1, det_j));
        Me += rIntegrationWeight * det_j * outer_prod(n1, n1);
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

/** \brief DualLagrangeMultiplierOperatorsWithDerivatives
 * This is the definition dual lagrange multiplier operators with derivatives
 */

template< const unsigned int TDim, const unsigned int TNumNodes, bool TFrictional>
class DualLagrangeMultiplierOperatorsWithDerivatives : public DualLagrangeMultiplierOperators<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
        
    typedef DualLagrangeMultiplierOperators<TNumNodes> BaseClassType;  
    
    typedef MortarKinematicVariablesWithDerivatives<TDim, TNumNodes> KinematicVariables;
    
    typedef DerivativeData<TDim, TNumNodes, TFrictional> DerivativeDataType;
    
    /// Counted pointer of DualLagrangeMultiplierOperatorsWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( DualLagrangeMultiplierOperatorsWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    DualLagrangeMultiplierOperatorsWithDerivatives(){}
    
    ~DualLagrangeMultiplierOperatorsWithDerivatives() override{}
    
    // Auxiliar sizes
    static const unsigned int size_1 = (TNumNodes * TDim);
    
    // Derivatives matrices
    array_1d<boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>, size_1> DeltaMe;
    array_1d<boost::numeric::ublas::bounded_matrix<double, TNumNodes, TNumNodes>, size_1> DeltaDe;
        
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
        const Vector n1 = rKinematicVariables.NSlave;
        
        BaseClassType::CalculateAeComponents(rKinematicVariables, rIntegrationWeight);
        
        for (unsigned int i = 0; i < TDim * TNumNodes; i++)
        {
            const double delta_det_j = rDerivativeData.DeltaDetjSlave[i];
            
            DeltaDe[i] += rIntegrationWeight * this->ComputeDe( n1, delta_det_j );
            DeltaMe[i] += rIntegrationWeight * delta_det_j * outer_prod(n1, n1);
        }
    }
 
    /**
     * Calculates the matrix DeltaAe
     */
    bool CalculateDeltaAe(DerivativeDataType& rDerivativeData)
    {        
        double aux_det;
        const double tolerance = std::numeric_limits<double>::epsilon();
        
        // We compute the norm
        const double norm_Me = norm_frobenius(BaseClassType::Me);
        
        // Now we normalize the matrix
        const bounded_matrix<double, TNumNodes, TNumNodes> normalized_Me = BaseClassType::Me/norm_Me;
        
        // We compute the normalized inverse
        aux_det = MathUtils<double>::DetMat<TNumNodes>(normalized_Me);
        if (std::abs(aux_det) < tolerance)
        {
            return false;
        }
        
        const bounded_matrix<double, TNumNodes, TNumNodes> normalized_inv_Me = MathUtils<double>::InvertMatrix<TNumNodes>(normalized_Me, aux_det, tolerance); 
        
        // Now we compute the inverse
        const bounded_matrix<double, TNumNodes, TNumNodes> inv_Me = normalized_inv_Me/norm_Me;
        
        noalias(rDerivativeData.Ae) = prod(BaseClassType::De, inv_Me);
        
        static const unsigned int size_1 = (TNumNodes * TDim);
        array_1d<bounded_matrix<double, TNumNodes, TNumNodes> , size_1>& DeltaAe = rDerivativeData.DeltaAe;
        
        for (unsigned int i = 0; i < TDim * TNumNodes; i++)
        {
            DeltaAe[i] = DeltaDe[i] - prod(rDerivativeData.Ae, DeltaMe[i]);
            DeltaAe[i] = prod(rDerivativeData.DeltaAe[i], inv_Me);
    //         DeltaAe[i] = ZeroMatrix(TNumNodes, TNumNodes); // NOTE: Test with zero derivative
        }
        
        return true;
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

}; // Class DualLagrangeMultiplierOperatorsWithDerivatives

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_DUAL_LM_OPERATORS  defined 
