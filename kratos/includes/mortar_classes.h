//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_MORTAR_CLASSES )
#define  KRATOS_MORTAR_CLASSES

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "utilities/mortar_utilities.h"

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
    
    #if !defined(SHARED_POINTER_HASHER)
    #define SHARED_POINTER_HASHER
    template<class TSharedPointer>
    struct SharedPointerHasher
    {
        size_t operator()(const TSharedPointer& pCond) const
        {
            return reinterpret_cast<size_t>(pCond.get());
        }
    };
    #endif
    #if !defined(SHARED_POINTER_COMPARATOR)
    #define SHARED_POINTER_COMPARATOR
    template<class TSharedPointer>
    struct SharedPointerComparator
    {
        bool operator()(const TSharedPointer& first, const TSharedPointer& second) const
        {
            return first.get() == second.get();
        }
    };
    #endif
    
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
    Vector NMaster, NSlave, PhiLagrangeMultipliers;

    // Determinant of slave cell's jacobian
    double DetjSlave;
        
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
        DetjSlave = 0.0;
    }
    
    /**
     * This method prints the current operators
     */
    void print( )
    {
        KRATOS_WATCH( NSlave );
        KRATOS_WATCH( NMaster );
        KRATOS_WATCH( PhiLagrangeMultipliers );
        KRATOS_WATCH( DetjSlave );
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

template< const unsigned int TDim, const unsigned int TNumNodes>
class MortarKinematicVariablesWithDerivatives : public MortarKinematicVariables<TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef MortarKinematicVariables<TNumNodes> BaseClassType;    
    
    /// Counted pointer of MortarKinematicVariables
    KRATOS_CLASS_POINTER_DEFINITION( MortarKinematicVariablesWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarKinematicVariablesWithDerivatives(){}
    
    ~MortarKinematicVariablesWithDerivatives() override{}
  
    // Shape functions local derivatives for contact pair
    Matrix DNDeMaster, DNDeSlave;
    
    /*
    * Jacobians in current configuration on all integration points of slave segment
    * Only those two variables contain info on all GP
    * other variables contain info only on the currently-calculated GP
    */
    Matrix jSlave;
        
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
        
        // Shape functions local derivatives
        DNDeMaster = ZeroMatrix(TNumNodes, TDim - 1);
        DNDeSlave  = ZeroMatrix(TNumNodes, TDim - 1);
        
        // Jacobians on all integration points
        jSlave = ZeroMatrix(TDim, TDim - 1);
    }
    
    /**
     * This method prints the current operators
     */
    void print( )
    {
        BaseClassType::print();
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

}; // Class MortarKinematicVariablesWithDerivatives

/** \brief DerivativeData
 * This data will be used to compute the derivatives
 */
template< unsigned int TDim, unsigned int TNumNodes, bool TFrictional>
class DerivativeData
{
public:
    ///@name Type Definitions
    ///@{
    
    // Auxiliar types
//     typedef int zero[0]; // NOTE: Problems in Windows
    typedef array_1d<double, TNumNodes> type_1;
    typedef bounded_matrix<double, TNumNodes, TDim> type_2;
    typedef bounded_matrix<double, TNumNodes, TNumNodes> type_3;
    
    // Auxiliar sizes
    static const unsigned int size_1 =     (TNumNodes * TDim);
    static const unsigned int size_2 = 2 * (TNumNodes * TDim);
    
    ///@}
    ///@name Life Cycle
    ///@{

    DerivativeData(){}
    
    virtual ~DerivativeData(){}
    
    // The ALM parameters
    array_1d<double, TNumNodes> PenaltyParameter;
    double ScaleFactor;
    
    typename std::conditional< TFrictional,double,int >::type TangentFactor;
    
    // The normals of the nodes
    type_2 NormalMaster, NormalSlave;
    
    // Displacements and velocities
    type_2 X1, X2, u1, u2;
    
    typename std::conditional< TFrictional,type_2,int >::type u1old, u2old;
    
    // Derivatives    
    array_1d<double, size_1> DeltaDetjSlave;
    array_1d<type_1, size_1> DeltaPhi;
    array_1d<type_1, size_2> DeltaN1, DeltaN2;
    array_1d<type_2, size_1> DeltaNormalSlave, DeltaNormalMaster;
    
    // Ae
    type_3 Ae;
    
    // Derivatives Ae
    array_1d<type_3, size_1> DeltaAe;
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Initializer method 
     * @param SlaveGeometry: The geometry of the slave 
     * @param rCurrentProcessInfo: The process info from the system
     */
    
    void Initialize(
        const GeometryType& SlaveGeometry,
        const ProcessInfo& rCurrentProcessInfo
        )
    {        
        // The normals of the nodes
        NormalSlave = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry,  NORMAL);
        
        // Displacements and velocities of the slave       
        u1 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 0) 
           - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1);
        X1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(SlaveGeometry, false, 1);
        
        // We get the ALM variables
//         const double penalty_parameter = rCurrentProcessInfo[PENALTY];
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
//             PenaltyParameter[i] = penalty_parameter;
            PenaltyParameter[i] = SlaveGeometry[i].GetValue(PENALTY);
        }
        ScaleFactor = rCurrentProcessInfo[SCALE_FACTOR];
        
        // Derivatives 
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaPhi[i] = ZeroVector(TNumNodes);
            DeltaN1[i] = ZeroVector(TNumNodes);
            DeltaN1[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
            DeltaN2[i] = ZeroVector(TNumNodes);
            DeltaN2[i + TNumNodes * TDim] = ZeroVector(TNumNodes);
            DeltaNormalSlave[i]      = ZeroMatrix(TNumNodes, TDim);
        }
        
        #if (TFrictional == true)
            TangentFactor = rCurrentProcessInfo[TANGENT_FACTOR];
            u1old = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1) 
                  - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 2);
        #endif
    }
    
    /**
     * Initialize the DeltaAe components
     */
    
    void InitializeDeltaAeComponents()
    {
        // Ae
        Ae = ZeroMatrix(TNumNodes, TNumNodes);
        
        // Derivatives Ae
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaAe[i] = ZeroMatrix(TNumNodes, TNumNodes);
        }
    }

    /**
     * Updating the Master pair
     * @param  pCond: The pointer of the current master
     */
    
    virtual void UpdateMasterPair(const Condition::Pointer& pCond)
    {
        GeometryType MasterGeometry =  pCond->GetGeometry();
        
        NormalMaster = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry,  NORMAL);
        
        // Displacements, coordinates and normals of the master
        u2 = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 0)
           - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1);
        X2 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(MasterGeometry, false, 1);

        // Derivative of master's normal
        for (unsigned int i = 0; i < TNumNodes * TDim; i++)
        {
            DeltaNormalMaster[i] = ZeroMatrix(TNumNodes, TDim);
        }
        
        #if (TFrictional == true)
            u2old = MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1)
                  - MortarUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 2);
        #endif
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
    
};  // Class DerivativeData

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
    
    typedef DerivativeData<TDim, TNumNodes, TFrictional> DerivativeDataType;
    
    /// Counted pointer of MortarOperatorWithDerivatives
    KRATOS_CLASS_POINTER_DEFINITION( MortarOperatorWithDerivatives );
         
    ///@}
    ///@name Life Cycle
    ///@{

    MortarOperatorWithDerivatives(){}
    
    ~MortarOperatorWithDerivatives() override{}
    
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
    
/** @brief Custom Point container to be used by the mapper
 */
class ConditionMap : public std::unordered_map<Condition::Pointer, bool, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> >
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of ConditionMap
    KRATOS_CLASS_POINTER_DEFINITION( ConditionMap );

    typedef std::unordered_map<Condition::Pointer, bool, SharedPointerHasher<Condition::Pointer>, SharedPointerComparator<Condition::Pointer> > BaseType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ConditionMap(){}

    /// Destructor
    virtual ~ConditionMap(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * It removes one particular condition from the map
     * @param pCond: The condition to remove
     */     
    void RemoveCondition(Condition::Pointer pCond)
    {
        BaseType::iterator set = find(pCond);
        if(set != end())
        {
            erase(set);
        }
    }
    
    /**
     * It adds one new condition
     * @param pCond: The condition to set
     */
    void AddNewCondition(Condition::Pointer pCond)
    {
        insert({pCond, true}); // True by default when adding a new one
    }
    
    /**
     * It adds one new condition, as active
     * @param pCond: The condition to set
     */
    void AddNewActiveCondition(Condition::Pointer pCond)
    {
        insert({pCond, true});
    }
    
    /**
     * It adds one new condition, as inactive
     * @param pCond: The condition to set
     */
    void AddNewInactiveCondition(Condition::Pointer pCond)
    {
        insert({pCond, false});
    }
    
    /**
     * It sets one particular condition as active or not
     * @param pCond: The condition to set
     * @param Active: The flag, true if active, false otherwise
     */
    void SetActive(Condition::Pointer pCond, const bool Active = true)
    {
        BaseType::iterator set = find(pCond);
        if(set != end())
        {
            set->second = Active;
        }
    }
    
    /**
     * It checks if one particular condition is active
     * @param pCond: The condition to check
     * @return True if it is active, false otherwise
     */
    bool IsActive(Condition::Pointer pCond) const 
    {
        BaseType::const_iterator set = find(pCond);
        return (set->second);
    }
    
    /**
     * It checks if at least one pair is active
     * @return True if at least one pair is active, false otherwise
     */
    bool AtLeastOnePairActive()
    {
        for ( auto it = begin(); it != end(); ++it )
        {
            if (it->second == true)
            {
                return true;
            }
        }
        
        return false;
    }
    
    /**
     * Print the map information
     */
    void print()
    {
        for ( auto it = begin(); it != end(); ++it )
        {
            std::cout << "The condition " << (it->first)->Id() << " is ACTIVE: " << it->second;
            
            KRATOS_WATCH((it->first)->GetGeometry());
        }
    }

    void save( Serializer& rSerializer ) const
    {
        // TODO: Fill if necessary
    }

    void load( Serializer& rSerializer )
    {
        // TODO: Fill if necessary
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
}; // Class ConditionMap 

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MORTAR_CLASSES  defined 
