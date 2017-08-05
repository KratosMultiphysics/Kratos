 
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
// #include "utilities/math_utils.h"
#include "custom_utilities/contact_utilities.h"

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

    MortarKinematicVariables()= default;
    
    virtual ~MortarKinematicVariables()= default;
    
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

    MortarKinematicVariablesWithDerivatives()= default;
    
    ~MortarKinematicVariablesWithDerivatives() override= default;
  
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
template< unsigned int TDim, unsigned int TNumNodes>
class DerivativeData
{
public:
    ///@name Type Definitions
    ///@{
    
    // Auxiliar types
    typedef array_1d<double, TNumNodes>                  type_1;
    typedef bounded_matrix<double, TNumNodes, TDim>      type_2;
    typedef bounded_matrix<double, TNumNodes, TNumNodes> type_3;
    typedef bounded_matrix<double, 3, 3>                 type_4;
    
    // Auxiliar sizes
    static const unsigned int size_1 =     (TNumNodes * TDim);
    static const unsigned int size_2 = 2 * (TNumNodes * TDim);
    
    ///@}
    ///@name Life Cycle
    ///@{

    DerivativeData()= default;
    
    virtual ~DerivativeData()= default;
    
    // The ALM parameters
    array_1d<double, TNumNodes> PenaltyParameter;
    double ScaleFactor;
    
    // The normals of the nodes
    type_2 NormalMaster, NormalSlave;
    
    // Displacements and velocities
    type_2 X1, X2, u1, u2;
    
    // Derivatives    
    array_1d<double, size_1> DeltaDetjSlave;
    array_1d<type_1, size_1> DeltaPhi;
    array_1d<type_1, size_2> DeltaN1, DeltaN2;
    array_1d<type_2, size_1> DeltaNormalSlave, DeltaNormalMaster;
    array_1d<type_4, size_2> DeltaCellVertex;
    
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
    
    virtual void Initialize(
        const GeometryType& SlaveGeometry,
        const ProcessInfo& rCurrentProcessInfo
        )
    {        
        // The normals of the nodes
        NormalSlave = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry,  NORMAL);
        
        // Displacements and velocities of the slave       
        u1 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 0) 
           - ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1);
        X1 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(SlaveGeometry, false, 1);
        
        // We get the ALM variables
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            PenaltyParameter[i] = SlaveGeometry[i].GetValue(PENALTY_PARAMETER);
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
            DeltaNormalSlave[i] = ZeroMatrix(TNumNodes, TDim);
        }
    
        if (TDim == 3)
        {
            for (unsigned int i = 0; i < 2 * TNumNodes * TDim; i++)
            {
                DeltaCellVertex[i] = ZeroMatrix(3, 3);
            }
        }
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
        
        NormalMaster = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry,  NORMAL);
        
        // Displacements, coordinates and normals of the master
        u2 = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 0)
           - ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1);
        X2 = ContactUtilities::GetCoordinates<TDim,TNumNodes>(MasterGeometry, false, 1);

        // Derivative of master's normal
        if (TDim == 2)
        {
            for (unsigned int i = 0; i < TNumNodes * TDim; i++)
            {
                DeltaNormalMaster[i] = ZeroMatrix(TNumNodes, TDim);
            }
        }
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

/** \brief DerivativeData
 * This data will be used to compute the derivatives
 */
template< unsigned int TDim, unsigned int TNumNodes>
class DerivativeDataFrictional : public DerivativeData<TDim, TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    
    // Auxiliar types
    typedef DerivativeData<TDim, TNumNodes>           BaseClass;
    typedef array_1d<double, TNumNodes>                  type_1;
    typedef bounded_matrix<double, TNumNodes, TDim>      type_2;
    typedef bounded_matrix<double, TNumNodes, TNumNodes> type_3;
    typedef bounded_matrix<double, 3, 3>                 type_4;
    
    // Auxiliar sizes
    static const unsigned int size_1 =     (TNumNodes * TDim);
    static const unsigned int size_2 = 2 * (TNumNodes * TDim);
    
    ///@}
    ///@name Life Cycle
    ///@{

    DerivativeDataFrictional()= default;
    
    virtual ~DerivativeDataFrictional()= default;
    
    // The ALM parameters
    double TangentFactor;
    
    // Displacements and velocities
    type_2 u1old, u2old;
    
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
        ) override
    {        
        BaseClass::Initialize(SlaveGeometry, rCurrentProcessInfo);
        
        TangentFactor = rCurrentProcessInfo[TANGENT_FACTOR];
        u1old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 1) 
                - ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(SlaveGeometry, DISPLACEMENT, 2);
    }
    
    /**
     * Updating the Master pair
     * @param  pCond: The pointer of the current master
     */
    
    void UpdateMasterPair(const Condition::Pointer& pCond) override
    {
        BaseClass::UpdateMasterPair(pCond);
        
        GeometryType MasterGeometry =  pCond->GetGeometry();
        
        u2old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 1)
                  - ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(MasterGeometry, DISPLACEMENT, 2);
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
    
};  // Class DerivativeDataFrictional

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_MORTAR_KINEMATIC_VARIABLES  defined 
