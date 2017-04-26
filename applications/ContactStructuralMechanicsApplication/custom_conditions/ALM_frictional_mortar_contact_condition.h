// KRATOS  ___|  |       |       |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//           | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License: BSD License
//   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED )
#define  KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/ALM_mortar_contact_condition.h"

namespace Kratos 
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    
    typedef Point<3>                                  PointType;
    typedef Node<3>                                    NodeType;
    typedef Geometry<NodeType>                     GeometryType;
    typedef Geometry<PointType>               GeometryPointType;
    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod   IntegrationMethod;
    
///@}
///@name  Enum's
///@{
    
///@}
///@name  Functions
///@{
    
///@}
///@name Kratos Classes
///@{

/** \brief DerivativeData
 * This data will be used to compute the derivatives
 */
template< unsigned int TDim, unsigned int TNumNodes>
class FrictionalDerivativeData : public DerivativeData<TDim, TNumNodes>
{
public:
    
    // Auxiliar types
    typedef DerivativeData<TDim, TNumNodes>      BaseType;
    typedef bounded_matrix<double, TNumNodes, TDim> Type2;
    
    // Displacements and velocities
    Type2 u1old;
    Type2 u2old;
    
    // Default destructor
    ~FrictionalDerivativeData(){}
    
    /**
     * Initializer method 
     * @param  GeometryInput: The geometry of the slave 
     */
    
    void Initialize(const GeometryType& GeometryInput) override
    {        
        BaseType::Initialize(GeometryInput);
        
        u1old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(this->GetGeometry(), DISPLACEMENT, 1);
    }

    /**
     * Updating the Master pair
     * @param  pCond: The pointer of the current master
     */
    
    void UpdateMasterPair(const Condition::Pointer& pCond) override
    {
        BaseType::UpdateMasterPair(pCond);
        
        u2old = ContactUtilities::GetVariableMatrix<TDim,TNumNodes>(pCond->GetGeometry(), DISPLACEMENT, 1);
    }
};  // Class DerivativeData
    
/** \brief AugmentedLagrangianMethodFrictionalMortarContactCondition
 * This is a contact condition which employes the mortar method with dual lagrange multiplier 
 * The method has been taken from the Alexander Popps thesis:
 * Popp, Alexander: Mortar Methods for Computational Contact Mechanics and General Interface Problems, Technische Universität München, jul 2012
 */
template< unsigned int TDim, unsigned int TNumNodes>
class AugmentedLagrangianMethodFrictionalMortarContactCondition: public AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, true> 
{
public:
    ///@name Type Definitions
    ///@{
        
    /// Counted pointer of AugmentedLagrangianMethodFrictionalMortarContactCondition
    KRATOS_CLASS_POINTER_DEFINITION( AugmentedLagrangianMethodFrictionalMortarContactCondition );

    typedef AugmentedLagrangianMethodMortarContactCondition<TDim, TNumNodes, true>       BaseType;
    
    typedef typename BaseType::MortarConditionMatrices                    MortarConditionMatrices;
    
    typedef Condition                                                           ConditionBaseType;
    
    typedef typename ConditionBaseType::VectorType                                     VectorType;

    typedef typename ConditionBaseType::MatrixType                                     MatrixType;

    typedef typename ConditionBaseType::IndexType                                       IndexType;

    typedef typename ConditionBaseType::GeometryType::Pointer                 GeometryPointerType;

    typedef typename ConditionBaseType::NodesArrayType                             NodesArrayType;

    typedef typename ConditionBaseType::PropertiesType::Pointer             PropertiesPointerType;
    
    typedef typename ConditionBaseType::EquationIdVectorType                 EquationIdVectorType;
    
    typedef typename ConditionBaseType::DofsVectorType                             DofsVectorType;
    
    typedef typename std::vector<array_1d<PointType,TDim>>                 ConditionArrayListType;
    
    typedef Line2D2<Point<3>>                                                            LineType;
    
    typedef Triangle3D3<Point<3>>                                                    TriangleType;
    
    typedef typename std::conditional<TDim == 2, LineType, TriangleType >::type DecompositionType;
    
    typedef DerivativeData<TDim, TNumNodes>                                    DerivativeDataType;
    
    static constexpr unsigned int MatrixSize = TDim * (TNumNodes + TNumNodes + TNumNodes);
         
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    AugmentedLagrangianMethodFrictionalMortarContactCondition(): BaseType() 
    {
    }
    
    // Constructor 1
    AugmentedLagrangianMethodFrictionalMortarContactCondition(IndexType NewId, GeometryPointerType pGeometry):BaseType(NewId, pGeometry)
    {
    }
    
    // Constructor 2
    AugmentedLagrangianMethodFrictionalMortarContactCondition(IndexType NewId, GeometryPointerType pGeometry, PropertiesPointerType pProperties):BaseType( NewId, pGeometry, pProperties )
    {
    }

    ///Copy constructor
    AugmentedLagrangianMethodFrictionalMortarContactCondition( AugmentedLagrangianMethodFrictionalMortarContactCondition const& rOther)
    {
    }

    /// Destructor.
   virtual  ~AugmentedLagrangianMethodFrictionalMortarContactCondition();

    /**
     * Creates a new element pointer from an arry of nodes
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    
    Condition::Pointer Create( 
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesPointerType pProperties 
        ) const override;
    
    /**
     * Creates a new element pointer from an existing geometry
     * @param NewId: the ID of the new element
     * @param pGeom: the  geometry taken to create the condition
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    
    Condition::Pointer Create(
        IndexType NewId,
        GeometryPointerType pGeom,
        PropertiesPointerType pProperties
        ) const override;
        
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/

    /**
     * Sets on rResult the ID's of the element degrees of freedom
     * @return rResult: The result vector with the ID's of the DOF
     * @param rCurrentProcessInfo: the current process info instance
     */
    
    void EquationIdVector( 
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @return rConditionalDofList
     * @param rCurrentProcessInfo: the current process info instance
     */
    
    void GetDofList( 
        DofsVectorType& rConditionalDofList,
        ProcessInfo& rCurrentProcessInfo 
        ) override;

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
    
    // TODO: Define the "CL" or friction law to compute this
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /********************************************************************************/
    /**************** METHODS TO CALCULATE MORTAR CONDITION MATRICES ****************/
    /********************************************************************************/

    /*
     * Calculates the local contibution of the LHS
     */
    
    bounded_matrix<double, MatrixSize, MatrixSize> CalculateLocalLHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        ) override;
    
    /*
     * Calculates the local contibution of the LHS
     */
    
    array_1d<double, MatrixSize> CalculateLocalRHS(
        const MortarConditionMatrices& rMortarConditionMatrices,
        const DerivativeDataType& rDerivativeData,
        const double& rPenaltyFactor,
        const double& rScaleFactor,
        const unsigned int& rActiveInactive
        ) override;
        
    /******************************************************************/
    /********** AUXILLIARY METHODS FOR GENERAL CALCULATIONS ***********/
    /******************************************************************/
    
    /*
     * Returns a value depending of the active/inactive set
     */
    
    unsigned int GetActiveInactiveValue(GeometryType& CurrentGeometry) const override
    {
        unsigned int value = 0;
        
        for (unsigned int i_node = 0; i_node < CurrentGeometry.size(); i_node++)
        {
            if (CurrentGeometry[i_node].Is(ACTIVE) == true)
            {
                if (CurrentGeometry[i_node].Is(SLIP) == true)
                {
                    value += std::pow(3, i_node);
                }
                else
                {
                    value += 2 * std::pow(3, i_node);
                }
            }
        }
        
        return value;
    }
    
    /*
     * Returns a value depending of the active/inactive set
     */
    
    array_1d<double, TNumNodes> GetFrictionCoefficient()
    {
        // The friction coefficient
        const double mu = this->GetProperties().GetValue(FRICTION_COEFFICIENT);
        
        array_1d<double, TNumNodes> FrictionCoefficientVector(TNumNodes, mu);
        
        // TODO: Define the "CL" or friction law to compute this
        
        return FrictionCoefficientVector;
    }
    
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

    // Serialization 
    
    friend class Serializer;
    
    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
    }

    ///@}

}; // Class AugmentedLagrangianMethodFrictionalMortarContactCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_ALM_FRICTIONAL_MORTAR_CONTACT_CONDITION_H_INCLUDED  defined 
