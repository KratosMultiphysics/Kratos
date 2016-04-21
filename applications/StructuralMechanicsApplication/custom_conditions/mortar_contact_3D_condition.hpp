// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_MORTAR_CONTACT_3D_CONDITION_H_INCLUDED )
#define  KRATOS_MORTAR_CONTACT_3D_CONDITION_H_INCLUDED

// System includes

// External includes
#include "boost/smart_ptr.hpp"
#include <vector>

// Project includes
#include "includes/serializer.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "utilities/math_utils.h"
#include "custom_utilities/projection.h"

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class MortarContact3DCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of MortarContact3DCondition
    KRATOS_CLASS_POINTER_DEFINITION( MortarContact3DCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MortarContact3DCondition(IndexType NewId, GeometryType::Pointer pGeometry);
    MortarContact3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

    ///Copy constructor
    MortarContact3DCondition( MortarContact3DCondition const& rOther);

    /// Destructor.
    virtual ~MortarContact3DCondition();

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //************* STARTING - ENDING  METHODS

    /**
     * Called at the beginning of each solution step
     */
    void Initialize();

    /**
     * Called at the beginning of each iteration
     */
    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    /**
     * Creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

    /**
     * This is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental leTransverseGradientFt hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

   /**
    * This function provides a more general interface to the element.
    * it is designed so that rRHSvariables are passed TO the element
    * thus telling what is the desired output
    * @param rRightHandSideVectors: container for the desired RHS output
    * @param rRHSVariables: parameter describing the expected RHSs
    */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

   /**
    * Sets on rResult the ID's of the element degrees of freedom
    * @return rResult: The result vector with the ID's of the DOF
    * @param rCurrentProcessInfo: the current process info instance
    */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

    /**
     * Sets on ConditionalDofList the degrees of freedom of the considered element geometry
     * @return rConditionalDofList
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rConditionalDofList,ProcessInfo& rCurrentProcessInfo);

    /**
     * Get on rVariable a double Value
     */
    void GetValueOnIntegrationPoints( const Variable<double>& rVariable, 
				      std::vector<double>& rValues, 
				      const ProcessInfo& rCurrentProcessInfo );

    /**
     * Calculate a double Variable
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, 
				      std::vector<double>& rOutput, 
				      const ProcessInfo& rCurrentProcessInfo);
    
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

    /**
     * Energy variable for loads
     */
    double mEnergy; 

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

    friend class Serializer;

    // A private default constructor necessary for serialization
    MortarContact3DCondition() {};

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);

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

}; // Class MortarContact3DCondition

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_MORTAR_CONTACT_3D_CONDITION_H_INCLUDED  defined 