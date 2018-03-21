//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on the work of Massimo Petracca and Peter Wilson)
//

#if !defined(KRATOS_BASE_SHELL_ELEMENT_H_INCLUDED)
#define KRATOS_BASE_SHELL_ELEMENT_H_INCLUDED


// System includes

// External includes


// Project includes
#include "includes/element.h"
#include "utilities/quaternion.h"
#include "custom_utilities/shell_cross_section.hpp"
#include "structural_mechanics_application_variables.h"


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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) BaseShellElement
    : public Element
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of BaseShellElement
    KRATOS_CLASS_POINTER_DEFINITION(BaseShellElement);

    typedef std::vector< ShellCrossSection::Pointer > CrossSectionContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Constructor using Geometry
    */
    BaseShellElement(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
    * Constructor using Properties
    */
    BaseShellElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
    * Destructor
    */
    ~BaseShellElement() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * ELEMENTS inherited from this class have to implement next
    * Create and Clone methods: MANDATORY
    */

    /**
    * this determines the elemental equation ID vector for all elemental
    * DOFs
    * @param rResult: the elemental equation ID vector
    * @param rCurrentProcessInfo: the current process info instance
    */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override;

    /**
    * determines the elemental list of DOFs
    * @param ElementalDofList: the list of DOFs
    * @param rCurrentProcessInfo: the current process info instance
    */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;


    void GetValuesVector(Vector& rValues, int Step = 0) override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

    void ResetConstitutiveLaw() override;

    /**
    * ELEMENTS inherited from this class have to implement next
    * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
    * they can be managed internally with a private method to do the same calculations
    * only once: MANDATORY
    */

//   /**
//    * this is called during the assembling process in order
//    * to calculate all elemental contributions to the global system
//    * matrix and the right hand side
//    * @param rLeftHandSideMatrix: the elemental left hand side matrix
//    * @param rRightHandSideVector: the elemental right hand side
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateLocalSystem(
//       MatrixType& rLeftHandSideMatrix,
//       VectorType& rRightHandSideVector,
//       ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental left hand side matrix only
//    * @param rLeftHandSideMatrix: the elemental left hand side matrix
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental right hand side vector only
//    * @param rRightHandSideVector: the elemental right hand side vector
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the first derivatives contributions for the LHS and RHS
//    * @param rLeftHandSideMatrix: the elemental left hand side matrix
//    * @param rRightHandSideVector: the elemental right hand side
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateFirstDerivativesContributions(
//       MatrixType& rLeftHandSideMatrix,
//       VectorType& rRightHandSideVector,
//       ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental left hand side matrix for the first derivatives constributions
//    * @param rLeftHandSideMatrix: the elemental left hand side matrix
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental right hand side vector for the first derivatives constributions
//    * @param rRightHandSideVector: the elemental right hand side vector
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);

//   /**
//    * ELEMENTS inherited from this class must implement this methods
//    * if they need to add dynamic element contributions
//    * note: second derivatives means the accelerations if the displacements are the dof of the analysis
//    * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
//    * CalculateSecondDerivativesContributions,
//    * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
//    */


//  /**
//    * this is called during the assembling process in order
//    * to calculate the second derivative contributions for the LHS and RHS
//    * @param rLeftHandSideMatrix: the elemental left hand side matrix
//    * @param rRightHandSideVector: the elemental right hand side
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateSecondDerivativesContributions(
//       MatrixType& rLeftHandSideMatrix,
//       VectorType& rRightHandSideVector,
//       ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental left hand side matrix for the second derivatives constributions
//    * @param rLeftHandSideMatrix: the elemental left hand side matrix
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateSecondDerivativesLHS(
//       MatrixType& rLeftHandSideMatrix,
//       ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental right hand side vector for the second derivatives constributions
//    * @param rRightHandSideVector: the elemental right hand side vector
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateSecondDerivativesRHS(
//       VectorType& rRightHandSideVector,
//       ProcessInfo& rCurrentProcessInfo);

//   /**
//    * this is called during the assembling process in order
//    * to calculate the elemental mass matrix
//    * @param rMassMatrix: the elemental mass matrix
//    * @param rCurrentProcessInfo: the current process info instance
//    */
//   virtual void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo);

  /**
   * this is called during the assembling process in order
   * to calculate the elemental damping matrix
   * @param rDampingMatrix: the elemental damping matrix
   * @param rCurrentProcessInfo: the current process info instance
   */
  void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;

//   /**
//    * This method provides the place to perform checks on the completeness of the input
//    * and the compatibility with the problem options as well as the contitutive laws selected
//    * It is designed to be called only once (or anyway, not often) typically at the beginning
//    * of the calculations, so to verify that nothing is missing from the input
//    * or that no common error is found.
//    * @param rCurrentProcessInfo
//    * this method is: MANDATORY
//    */
//   virtual int Check(const ProcessInfo& rCurrentProcessInfo);

    /**
    * returns the used integration method. In the general case this is the
    * default integration method of the used geometry. I an other integration
    * method is used the method has to be overwritten within the element
    * @return default integration method of the used Geometry
    */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mIntegrationMethod;
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

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

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

    SizeType mNumDofs;
    SizeType mNumGPs;
    IntegrationMethod mIntegrationMethod = GeometryData::GI_GAUSS_2;

    CrossSectionContainerType mSections; /*!< Container for cross section associated to each integration point */

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    /**
    * Protected empty constructor
    */
    BaseShellElement() : Element()
    {
    }

    void SetBaseMembers();

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    virtual void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    );

    void BaseInitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    void BaseFinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo);

    void BaseInitializeSolutionStep(ProcessInfo& rCurrentProcessInfo);

    void BaseFinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

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

}; // Class BaseShellElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_BASE_SHELL_ELEMENT_H_INCLUDED  defined
