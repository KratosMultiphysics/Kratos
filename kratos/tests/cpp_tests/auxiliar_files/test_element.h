//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_TEST_ELEMENT_H_INCLUDED )
#define  KRATOS_TEST_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{
    namespace Testing
    {
    ///@name Testing Globals
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
    ///@name Testing Classes
    ///@{

    /**
    * @class TestElement
    * @ingroup KratosCore
    * @brief This is simple test element
    * @details It is designed to create a simple LHS and RHS in order to test strategies, processes. etc.... This way the common interface of the elements/conditions can be used to minimize the difference between the actual implementation and the test
    * @author Vicente Mataix Ferrandiz
    */
    class TestElement
        : public Element
    {
    public:

        ///@name Type Definitions
        ///@{
        /// Counted pointer of TestElement
        KRATOS_CLASS_POINTER_DEFINITION( TestElement);
        ///@}

        ///@name  Enum's
        ///@{

        /**
        * @brief This enum is used in order of differentiante
        * @details If more implementations are added, add the corresponding enum
        */
        enum class ResidualType {LINEAR = 0, NON_LINEAR = 1, ARC_LENGTH = 2};

        ///@}
        ///@name Life Cycle
        ///@{

        /**
        * @brief Default constructor
        * @param NewId the ID of the new element
        * @param pGeometry the nodes of the new element
        * @param TheResidualType The problem to be solved (linear, non-linear, arc-length ...)
        */
        TestElement(IndexType NewId, GeometryType::Pointer pGeometry, const ResidualType TheResidualType = ResidualType::LINEAR);

        /**
        * @brief Default constructor
        * @param NewId The ID of the new element
        * @param pGeometry The nodes of the new element
        * @param pProperties The properties assigned to the new element
        * @param TheResidualType The problem to be solved (linear, non-linear, arc-length ...)
        */
        TestElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, const ResidualType TheResidualType = ResidualType::LINEAR);

        ///Copy constructor
        TestElement(TestElement const& rOther);

        /// Destructor.
        ~TestElement() override;

        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.
        TestElement& operator=(TestElement const& rOther);

        ///@}
        ///@name Operations
        ///@{

        /**
        * @brief Creates a new total lagrangian updated element pointer
        * @param NewId the ID of the new element
        * @param ThisNodes the nodes of the new element
        * @param pProperties the properties assigned to the new element
        * @return a Pointer to the new element
        */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /**
        * @brief Clones the selected element variables, creating a new one
        * @param NewId the ID of the new element
        * @param ThisNodes the nodes of the new element
        * @return a Pointer to the new element
        */
        Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

        //************* GETTING METHODS

        /**
        * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the list of the DoFs
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        /**
        * @brief Sets on rResult the ID's of the element degrees of freedom
        * @param rResult The vector containing the ids of the DoFs
        * @param rCurrentProcessInfo The current process info instance
        */
        void EquationIdVector(
            EquationIdVectorType& rResult,
            const ProcessInfo& rCurrentProcessInfo
            ) const override;

        /**
        * @brief Sets on rValues the nodal displacements
        * @param rValues The vector containing the dofs values
        * @param Step The time step computed (must be in the buffer)
        */
        void GetValuesVector(Vector& rValues, int Step = 0) const override;

        /**
        * @brief Sets on rValues the nodal velocities
        * @param rValues The vector containing the first derivatives
        * @param Step The time step computed (must be in the buffer)
        */
        void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

        /**
        * @brief Sets on rValues the nodal accelerations
        * @param rValues The vector containing the second derivatives
        * @param Step The time step computed (must be in the buffer)
        */
        void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

        //************* COMPUTING  METHODS

        /**
        * @brief This is called during the assembling process in order to calculate all elemental contributions to the global system matrix and the right hand side
        * @param rLeftHandSideMatrix The elemental left hand side matrix
        * @param rRightHandSideVector The elemental right hand side
        * @param rCurrentProcessInfo The current process info instance
        */

        void CalculateLocalSystem(
            MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        /**
        * @brief This calculates just the RHS
        * @param rLeftHandSideMatrix The elemental left hand side matrix
        * @param rRightHandSideVector The elemental right hand side
        * @param rCurrentProcessInfo The current process info instance
        */

        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        /**
        * @brief This calculates just the LHS
        * @param rLeftHandSideMatrix The elemental left hand side matrix
        * @param rRightHandSideVector The elemental right hand side
        * @param rCurrentProcessInfo The current process info instance
        */

        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        /**
        * @brief This is called during the assembling process in order to calculate the elemental mass matrix
        * @param rMassMatrix the elemental mass matrix
        * @param rCurrentProcessInfo the current process info instance
        */
        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        /**
        * @brief This is called during the assembling process in order
        * to calculate the elemental damping matrix
        * @param rDampingMatrix The elemental damping matrix
        * @param rCurrentProcessInfo The current process info instance
        */
        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            const ProcessInfo& rCurrentProcessInfo
            ) override;

        /**
        * @brief This function provides the place to perform checks on the completeness of the input.
        * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
        * @param rCurrentProcessInfo The current process info instance
        */
        int Check(const ProcessInfo& rCurrentProcessInfo) override;

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

        ResidualType mResidualType;

        ///@name Protected Operators
        ///@{

        TestElement() : Element()
        {
        }

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

        friend class Serializer;

        // A private default constructor necessary for serialization

        void save(Serializer& rSerializer) const override;

        void load(Serializer& rSerializer) override;

        ///@name Private Inquiry
        ///@{
        ///@}
        ///@name Un accessible methods
        ///@{
        ///@}

    }; // Class TestElement

    ///@}
    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    } // namespace Testing.
} // namespace Kratos.
#endif // KRATOS_TEST_ELEMENT_H_INCLUDED  defined
