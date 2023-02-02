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

#if !defined(KRATOS_TEST_BAR_ELEMENT_H_INCLUDED )
#define  KRATOS_TEST_BAR_ELEMENT_H_INCLUDED

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
    * @class TestBarElement
    * @ingroup KratosCore
    * @brief This is simple test element
    * @details It is designed to create a simple LHS and RHS in order to test builder and solvers (a minimal connectivity of two nodes). This way the common interface of the elements/conditions can be used to minimize the difference between the actual implementation and the test
    * @author Vicente Mataix Ferrandiz
    */
    class TestBarElement
        : public Element
    {
    public:

        ///@name Type Definitions
        ///@{

        /// Counted pointer of TestBarElement
        KRATOS_CLASS_POINTER_DEFINITION( TestBarElement);

        ///@}

        ///@name  Enum's
        ///@{

        ///@}
        ///@name Life Cycle
        ///@{

        /**
        * @brief Default constructor
        * @param NewId the ID of the new element
        * @param pGeometry the nodes of the new element
        * @param TheResidualType The problem to be solved (linear, non-linear, arc-length ...)
        */
        TestBarElement(IndexType NewId, GeometryType::Pointer pGeometry)
            : Element( NewId, pGeometry )
        {
        }

        /**
        * @brief Default constructor
        * @param NewId The ID of the new element
        * @param pGeometry The nodes of the new element
        * @param pProperties The properties assigned to the new element
        * @param TheResidualType The problem to be solved (linear, non-linear, arc-length ...)
        */
        TestBarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
            : Element( NewId, pGeometry, pProperties )
        {
        }

        ///Copy constructor
        TestBarElement(TestBarElement const& rOther)
            :Element(rOther)
        {
        }

        /// Destructor.
        ~TestBarElement() override
        {
        }

        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.
        TestBarElement& operator=(TestBarElement const& rOther)
        {
            Element::operator=(rOther);

            return *this;
        }

        ///@}
        ///@name Operations
        ///@{

        /**
        * @brief Creates a new total lagrangian updated element pointer
        * @param NewId the ID of the new element
        * @param rThisNodes the nodes of the new element
        * @param pProperties the properties assigned to the new element
        * @return a Pointer to the new element
        */
        Element::Pointer Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const override
        {
            return Kratos::make_intrusive<TestBarElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
        }

        /**
        * @brief Clones the selected element variables, creating a new one
        * @param NewId the ID of the new element
        * @param rThisNodes the nodes of the new element
        * @return a Pointer to the new element
        */
        Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
        {
            TestBarElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

            return Kratos::make_intrusive<TestBarElement>(new_element);
        }

        //************* GETTING METHODS

        /**
        * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
        * @param rElementalDofList The vector containing the list of the DoFs
        * @param rCurrentProcessInfo The current process info instance
        */
        void GetDofList(
            DofsVectorType& rElementalDofList,
            const ProcessInfo& rCurrentProcessInfo
            ) const override
        {
            rElementalDofList.resize( 0 );

            rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
            rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
            rElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_X ) );
            rElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_Y ) );
        }

        /**
        * @brief Sets on rResult the ID's of the element degrees of freedom
        * @param rResult The vector containing the ids of the DoFs
        * @param rCurrentProcessInfo The current process info instance
        */
        void EquationIdVector(
            EquationIdVectorType& rResult,
            const ProcessInfo& rCurrentProcessInfo
            ) const override
        {
            if ( rResult.size() != 4 )
                rResult.resize( 4, false );

            rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[2] = GetGeometry()[1].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[3] = GetGeometry()[1].GetDof( DISPLACEMENT_Y ).EquationId();
        }

        /**
        * @brief Sets on rValues the nodal displacements
        * @param rValues The vector containing the dofs values
        * @param Step The time step computed (must be in the buffer)
        */
        void GetValuesVector(Vector& rValues, int Step = 0) const override
        {
            if ( rValues.size() != 4 )
                rValues.resize( 4, false );

            rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );
            rValues[2] = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[3] = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Y, Step );
        }

        /**
        * @brief Sets on rValues the nodal velocities
        * @param rValues The vector containing the first derivatives
        * @param Step The time step computed (must be in the buffer)
        */
        void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override
        {
            if ( rValues.size() != 4 )
                rValues.resize( 4, false );

            rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
            rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );
            rValues[2] = GetGeometry()[1].GetSolutionStepValue( VELOCITY_X, Step );
            rValues[3] = GetGeometry()[1].GetSolutionStepValue( VELOCITY_Y, Step );
        }

        /**
        * @brief Sets on rValues the nodal accelerations
        * @param rValues The vector containing the second derivatives
        * @param Step The time step computed (must be in the buffer)
        */
        void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override
        {
            if ( rValues.size() != 4 )
                rValues.resize( 4, false );

            rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
            rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );
            rValues[2] = GetGeometry()[1].GetSolutionStepValue( ACCELERATION_X, Step );
            rValues[3] = GetGeometry()[1].GetSolutionStepValue( ACCELERATION_Y, Step );
        }

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
            ) override
        {
            KRATOS_TRY;

            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

            // Compute LHS
            this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

            KRATOS_CATCH( "" );
        }

        /**
        * @brief This calculates just the RHS
        * @param rLeftHandSideMatrix The elemental left hand side matrix
        * @param rRightHandSideVector The elemental right hand side
        * @param rCurrentProcessInfo The current process info instance
        */

        void CalculateRightHandSide(
            VectorType& rRightHandSideVector,
            const ProcessInfo& rCurrentProcessInfo
            ) override
        {
            const unsigned int system_size = 4;

            if ( rRightHandSideVector.size() != system_size )
                rRightHandSideVector.resize( system_size, false );

            rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

            // We don't care about filling the RHS, we are just interested in the LHS for the test
        }

        /**
        * @brief This calculates just the LHS
        * @param rLeftHandSideMatrix The elemental left hand side matrix
        * @param rRightHandSideVector The elemental right hand side
        * @param rCurrentProcessInfo The current process info instance
        */

        void CalculateLeftHandSide(
            MatrixType& rLeftHandSideMatrix,
            const ProcessInfo& rCurrentProcessInfo
            ) override
        {
            const unsigned int system_size = 4;

            if ( rLeftHandSideMatrix.size1() != system_size )
                rLeftHandSideMatrix.resize( system_size, system_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

            const double E = this->GetProperties()[YOUNG_MODULUS];
            const double A = this->GetProperties()[NODAL_AREA];
            const double L = this->GetGeometry().Length();
            const double k = E*A/L;

            const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
            const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
            const double theta = std::atan2(dy, dx);
            const double l = std::cos(theta);
            const double m = std::sin(theta);

            rLeftHandSideMatrix(0,0) = k * l * l;
            rLeftHandSideMatrix(0,1) = k * m * l;
            rLeftHandSideMatrix(0,2) = -k * l * l;
            rLeftHandSideMatrix(0,3) = -k * m * l;

            rLeftHandSideMatrix(1,0) = k * l * m;
            rLeftHandSideMatrix(1,1) = k * m * m;
            rLeftHandSideMatrix(1,2) = -k * m * l;
            rLeftHandSideMatrix(1,3) = -k * m * m;

            rLeftHandSideMatrix(2,0) = -k * l * l;
            rLeftHandSideMatrix(2,1) = -k * m * l;
            rLeftHandSideMatrix(2,2) = k * l * l;
            rLeftHandSideMatrix(2,3) = k * m * l;

            rLeftHandSideMatrix(3,0) = -k * l * m;
            rLeftHandSideMatrix(3,1) = -k * m * m;
            rLeftHandSideMatrix(3,2) = k * m * l;
            rLeftHandSideMatrix(3,3) = k * m * m;
        }

        /**
        * @brief This is called during the assembling process in order to calculate the elemental mass matrix
        * @param rMassMatrix the elemental mass matrix
        * @param rCurrentProcessInfo the current process info instance
        */
        void CalculateMassMatrix(
            MatrixType& rMassMatrix,
            const ProcessInfo& rCurrentProcessInfo
            ) override
        {
            KRATOS_TRY

            unsigned int system_size = 4;
            if ( rMassMatrix.size1() != system_size || rMassMatrix.size2() != system_size) {
                rMassMatrix.resize( system_size, system_size, false );
            }

            rMassMatrix = ZeroMatrix( system_size, system_size );
            const double rho = this->GetProperties()[DENSITY];
            const double A = this->GetProperties()[NODAL_AREA];
            const double L = this->GetGeometry().Length();
            const double aux_val = A * rho * L / 6.0;

            const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
            const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
            const double l = dx / L;
            const double m = dy / L;

            rMassMatrix(0,0) = 2.0 * aux_val * l * l;
            rMassMatrix(0,1) = 2.0 * aux_val * l * m;
            rMassMatrix(0,2) = aux_val * l * l;
            rMassMatrix(0,3) = aux_val * l * m;

            rMassMatrix(1,0) = 2.0 * aux_val * l * m;
            rMassMatrix(1,1) = 2.0 * aux_val * m * m;
            rMassMatrix(1,2) = aux_val * l * m;
            rMassMatrix(1,3) = aux_val * m * m;

            rMassMatrix(2,0) = aux_val * l * l;
            rMassMatrix(2,1) = aux_val * l * m;
            rMassMatrix(2,2) = 2.0 * aux_val * l * l;
            rMassMatrix(2,3) = 2.0 * aux_val * l * m;

            rMassMatrix(3,0) = aux_val * l * m;
            rMassMatrix(3,1) = aux_val * m * m;
            rMassMatrix(3,2) = 2.0 * aux_val * l * m;
            rMassMatrix(3,3) = 2.0 * aux_val * m * m;

            KRATOS_CATCH( "" );
        }

        /**
        * @brief This is called during the assembling process in order
        * to calculate the elemental damping matrix
        * @param rDampingMatrix The elemental damping matrix
        * @param rCurrentProcessInfo The current process info instance
        */
        void CalculateDampingMatrix(
            MatrixType& rDampingMatrix,
            const ProcessInfo& rCurrentProcessInfo
            ) override
        {
            KRATOS_TRY;

            // Resizing as needed the LHS
            const unsigned int system_size = 4;

            rDampingMatrix = ZeroMatrix( system_size, system_size );

            KRATOS_CATCH( "" );
        }

        /**
        * this is called during the initialize of the builder
        * to calculate the lumped mass vector
        * @param rLumpedMassVector the elemental lumped mass vector
        * @param rCurrentProcessInfo the current process info instance
        */
        void CalculateLumpedMassVector(
            VectorType& rLumpedMassVector,
            const ProcessInfo& rCurrentProcessInfo
            ) const override
        {
            KRATOS_TRY;

            unsigned int system_size = 4;
            // Initialize the lumped mass vector
            if (rLumpedMassVector.size() != system_size) {
                rLumpedMassVector.resize(system_size, false);
            }

            const double rho = this->GetProperties()[DENSITY];
            const double A = this->GetProperties()[NODAL_AREA];
            const double L = this->GetGeometry().Length();
            const double aux_val = A * rho * L / 6.0;

            const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
            const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
            const double l = dx / L;
            const double m = dy / L;

            // Fill the lumped mass vector
            rLumpedMassVector(0) = 3.0 * aux_val * l * l + 3.0 * aux_val * l * m;
            rLumpedMassVector(1) = 3.0 * aux_val * l * m + 3.0 * aux_val * m * m;
            rLumpedMassVector(2) = 3.0 * aux_val * l * l + 3.0 * aux_val * l * m;
            rLumpedMassVector(3) = 3.0 * aux_val * l * m + 3.0 * aux_val * m * m;

            KRATOS_CATCH( "" );
        }

        /**
        * @brief This function provides the place to perform checks on the completeness of the input.
        * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
        * @param rCurrentProcessInfo The current process info instance
        */
        int Check(const ProcessInfo& rCurrentProcessInfo) const override
        {
            KRATOS_TRY

            KRATOS_ERROR_IF_NOT(GetProperties().Has(YOUNG_MODULUS)) << "YOUNG_MODULUS not defined" << std::endl;
            KRATOS_ERROR_IF_NOT(GetProperties().Has(NODAL_AREA)) << "NODAL_AREA not defined" << std::endl;

            // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
            for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
                const Node<3>& rnode = this->GetGeometry()[i];

                KRATOS_EXPECT_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
                KRATOS_EXPECT_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
                KRATOS_EXPECT_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

                KRATOS_EXPECT_DOF_IN_NODE(DISPLACEMENT_X,rnode)
                KRATOS_EXPECT_DOF_IN_NODE(DISPLACEMENT_Y,rnode)
                KRATOS_EXPECT_DOF_IN_NODE(DISPLACEMENT_Z,rnode)
            }

            return 0;

            KRATOS_CATCH( "Problem in the Check in the TestBarElement" )
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

    }; // Class TestBarElement

    void TestBarElement::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void TestBarElement::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    ///@}
    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Input and output
    ///@{
    ///@}
    } // namespace Testing.
} // namespace Kratos.
#endif // KRATOS_TEST_BAR_ELEMENT_H_INCLUDED  defined
