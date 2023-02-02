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
        TestElement(IndexType NewId, GeometryType::Pointer pGeometry, const ResidualType TheResidualType = ResidualType::LINEAR)
            : Element( NewId, pGeometry )
            , mResidualType( TheResidualType )
        {

        }

        /**
        * @brief Default constructor
        * @param NewId The ID of the new element
        * @param pGeometry The nodes of the new element
        * @param pProperties The properties assigned to the new element
        * @param TheResidualType The problem to be solved (linear, non-linear, arc-length ...)
        */
        TestElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, const ResidualType TheResidualType = ResidualType::LINEAR)
            : Element( NewId, pGeometry, pProperties )
            , mResidualType( TheResidualType )
        {

        }

        ///Copy constructor
        TestElement(TestElement const& rOther)
            :Element(rOther)
            ,mResidualType(rOther.mResidualType)
        {

        }

        /// Destructor.
        ~TestElement() override 
        {

        }

        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator.
        TestElement& operator=(TestElement const& rOther)
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
            return Kratos::make_intrusive<TestElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mResidualType );
        }

        /**
        * @brief Clones the selected element variables, creating a new one
        * @param NewId the ID of the new element
        * @param rThisNodes the nodes of the new element
        * @return a Pointer to the new element
        */
        Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override
        {
            TestElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mResidualType );

            return Kratos::make_intrusive<TestElement>(new_element);
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
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            rElementalDofList.resize( 0 );

            rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
            rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
            if( dimension == 3 )
                rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );
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
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            if ( rResult.size() != dimension )
                rResult.resize( dimension, false );

            rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
            if( dimension == 3)
                rResult[2] = GetGeometry()[0].GetDof( DISPLACEMENT_Z ).EquationId();
        }

        /**
        * @brief Sets on rValues the nodal displacements
        * @param rValues The vector containing the dofs values
        * @param Step The time step computed (must be in the buffer)
        */
        void GetValuesVector(Vector& rValues, int Step = 0) const override
        {
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            if ( rValues.size() != dimension )
                rValues.resize( dimension, false );

            rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
            rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dimension == 3 )
                rValues[2] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }

        /**
        * @brief Sets on rValues the nodal velocities
        * @param rValues The vector containing the first derivatives
        * @param Step The time step computed (must be in the buffer)
        */
        void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override
        {
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            if ( rValues.size() != dimension )
                rValues.resize( dimension, false );

            rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
            rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );

            if ( dimension == 3 )
                rValues[2] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Z, Step );
        }

        /**
        * @brief Sets on rValues the nodal accelerations
        * @param rValues The vector containing the second derivatives
        * @param Step The time step computed (must be in the buffer)
        */
        void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override
        {
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            if ( rValues.size() != dimension )
                rValues.resize( dimension, false );

            rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
            rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dimension == 3 )
                rValues[2] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Z, Step );
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

            // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
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
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            // Resizing as needed the RHS
            const unsigned int system_size = dimension;

            if ( rRightHandSideVector.size() != system_size )
                rRightHandSideVector.resize( system_size, false );

            rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

            const array_1d<double, 3 >& delta_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

            switch ( mResidualType )
            {
                case ResidualType::LINEAR:
                    for ( unsigned int j = 0; j < dimension; ++j )
                        rRightHandSideVector[j] -= delta_displacement[j] - 1.0;
                    break;
                case ResidualType::NON_LINEAR:
                    for ( unsigned int j = 0; j < dimension; ++j )
                        rRightHandSideVector[j] -= std::pow(delta_displacement[j], 2) - 1.0;
                    break;
                default:
                    KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
            }
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
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            // Resizing as needed the LHS
            const unsigned int system_size = dimension;

            if ( rLeftHandSideMatrix.size1() != system_size )
                rLeftHandSideMatrix.resize( system_size, system_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

            const array_1d<double, 3 >& delta_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

            switch ( mResidualType )
            {
                case ResidualType::LINEAR:
                    for ( unsigned int j = 0; j < dimension; ++j )
                        rLeftHandSideMatrix(j, j) += 1.0;
                    break;
                case ResidualType::NON_LINEAR:
                    for ( unsigned int j = 0; j < dimension; ++j )
                        rLeftHandSideMatrix(j, j) += delta_displacement[j] * 2;
                    break;
                default:
                    KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
            }
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

            //lumped
            unsigned int dimension = GetGeometry().WorkingSpaceDimension();
            unsigned int system_size = dimension;

            if ( rMassMatrix.size1() != system_size )
                rMassMatrix.resize( system_size, system_size, false );

            rMassMatrix = ZeroMatrix( system_size, system_size );

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

            //0.-Initialize the DampingMatrix:
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

            // Resizing as needed the LHS
            const unsigned int system_size = dimension;

            rDampingMatrix = ZeroMatrix( system_size, system_size );

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

            KRATOS_CATCH( "Problem in the Check in the TestElement" )
        }

        /**
         * @brief Get on rVariable Constitutive Law from the element
         * @param rVariable The variable we want to get
         * @param rValues The results in the integration points
         * @param rCurrentProcessInfo the current process info instance
         */
        void CalculateOnIntegrationPoints(
            const Variable<ConstitutiveLaw::Pointer>& rVariable,
            std::vector<ConstitutiveLaw::Pointer>& rValues,
            const ProcessInfo& rCurrentProcessInfo
            ) override
        {
            if (rVariable == CONSTITUTIVE_LAW) {
                const SizeType integration_points_number = mConstitutiveLawVector.size();
                if (rValues.size() != integration_points_number) {
                    rValues.resize(integration_points_number);
                }
                for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
                    rValues[point_number] = mConstitutiveLawVector[point_number];
                }
            }
        }

        void Initialize(const ProcessInfo& rCurrentProcessInfo) override
        {
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
            const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

            //Constitutive Law initialisation
            if ( mConstitutiveLawVector.size() != integration_points.size() )
                mConstitutiveLawVector.resize( integration_points.size() );

            const GeometryType& r_geometry = GetGeometry();
            const Properties& r_properties = GetProperties();
            const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
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

        IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods
        std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

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

    void TestElement::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void TestElement::load( Serializer& rSerializer )
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
#endif // KRATOS_TEST_ELEMENT_H_INCLUDED  defined
