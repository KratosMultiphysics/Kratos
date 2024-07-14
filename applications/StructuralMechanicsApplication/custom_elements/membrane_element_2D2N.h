// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/integration_utilities.h"

// Application includes
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

/**
 * @class MembraneElement2D2N
 * @ingroup StructuralMechanicsApplication
 * @brief Simplified two-dimensional membrane element
 * @details This elements implements a simplified linear elastic membrane model to be used in two-dimensional problems (see Gordnier (2009) and Smith and Shyy (1995))
 * @author Ruben Zorrilla
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) MembraneElement2D2N
    : public Element
{
protected:
    /**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        Vector  N;
        Matrix  B;
        double  detF;
        Matrix  F;
        double  detJ0;
        Matrix  J0;
        Matrix  InvJ0;
        Matrix  DN_DX;
        Vector Displacements;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         * @param Dimension The problem dimension: 2D or 3D
         * @param NumberOfNodes The number of nodes in the element
         */
        KinematicVariables(
            const SizeType StrainSize,
            const SizeType Dimension,
            const SizeType NumberOfNodes
            )
        {
            detF = 1.0;
            detJ0 = 1.0;
            N = ZeroVector(NumberOfNodes);
            B = ZeroMatrix(StrainSize, Dimension * NumberOfNodes);
            F = IdentityMatrix(Dimension);
            DN_DX = ZeroMatrix(NumberOfNodes, Dimension);
            J0 = ZeroMatrix(Dimension, Dimension);
            InvJ0 = ZeroMatrix(Dimension, Dimension);
            Displacements = ZeroVector(Dimension * NumberOfNodes);
        }
    };

    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        ConstitutiveLaw::StrainVectorType StrainVector;
        ConstitutiveLaw::StressVectorType StressVector;
        ConstitutiveLaw::VoigtSizeMatrixType D;

        /**
         * The default constructor
         * @param StrainSize The size of the strain vector in Voigt notation
         */
        ConstitutiveVariables(const SizeType StrainSize)
        {
            if (StrainVector.size() != StrainSize)
                StrainVector.resize(StrainSize);

            if (StressVector.size() != StrainSize)
                StressVector.resize(StrainSize);

            if (D.size1() != StrainSize || D.size2() != StrainSize)
                D.resize(StrainSize, StrainSize);

            noalias(StrainVector) = ZeroVector(StrainSize);
            noalias(StressVector) = ZeroVector(StrainSize);
            noalias(D)            = ZeroMatrix(StrainSize, StrainSize);
        }
    };
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///StressMeasure from constitutive laws
    typedef ConstitutiveLawType::StressMeasure StressMeasureType;

    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// This is the definition of the node.
    typedef Node NodeType;

    /// The base element type
    typedef Element BaseType;

    // Counted pointer of MembraneElement2D2N
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( MembraneElement2D2N );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    MembraneElement2D2N()
    {};

    // Constructor using an array of nodes
    MembraneElement2D2N(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    MembraneElement2D2N(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
        // This is needed to prevent uninitialised integration method in inactive elements
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    };

    // Copy constructor
    MembraneElement2D2N(MembraneElement2D2N const& rOther)
        : BaseType(rOther)
        , mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    {};

    // Destructor
    ~MembraneElement2D2N() override
    {};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Called to initialize the element.
     * @warning Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of eahc solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const &ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes) const override;

    /**
     * @brief Sets on rResult the ID's of the element degrees of freedom
     * @param rResult The vector containing the equation id
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Sets on rElementalDofList the degrees of freedom of the considered element geometry
     * @param rElementalDofList The vector containing the dof of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Returns the used integration method
     * @return default integration method of the used Geometry
     */
    IntegrationMethod GetIntegrationMethod() const override
    {
        return mThisIntegrationMethod;
    }

    /**
    * element can be integrated using the GP provided by the geometry or custom ones
    * by default, the base element will use the standard integration provided by the geom
    * @return bool to select if use/not use GPs given by the geometry
    */
    bool virtual UseGeometryIntegrationMethod() const
    {
        return true;
    }

    const virtual GeometryType::IntegrationPointsArrayType IntegrationPoints() const
    {
        return GetGeometry().IntegrationPoints();
    }

    const virtual GeometryType::IntegrationPointsArrayType IntegrationPoints(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().IntegrationPoints(ThisMethod);
    }

    const virtual Matrix& ShapeFunctionsValues(IntegrationMethod ThisMethod) const
    {
        return GetGeometry().ShapeFunctionsValues(ThisMethod);
    }

    /**
     * @brief Sets on rValues the nodal velocities
     * @param rValues The values of velocities
     * @param Step The step to be computed
     */
    void GetFirstDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /**
     * @brief Sets on rValues the nodal accelerations
     * @param rValues The values of accelerations
     * @param Step The step to be computed
     */
    void GetSecondDerivativesVector(
        Vector& rValues,
        int Step = 0) const override;

    /**
     * @brief This function provides a more general interface to the element.
     * @details It is designed so that rLHSvariables and rRHSvariables are passed to the element thus telling what is the desired output
     * @param rLeftHandSideMatrix container with the output Left Hand Side matrix
     * @param rRightHandSideVector container for the desired RHS output
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief This is called during the assembling process in order to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental right hand side vector only
      * @param rRightHandSideVector the elemental right hand side vector
      * @param rCurrentProcessInfo the current process info instance
      */
    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental mass matrix
      * @param rMassMatrix The elemental mass matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
      * @brief This is called during the assembling process in order to calculate the elemental damping matrix
      * @param rDampingMatrix The elemental damping matrix
      * @param rCurrentProcessInfo The current process info instance
      */
    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    const Parameters GetSpecifications() const override;

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MembraneElement2D2N #";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MembraneElement2D2N #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Sets the used integration method
     * @param ThisIntegrationMethod Integration method used
     */
    void SetIntegrationMethod(const IntegrationMethod &ThisIntegrationMethod)
    {
        mThisIntegrationMethod = ThisIntegrationMethod;
    }

    double GetMembranePrestress() const;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override;

    void load( Serializer& rSerializer ) override;

}; // class MembraneElement2D2N.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
