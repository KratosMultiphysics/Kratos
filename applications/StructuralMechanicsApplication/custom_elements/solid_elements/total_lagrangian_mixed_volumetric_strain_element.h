// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//                   Guglielmo Scovazzi
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
 * @class TotalLagrangianMixedVolumetricStrainElement
 * @ingroup StructuralMechanicsApplication
 * @brief Total Lagrangian mixed displacement volumetric strain formulation element
 * @details This implements an ASGS stabilized total Lagrangian element mixed formulation.
 * The nodal DOFs are the displacement and the volumetric strain, which is defined as 1 minus
 * the Jacobian determinant. This is a strain-driven formulation, meaning that it can be used
 * in combination with any constitutive model, including (almost) incompressible materials.
 * Note that the ASGS stabilization is specifically designed for linear triangle and tetrahedra.
 * Also note that the volumetric strain definition in the formulation implies that the plane
 * stress case cannot be considered in 2D.
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 * @author Guglielmo Scovazzi
 */
template<std::size_t TDim>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) TotalLagrangianMixedVolumetricStrainElement
    : public Element
{

    static constexpr std::size_t NumNodes = TDim + 1;
    static constexpr std::size_t StrainSize = TDim == 2 ? 3 : 6;
    static constexpr std::size_t BlockSize = TDim +1;
    static constexpr std::size_t LocalSize = NumNodes*BlockSize;

protected:

    /**
     * Internal variables used in the kinematic calculations
     */
    struct KinematicVariables
    {
        double detF;
        Matrix F;
        double detJ0;
        Matrix J0;
        Matrix InvJ0;
        Vector N;
        Matrix DN_DX;
        BoundedMatrix<double,NumNodes, TDim> Displacements;
        BoundedVector<double,NumNodes> JacobianDeterminant;
        Vector EquivalentStrain;

        KinematicVariables()
        {
            detF = 1.0;
            detJ0 = 1.0;
            F = IdentityMatrix(TDim);
            J0 = ZeroMatrix(TDim, TDim);
            InvJ0 = ZeroMatrix(TDim, TDim);
            N = ZeroVector(NumNodes);
            DN_DX = ZeroMatrix(NumNodes, TDim);
            Displacements = ZeroMatrix(NumNodes, TDim);
            JacobianDeterminant = ZeroVector(NumNodes);
            EquivalentStrain = ZeroVector(StrainSize);
        }
    };

    /**
     * Internal variables used in the kinematic calculations
     */
    struct ConstitutiveVariables
    {
        Vector StrainVector;
        Vector StressVector;
        Matrix ConstitutiveMatrix;

        /**
         * The default constructor
         */
        ConstitutiveVariables()
        {
            StrainVector = ZeroVector(StrainSize);
            StressVector = ZeroVector(StrainSize);
            ConstitutiveMatrix = ZeroMatrix(StrainSize, StrainSize);
        }
    };

public:

    ///@name Type Definitions
    ///@{

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

    /// The definition of nodes container
    using NodesArrayType = BaseType::NodesArrayType;

    /// The definition of the properties pointer
    using PropertiesType = BaseType::PropertiesType;

    /// The definition of the index type
    using IndexType = BaseType::IndexType;

    /// The definition of the sizetype
    using SizeType = BaseType::SizeType;

    // Counted pointer of TotalLagrangianMixedVolumetricStrainElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( TotalLagrangianMixedVolumetricStrainElement );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor void
    TotalLagrangianMixedVolumetricStrainElement()
    {};

    // Constructor using an array of nodes
    TotalLagrangianMixedVolumetricStrainElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {};

    // Constructor using an array of nodes with properties
    TotalLagrangianMixedVolumetricStrainElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {};

    // Copy constructor
    TotalLagrangianMixedVolumetricStrainElement(TotalLagrangianMixedVolumetricStrainElement const& rOther)
        : BaseType(rOther)
        , mThisIntegrationMethod(rOther.mThisIntegrationMethod)
        , mConstitutiveLawVector(rOther.mConstitutiveLawVector)
    {};

    // Destructor
    ~TotalLagrangianMixedVolumetricStrainElement() override
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
    void Initialize(const ProcessInfo &rCurrentProcessInfo) override;

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of each solution step
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
        NodesArrayType const& ThisNodes,
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
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo the current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate a Vector Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Small Displacement Mixed Strain Element #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()) {
          buffer << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          buffer << "\nNo constitutive law.";
        }
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Small Displacement Mixed Strain Element #" << Id();
        if (!BaseType::mConstitutiveLawVector.empty()) {
          rOStream << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        } else {
          rOStream << "\nNo constitutive law.";
        }
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    Vector mMinShearModulusVector; /// Gauss pt. minimum historical shear modulus

    IntegrationMethod mThisIntegrationMethod; /// Currently selected integration methods

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector; /// The vector containing the constitutive laws

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
    void SetIntegrationMethod(const IntegrationMethod& ThisIntegrationMethod)
    {
        mThisIntegrationMethod = ThisIntegrationMethod;
    }

    /**
     * @brief Sets the used constitutive laws
     * @param ThisConstitutiveLawVector Constitutive laws used
     */
    void SetConstitutiveLawVector(const std::vector<ConstitutiveLaw::Pointer>& ThisConstitutiveLawVector)
    {
        mConstitutiveLawVector = ThisConstitutiveLawVector;
    }

    /**
     * @brief It initializes the material
     */
    virtual void InitializeMaterial();

    /**
     * @brief This method returns if the element provides the strain
     */
    virtual bool UseElementProvidedStrain() const;

    /**
     * @brief This functions updates the data structure passed to the CL
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     */
    virtual void SetConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const;

    /**
     * @brief This functions updates the constitutive variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param rThisConstitutiveVariables The constitutive variables
     * @param rValues The CL parameters
     * @param PointNumber The integration point considered
     * @param IntegrationPoints The list of integration points
     * @param ThisStressMeasure The stress measure considered
     */
    virtual void CalculateConstitutiveVariables(
        KinematicVariables& rThisKinematicVariables,
        ConstitutiveVariables& rThisConstitutiveVariables,
        ConstitutiveLaw::Parameters& rValues,
        const IndexType PointNumber,
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
        const ConstitutiveLaw::StressMeasure ThisStressMeasure = ConstitutiveLaw::StressMeasure_PK2) const;

    /**
     * @brief Calculate the kinematics
     * @details This method calculates the kinematics of the element for a given integration point
     * @param rThisKinematicVariables Integration point kinematics container
     * @param PointNumber Integration point index
     * @param rIntegrationMethod Integration rule
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod) const;

    /**
     * @brief Calculate the equivalent strain
     * This function computes the equivalent strain vector.
     * The equivalent strain is defined from a multiplicative decomposition of
     * the strain. The deviatoric part comes from the nodal interpolation of the
     * displacements while the volumetric part comes from the nodal interpolation
     * of the nodal Jacobian determinant.
     * @param rThisKinematicVariables Kinematic variables container
     */
    void CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const;

    /**
     * @brief Calculation of the deformation gradient F
     * The equivalent deformation gradient is defined from a multiplicative decomposition of
     * the strain. The deviatoric part comes from the nodal interpolation of the
     * displacements while the volumetric part comes from the nodal interpolation
     * of the nodal Jacobian determinant.
     * @param rThisKinematicVariables Kinematic variables container
     */
    void CalculateEquivalentF(KinematicVariables& rThisKinematicVariables) const;

    /**
     * @brief Returns the shear modulus from the given C matrix
     * This method calculates the shear modulus (2nd Lame parameter) for the given constitutive
     * matrix in Voigt notation.
     * @param rC The constitutive matrix in Voigt notation
     * @return double The shear modulus
     */
    double CalculateShearModulus(const Matrix &rC) const;

    /**
     * @brief Returns the bulk modulus from the given C matrix
     * This method calculates the bulk modulus for the given constitutive matrix in Voigt notation.
     * @param rC The constitutive matrix in Voigt notation
     * @return double The bulk modulus
     */
    double CalculateBulkModulus(const Matrix &rC) const;

    /**
     * @brief This method gets a value directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @tparam TType The type considered
     */
    template<class TType>
    void GetValueOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput)
    {
        const auto& r_geometry = GetGeometry();
        const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());

        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            mConstitutiveLawVector[i_gauss]->GetValue(rVariable, rOutput[i_gauss]);
        }
    }

    /**
     * @brief This method computes directly in the CL
     * @details Avoids code repetition
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo the current process info instance
     * @tparam TType The type considered
     */
    template<class TType>
    void CalculateOnConstitutiveLaw(
        const Variable<TType>& rVariable,
        std::vector<TType>& rOutput,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const auto& r_geometry = GetGeometry();
        const SizeType n_nodes = r_geometry.size();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
        const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables;
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node, d) = r_disp[d];
            }
            kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Create the constitutive variables and values containers
        ConstitutiveVariables constitutive_variables;
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Call the initialize material response
        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Recompute the kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Set the constitutive variables
            SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

            // Calculate the output value
            rOutput[i_gauss] = mConstitutiveLawVector[i_gauss]->CalculateValue(cons_law_values, rVariable, rOutput[i_gauss] );
        }
    }

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

}; // class TotalLagrangianMixedVolumetricStrainElement.

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.
