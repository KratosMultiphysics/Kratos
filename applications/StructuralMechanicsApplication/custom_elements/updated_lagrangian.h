// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED


// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/base_solid_element.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

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
 * @class UpdatedLagrangian
 * @ingroup StructuralMechanicsApplication
 * @brief Updated Lagrangian element for 2D and 3D geometries.
 * @details Implements an Updated Lagrangian definition for structural analysis. This works for arbitrary geometries in 2D and 3D
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) UpdatedLagrangian
    : public BaseSolidElement
{
public:
    ///@name Type Definitions
    ///@{

    ///Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    ///Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    ///Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The base element type
    typedef BaseSolidElement BaseType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the sizetype
    typedef std::size_t SizeType;

    /// Counted pointer of UpdatedLagrangian
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UpdatedLagrangian);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry);
    UpdatedLagrangian(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    // Copy constructor
    UpdatedLagrangian(UpdatedLagrangian const& rOther)
        :BaseType(rOther)
        ,mF0Computed(rOther.mF0Computed)
        ,mDetF0(rOther.mDetF0)
        ,mF0(rOther.mF0)
    {};

    /// Destructor.
    ~UpdatedLagrangian() override;

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Called to initialize the element.
     * Must be called before any calculation is done
     */
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the beginning of each solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Called at the end of eahc solution step
     * @param rCurrentProcessInfo the current process info instance
     */
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(
        const Variable<Matrix >& rVariable,
        std::vector< Matrix >& rOutput,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a double Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValuesOnIntegrationPoints(
        const Variable<double>& rVariable,
        const std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     /**
      * @brief Set a Matrix Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValuesOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        const std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Updated Lagrangian Solid Element #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Updated Lagrangian Solid Element #" << Id() << "\nConstitutive law: " << BaseType::mConstitutiveLawVector[0]->Info();
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

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    /* Historical total elastic deformation measure */
    bool mF0Computed;           // To avoid computing more than once the historical total elastic deformation measure
    std::vector<double> mDetF0; // The historical total elastic deformation measure determinant
    std::vector<Matrix> mF0;    // The historical total elastic deformation measure

    ///@}
    ///@name Protected Operators
    ///@{

    UpdatedLagrangian() : BaseSolidElement()
    {
    }

    /**
     * @brief This method clones the element database
     * @param rF0Computed To avoid computing more than once the historical total elastic deformation measure
     * @param rDetF0 The historical total elastic deformation measure determinant
     * @param rF0 The historical total elastic deformation measure
     */
    void CloneUpdatedLagrangianDatabase(
        const bool rF0Computed,
        const std::vector<double>& rDetF0,
        const std::vector<Matrix>& rF0
        )
    {
        mF0Computed = rF0Computed;
        mDetF0 = rDetF0;
        mF0 = rF0;
    }

    /**
     * Gives the StressMeasure used
     */
    ConstitutiveLaw::StressMeasure GetStressMeasure() const override;

    /**
     * @brief It updates the historical database
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     */
    void UpdateHistoricalDatabase(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber
        );

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    void CalculateKinematicVariables(
        KinematicVariables& rThisKinematicVariables,
        const IndexType PointNumber,
        const GeometryType::IntegrationMethod& rIntegrationMethod
        ) override;

    /**
     * @brief This functions calculate the derivatives in the reference frame
     * @param J0 The jacobian in the reference configuration
     * @param InvJ0 The inverse of the jacobian in the reference configuration
     * @param DN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the reference configuration
     */
    double CalculateDerivativesOnReferenceConfiguration(
        Matrix& J0,
        Matrix& InvJ0,
        Matrix& DN_DX,
        const IndexType PointNumber,
        IntegrationMethod ThisIntegrationMethod
        ) const override;

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

    /**
     * @brief This method computes the deformation matrix B
     * @param rB The deformation matrix
     * @param rDN_DX The gradient derivative of the shape function
     * @param StrainSize The size of the Voigt notation stress vector
     * @param PointNumber The integration point considered
     */
    void CalculateB(
        Matrix& rB,
        const Matrix& rDN_DX,
        const SizeType StrainSize,
        const IndexType PointNumber
        );

    /**
     * It returns the reference configuration deformation gradient determinant
     * @param PointNumber The integration point considered
     * @return The reference configuration deformation gradient determinant
     */
    double ReferenceConfigurationDeformationGradientDeterminant(const IndexType PointNumber) const;

    /**
     * It returns the reference configuration deformation gradient
     * @param PointNumber The integration point considered
     * @return The reference configuration deformation gradient
     */
    Matrix ReferenceConfigurationDeformationGradient(const IndexType PointNumber) const;

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
    /// Assignment operator.
    //UpdatedLagrangian& operator=(const UpdatedLagrangian& rOther);
    /// Copy constructor.
    //UpdatedLagrangian(const UpdatedLagrangian& rOther);
    ///@}

}; // Class UpdatedLagrangian

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_UPDATED_LAGRANGIAN_H_INCLUDED  defined
