// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz,
//                   Vahid Galavi
//


#if !defined(KRATOS_GEO_U_PW_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED)
#define  KRATOS_GEO_U_PW_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED


// System includes


// External includes

// Project includes
#include "custom_elements/U_Pw_base_element.hpp"
#include "custom_elements/U_Pw_small_strain_element.hpp"
#include "custom_utilities/comparison_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "geo_mechanics_application_variables.h"


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
 * @class UPwUpdatedLagrangianElement
 * @brief Updated Lagrangian element for 2D and 3D geometries.
 * @details Implements an Updated Lagrangian definition for U-P elements. This works for arbitrary geometries in 2D and 3D
 * @author Vicente Mataix Ferrandiz (StructuralMechanics)
 * @author Vahid Galavi (Geomechanics)
 */
template< unsigned int TDim, unsigned int TNumNodes >
class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwUpdatedLagrangianElement
    : public UPwSmallStrainElement<TDim,TNumNodes>
{
public:
    ///@name Type Definitions
    ///@{
    typedef std::size_t IndexType;
    typedef Properties PropertiesType;
    typedef Node <3> NodeType;
    typedef Geometry<NodeType> GeometryType;
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
    typedef Vector VectorType;
    typedef Matrix MatrixType;

    /// Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /// The definition of the sizetype
    typedef std::size_t SizeType;
    using UPwBaseElement<TDim,TNumNodes>::mConstitutiveLawVector;
    using UPwBaseElement<TDim,TNumNodes>::mRetentionLawVector;
    using UPwBaseElement<TDim,TNumNodes>::mStressVector;
    using UPwBaseElement<TDim,TNumNodes>::mStateVariablesFinalized;
    using UPwBaseElement<TDim,TNumNodes>::CalculateDerivativesOnInitialConfiguration;

    typedef typename UPwSmallStrainElement<TDim,TNumNodes>::ElementVariables ElementVariables;
    using UPwSmallStrainElement<TDim,TNumNodes>::CalculateBulkModulus;

    /// Counted pointer of UPwUpdatedLagrangianElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwUpdatedLagrangianElement);

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    UPwUpdatedLagrangianElement(IndexType NewId = 0) : UPwSmallStrainElement<TDim,TNumNodes>( NewId ) {}

    /// Constructor using an array of nodes
    UPwUpdatedLagrangianElement(IndexType NewId,
                                const NodesArrayType& ThisNodes)
                                : UPwSmallStrainElement<TDim,TNumNodes>(NewId, ThisNodes) {}

    /// Constructor using Geometry
    UPwUpdatedLagrangianElement(IndexType NewId,
                                GeometryType::Pointer pGeometry)
                                : UPwSmallStrainElement<TDim,TNumNodes>(NewId, pGeometry) {}

    /// Constructor using Properties
    UPwUpdatedLagrangianElement(IndexType NewId,
                                GeometryType::Pointer pGeometry,
                                PropertiesType::Pointer pProperties)
                                : UPwSmallStrainElement<TDim,TNumNodes>( NewId, pGeometry, pProperties ) {}

    /// Destructor
    ~UPwUpdatedLagrangianElement() override {}


    int Check(const ProcessInfo& rCurrentProcessInfo) const override;


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
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param pGeom The pointer to the geometry of the element
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Creates a new element
     * @param NewId The Id of the new created element
     * @param ThisNodes The array containing nodes
     * @param pProperties The pointer to property
     * @return The pointer to the created element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    /**
     * @brief Calculate a double Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * @brief Calculate a Matrix Variable on the Element Constitutive Law
     * @param rVariable The variable we want to get
     * @param rOutput The values obtained int the integration points
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
                                      std::vector< Matrix >& rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

     /**
      * @brief Set a double Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
                                      const std::vector<double>& rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

     /**
      * @brief Set a Matrix Value on the Element Constitutive Law
      * @param rVariable The variable we want to set
      * @param rValues The values to set in the integration points
      * @param rCurrentProcessInfo the current process info instance
      */
    void SetValuesOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      const std::vector<Matrix>& rValues,
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Updated Lagrangian U-Pw Element #" << this->Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Updated Lagrangian U-Pw Element #" << this->Id() << "\nConstitutive law: " << mConstitutiveLawVector[0]->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        this->pGetGeometry()->PrintData(rOStream);
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
    // To avoid computing more than once the historical total elastic deformation measure
    bool mF0Computed;

    // The historical total elastic deformation measure determinant
    std::vector<double> mDetF0;

    // The historical total elastic deformation measure
    std::vector<Matrix> mF0;

    ///@}
    ///@name Protected Operators
    ///@{

    Matrix& CalculateDeltaDisplacement(Matrix& DeltaDisplacement) const;

    /**
     * @brief This method clones the element database
     * @param rF0Computed To avoid computing more than once the historical total elastic deformation measure
     * @param rDetF0 The historical total elastic deformation measure determinant
     * @param rF0 The historical total elastic deformation measure
     */
    void CloneUpdatedLagrangianDatabase(const bool rF0Computed,
                                        const std::vector<double>& rDetF0,
                                        const std::vector<Matrix>& rF0)
    {
        mF0Computed = rF0Computed;
        mDetF0 = rDetF0;
        mF0 = rF0;
    }

    /**
     * @brief It updates the historical database
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     */
    void UpdateHistoricalDatabase(ElementVariables& rThisKinematicVariables,
                                  const IndexType PointNumber);

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix The LHS
     * @param rRightHandSideVector The RHS
     * @param rCurrentProcessInfo The current process info instance
     * @param CalculateStiffnessMatrixFlag The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag The flag to set if compute the RHS
     */
    void CalculateAll(MatrixType& rLeftHandSideMatrix,
                      VectorType& rRightHandSideVector,
                      const ProcessInfo& rCurrentProcessInfo,
                      const bool CalculateStiffnessMatrixFlag,
                      const bool CalculateResidualVectorFlag) override;

    /**
     * @brief This functions updates the kinematics variables
     * @param rThisKinematicVariables The kinematic variables to be calculated
     * @param PointNumber The integration point considered
     * @param rIntegrationMethod The integration method considered
     */
    void CalculateKinematics( ElementVariables& rVariables, const unsigned int &PointNumber ) override;

    /**
     * @brief This functions calculate the derivatives in the reference frame
     * @param J0 The jacobian in the reference configuration
     * @param InvJ0 The inverse of the jacobian in the reference configuration
     * @param DN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the reference configuration
     */
    double CalculateDerivativesOnReferenceConfiguration(Matrix& J0,
                                                        Matrix& InvJ0,
                                                        Matrix& DN_DX,
                                                        const IndexType &PointNumber,
                                                        IntegrationMethod ThisIntegrationMethod) const;

    /**
     * @brief This functions calculate the derivatives in the current frame
     * @param rJ The jacobian in the current configuration
     * @param rInvJ The inverse of the jacobian in the current configuration
     * @param rDN_DX The gradient derivative of the shape function
     * @param PointNumber The id of the integration point considered
     * @param ThisIntegrationMethod The integration method considered
     * @return The determinant of the jacobian in the current configuration
     */
    double CalculateDerivativesOnCurrentConfiguration(Matrix& rJ,
                                                      Matrix& rInvJ,
                                                      Matrix& rDN_DX,
                                                      const IndexType &PointNumber,
                                                      IntegrationMethod ThisIntegrationMethod) const;

    void CalculateAndAddGeometricStiffnessMatrix( MatrixType& rLeftHandSideMatrix,
                                                  ElementVariables& rVariables,
                                                  unsigned int GPoint );

    void CalculateStrain( ElementVariables& rVariables ) override;

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


    // Copy constructor
    UPwUpdatedLagrangianElement(UPwUpdatedLagrangianElement const& rOther);
        // : UPwSmallStrainElement<TDim,TNumNodes>(rOther),
        // mF0Computed(rOther.mF0Computed),
        // mDetF0(rOther.mDetF0),
        // mF0(rOther.mF0) {}

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

    void save(Serializer& rSerializer) const override
    {
        typedef UPwSmallStrainElement<TDim,TNumNodes> BaseClass;
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseClass );
        rSerializer.save("F0Computed", mF0Computed);
        rSerializer.save("DetF0", mDetF0);
        rSerializer.save("F0", mF0);
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPwSmallStrainElement<TDim,TNumNodes> BaseClass;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseClass );
        rSerializer.load("F0Computed", mF0Computed);
        rSerializer.load("DetF0", mDetF0);
        rSerializer.load("F0", mF0);
    }


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    /// Assignment operator.
    //UPwUpdatedLagrangianElement& operator=(const UPwUpdatedLagrangianElement& rOther);
    /// Copy constructor.
    //UPwUpdatedLagrangianElement(const UPwUpdatedLagrangianElement& rOther);
    ///@}

}; // Class UPwUpdatedLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_GEO_U_PW_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED  defined
