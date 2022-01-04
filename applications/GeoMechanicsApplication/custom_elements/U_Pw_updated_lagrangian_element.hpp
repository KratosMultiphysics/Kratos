// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
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
    using UPwBaseElement<TDim,TNumNodes>::mThisIntegrationMethod;

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

    ///@}
    ///@name Protected Operators
    ///@{


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


    // Copy constructor
    UPwUpdatedLagrangianElement(UPwUpdatedLagrangianElement const& rOther);

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
    }

    void load(Serializer& rSerializer) override
    {
        typedef UPwSmallStrainElement<TDim,TNumNodes> BaseClass;
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseClass );
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
