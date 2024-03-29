//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                     PNavas $
//   Last modified by:    $Co-Author:               LMonforte $
//   Date:                $Date:                 October 2017 $
//   Revision:            $Revision:                     -0.1 $
//
//

#if !defined(KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_J_W_WP_DME_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_J_W_WP_DME_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_U_J_W_wP_element.hpp"

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

/// Updated Lagrangian Large Displacement Lagrangian U-W Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)


class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) AxisymUpdatedLagrangianUJWwPDMEElement
    : public AxisymUpdatedLagrangianUJWwPElement
{
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

    /// Counted pointer of LargeDisplacementUPElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( AxisymUpdatedLagrangianUJWwPDMEElement );
    ///@}

    ///@name Life Cycle
    ///@{

    /// Empty constructor needed for serialization
    AxisymUpdatedLagrangianUJWwPDMEElement();

    /// Default constructors
    AxisymUpdatedLagrangianUJWwPDMEElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AxisymUpdatedLagrangianUJWwPDMEElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    AxisymUpdatedLagrangianUJWwPDMEElement(AxisymUpdatedLagrangianUJWwPDMEElement const& rOther);


    /// Destructor.
    virtual ~AxisymUpdatedLagrangianUJWwPDMEElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    AxisymUpdatedLagrangianUJWwPDMEElement& operator=(AxisymUpdatedLagrangianUJWwPDMEElement const& rOther);


    ///@}
    ///@name Operations
    ///@{
    /**
     * Returns the currently selected integration method
     * @return current integration method selected
     */
    /**
     * creates a new total lagrangian updated element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

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


    /**
     * Part of the mass matrix due to the stabilization
     */
    virtual void CalculateAndAddMassStabilizationMatrix(MatrixType& rMassMatrix,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          ) override; 

    /**
     * Part of the damping matrix due to the stabilization
     */
    virtual void CalculateAndAddDampingStabilizationMatrix(MatrixType& rDampingMatrix,
          ElementDataType & rVariables,
          double& rIntegrationWeight
          ) override;

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

    virtual void save(Serializer& rSerializer) const override;

    virtual void load(Serializer& rSerializer) override;


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}


}; // Class UpdatedLagrangianUJWwPPastorLo Element



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_W_wP_PASTOR_ELEMENT_H_INCLUDED

