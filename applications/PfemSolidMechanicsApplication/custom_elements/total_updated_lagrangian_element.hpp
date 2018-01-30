//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_TOTAL_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_TOTAL_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/updated_lagrangian_element.hpp"


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

/// Total Updated Lagrangian Element for 3D and 2D geometries.

/**
 * Implements a updated Lagrangian definition for structural analysis.
 * This works for arbitrary geometries in 3D and 2D
 */

class TotalUpdatedLagrangianElement
    : public UpdatedLagrangianElement
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

    /// Counted pointer of TotalUpdatedLagrangianElement
    KRATOS_CLASS_POINTER_DEFINITION( TotalUpdatedLagrangianElement );
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    TotalUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry);

    TotalUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ///Copy constructor
    TotalUpdatedLagrangianElement(TotalUpdatedLagrangianElement const& rOther);

    /// Destructor.
    virtual ~TotalUpdatedLagrangianElement();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    TotalUpdatedLagrangianElement& operator=(TotalUpdatedLagrangianElement const& rOther);

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
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;


    /**
     * clones the selected element variables, creating a new one
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

    //************************************************************************************
    //************************************************************************************
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */
    //int Check(const ProcessInfo& rCurrentProcessInfo);


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
    TotalUpdatedLagrangianElement() : UpdatedLagrangianElement()
    {
    }


    /**
     * Initialize Element General Variables
     */
    void InitializeElementVariables(ElementVariables& rVariables, const ProcessInfo& rCurrentProcessInfo);


    /**
     * Transform Element General Variables
     */
    void TransformElementVariables(ElementVariables & rVariables,
				   const double& rPointNumber);


    /**
     * Calculate Element Kinematics
     */
    void CalculateKinematics(ElementVariables& rVariables,
                             const double& rPointNumber);


    /**
     * Get the Historical Deformation Gradient to calculate after finalize the step
     */
    void GetHistoricalVariables( ElementVariables& rVariables, 
				 const double& rPointNumber );

    /**
     * Calculation of the Deformation Matrix  BL
     */
    void CalculateDeformationMatrix(Matrix& rB,
                                    Matrix& rF,
                                    Matrix& rDN_DX);

    /**
     * Calculation of the Volume Change of the Element
     */
    virtual double& CalculateVolumeChange(double& rVolumeChange, ElementVariables& rVariables);

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

    virtual void save(Serializer& rSerializer) const;

    virtual void load(Serializer& rSerializer);


    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class TotalUpdatedLagrangianElement

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_TOTAL_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED  defined 
